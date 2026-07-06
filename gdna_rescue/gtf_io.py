"""GTF parsing, per-chromosome interval indexing and annotation masking.

We use a small hand-rolled GTF parser rather than gffutils so the tool has no
mandatory database-build step and stays fast for the one thing we need:
building per-chromosome masks and a gene index for context labelling.

All coordinates are converted to 0-based half-open on the way in.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from .utils import merge_runs_with_gap

_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


@dataclass
class Gene:
    """Minimal gene record used for nearest-feature / context labelling."""

    gene_id: str
    gene_name: str
    chrom: str
    start: int          # 0-based
    end: int            # half-open
    strand: str


@dataclass
class Annotation:
    """Parsed annotation: per-chrom merged mask intervals + a gene index."""

    # chrom -> list of merged [start, end) intervals used for masking.
    mask_intervals: Dict[str, List[Tuple[int, int]]] = field(default_factory=dict)
    # Per-strand masks for stranded libraries (antisense-over-feature stays
    # discoverable because only the feature's own strand is masked).
    mask_intervals_plus: Dict[str, List[Tuple[int, int]]] = field(default_factory=dict)
    mask_intervals_minus: Dict[str, List[Tuple[int, int]]] = field(default_factory=dict)
    # chrom -> merged exon intervals (used to distinguish intronic vs exonic).
    exon_intervals: Dict[str, List[Tuple[int, int]]] = field(default_factory=dict)
    # chrom -> list of (start, end, strand) exons, kept for strandedness
    # inference (needs per-exon strand, which merging discards).
    stranded_exons: Dict[str, List[Tuple[int, int, str]]] = field(default_factory=dict)
    # chrom -> list[Gene]
    genes: Dict[str, List[Gene]] = field(default_factory=dict)
    # chrom -> (sorted start array, sorted end array) for fast nearest lookup.
    _gene_starts: Dict[str, np.ndarray] = field(default_factory=dict)
    _gene_ends: Dict[str, np.ndarray] = field(default_factory=dict)

    def build_indexes(self) -> None:
        """Precompute sorted arrays for nearest-feature queries."""
        for chrom, genes in self.genes.items():
            genes.sort(key=lambda g: g.start)
            self._gene_starts[chrom] = np.array([g.start for g in genes], dtype=np.int64)
            self._gene_ends[chrom] = np.array([g.end for g in genes], dtype=np.int64)

    def mask_array(self, chrom: str, length: int) -> np.ndarray:
        """Return a boolean array of ``length`` marking annotated positions True."""
        mask = np.zeros(length, dtype=bool)
        for start, end in self.mask_intervals.get(chrom, []):
            s = max(0, start)
            e = min(length, end)
            if e > s:
                mask[s:e] = True
        return mask

    def exon_mask_array(self, chrom: str, length: int) -> np.ndarray:
        mask = np.zeros(length, dtype=bool)
        for start, end in self.exon_intervals.get(chrom, []):
            s = max(0, start)
            e = min(length, end)
            if e > s:
                mask[s:e] = True
        return mask

    def stranded_mask_arrays(self, chrom: str, length: int):
        """Return (plus_mask, minus_mask) boolean arrays for stranded masking.

        ``plus_mask`` marks positions annotated by a + (or strandless) feature,
        ``minus_mask`` positions annotated by a - (or strandless) feature.
        """
        def build(intervals):
            m = np.zeros(length, dtype=bool)
            for start, end in intervals:
                s = max(0, start)
                e = min(length, end)
                if e > s:
                    m[s:e] = True
            return m

        return (
            build(self.mask_intervals_plus.get(chrom, [])),
            build(self.mask_intervals_minus.get(chrom, [])),
        )

    def nearest_gene(
        self, chrom: str, start: int, end: int
    ) -> Tuple[Optional[Gene], int, List[Gene]]:
        """Return (nearest_gene, distance, overlapping_genes) for [start, end).

        ``distance`` is 0 when the region overlaps a gene span. Overlapping
        genes are returned separately so callers can decide sense/antisense.
        """
        genes = self.genes.get(chrom)
        if not genes:
            return None, -1, []

        starts = self._gene_starts[chrom]
        ends = self._gene_ends[chrom]

        # Overlap: gene.start < end AND gene.end > start.
        overlap_idx = np.flatnonzero((starts < end) & (ends > start))
        overlapping = [genes[i] for i in overlap_idx.tolist()]
        if overlapping:
            return overlapping[0], 0, overlapping

        # No overlap -> nearest by gap distance.
        # Distance to genes ending before the region start:
        left_gap = start - ends            # positive when gene is to the left
        right_gap = starts - end           # positive when gene is to the right
        dist = np.where(left_gap >= 0, left_gap,
                        np.where(right_gap >= 0, right_gap, 0))
        # For genes that overlap we already returned; remaining are all disjoint.
        nearest_i = int(np.argmin(dist))
        return genes[nearest_i], int(dist[nearest_i]), []


def _open(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def _parse_attributes(attr_field: str) -> Dict[str, str]:
    return {k: v for k, v in _ATTR_RE.findall(attr_field)}


def parse_gtf(path: str, annotation_mode: str) -> Annotation:
    """Parse a GTF into an :class:`Annotation`.

    ``annotation_mode`` selects what is masked:
      * exon        -> exon features
      * transcript  -> transcript spans (min..max of a transcript's features)
      * gene        -> gene spans
      * all         -> every feature line

    Genes and exons are always collected regardless of mode because they are
    needed for context labelling (intronic / antisense / intergenic).
    """
    ann = Annotation()

    # Raw interval accumulators before merging.
    mask_raw: Dict[str, List[Tuple[int, int]]] = {}          # positional (any strand)
    mask_raw_plus: Dict[str, List[Tuple[int, int]]] = {}     # + features only
    mask_raw_minus: Dict[str, List[Tuple[int, int]]] = {}    # - features only
    exon_raw: Dict[str, List[Tuple[int, int]]] = {}
    genes_seen: Dict[str, Gene] = {}
    # For transcript-span mode when explicit 'transcript' lines are absent.
    tx_span: Dict[str, List] = {}  # tx_id -> [chrom, start, end, strand]

    def add(dic, chrom, start, end):
        dic.setdefault(chrom, []).append((start, end))

    def add_mask(chrom, start, end, strand):
        """Record a masked interval both positionally and per-strand.

        Strandless features ('.') are masked on both strands (they block
        discovery on either strand)."""
        add(mask_raw, chrom, start, end)
        if strand == "+":
            add(mask_raw_plus, chrom, start, end)
        elif strand == "-":
            add(mask_raw_minus, chrom, start, end)
        else:
            add(mask_raw_plus, chrom, start, end)
            add(mask_raw_minus, chrom, start, end)

    with _open(path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start_s, end_s, score, strand, frame, attrs = parts[:9]
            try:
                start = int(start_s) - 1          # -> 0-based
                end = int(end_s)                  # half-open == 1-based inclusive end
            except ValueError:
                continue
            attr = _parse_attributes(attrs)

            if feature == "exon":
                add(exon_raw, chrom, start, end)
                ann.stranded_exons.setdefault(chrom, []).append((start, end, strand))

            # Collect gene records.
            if feature == "gene":
                gid = attr.get("gene_id", f"gene@{chrom}:{start}")
                genes_seen[gid] = Gene(
                    gene_id=gid,
                    gene_name=attr.get("gene_name", gid),
                    chrom=chrom,
                    start=start,
                    end=end,
                    strand=strand,
                )

            # Track transcript spans (from explicit lines or by aggregation).
            tid = attr.get("transcript_id")
            if tid:
                span = tx_span.get(tid)
                if span is None:
                    tx_span[tid] = [chrom, start, end, strand]
                else:
                    span[1] = min(span[1], start)
                    span[2] = max(span[2], end)

            # Build the mask according to mode.
            if annotation_mode == "all":
                add_mask(chrom, start, end, strand)
            elif annotation_mode == "exon" and feature == "exon":
                add_mask(chrom, start, end, strand)
            elif annotation_mode == "gene" and feature == "gene":
                add_mask(chrom, start, end, strand)
            elif annotation_mode == "transcript" and feature == "transcript":
                add_mask(chrom, start, end, strand)

    # If we asked for transcript/gene masking but the GTF lacked those explicit
    # feature lines, fall back to aggregated transcript spans.
    if annotation_mode == "transcript" and not any(mask_raw.values()):
        for tid, (chrom, start, end, strand) in tx_span.items():
            add_mask(chrom, start, end, strand)
    if annotation_mode == "gene" and not genes_seen:
        # No gene lines: aggregate transcript spans as a proxy for genes.
        for tid, (chrom, start, end, strand) in tx_span.items():
            add_mask(chrom, start, end, strand)

    # If no explicit gene lines exist, synthesise genes from transcript spans so
    # context labelling still works.
    if not genes_seen and tx_span:
        for tid, (chrom, start, end, strand) in tx_span.items():
            genes_seen[tid] = Gene(tid, tid, chrom, start, end, strand)

    # Merge raw intervals per chrom (sorted, gap 0 = touching allowed).
    for chrom, ivals in mask_raw.items():
        ann.mask_intervals[chrom] = merge_runs_with_gap(ivals, max_gap=0)
    for chrom, ivals in mask_raw_plus.items():
        ann.mask_intervals_plus[chrom] = merge_runs_with_gap(ivals, max_gap=0)
    for chrom, ivals in mask_raw_minus.items():
        ann.mask_intervals_minus[chrom] = merge_runs_with_gap(ivals, max_gap=0)
    for chrom, ivals in exon_raw.items():
        ann.exon_intervals[chrom] = merge_runs_with_gap(ivals, max_gap=0)

    for g in genes_seen.values():
        ann.genes.setdefault(g.chrom, []).append(g)

    ann.build_indexes()
    return ann


def read_chrom_sizes(fai_or_sizes: str) -> Dict[str, int]:
    """Read chromosome sizes from a ``.fai`` index or a two-column sizes file."""
    sizes: Dict[str, int] = {}
    with _open(fai_or_sizes) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                try:
                    sizes[parts[0]] = int(parts[1])
                except ValueError:
                    continue
    return sizes
