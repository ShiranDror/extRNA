"""External annotation overlay for consensus regions.

The cross-sample step produces a list of reproducible coordinates. Those
coordinates often fall on features that are already known — regulatory elements,
repeats, CAGE TSS peaks — and knowing that *upfront*, in the consensus table, is
far more useful than discovering it later in a browser.

This module loads one or more external annotation files (Ensembl GFF3 / GFF /
GTF) into a per-chromosome interval index and attaches overlap information to
each consensus region.

Design decisions
----------------
  * **Annotate, never classify.** Overlay data is attached AFTER the classifier
    and the consensus vote have run. It adds columns; it never changes a call.
    That keeps the classifier coverage-only and auditable, and it keeps the
    overlay independent of the calls — which is what makes a downstream
    "novel transcripts are enriched at enhancers" statement a real finding
    rather than a circular one.
  * **Strand-agnostic by default.** Ensembl regulatory features are mostly
    strandless ('.'), and where a strand exists (CTCF_binding_site) it denotes
    motif orientation, not transcription. Requiring a strand match would drop
    real hits. ``stranded=True`` is available for sources where strand does mean
    transcription (e.g. CAGE TSS peaks).
  * **Chromosome naming is normalised, and a total mismatch is an error.**
    Ensembl regulation uses '1'/'MT'; UCSC-style BAMs give 'chr1'/'chrM'. A
    silent zero-overlap result is the worst outcome here (it looks exactly like
    "nothing is regulatory"), so if no chromosome name can possibly match we
    raise instead of writing a table full of NA.

Pure standard library (bisect + gzip). No pysam, no numpy, no polars, so this
runs natively anywhere the merge step runs.
"""

from __future__ import annotations

import gzip
import os
import re
from bisect import bisect_left, bisect_right
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from urllib.parse import unquote

# GTF style:  gene_id "ENSG..."; transcript_id "ENST...";
_GTF_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')

# Attribute keys tried in order when naming a feature.
_ID_KEYS = ("ID", "Id", "id", "Name", "name", "gene_id", "transcript_id")
# Attribute keys carrying linked gene labels (Ensembl promoters carry these).
_GENE_KEYS = ("gene_name", "gene_id")

# Filenames -> short, readable column prefixes. Checked case-insensitively as
# substrings of the basename, first match wins.
_LABEL_KEYWORDS = (
    ("regulatory_features", "regulatory"),
    ("regulatory", "regulatory"),
    ("emar", "emar"),
    ("motif_features", "motif"),
    ("motif", "motif"),
    ("cage", "cage"),
    ("repeat", "repeat"),
    ("rmsk", "repeat"),
    ("enhancer", "enhancer"),
    ("promoter", "promoter"),
    ("lncrna", "lncrna"),
)

_SEP = ","          # inside a single column value (types, genes)
_ID_SEP = ";"       # between feature ids
_MAX_IDS = 10       # cap the id list so the column stays readable
_NA = "NA"


# --------------------------------------------------------------------------- #
# Chromosome naming
# --------------------------------------------------------------------------- #

def normalise_chrom(name: str) -> str:
    """Normalise a chromosome name for cross-source matching.

    'chr1' -> '1', 'chrM'/'M'/'chrMT' -> 'MT', 'chrX' -> 'X'. Anything we do not
    recognise is returned with only a 'chr' prefix stripped, so scaffold names
    still match when both sources use the same convention.
    """
    n = name.strip()
    if n[:3].lower() == "chr":
        n = n[3:]
    if n.upper() in ("M", "MT"):
        return "MT"
    return n


# --------------------------------------------------------------------------- #
# Attribute parsing (GFF3 and GTF)
# --------------------------------------------------------------------------- #

def parse_attributes(attr_field: str) -> Dict[str, str]:
    """Parse a GFF3 (``key=value;``) or GTF (``key "value";``) attribute column.

    The format is detected per line: GTF quoting is unambiguous, so if the GTF
    pattern matches anything we treat the line as GTF, otherwise as GFF3. GFF3
    values are percent-decoded per the spec.
    """
    gtf = _GTF_ATTR_RE.findall(attr_field)
    if gtf:
        return {k: v for k, v in gtf}

    out: Dict[str, str] = {}
    for chunk in attr_field.split(";"):
        chunk = chunk.strip()
        if not chunk or "=" not in chunk:
            continue
        key, _, value = chunk.partition("=")
        key = key.strip()
        if key:
            out[key] = unquote(value.strip())
    return out


def derive_label(path: str, used: Sequence[str] = ()) -> str:
    """Pick a short, readable column prefix for an annotation file.

    Tries a keyword match against the filename first (so the Ensembl
    ``...regulatory_features.v116.gff3.gz`` becomes ``regulatory``), then falls
    back to the longest dot-separated token of the basename. Uniquified against
    ``used``.
    """
    base = os.path.basename(path)
    lower = base.lower()

    label = None
    for needle, canonical in _LABEL_KEYWORDS:
        if needle in lower:
            label = canonical
            break

    if label is None:
        stem = base
        for ext in (".gz", ".bgz"):
            if stem.lower().endswith(ext):
                stem = stem[: -len(ext)]
        for ext in (".gff3", ".gff", ".gtf", ".bed", ".txt", ".tsv"):
            if stem.lower().endswith(ext):
                stem = stem[: -len(ext)]
        tokens = [t for t in re.split(r"[.]+", stem) if t]
        label = max(tokens, key=len) if tokens else "annotation"

    label = re.sub(r"\W+", "_", label).strip("_").lower() or "annotation"
    if label not in used:
        return label
    i = 2
    while f"{label}{i}" in used:
        i += 1
    return f"{label}{i}"


# --------------------------------------------------------------------------- #
# Features and index
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class OverlayFeature:
    """One external annotation feature, 0-based half-open."""

    start: int
    end: int
    ftype: str
    fid: str
    strand: str
    genes: Tuple[str, ...] = ()

    @property
    def label(self) -> str:
        return f"{self.ftype}:{self.fid}" if self.fid != _NA else self.ftype


@dataclass
class OverlayIndex:
    """Per-chromosome sorted interval index for one annotation source."""

    label: str
    path: str
    stranded: bool = False
    # normalised chrom -> features sorted by start
    features: Dict[str, List[OverlayFeature]] = field(default_factory=dict)
    # parallel sorted start lists + the longest feature per chrom, which bounds
    # how far back a query has to look.
    _starts: Dict[str, List[int]] = field(default_factory=dict)
    _max_len: Dict[str, int] = field(default_factory=dict)
    n_features: int = 0
    n_skipped: int = 0
    type_counts: Dict[str, int] = field(default_factory=dict)
    # Union bp occupied by each feature type across all chromosomes. Needed to
    # judge an overlap count against chance: with ~237k enhancers covering a
    # sizeable slice of the genome, a raw "N regions overlap an enhancer" is
    # uninterpretable without knowing the footprint being hit.
    covered_bp_by_type: Dict[str, int] = field(default_factory=dict)

    def finalise(self) -> None:
        for chrom, feats in self.features.items():
            feats.sort(key=lambda f: (f.start, f.end))
            self._starts[chrom] = [f.start for f in feats]
            self._max_len[chrom] = max((f.end - f.start) for f in feats) if feats else 0
        self._compute_footprint()

    def _compute_footprint(self) -> None:
        """Union bp per feature type (features of one type may overlap)."""
        totals: Dict[str, int] = {}
        for feats in self.features.values():
            by_type: Dict[str, List[Tuple[int, int]]] = {}
            for f in feats:
                by_type.setdefault(f.ftype, []).append((f.start, f.end))
            for ftype, spans in by_type.items():
                spans.sort()
                cur_s, cur_e = spans[0]
                acc = 0
                for s, e in spans[1:]:
                    if s <= cur_e:
                        cur_e = max(cur_e, e)
                    else:
                        acc += cur_e - cur_s
                        cur_s, cur_e = s, e
                acc += cur_e - cur_s
                totals[ftype] = totals.get(ftype, 0) + acc
        self.covered_bp_by_type = dict(sorted(totals.items()))

    @property
    def chroms(self):
        return set(self.features)

    def overlap(
        self, chrom: str, start: int, end: int, strand: str = "."
    ) -> List[OverlayFeature]:
        """Features overlapping the 0-based half-open interval [start, end)."""
        return self._scan(chrom, start, end, strand, require_overlap=True)

    def within(
        self, chrom: str, start: int, end: int, window: int, strand: str = "."
    ) -> List[OverlayFeature]:
        """Features overlapping [start, end) expanded by ``window`` on each side."""
        lo = max(0, start - window)
        return self._scan(chrom, lo, end + window, strand, require_overlap=True)

    def _scan(
        self, chrom: str, start: int, end: int, strand: str, require_overlap: bool
    ) -> List[OverlayFeature]:
        key = normalise_chrom(chrom)
        feats = self.features.get(key)
        if not feats or end <= start:
            return []

        starts = self._starts[key]
        max_len = self._max_len[key]
        # A feature starting before start-max_len cannot reach start.
        lo = bisect_left(starts, start - max_len)
        # A feature starting at or after end cannot overlap.
        hi = bisect_right(starts, end - 1)

        hits = []
        for i in range(lo, hi):
            f = feats[i]
            if f.end <= start:
                continue
            if self.stranded and strand in ("+", "-") and f.strand in ("+", "-"):
                if f.strand != strand:
                    continue
            hits.append(f)
        return hits


def parse_overlay(
    path: str,
    label: Optional[str] = None,
    feature_types: Optional[Iterable[str]] = None,
    stranded: bool = False,
) -> OverlayIndex:
    """Load a GFF3/GFF/GTF annotation file into an :class:`OverlayIndex`.

    ``feature_types`` optionally restricts which column-3 values are kept
    (case-insensitive), e.g. ``["enhancer", "promoter"]``.
    Input coordinates are 1-based inclusive and converted to 0-based half-open.
    """
    wanted = {t.lower() for t in feature_types} if feature_types else None
    idx = OverlayIndex(
        label=label or derive_label(path), path=path, stranded=stranded
    )

    opener = gzip.open if path.endswith((".gz", ".bgz")) else open
    with opener(path, "rt") as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            parts = line.rstrip("\n").rstrip("\r").split("\t")
            if len(parts) < 9:
                continue
            chrom, _source, ftype, start_s, end_s, _score, strand = parts[:7]
            attrs_field = parts[8]

            if wanted is not None and ftype.lower() not in wanted:
                continue
            try:
                start = int(start_s) - 1     # -> 0-based
                end = int(end_s)             # half-open == 1-based inclusive end
            except ValueError:
                idx.n_skipped += 1
                continue
            if end <= start:
                idx.n_skipped += 1
                continue

            attrs = parse_attributes(attrs_field)
            fid = _NA
            for k in _ID_KEYS:
                if attrs.get(k):
                    fid = attrs[k]
                    break

            genes: Tuple[str, ...] = ()
            for k in _GENE_KEYS:
                raw = attrs.get(k)
                if raw:
                    # Ensembl promoters can list several genes, comma-separated.
                    genes = tuple(g for g in (v.strip() for v in raw.split(",")) if g)
                    break

            key = normalise_chrom(chrom)
            idx.features.setdefault(key, []).append(
                OverlayFeature(start, end, ftype, fid, strand, genes)
            )
            idx.n_features += 1
            idx.type_counts[ftype] = idx.type_counts.get(ftype, 0) + 1

    idx.finalise()
    return idx


# --------------------------------------------------------------------------- #
# Annotation of regions
# --------------------------------------------------------------------------- #

def sanitise_name_token(token: str) -> str:
    """Make an arbitrary feature type safe to embed in a transcript identifier.

    Feature types come from column 3 of a user-supplied file, so they cannot be
    assumed to be identifier-safe. '-' is the separator used to attach the token
    to the transcript name, so it is collapsed too.
    """
    cleaned = re.sub(r"[^A-Za-z0-9_.]+", "_", token.strip()).strip("_.")
    return cleaned or "annotated"


def annotation_name_token(
    overlay_types: Dict[str, Dict[str, int]], labels: Sequence[str]
) -> Optional[str]:
    """The feature type to attach to a consensus transcript name, or None.

    Deliberately **vocabulary-agnostic**: there is no built-in list of known
    feature types and no priority table, so adding an annotation source never
    requires a code change here.

    "First" means the first type of the first ``--annotate`` source that overlaps
    this region, using the same ordering as the ``<label>_types`` column — so the
    suffix can always be read straight off that column. Source precedence is the
    order the user passed ``--annotate`` in, which keeps the choice in their hands
    rather than hardcoded.
    """
    for label in labels:
        counts = overlay_types.get(label) or {}
        if counts:
            return sanitise_name_token(sorted(counts)[0])
    return None


def overlay_columns(label: str) -> List[str]:
    """Column names contributed by one annotation source, in output order."""
    return [
        f"{label}_n",
        f"{label}_types",
        f"{label}_ids",
        f"{label}_overlap_bp",
        f"{label}_overlap_frac",
        f"{label}_genes",
        f"{label}_nearest",
        f"{label}_nearest_distance",
    ]


def _merged_overlap_bp(
    feats: Sequence[OverlayFeature], start: int, end: int
) -> int:
    """Bases of [start, end) covered by at least one feature (union, not sum)."""
    clipped = sorted(
        (max(f.start, start), min(f.end, end))
        for f in feats
        if min(f.end, end) > max(f.start, start)
    )
    total = 0
    cur_s = cur_e = None
    for s, e in clipped:
        if cur_e is None:
            cur_s, cur_e = s, e
        elif s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            total += cur_e - cur_s
            cur_s, cur_e = s, e
    if cur_e is not None:
        total += cur_e - cur_s
    return total


def _distance(f: OverlayFeature, start: int, end: int) -> int:
    """Gap in bp between a feature and [start, end); 0 when they overlap."""
    if f.end > start and f.start < end:
        return 0
    return f.start - end if f.start >= end else start - f.end


@dataclass
class RegionAnnotation:
    """One source's overlay result for one region.

    ``columns`` goes straight into the output table; ``type_counts`` is kept
    separately so the summary can cross-tabulate feature type against consensus
    class without re-parsing the formatted column.
    """

    columns: Dict[str, object] = field(default_factory=dict)
    type_counts: Dict[str, int] = field(default_factory=dict)


def annotate_region(
    idx: OverlayIndex,
    chrom: str,
    start: int,
    end: int,
    strand: str,
    nearest_window: int,
) -> RegionAnnotation:
    """Compute one source's overlay fields for a single 0-based region."""
    label = idx.label
    hits = idx.overlap(chrom, start, end, strand)
    length = max(1, end - start)

    if hits:
        counts: Dict[str, int] = {}
        for f in hits:
            counts[f.ftype] = counts.get(f.ftype, 0) + 1
        types = _SEP.join(f"{t}:{counts[t]}" for t in sorted(counts))

        ids = [f.fid for f in hits if f.fid != _NA]
        ids_str = _ID_SEP.join(ids[:_MAX_IDS]) if ids else _NA
        if len(ids) > _MAX_IDS:
            ids_str += f"{_ID_SEP}+{len(ids) - _MAX_IDS}_more"

        genes = sorted({g for f in hits for g in f.genes})
        overlap_bp = _merged_overlap_bp(hits, start, end)
        return RegionAnnotation(
            columns={
                f"{label}_n": len(hits),
                f"{label}_types": types,
                f"{label}_ids": ids_str,
                f"{label}_overlap_bp": overlap_bp,
                f"{label}_overlap_frac": round(overlap_bp / length, 4),
                f"{label}_genes": _SEP.join(genes) if genes else _NA,
                # Overlapping, so the nearest feature is the region itself.
                f"{label}_nearest": _NA,
                f"{label}_nearest_distance": 0,
            },
            type_counts=counts,
        )

    # No overlap: report the nearest feature inside the window, if any.
    nearest_label, nearest_dist = _NA, -1
    if nearest_window > 0:
        near = idx.within(chrom, start, end, nearest_window, strand)
        best = None
        best_d = None
        for f in near:
            d = _distance(f, start, end)
            if best_d is None or d < best_d:
                best, best_d = f, d
        if best is not None:
            nearest_label, nearest_dist = best.label, best_d

    return RegionAnnotation(
        columns={
            f"{label}_n": 0,
            f"{label}_types": _NA,
            f"{label}_ids": _NA,
            f"{label}_overlap_bp": 0,
            f"{label}_overlap_frac": 0.0,
            f"{label}_genes": _NA,
            f"{label}_nearest": nearest_label,
            f"{label}_nearest_distance": nearest_dist,
        },
        type_counts={},
    )


def check_chrom_compatibility(
    idx: OverlayIndex, region_chroms: Iterable[str]
) -> None:
    """Raise if no chromosome name can possibly match between the two sources.

    A naming mismatch (Ensembl '1' vs UCSC 'chr1') otherwise produces a table of
    NA that reads exactly like a genuine "nothing overlaps" result.
    """
    region_norm = {normalise_chrom(c) for c in region_chroms}
    if not region_norm or not idx.chroms:
        return
    if region_norm & idx.chroms:
        return
    raise ValueError(
        f"annotation source {idx.label!r} ({os.path.basename(idx.path)}) shares no "
        f"chromosome names with the candidate regions, so nothing could ever "
        f"overlap. Candidate chromosomes look like "
        f"{sorted(region_norm)[:5]}; annotation chromosomes look like "
        f"{sorted(idx.chroms)[:5]}. Check that both use the same assembly and "
        f"naming convention (Ensembl '1' vs UCSC 'chr1')."
    )


def load_sources(
    specs: Sequence[str],
    labels: Optional[Sequence[str]] = None,
    feature_types: Optional[Iterable[str]] = None,
    stranded: bool = False,
) -> List[OverlayIndex]:
    """Load annotation sources from ``--annotate`` specs.

    Each spec is either ``path`` or ``label=path``. Explicit ``labels`` (from
    ``--annotate-labels``) take precedence over both.
    """
    if labels is not None and len(labels) != len(specs):
        raise ValueError("--annotate-labels count must match --annotate count")

    parsed: List[Tuple[Optional[str], str]] = []
    for spec in specs:
        # Only treat '=' as a label separator when the left side looks like a
        # label, so Windows paths and '=' inside filenames are safe.
        if "=" in spec:
            head, _, tail = spec.partition("=")
            if head and tail and re.fullmatch(r"\w+", head):
                parsed.append((head, tail))
                continue
        parsed.append((None, spec))

    indexes: List[OverlayIndex] = []
    used: List[str] = []
    for i, (spec_label, path) in enumerate(parsed):
        if not os.path.exists(path):
            raise FileNotFoundError(f"--annotate file not found: {path!r}")
        label = (labels[i] if labels is not None else spec_label) or derive_label(
            path, used
        )
        label = re.sub(r"\W+", "_", label).strip("_").lower()
        if label in used:
            raise ValueError(f"duplicate annotation label: {label!r}")
        used.append(label)
        indexes.append(
            parse_overlay(
                path, label=label, feature_types=feature_types, stranded=stranded
            )
        )
    return indexes
