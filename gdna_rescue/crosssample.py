"""Cross-sample consensus of candidate regions.

Given the per-sample ``*.candidate_regions.tsv`` files from N samples, cluster
overlapping candidates across samples and keep the loci reproduced in at least
``--min-samples`` samples. Genuine novel transcripts should recur across
biological replicates; sample-specific candidates are more likely to be noise or
alignment artifacts.

Design decisions (see README):
  * Candidates never share exact coordinates across samples, so loci are matched
    by RECIPROCAL OVERLAP (each must cover >= --min-reciprocal-overlap of the
    other). The default (0.85) is deliberately high: two intervals that are
    "the same transcript" should be nearly co-extensive, not merely touching.
  * Matching is STRAND-AWARE by default (this is RNA): a + and a - locus at the
    same position are not the same feature.
  * Every candidate contributes to a cluster's sample count, but the class is
    decided by a majority VOTE across samples, with ties broken toward the more
    conservative (reject) call.
  * Only reproducible NOVEL transcripts are written to the consensus GTF.
    Recurrent gDNA / multimapper loci get their own consensus classes and are
    reported in the table but NOT added to the GTF: recurrent contamination is
    worth investigating (and confirming biochemically) but must not be added to
    an annotation. Losing a true annotation is preferable to adding a bad one.

Pure polars + Python (no pysam), so this runs natively anywhere.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import polars as pl

from .classify import (
    LIKELY_GDNA,
    LIKELY_MULTIMAPPER,
    LIKELY_NOVEL,
    POSSIBLE_BIDIRECTIONAL,
)

# Consensus-level classes.
CONS_NOVEL = "reproducible_novel"
CONS_BIDIR = "reproducible_bidirectional"
CONS_GDNA = "recurrent_gDNA"
CONS_MULTI = "recurrent_multimapper_artifact"

# Most-conservative-first ordering, used to break class-vote ties.
_CONSERVATIVE_ORDER = [
    LIKELY_GDNA,
    LIKELY_MULTIMAPPER,
    POSSIBLE_BIDIRECTIONAL,
    LIKELY_NOVEL,
]

_PER_SAMPLE_TO_CONSENSUS = {
    LIKELY_NOVEL: CONS_NOVEL,
    POSSIBLE_BIDIRECTIONAL: CONS_BIDIR,
    LIKELY_GDNA: CONS_GDNA,
    LIKELY_MULTIMAPPER: CONS_MULTI,
}

GTF_SOURCE = "gdna_rescue_consensus"

# Columns we rely on from the per-sample TSV.
_REQUIRED_COLS = ("chrom", "start", "end", "class")


@dataclass
class ConsensusConfig:
    tsvs: List[str] = field(default_factory=list)
    sample_names: Optional[List[str]] = None
    out_prefix: str = "cohort"
    min_samples: int = 2
    min_reciprocal_overlap: float = 0.85
    strand_aware: bool = True
    include_bidirectional: bool = False  # add reproducible bidirectional to GTF
    verbose: bool = False

    def validate(self) -> None:
        if not self.tsvs:
            raise ValueError("at least one --tsv is required")
        for p in self.tsvs:
            if not os.path.exists(p):
                raise FileNotFoundError(f"candidate TSV not found: {p!r}")
        if self.sample_names and len(self.sample_names) != len(self.tsvs):
            raise ValueError("--sample-names count must match --tsv count")
        if self.min_samples < 1:
            raise ValueError("--min-samples must be >= 1")
        if not (0.0 < self.min_reciprocal_overlap <= 1.0):
            raise ValueError("--min-reciprocal-overlap must be in (0, 1]")


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #

def _default_sample_name(path: str) -> str:
    base = os.path.basename(path)
    for suffix in (".candidate_regions.tsv", ".tsv"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def load_candidates(cfg: ConsensusConfig) -> pl.DataFrame:
    """Load and vertically concatenate per-sample candidate TSVs (tagged)."""
    names = cfg.sample_names or [_default_sample_name(p) for p in cfg.tsvs]
    if len(set(names)) != len(names):
        raise ValueError(f"sample names are not unique: {names}")

    frames = []
    for path, name in zip(cfg.tsvs, names):
        df = pl.read_csv(path, separator="\t", infer_schema_length=10000)
        missing = [c for c in _REQUIRED_COLS if c not in df.columns]
        if missing:
            raise ValueError(f"{path!r} is missing required columns: {missing}")
        df = df.with_columns(pl.lit(name).alias("sample"))
        frames.append(df)
    return pl.concat(frames, how="diagonal_relaxed")


# --------------------------------------------------------------------------- #
# Clustering by reciprocal overlap (strand-aware)
# --------------------------------------------------------------------------- #

def _reciprocal_overlap(s1: int, e1: int, s2: int, e2: int) -> float:
    """Return min(overlap/len1, overlap/len2) for two half-open intervals."""
    inter = min(e1, e2) - max(s1, s2)
    if inter <= 0:
        return 0.0
    len1 = e1 - s1
    len2 = e2 - s2
    if len1 <= 0 or len2 <= 0:
        return 0.0
    return min(inter / len1, inter / len2)


class _UnionFind:
    def __init__(self, n: int):
        self.parent = list(range(n))

    def find(self, x: int) -> int:
        root = x
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[x] != root:  # path compression
            self.parent[x], x = root, self.parent[x]
        return root

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def cluster_rows(rows: List[dict], cfg: ConsensusConfig) -> List[List[int]]:
    """Cluster candidate rows by reciprocal overlap (single linkage).

    ``rows`` must carry integer 0-based half-open ``_s`` / ``_e`` and a group key
    ``_group`` (chrom or chrom+strand). Returns a list of clusters (each a list
    of row indices).
    """
    uf = _UnionFind(len(rows))

    # Group indices, then sweep within each group so we only compare intervals
    # that can overlap.
    groups: Dict[Tuple, List[int]] = {}
    for i, r in enumerate(rows):
        groups.setdefault(r["_group"], []).append(i)

    for idxs in groups.values():
        idxs.sort(key=lambda i: rows[i]["_s"])
        active: List[int] = []
        for i in idxs:
            si, ei = rows[i]["_s"], rows[i]["_e"]
            active = [j for j in active if rows[j]["_e"] > si]  # prune finished
            for j in active:
                if _reciprocal_overlap(
                    si, ei, rows[j]["_s"], rows[j]["_e"]
                ) >= cfg.min_reciprocal_overlap:
                    uf.union(i, j)
            active.append(i)

    clusters: Dict[int, List[int]] = {}
    for i in range(len(rows)):
        clusters.setdefault(uf.find(i), []).append(i)
    return list(clusters.values())


# --------------------------------------------------------------------------- #
# Consensus building
# --------------------------------------------------------------------------- #

def _vote_class(member_classes: List[str]) -> str:
    """Majority vote; ties broken toward the most conservative (reject) class."""
    counts: Dict[str, int] = {}
    for c in member_classes:
        counts[c] = counts.get(c, 0) + 1
    best = max(counts.values())
    tied = [c for c, n in counts.items() if n == best]
    if len(tied) == 1:
        return tied[0]
    for c in _CONSERVATIVE_ORDER:  # earliest = most conservative
        if c in tied:
            return c
    return sorted(tied)[0]


@dataclass
class ConsensusRegion:
    consensus_id: str
    chrom: str
    start: int              # 0-based half-open
    end: int
    strand: str
    n_samples: int
    n_members: int
    samples: List[str]
    member_region_ids: List[str]
    member_classes: List[str]
    majority_class: str
    consensus_class: str
    class_agreement: float
    mean_unique_fraction: Optional[float]
    mean_dual_strand_fraction: Optional[float]
    mean_profile_correlation: Optional[float]
    mean_avg_depth: Optional[float]
    passes_min_samples: bool
    in_consensus_gtf: bool
    consensus_transcript_name: Optional[str] = None


def _mean_opt(rows: List[dict], key: str) -> Optional[float]:
    vals = [r[key] for r in rows if r.get(key) is not None]
    vals = [v for v in vals if isinstance(v, (int, float))]
    return round(sum(vals) / len(vals), 4) if vals else None


def build_consensus(df: pl.DataFrame, cfg: ConsensusConfig) -> List[ConsensusRegion]:
    """Cluster candidates and derive consensus regions."""
    has_strand = "dominant_strand" in df.columns
    rows = df.to_dicts()

    # Attach clustering helpers (TSV coords are 1-based inclusive -> half-open).
    for r in rows:
        strand = r.get("dominant_strand", ".") if has_strand else "."
        r["_s"] = int(r["start"]) - 1
        r["_e"] = int(r["end"])
        r["_group"] = (r["chrom"], strand) if cfg.strand_aware else (r["chrom"],)

    clusters = cluster_rows(rows, cfg)

    consensus: List[ConsensusRegion] = []
    for members in clusters:
        mrows = [rows[i] for i in members]
        samples = sorted({r["sample"] for r in mrows})
        n_samples = len(samples)
        classes = [str(r["class"]) for r in mrows]
        majority = _vote_class(classes)
        cons_class = _PER_SAMPLE_TO_CONSENSUS.get(majority, majority)
        agreement = round(classes.count(majority) / len(classes), 4)

        start = min(r["_s"] for r in mrows)
        end = max(r["_e"] for r in mrows)

        # Consensus strand: fixed by the group when strand-aware, else majority.
        if cfg.strand_aware:
            strand = mrows[0].get("dominant_strand", ".")
        else:
            strand = _vote_class([r.get("dominant_strand", ".") for r in mrows])

        passes = n_samples >= cfg.min_samples
        in_gtf = passes and (
            cons_class == CONS_NOVEL
            or (cons_class == CONS_BIDIR and cfg.include_bidirectional)
        )

        consensus.append(
            ConsensusRegion(
                consensus_id=f"{mrows[0]['chrom']}:{start}-{end}:{strand}",
                chrom=mrows[0]["chrom"],
                start=start,
                end=end,
                strand=strand if strand in ("+", "-") else ".",
                n_samples=n_samples,
                n_members=len(mrows),
                samples=samples,
                member_region_ids=[str(r.get("region_id", "NA")) for r in mrows],
                member_classes=classes,
                majority_class=majority,
                consensus_class=cons_class,
                class_agreement=agreement,
                mean_unique_fraction=_mean_opt(mrows, "unique_fraction"),
                mean_dual_strand_fraction=_mean_opt(mrows, "dual_strand_fraction"),
                mean_profile_correlation=_mean_opt(mrows, "profile_correlation"),
                mean_avg_depth=_mean_opt(mrows, "avg_depth"),
                passes_min_samples=passes,
                in_consensus_gtf=in_gtf,
            )
        )

    # Deterministic order: genomic position.
    consensus.sort(key=lambda c: (c.chrom, c.start, c.end, c.strand))
    counter = 0
    for c in consensus:
        if c.in_consensus_gtf:
            counter += 1
            c.consensus_transcript_name = f"consensus_transcript_{counter}"
    return consensus


# --------------------------------------------------------------------------- #
# Writers
# --------------------------------------------------------------------------- #

def consensus_to_dataframe(regions: List[ConsensusRegion]) -> pl.DataFrame:
    rows = []
    for c in regions:
        rows.append(
            {
                "consensus_id": c.consensus_id,
                "chrom": c.chrom,
                "start": c.start + 1,   # back to 1-based inclusive for the table
                "end": c.end,
                "length": c.end - c.start,
                "strand": c.strand,
                "consensus_class": c.consensus_class,
                "majority_class": c.majority_class,
                "class_agreement": c.class_agreement,
                "n_samples": c.n_samples,
                "n_members": c.n_members,
                "samples": ",".join(c.samples),
                "member_classes": ",".join(c.member_classes),
                "mean_unique_fraction": c.mean_unique_fraction,
                "mean_dual_strand_fraction": c.mean_dual_strand_fraction,
                "mean_profile_correlation": c.mean_profile_correlation,
                "mean_avg_depth": c.mean_avg_depth,
                "passes_min_samples": "yes" if c.passes_min_samples else "no",
                "in_consensus_gtf": "yes" if c.in_consensus_gtf else "no",
                "consensus_transcript_name": c.consensus_transcript_name or "NA",
                "member_region_ids": ";".join(c.member_region_ids),
            }
        )
    schema_hint = rows if rows else [{
        "consensus_id": None, "chrom": None, "start": None, "end": None,
        "length": None, "strand": None, "consensus_class": None,
        "majority_class": None, "class_agreement": None, "n_samples": None,
        "n_members": None, "samples": None, "member_classes": None,
        "mean_unique_fraction": None, "mean_dual_strand_fraction": None,
        "mean_profile_correlation": None, "mean_avg_depth": None,
        "passes_min_samples": None, "in_consensus_gtf": None,
        "consensus_transcript_name": None, "member_region_ids": None,
    }]
    df = pl.DataFrame(schema_hint)
    if not rows:
        df = df.clear()
    return df


def _gtf_attr(d: Dict[str, str]) -> str:
    return " ".join(f'{k} "{v}";' for k, v in d.items())


def write_consensus_gtf(regions: List[ConsensusRegion], path: str) -> int:
    kept = [c for c in regions if c.in_consensus_gtf and c.consensus_transcript_name]
    with open(path, "w") as fh:
        fh.write("##description: reproducible novel transcripts (cross-sample consensus)\n")
        for c in kept:
            name = c.consensus_transcript_name
            start1 = c.start + 1
            common = {
                "gene_id": f"{name}_gene",
                "transcript_id": name,
                "gene_name": name,
                "source": GTF_SOURCE,
                "consensus_class": c.consensus_class,
                "n_samples": str(c.n_samples),
                "samples": ",".join(c.samples),
                "member_region_ids": ";".join(c.member_region_ids),
            }
            tx = _gtf_attr(common)
            exon = _gtf_attr({**common, "exon_number": "1"})
            fh.write(
                f"{c.chrom}\t{GTF_SOURCE}\ttranscript\t{start1}\t{c.end}\t.\t"
                f"{c.strand}\t.\t{tx}\n"
            )
            fh.write(
                f"{c.chrom}\t{GTF_SOURCE}\texon\t{start1}\t{c.end}\t.\t"
                f"{c.strand}\t.\t{exon}\n"
            )
    return len(kept)


def write_consensus_summary(
    cfg: ConsensusConfig, regions: List[ConsensusRegion], path: str
) -> dict:
    def count(cls):
        return sum(1 for c in regions if c.consensus_class == cls and c.passes_min_samples)

    passing = [c for c in regions if c.passes_min_samples]
    summary = {
        "parameters": {
            "tsvs": cfg.tsvs,
            "sample_names": cfg.sample_names,
            "n_samples_input": len(cfg.tsvs),
            "min_samples": cfg.min_samples,
            "min_reciprocal_overlap": cfg.min_reciprocal_overlap,
            "strand_aware": cfg.strand_aware,
            "include_bidirectional": cfg.include_bidirectional,
        },
        "n_clusters_total": len(regions),
        "n_clusters_passing_min_samples": len(passing),
        "n_reproducible_novel": count(CONS_NOVEL),
        "n_reproducible_bidirectional": count(CONS_BIDIR),
        "n_recurrent_gDNA": count(CONS_GDNA),
        "n_recurrent_multimapper_artifact": count(CONS_MULTI),
        "n_written_to_consensus_gtf": sum(1 for c in regions if c.in_consensus_gtf),
    }
    with open(path, "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    return summary


def run(cfg: ConsensusConfig) -> dict:
    """Run the cross-sample consensus and write outputs. Returns the summary."""
    from .utils import get_logger

    logger = get_logger(cfg.verbose)
    cfg.validate()

    df = load_candidates(cfg)
    logger.info("Loaded %d candidate rows from %d samples.", df.height, len(cfg.tsvs))

    regions = build_consensus(df, cfg)
    passing = [c for c in regions if c.passes_min_samples]

    table = f"{cfg.out_prefix}.consensus_regions.tsv"
    gtf = f"{cfg.out_prefix}.consensus_transcripts.gtf"
    summ = f"{cfg.out_prefix}.consensus_summary.json"

    # Write only the passing clusters to the main table (that IS the filter),
    # keeping every consensus class so recurrent gDNA/artifacts are visible.
    consensus_to_dataframe(passing).write_csv(table, separator="\t")
    n_gtf = write_consensus_gtf(regions, gtf)
    summary = write_consensus_summary(cfg, regions, summ)

    logger.info(
        "Consensus: %d clusters, %d reproducible (>=%d samples) | "
        "%d novel | %d bidirectional | %d recurrent_gDNA | %d multimapper | "
        "%d written to GTF.",
        summary["n_clusters_total"],
        summary["n_clusters_passing_min_samples"],
        cfg.min_samples,
        summary["n_reproducible_novel"],
        summary["n_reproducible_bidirectional"],
        summary["n_recurrent_gDNA"],
        summary["n_recurrent_multimapper_artifact"],
        n_gtf,
    )
    logger.info("Consensus outputs written with prefix: %s", cfg.out_prefix)
    return summary
