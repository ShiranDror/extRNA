"""Output writers: candidate TSV, unknown-transcript GTF, merged GTF, summary JSON,
optional BED and bedGraph.

Coordinate conventions on output:
  * TSV  : 1-based inclusive start/end (matches GTF / genome browsers).
  * GTF  : 1-based inclusive (spec).
  * BED / bedGraph : 0-based half-open (spec).
"""

from __future__ import annotations

import gzip
import json
import os
from typing import Dict, List

import pandas as pd

from .classify import (
    LIKELY_GDNA,
    LIKELY_NOVEL,
    LIKELY_MULTIMAPPER,
    POSSIBLE_BIDIRECTIONAL,
)
from .config import Config
from .discovery import Candidate

GTF_SOURCE = "gdna_rescue"


def candidates_to_dataframe(candidates: List[Candidate]) -> pd.DataFrame:
    """Flatten candidates + metrics into the QC table."""
    rows = []
    for c in candidates:
        m = c.metrics
        rows.append(
            {
                "region_id": c.region_id,
                "chrom": c.chrom,
                "start": c.start + 1,          # 1-based inclusive
                "end": c.end,
                "length": c.length,
                "class": c.label,
                "total_coverage": round(m.total_coverage, 3),
                "covered_bases": m.covered_bases,
                "covered_fraction": round(m.covered_fraction, 4),
                "avg_depth": round(m.avg_depth, 3),
                "max_depth": round(m.max_depth, 3),
                "unique_total": round(m.unique_total, 1),
                "multi_total": round(m.multi_total, 1),
                "unique_fraction": round(m.unique_fraction, 4),
                "plus_covered_len": m.plus_covered_len,
                "minus_covered_len": m.minus_covered_len,
                "plus_mean_depth": round(m.plus_mean_depth, 3),
                "minus_mean_depth": round(m.minus_mean_depth, 3),
                "pm_depth_ratio": round(m.pm_depth_ratio, 4),
                "covered_len_ratio": round(m.covered_len_ratio, 4),
                "strand_length_ratio_diff": round(m.strand_length_ratio_diff, 4),
                "dominant_strand": m.dominant_strand,
                "dominant_strand_fraction": round(m.dominant_strand_fraction, 4),
                "dual_strand_fraction": round(m.dual_strand_fraction, 4),
                "strand_overlap_jaccard": round(m.strand_overlap_jaccard, 4),
                "profile_correlation": round(m.profile_correlation, 4),
                "plus_cv": round(m.plus_cv, 4),
                "minus_cv": round(m.minus_cv, 4),
                "strand_balance": round(m.strand_entropy, 4),
                "nearest_feature_id": c.nearest_feature_id or "NA",
                "nearest_feature_distance": c.nearest_feature_distance,
                "context_label": c.context_label,
                "kept_as_unknown_transcript": "yes" if c.kept else "no",
                "unknown_transcript_name": c.unknown_transcript_name or "NA",
                "reason_for_classification": c.reason,
            }
        )
    return pd.DataFrame(rows)


def write_tsv(candidates: List[Candidate], path: str) -> None:
    df = candidates_to_dataframe(candidates)
    df.to_csv(path, sep="\t", index=False)


def _gtf_attr(d: Dict[str, str]) -> str:
    return " ".join(f'{k} "{v}";' for k, v in d.items())


def _gtf_lines_for_candidate(c: Candidate) -> List[str]:
    """Emit transcript + single-exon lines (1-based inclusive) for a kept region."""
    name = c.unknown_transcript_name
    strand = c.metrics.dominant_strand
    if strand not in ("+", "-"):
        strand = "."
    start1 = c.start + 1
    end1 = c.end

    common = {
        "gene_id": f"{name}_gene",
        "transcript_id": name,
        "gene_name": name,
        "source": GTF_SOURCE,
        "classification": c.label,
        "context": c.context_label,
        "original_region_id": c.region_id,
    }
    tx_attr = _gtf_attr(common)
    exon_attr = _gtf_attr({**common, "exon_number": "1"})

    tx = f"{c.chrom}\t{GTF_SOURCE}\ttranscript\t{start1}\t{end1}\t.\t{strand}\t.\t{tx_attr}"
    exon = f"{c.chrom}\t{GTF_SOURCE}\texon\t{start1}\t{end1}\t.\t{strand}\t.\t{exon_attr}"
    return [tx, exon]


def write_unknown_gtf(candidates: List[Candidate], path: str) -> int:
    """Write only the rescued unknown transcripts. Returns count written."""
    kept = [c for c in candidates if c.kept and c.unknown_transcript_name]
    with open(path, "w") as fh:
        fh.write("##description: novel unannotated transcripts rescued by gdna_rescue\n")
        for c in kept:
            for line in _gtf_lines_for_candidate(c):
                fh.write(line + "\n")
    return len(kept)


def _open_maybe_gz(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def write_merged_gtf(
    original_gtf: str, candidates: List[Candidate], path: str
) -> None:
    """Write original annotation verbatim followed by rescued transcripts."""
    kept = [c for c in candidates if c.kept and c.unknown_transcript_name]
    with open(path, "w") as out:
        with _open_maybe_gz(original_gtf) as src:
            for line in src:
                out.write(line if line.endswith("\n") else line + "\n")
        out.write("##gdna_rescue: appended rescued unknown transcripts below\n")
        for c in kept:
            for line in _gtf_lines_for_candidate(c):
                out.write(line + "\n")


def write_bed(candidates: List[Candidate], path: str) -> None:
    """BED6 of all candidate regions (0-based half-open), name=class."""
    with open(path, "w") as fh:
        for c in candidates:
            strand = c.metrics.dominant_strand
            if strand not in ("+", "-"):
                strand = "."
            score = min(1000, int(round(c.metrics.avg_depth)))
            fh.write(
                f"{c.chrom}\t{c.start}\t{c.end}\t{c.label}\t{score}\t{strand}\n"
            )


def write_bedgraph(candidates: List[Candidate], prefix: str) -> None:
    """Write per-strand mean-depth bedGraph over candidate regions.

    Region-level (not per-base) summary — a lightweight, browser-loadable view.
    """
    for strand_name, attr in (("plus", "plus_mean_depth"), ("minus", "minus_mean_depth")):
        with open(f"{prefix}.candidates.{strand_name}.bedgraph", "w") as fh:
            fh.write(f'track type=bedGraph name="candidates_{strand_name}"\n')
            for c in candidates:
                val = round(getattr(c.metrics, attr), 3)
                fh.write(f"{c.chrom}\t{c.start}\t{c.end}\t{val}\n")


def compute_gdna_qc(
    candidates: List[Candidate],
    genome_unique_coverage: int,
    genome_multi_coverage: int,
) -> Dict:
    """Compute gDNA contamination QC metrics.

    ``genome_*_coverage`` are genome-wide summed per-base depths (aligned bases)
    from unique / multimapped reads. The headline metric,
    ``pct_gDNA_of_mapped_coverage``, is the share of ALL uniquely-mapped signal
    that falls in gDNA-flagged unannotated regions — a library-level
    contamination proxy. We also report the share relative to candidate regions
    only, since gDNA is only tested among unannotated candidates.
    """
    def pct(a, b):
        return round(100.0 * a / b, 4) if b > 0 else 0.0

    cand_cov = sum(c.metrics.total_coverage for c in candidates)
    cand_bases = sum(c.length for c in candidates)
    gdna = [c for c in candidates if c.label == LIKELY_GDNA]
    gdna_cov = sum(c.metrics.total_coverage for c in gdna)
    gdna_bases = sum(c.length for c in gdna)

    return {
        "genome_mapped_unique_coverage": genome_unique_coverage,
        "genome_mapped_multi_coverage": genome_multi_coverage,
        "candidate_unique_coverage": cand_cov,
        "candidate_bases": cand_bases,
        "gDNA_unique_coverage": gdna_cov,
        "gDNA_bases": gdna_bases,
        "n_gDNA_regions": len(gdna),
        "pct_gDNA_of_mapped_coverage": pct(gdna_cov, genome_unique_coverage),
        "pct_gDNA_of_candidate_coverage": pct(gdna_cov, cand_cov),
        "pct_gDNA_of_candidate_bases": pct(gdna_bases, cand_bases),
    }


def write_multiqc_tsv(
    cfg: Config, candidates: List[Candidate], gdna_qc: Dict, path: str
) -> None:
    """Write a MultiQC custom-content TSV that feeds the General Statistics table.

    MultiQC auto-detects files matching ``*_mqc.tsv``. The commented YAML header
    configures the columns; the single data row is keyed by the sample name, so
    MultiQC merges rows across all samples in a run.
    """
    sample = cfg.sample_name or os.path.basename(cfg.out_prefix) or "sample"
    n_rescued = sum(1 for c in candidates if c.kept)

    header = [
        "# id: 'extrna_gdna'",
        "# section_name: 'extRNA gDNA contamination'",
        "# description: 'gDNA-like signal in unannotated regions and rescued novel transcripts (extRNA).'",
        "# plot_type: 'generalstats'",
        "# pconfig:",
        "#     - extrna_pct_gDNA:",
        "#         title: '% gDNA'",
        "#         description: 'Percent of uniquely-mapped coverage in gDNA-flagged unannotated regions'",
        "#         min: 0",
        "#         max: 100",
        "#         suffix: '%'",
        "#         scale: 'OrRd'",
        "#         format: '{:,.2f}'",
        "#     - extrna_gDNA_regions:",
        "#         title: 'gDNA regions'",
        "#         description: 'Number of unannotated regions classified as likely gDNA'",
        "#         format: '{:,.0f}'",
        "#     - extrna_novel_transcripts:",
        "#         title: 'novel transcripts'",
        "#         description: 'Rescued candidate novel transcripts'",
        "#         format: '{:,.0f}'",
    ]
    cols = ["Sample", "extrna_pct_gDNA", "extrna_gDNA_regions", "extrna_novel_transcripts"]
    row = [
        sample,
        f"{gdna_qc['pct_gDNA_of_mapped_coverage']}",
        f"{gdna_qc['n_gDNA_regions']}",
        f"{n_rescued}",
    ]
    with open(path, "w") as fh:
        fh.write("\n".join(header) + "\n")
        fh.write("\t".join(cols) + "\n")
        fh.write("\t".join(row) + "\n")


def write_summary_json(
    cfg: Config,
    candidates: List[Candidate],
    strandedness_metrics: Dict,
    path: str,
    gdna_qc: Dict | None = None,
) -> Dict:
    """Write the run summary and return it."""
    n_gdna = sum(1 for c in candidates if c.label == LIKELY_GDNA)
    n_bidir = sum(1 for c in candidates if c.label == POSSIBLE_BIDIRECTIONAL)
    n_novel = sum(1 for c in candidates if c.label == LIKELY_NOVEL)
    n_multi = sum(1 for c in candidates if c.label == LIKELY_MULTIMAPPER)
    rescued = [c for c in candidates if c.kept]
    total_rescued_bases = sum(c.length for c in rescued)

    summary = {
        "run_parameters": cfg.to_serialisable(),
        "inferred_library_strandedness": cfg.inferred_strandedness,
        "strandedness_metrics": strandedness_metrics,
        "n_candidate_regions": len(candidates),
        "n_likely_gDNA": n_gdna,
        "n_possible_bidirectional_RNA": n_bidir,
        "n_likely_novel_transcript": n_novel,
        "n_likely_multimapper_artifact": n_multi,
        "n_rescued_unknown_transcripts": len(rescued),
        "total_bases_in_rescued_unknown_transcripts": total_rescued_bases,
        "gdna_contamination_qc": gdna_qc or {},
    }
    with open(path, "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    return summary
