"""End-to-end orchestration: BAM + GTF -> candidate table, GTFs, summary.

Chromosome-wise processing keeps memory bounded (one chromosome's coverage
arrays at a time). With --threads > 1 chromosomes are processed in parallel
processes, each opening its own BAM handle.
"""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List

from .config import Config
from .discovery import Candidate, build_candidates_for_chrom
from .gtf_io import Annotation, parse_gtf, read_chrom_sizes
from .utils import get_logger
from . import writers


def _process_one_chrom(
    bam_path: str, chrom: str, length: int, cfg: Config,
    strandedness: str, annotation: Annotation,
):
    """Worker: build strand coverage for one chrom and discover candidates.

    Returns (candidates, total_unique_coverage, total_multi_coverage) where the
    totals are genome-wide mapped-base counts used for the library-level gDNA
    contamination percentage (summing the arrays is essentially free).
    """
    from .bam_io import strand_coverage_for_chrom  # local import: pysam only here

    plus, minus, multi = strand_coverage_for_chrom(
        bam_path, chrom, length, cfg, strandedness
    )
    total_unique = int(plus.sum()) + int(minus.sum())
    total_multi = int(multi.sum())
    candidates = build_candidates_for_chrom(
        chrom, plus, minus, multi, annotation, cfg
    )
    return candidates, total_unique, total_multi


def _assign_names(candidates: List[Candidate]) -> None:
    """Assign deterministic region ids and unknown_transcript_N names.

    Region ids are chrom:start-end already; we number kept transcripts in the
    order they appear (chrom order, then coordinate).
    """
    counter = 0
    for c in candidates:
        if c.kept:
            counter += 1
            c.unknown_transcript_name = f"unknown_transcript_{counter}"


def run(cfg: Config) -> Dict:
    """Run the whole pipeline and return the summary dict."""
    logger = get_logger(cfg.verbose)
    cfg.validate()

    for label, path in (("BAM", cfg.bam), ("GTF", cfg.gtf)):
        if not path or not os.path.exists(path):
            raise FileNotFoundError(f"{label} file not found: {path!r}")
    if cfg.fai and not os.path.exists(cfg.fai):
        raise FileNotFoundError(f"--fai file not found: {cfg.fai!r}")

    logger.info("Parsing GTF (%s, annotation-mode=%s) ...", cfg.gtf, cfg.annotation_mode)
    annotation = parse_gtf(cfg.gtf, cfg.annotation_mode)

    # Chromosome sizes from BAM header (authoritative for fetch bounds).
    from .bam_io import get_chrom_sizes_from_bam
    chrom_sizes = get_chrom_sizes_from_bam(cfg.bam)

    if cfg.fai:
        fai_sizes = read_chrom_sizes(cfg.fai)
        for chrom, blen in chrom_sizes.items():
            flen = fai_sizes.get(chrom)
            if flen is not None and flen != blen:
                logger.warning(
                    "Chromosome %s length mismatch: BAM=%d, fai=%d. Using BAM.",
                    chrom, blen, flen,
                )

    # Resolve strandedness (auto inference or explicit).
    from .strandedness import resolve_strandedness
    strandedness, strand_metrics = resolve_strandedness(cfg.bam, annotation, cfg)
    cfg.inferred_strandedness = strandedness
    cfg.strandedness_metrics = strand_metrics

    chroms = [c for c in chrom_sizes if chrom_sizes[c] > 0]
    logger.info("Processing %d chromosomes (threads=%d) ...", len(chroms), cfg.threads)

    candidates: List[Candidate] = []
    genome_unique_cov = 0
    genome_multi_cov = 0
    if cfg.threads > 1 and len(chroms) > 1:
        with ProcessPoolExecutor(max_workers=cfg.threads) as ex:
            futures = {
                ex.submit(
                    _process_one_chrom, cfg.bam, chrom, chrom_sizes[chrom],
                    cfg, strandedness, annotation,
                ): chrom
                for chrom in chroms
            }
            results: Dict[str, tuple] = {}
            for fut in as_completed(futures):
                chrom = futures[fut]
                results[chrom] = fut.result()
                logger.debug("Finished %s: %d candidates", chrom, len(results[chrom][0]))
        # Reassemble in header order for deterministic numbering.
        for chrom in chroms:
            cands, tu, tm = results.get(chrom, ([], 0, 0))
            candidates.extend(cands)
            genome_unique_cov += tu
            genome_multi_cov += tm
    else:
        for chrom in chroms:
            chrom_cands, tu, tm = _process_one_chrom(
                cfg.bam, chrom, chrom_sizes[chrom], cfg, strandedness, annotation
            )
            logger.debug("Finished %s: %d candidates", chrom, len(chrom_cands))
            candidates.extend(chrom_cands)
            genome_unique_cov += tu
            genome_multi_cov += tm

    # Sort within-chrom by start (header order already applied across chroms).
    candidates.sort(key=lambda c: (chroms.index(c.chrom), c.start, c.end))
    _assign_names(candidates)

    gdna_qc = writers.compute_gdna_qc(candidates, genome_unique_cov, genome_multi_cov)

    # --- write outputs ----------------------------------------------------
    tsv = f"{cfg.out_prefix}.candidate_regions.tsv"
    unknown_gtf = f"{cfg.out_prefix}.unknown_transcripts.gtf"
    merged_gtf = f"{cfg.out_prefix}.annotation_plus_unknowns.gtf"
    summary_json = f"{cfg.out_prefix}.summary.json"

    writers.write_tsv(candidates, tsv)
    n_unknown = writers.write_unknown_gtf(candidates, unknown_gtf)
    writers.write_merged_gtf(cfg.gtf, candidates, merged_gtf)
    if cfg.emit_bed:
        writers.write_bed(candidates, f"{cfg.out_prefix}.candidate_regions.bed")
    if cfg.emit_bedgraph:
        writers.write_bedgraph(candidates, cfg.out_prefix)
    summary = writers.write_summary_json(
        cfg, candidates, strand_metrics, summary_json, gdna_qc=gdna_qc
    )
    if cfg.emit_multiqc:
        writers.write_multiqc_tsv(
            cfg, candidates, gdna_qc, f"{cfg.out_prefix}.gdna_mqc.tsv"
        )

    logger.info(
        "Done. %d candidate regions | %d likely_gDNA | %d bidirectional | "
        "%d novel | %d multimapper_artifact | %d rescued transcripts.",
        summary["n_candidate_regions"],
        summary["n_likely_gDNA"],
        summary["n_possible_bidirectional_RNA"],
        summary["n_likely_novel_transcript"],
        summary["n_likely_multimapper_artifact"],
        n_unknown,
    )
    logger.info(
        "gDNA-like signal: %.2f%% of mapped unique coverage (%.2f%% of candidate "
        "coverage).",
        gdna_qc["pct_gDNA_of_mapped_coverage"],
        gdna_qc["pct_gDNA_of_candidate_coverage"],
    )
    logger.info("Outputs written with prefix: %s", cfg.out_prefix)
    return summary
