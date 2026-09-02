# Changelog

## 1.1.0 — 2026-09-02

### Changed (affects default results)

- **New pair-geometry filter, on by default** (`--pair-filter concordant`): a
  paired read whose mate is unmapped (half-mapped) or whose pair lacks the
  aligner's proper-pair flag (discordant / wrong orientation) is now classified
  as noise (`multi`), never as `unique` evidence, regardless of MAPQ. Such
  pairs are enriched for adapter chimeras, assembly-gap edges and contaminant
  fragments and must not define novel-transcription regions. Applied at the
  single `read_category()` chokepoint, so region discovery, strandedness
  inference, per-region read counting and the coverage store all inherit it
  consistently.
  - Decided from the **proper-pair FLAG**, never an insert-size/TLEN cutoff:
    splice-aware aligners (STAR/HISAT2) flag intron-spanning pairs as proper
    even when genomic TLEN is 10–100 kb, so a TLEN cutoff would bias against
    long-intron novel transcripts.
  - **gDNA detection is unaffected**: gDNA fragments are ordinary proper
    pairs, so gDNA evidence still flows into region classification.
  - Single-end reads are unaffected. `--pair-filter off` restores the pre-1.1
    behavior exactly.
- Candidate regions dominated by half-mapped/discordant pairs now fail the
  `--min-unique-fraction` test and are flagged `likely_multimapper_artifact`
  (or are not discovered at all when no concordant unique coverage remains).

### Added

- `summary.json` now carries a `read_totals` block with genome-wide read
  counts, including `n_half_mapped_reads` / `n_discordant_reads` — the
  would-otherwise-be-unique reads the pair filter reclassified — and the
  per-sample log reports the same counts. The MultiQC bargraph categories
  (`extrna_read_assignment`) are unchanged so reports merge across tool
  versions.

## 1.0.0

- Initial release.
