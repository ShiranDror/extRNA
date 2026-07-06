# detect_gdna_vs_novel

Detect likely **genomic-DNA (gDNA) contamination** in RNA-seq alignments and
**rescue candidate unannotated transcripts** into a new GTF.

The tool is deliberately **not** a splice-junction-centric RNA detector. It
respects the fact that genuine RNA can be continuous, unspliced, intronic,
antisense, intergenic, single-exon, low-complexity or bidirectional. Instead of
asking "does this look like a spliced mRNA?", it asks:

> Over this continuously-covered, unannotated interval, does the signal look like
> it came from **double-stranded genomic DNA** (both strands, symmetric, broad,
> correlated/flat) or from **transcription** (dominantly one strand, coherent)?

---

## Why this tool exists

As we process more **rRNA-depleted** RNA-seq libraries (total-RNA protocols,
rather than poly(A) selection), a growing share of the mapped signal falls into
**intronic and intergenic** regions. The reflexive interpretation of that signal
is *genomic-DNA contamination* — and standard QC that simply reports
"% intronic / % intergenic" reinforces that reading, treating anything outside
annotated exons as suspect.

In our data that interpretation is frequently **wrong**. When we actually look at
these unannotated regions, they often show hallmarks of genuine transcription
rather than gDNA:

- **Strandedness** — coverage resolves predominantly to one transcription
  strand, whereas double-stranded genomic DNA would contribute to both strands
  symmetrically.
- **Continuous coverage** — coherent, contiguous signal over a defined interval,
  consistent with a transcribed unit rather than the broad, uniform smear of
  randomly-fragmented DNA.

rRNA depletion is exactly the condition that *exposes* this material: it retains
non-poly(A), unspliced, intronic, antisense and other non-canonical RNA that
poly(A) selection discards, and much of it is simply **unannotated** — not
contamination. Discarding it wholesale as "intronic/intergenic gDNA" throws away
real biology.

So this tool replaces the coarse "intronic/intergenic ⇒ suspect" heuristic with
an explicit, per-region test of what the signal actually looks like. It uses
**strand asymmetry and continuous coverage as positive evidence of RNA**, flags
only the regions that genuinely look like double-stranded genomic DNA (symmetric,
broad, correlated-or-flat on both strands) or multimapper artifacts, and
**rescues the rest as candidate novel transcripts** instead of silently losing
them. Contamination annotation from a browser context label becomes a decision
you can audit per region.

---

## What it does

1. Reads a coordinate-sorted, indexed **BAM** (STAR output assumed) and a
   reference **GTF**.
2. Infers **library strandedness** (`auto`), or takes it from you.
3. Builds **strand-specific coverage** per chromosome (spliced reads contribute
   only their exonic blocks, so introns are not spuriously "covered").
4. Masks annotated features (`exon` / `transcript` / `gene` / `all`) and
   discovers **continuous unannotated covered intervals**.
5. Computes a transparent set of **per-region metrics** and applies a
   **rule-based classifier**:
   - `likely_multimapper_artifact` — dropped (repeat/alignment artifact),
   - `likely_gDNA` — dropped (contamination),
   - `likely_novel_transcript` — rescued,
   - `possible_bidirectional_RNA` — rescued (not gDNA-like).
6. Writes rescued regions as `unknown_transcript_N` features into a new GTF,
   plus a QC table, a merged GTF, a summary JSON, and (optionally) BED/bedGraph.
7. **(Optional, multi-sample)** `merge_candidates.py` takes the per-sample QC
   tables from N samples and keeps only loci **reproduced in ≥ x samples** —
   genuine novel transcripts should recur across replicates. See
   [Cross-sample consensus](#cross-sample-consensus).

---

## Installation

### Dependencies
- Python ≥ 3.9
- numpy, pandas, scipy — install anywhere.
- **pysam** — required for BAM reading. **pysam has no Windows wheels.**

### Recommended: conda (Linux / macOS / WSL)
```bash
conda env create -f environment.yml
conda activate gdna_rescue
```

### pip (Linux / macOS)
```bash
pip install -r requirements.txt
```

### Windows
`pysam` cannot be pip-installed on native Windows. Use **WSL** (Windows Subsystem
for Linux) or conda inside WSL. The pure-analysis modules and the unit tests
(everything except BAM reading) *do* run on native Windows:

```powershell
pip install numpy pandas scipy pytest
python -m pytest -q          # 22 pysam-free tests run; integration tests skip
```

---

## Usage

```bash
python detect_gdna_vs_novel.py \
  --bam sample.bam \
  --gtf reference.gtf \
  --out-prefix sample_analysis \
  --library-strandedness auto
```

The BAM must be coordinate-sorted and indexed:
```bash
samtools sort -o sample.bam aligned.bam
samtools index sample.bam
```

### Key options
```
--min-mapq 20                 # STAR: 255=unique, 3/1/0=multimapper; 20 keeps unique only
--min-baseq 0                 # >0 enables (slower) per-base quality filtering
--min-unique-fraction 0.5     # region kept only if >= this fraction of coverage is uniquely mapped
--no-count-secondary          # don't count secondary alignments in the multimapper track
--min-depth 10                # combined (both-strand) per-base depth for a base to count as covered
--strand-min-depth 3          # per-strand per-base depth for a base to be "supported"
--max-gap 50                  # merge covered runs across <= N unannotated bases
--min-region-length 200       # discard shorter candidate intervals
--min-covered-bases 100       # discard regions with fewer covered bases
--min-covered-fraction 0.7    # discard sparse/punctate regions
--annotation-mode exon        # exon | transcript | gene | all
--library-strandedness auto   # auto | forward | reverse | unstranded
--nearest-feature-window 10000
--threads 4                   # chromosome-level parallelism
--no-bed                      # skip the BED output
--emit-bedgraph               # write per-strand candidate bedGraph
--no-multiqc                  # skip the *.gdna_mqc.tsv MultiQC file
--sample-name S1              # sample label in the MultiQC row (default: out-prefix basename)
--verbose
```

### gDNA / novel classification thresholds (all tunable)
```
--gdna-min-dual-strand-fraction 0.60      # fraction of interval covered on BOTH strands
--gdna-max-strand-length-ratio-diff 0.25  # |plus_len - minus_len| / max
--gdna-min-profile-correlation 0.70       # Pearson r of per-base +/- profiles
--gdna-min-covered-fraction 0.50          # gDNA is broadly continuous
--gdna-min-depth-balance 0.50             # min/max strand mean depth
--gdna-flat-cv-threshold 0.40             # flat-on-both-strands alternative to correlation
--novel-min-dominant-strand-fraction 0.80 # one strand carries >= this => novel
```

---

## Outputs (`--out-prefix sample_analysis`)

| File | Contents |
|------|----------|
| `sample_analysis.candidate_regions.tsv` | Per-region metrics + class + reason |
| `sample_analysis.unknown_transcripts.gtf` | Rescued `unknown_transcript_N` features |
| `sample_analysis.annotation_plus_unknowns.gtf` | Original GTF + rescued features |
| `sample_analysis.summary.json` | Run parameters, strandedness, class counts |
| `sample_analysis.candidate_regions.bed` | BED6 of all candidates (optional) |
| `sample_analysis.candidates.{plus,minus}.bedgraph` | Per-strand depth (with `--emit-bedgraph`) |
| `sample_analysis.gdna_mqc.tsv` | **MultiQC** bargraph TSV — read counts per class (annotated / novel / gDNA / …), normalisable to % in MultiQC |

**Coordinate conventions:** TSV and GTF are 1-based inclusive; BED/bedGraph are
0-based half-open.

### MultiQC integration

Each per-sample run writes `sample_analysis.gdna_mqc.tsv`, a MultiQC
custom-content **bargraph** of **read assignment**. Point MultiQC at your output
directory and a stacked bar per sample appears in the report:

```bash
multiqc .    # discovers *_mqc.tsv from all samples automatically
```

The bar reports **raw uniquely-mapped read counts** per region class, so
MultiQC's built-in **counts / percentages toggle** does the normalisation:

| Category | Meaning |
|---|---|
| `annotated` | reads over existing annotated features |
| `novel_transcript` | reads in rescued `likely_novel_transcript` regions |
| `bidirectional_RNA` | reads in rescued `possible_bidirectional_RNA` regions |
| `gDNA` | reads in `likely_gDNA` regions |
| `multimapper_artifact` | reads in `likely_multimapper_artifact` regions |
| `other_unannotated` | unannotated reads that formed no candidate region |

The categories partition the uniquely-mapped reads (they sum to the sample
total). Counts are assigned by read midpoint: annotated reads and
candidate-region reads are disjoint because candidate regions contain no
annotated positions. The same counts, plus the coverage-based gDNA percentages,
are in `summary.json` under `read_assignment_counts` and `gdna_contamination_qc`.

> gDNA is only tested among *unannotated* candidate regions (annotated exons are
> not), so the `gDNA` bar reflects contamination surfacing as novel-looking
> signal, not total genomic DNA in the library.

Each rescued transcript carries attributes:
`gene_id "unknown_transcript_N_gene"`, `transcript_id "unknown_transcript_N"`,
`gene_name "unknown_transcript_N"`, `source "gdna_rescue"`,
`classification "..."`, `context "..."`, `original_region_id "..."`.

---

## Classification logic (transparent, rule-based)

For every unannotated candidate interval the tool computes per-base plus/minus
coverage and derives (all in the TSV): length, total/average/max depth, covered
bases and fraction, plus/minus covered lengths and mean depths, plus/minus depth
ratio, covered-length ratio and difference, dominant-strand fraction, strand
entropy (balance), strand-overlap Jaccard, **dual-strand fraction**, **per-base
profile correlation**, per-strand coefficient of variation, the **uniquely-mapped
fraction**, plus context (intergenic / intronic / antisense / near-gene) and
nearest annotated feature.

The decision tree (applied in order):

0. **Multimapper / repeat artifact → `likely_multimapper_artifact` (dropped).**
   Regions are *discovered* from uniquely-mapped reads (MAPQ ≥ `--min-mapq`;
   STAR unique = 255), but multimapped reads (MAPQ 3/1/0 and secondary
   alignments) are tracked separately as a noise signal. If the uniquely-mapped
   fraction of a region's coverage is below `--min-unique-fraction` (default
   0.50), the locus is swamped by multimappers and is almost certainly a repeat
   or alignment artifact — it is dropped regardless of its strand pattern. Local
   multimapper-only stretches are tolerated as long as the region as a whole is
   majority-unique. This is checked **first** so an artifact can never be rescued
   as a novel transcript. (The `unique_fraction` metric is in the TSV.)

1. **Single dominant strand → `likely_novel_transcript`.**
   If one strand carries ≥ `--novel-min-dominant-strand-fraction` (default 0.80)
   of the signal, the locus is transcribed from one template strand. This is the
   RNA signature and is checked first, so genuine single-strand RNA (even
   continuous/unspliced/antisense/intronic) is never mistaken for gDNA.

2. **Symmetric, balanced, broad dual-strand + consistent pattern → `likely_gDNA`.**
   Requires *all* of:
   - dual-strand fraction ≥ `--gdna-min-dual-strand-fraction`,
   - strand covered-length difference ≤ `--gdna-max-strand-length-ratio-diff`,
   - covered fraction ≥ `--gdna-min-covered-fraction` (broad, continuous),
   - strand depth balance ≥ `--gdna-min-depth-balance`,
   - **and** the per-base pattern is consistent with dsDNA, which means *either*
     the plus/minus profiles are correlated (≥ `--gdna-min-profile-correlation`)
     *or* both strands are **flat** (CV ≤ `--gdna-flat-cv-threshold`).

   > **Why the "flat OR correlated" rule matters:** randomly-fragmented genomic
   > DNA produces roughly **uniform** coverage on both strands. Two flat noisy
   > profiles have almost no variance to correlate, so a correlation-only rule
   > would *miss the most common gDNA case*. A flat, symmetric, balanced,
   > broad dual-strand profile is itself a strong gDNA signature. Correlation
   > catches the other case: gDNA whose depth varies with local mappability, so
   > both strands rise and fall together.

3. **Otherwise → `possible_bidirectional_RNA`.**
   Both strands contribute but the pattern is not gDNA-symmetric (asymmetric
   lengths, unbalanced depth, or offset/anti-correlated non-flat profiles). This
   is kept as a candidate transcript because real loci can be bidirectional.

**`likely_gDNA` and `likely_multimapper_artifact` are discarded.** Everything
else is rescued and numbered `unknown_transcript_1, 2, …`. The exact metric
values and which rule fired are written to `reason_for_classification` in the
TSV, so every call is auditable.

The classifier lives in `gdna_rescue/classify.py` and takes plain numpy arrays,
so a statistical/ML model can be dropped in later behind the same interface.

---

## STAR recommendations

- Default `--min-mapq 20` keeps STAR **uniquely-mapped** reads (MAPQ 255) and
  drops multimappers (MAPQ 3/1/0), which is usually what you want for
  contamination assessment.
- **2-pass mapping** (`--twopassMode Basic`) improves novel splice-junction
  detection and reduces spurious intronic coverage from misalignment, which in
  turn reduces false novel-transcript calls. Recommended when annotation is
  incomplete.
- Keep unsorted-vs-sorted straight: this tool needs a **coordinate-sorted,
  indexed** BAM (`--outSAMtype BAM SortedByCoordinate`, then `samtools index`).
- If you filtered the BAM upstream, make sure the header sort order is still
  `SO:coordinate`.

---

## Examples

Auto strandedness, exon-masking, 8 threads:
```bash
python detect_gdna_vs_novel.py --bam s.bam --gtf ref.gtf \
  --out-prefix s --library-strandedness auto --threads 8
```

Reverse (dUTP) library, mask whole gene spans, stricter discovery:
```bash
python detect_gdna_vs_novel.py --bam s.bam --gtf ref.gtf --out-prefix s \
  --library-strandedness reverse --annotation-mode gene \
  --min-depth 5 --min-region-length 200 --min-covered-fraction 0.7
```

More permissive gDNA calling (flag more contamination):
```bash
python detect_gdna_vs_novel.py --bam s.bam --gtf ref.gtf --out-prefix s \
  --gdna-min-dual-strand-fraction 0.5 --gdna-min-profile-correlation 0.6
```

---

## Cross-sample consensus

`merge_candidates.py` combines the per-sample `*.candidate_regions.tsv` files
from several samples and keeps loci reproduced in **at least `--min-samples`**
samples. Genuine novel transcripts should recur across biological replicates,
while contamination and one-off alignment artifacts tend to be sample-specific —
so this doubles as a noise filter and as independent evidence for the per-sample
calls.

```bash
python merge_candidates.py \
  --tsv A.candidate_regions.tsv B.candidate_regions.tsv \
        C.candidate_regions.tsv D.candidate_regions.tsv \
  --reference-gtf reference.gtf \
  --out-prefix cohort \
  --min-samples 2
```

It is pure polars/Python — **no pysam** — so it runs natively anywhere,
including Windows.

Passing `--reference-gtf` writes `cohort.reference_plus_consensus.gtf`: the
**analysis-ready annotation** (reference genes + reproducible novel transcripts,
with feature IDs consistent across all samples). Run `featureCounts` on the
**original STAR BAMs** against that single GTF — no manual concatenation needed.

### How matching works
- Candidates never share exact coordinates across samples, so loci are matched by
  **reciprocal overlap**: two candidates cluster only if each covers ≥
  `--min-reciprocal-overlap` (default **0.85**) of the other. The default is
  deliberately high — "the same transcript" should mean nearly co-extensive
  intervals, not a 50% touch.
- Matching is **strand-aware by default** (this is RNA); use `--ignore-strand`
  to disable.
- Each cluster's class is a **majority vote** across samples, ties broken toward
  the more conservative (reject) call.

### Consensus classes and what reaches the GTF
| Majority call | Consensus class | In consensus GTF? |
|---|---|---|
| `likely_novel_transcript` | `reproducible_novel` | **yes** |
| `possible_bidirectional_RNA` | `reproducible_bidirectional` | no (unless `--include-bidirectional`) |
| `likely_gDNA` | `recurrent_gDNA` | **no** — reported for manual review |
| `likely_multimapper_artifact` | `recurrent_multimapper_artifact` | **no** |

Recurrent gDNA is odd and worth investigating (ideally confirmed biochemically,
e.g. DNase treatment) — but it is **not** added to the annotation. The guiding
principle is conservative: **losing a true annotation is preferable to adding a
bad one**, so only reproducible novel transcripts are written to the consensus
GTF by default.

### Options
```
--tsv A.tsv B.tsv ...        # per-sample candidate_regions.tsv files (required)
--sample-names A B ...        # optional; defaults to filenames
--out-prefix cohort           # required
--min-samples 2               # keep loci present in >= this many samples
--min-reciprocal-overlap 0.85 # each candidate must cover >= this fraction of the other
--reference-gtf reference.gtf # if given, also write the analysis-ready reference+consensus GTF
--ignore-strand               # match regardless of strand (default: strand-aware)
--include-bidirectional       # also add reproducible bidirectional loci to the GTF
```

### Outputs (`--out-prefix cohort`)
| File | Contents |
|------|----------|
| `cohort.consensus_regions.tsv` | Reproducible clusters: consensus class, n_samples, per-sample classes, union coordinates, mean metrics, provenance |
| `cohort.consensus_transcripts.gtf` | Reproducible novel loci as `consensus_transcript_N` (union span; carries `n_samples`, `samples`, `member_region_ids`) |
| `cohort.reference_plus_consensus.gtf` | **(with `--reference-gtf`)** reference + consensus — the analysis-ready GTF for featureCounts on the original STAR BAMs |
| `cohort.consensus_summary.json` | Parameters and counts per consensus class |

Consensus coordinates use the **union span** of the clustered members.

> **Caveat:** reproducibility filtering will discard genuinely sample-specific
> biology (e.g. a transcript induced in a single condition). It is the right
> tool for finding *robust* novel transcripts, not a complete catalogue, and it
> assumes consistent chromosome naming across samples.

---

## Recommended end-to-end pipeline

```
FastQC / fastp (trim)
  -> STAR  (2-pass; KEEP multimappers; SortedByCoordinate + samtools index)
  -> extRNA  detect_gdna_vs_novel.py   (per sample)
  -> extRNA  merge_candidates.py --reference-gtf reference.gtf   (cohort consensus)
       => cohort.reference_plus_consensus.gtf
  -> featureCounts  (ORIGINAL STAR BAMs, matched strandedness)
  -> edgeR / DESeq2  (differential expression)

  (aggregate QC across samples with `multiqc .` — the per-sample
   *.gdna_mqc.tsv files add a stacked read-assignment bar: annotated vs
   novel vs gDNA, with a counts/percentage toggle)
```

Points that matter for correct results:

- **Keep multimappers in the STAR BAM.** extRNA's `likely_multimapper_artifact`
  detection needs them; if the BAM is pre-filtered to unique-only, that check
  silently does nothing. Let extRNA separate unique vs multi by MAPQ.
- **Run the consensus step before featureCounts.** extRNA runs per sample and
  each sample's `unknown_transcript_N` differ; `merge_candidates.py` collapses
  them to reproducible `consensus_transcript_N` with IDs consistent across
  samples. Use `--reference-gtf` to get the single analysis-ready GTF directly.
- **Match strandedness across tools.** Use extRNA's inferred strandedness (from
  `summary.json`) for featureCounts `-s` (`-s 1` forward, `-s 2` reverse). This
  is critical for the antisense/intronic novel features — a wrong `-s` miscounts
  exactly the loci this tool rescues.
- **Replicates + low-count filtering** at the DE step; treat novel features as
  exploratory (approximate single-exon models — sanity-check top hits in IGV).

---

## Testing

```bash
# Pure-numpy logic (runs anywhere, no pysam):
python -m pytest -q tests/test_classify.py tests/test_discovery.py tests/test_strandedness.py

# Full suite incl. end-to-end integration (needs pysam):
python -m pytest -q

# See the classifier decide on the three synthetic archetypes:
python -m tests.generate_test_data
```

The synthetic generator builds three archetypes — a symmetric gDNA region, a
single-strand novel transcript, and an asymmetric bidirectional region — and the
integration test asserts they are classified and rescued as expected.

---

## Limitations (please read)

- **Real bidirectional transcription exists** (e.g. promoters, eRNAs). Such loci
  can resemble gDNA. The tool labels ambiguous both-strand loci
  `possible_bidirectional_RNA` and keeps them rather than discarding them, but
  truly symmetric bidirectional RNA *can* be mislabelled `likely_gDNA`.
- **Unspliced / continuous RNA exists** and is fully supported — the tool does
  not require splice junctions.
- **Antisense and intronic transcription exist**; these are reported via the
  `context_label` column and are rescued when single-strand-dominant.
- **Incomplete GTF annotation inflates unannotated signal.** Genuine but
  unannotated genes will appear as candidates (that is partly the point), but a
  sparse annotation will produce many candidates; use a complete, matched
  annotation.
- **Unstranded libraries cannot be assessed for gDNA by strand symmetry.** In an
  unstranded library even genuine single-strand RNA maps to both strands ~50/50,
  so the core discriminator is uninformative. The tool warns loudly and you
  should not trust `likely_gDNA` calls on unstranded data.
- **Multimapper filtering is coverage-based, not locus-resolved.** A region is
  flagged when multimapped reads dominate its coverage; the tool does not resolve
  where those reads truly originate. Genuine transcripts from recently-duplicated
  gene families can therefore be flagged as artifacts — tune `--min-unique-fraction`
  if your biology involves such loci.
- **Cross-sample reproducibility discards sample-specific biology.** A transcript
  induced in only one condition/sample will not survive the `--min-samples`
  filter. Use the consensus step to find *robust* transcripts, not a full catalogue.
- This tool **estimates likely gDNA-like regions; it does not prove DNA origin.**
  Orthogonal evidence (e.g. RNase/DNase treatment, intron-retention patterns,
  qPCR) is needed for confirmation.
- Coverage arrays are built one chromosome at a time. Peak memory scales with the
  largest chromosome (≈ 2 × int32 × chrom length). For human chr1 that is roughly
  2 GB; reduce `--threads` if memory-constrained (each worker holds one
  chromosome).

---

## Nice-to-have / future work

- BED12 export of multi-exon models (currently single-exon models are emitted).
- Per-base bigWig strand coverage.
- Plotting of candidate strand profiles.
- Blacklist of problematic genomic regions (e.g. ENCODE blacklist) to pre-filter
  candidates.
- Pluggable statistical classifier behind `classify_region`.

---

## Package layout

```
detect_gdna_vs_novel.py     # thin CLI entry point (per-sample analysis)
merge_candidates.py         # thin CLI entry point (cross-sample consensus)
gdna_rescue/
  config.py        # all tunable thresholds (Config dataclass)
  utils.py         # logging + interval maths (no pysam)
  gtf_io.py        # GTF parsing / masking / gene index
  classify.py      # per-region metrics + rule-based classifier (numpy only)
  discovery.py     # unannotated region discovery + context labelling (numpy only)
  bam_io.py        # the ONLY pysam-dependent module
  strandedness.py  # automatic strandedness inference
  writers.py       # TSV / GTF / JSON / BED / bedGraph
  pipeline.py      # orchestration (chromosome-wise, optional multiprocessing)
  crosssample.py   # cross-sample consensus / reproducibility filter (polars only)
  cli.py           # argument parsing
tests/
  generate_test_data.py   # synthetic archetypes + synthetic BAM/GTF
  test_classify.py        # classifier unit tests (incl. multimapper artifact)
  test_discovery.py       # discovery/merging unit tests
  test_strandedness.py    # read -> strand mapping unit tests
  test_crosssample.py     # cross-sample consensus unit tests (polars)
  test_pipeline.py        # end-to-end integration (needs pysam)
```
