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

> **Scope: this tool is for STRANDED RNA-seq libraries** (forward /
> fr-secondstrand, or reverse / fr-firststrand / dUTP). The whole method rests on
> comparing the two transcription strands, which only carries information when
> the library preserves strand. **Unstranded libraries are out of scope** — in an
> unstranded library even genuine single-strand RNA maps to both strands ~50/50,
> so gDNA-vs-RNA cannot be told apart by strand symmetry. The tool will run on
> unstranded data and warn, but its calls are not meaningful there; don't use it
> for unstranded libraries.

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
python -m pytest -q          # 78 pysam-free tests run; integration tests skip
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
--annotation-mode exon        # what counts as "annotated" (see Annotation modes below)
--no-stranded-masking         # force positional masking (default: per-strand for stranded libs)
--library-strandedness auto   # auto | forward | reverse  (unstranded is out of scope)
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

### Annotation modes

`--annotation-mode` sets what is treated as "annotated" and therefore **masked
out** before discovery — i.e. what is *not* eligible to become a novel candidate:

| Mode | Masks | Introns discoverable? |
|------|-------|-----------------------|
| `exon` (default) | exon features only | **yes** — introns are unmasked, so intronic (e.g. pre-mRNA / retained-intron) signal can be rescued |
| `transcript` | full transcript spans (exons + introns) | no — the whole transcript body is masked |
| `gene` | full gene spans (exons + introns) | no — the whole gene body is masked |
| `all` | **every feature line in the GTF** (gene, transcript, exon, CDS, UTR, codons, …) | no — because gene/transcript lines span the locus, `all` masks the entire gene body, introns included |

So `all` is the **most aggressive** mask: any base mentioned by any feature is
excluded. Since gene/transcript records span introns, `all` behaves like `gene`
(introns masked) plus every sub-feature — use it when you only want to discover
signal in regions the GTF never mentions at all (strictly intergenic). Use `exon`
(the default) if you want to recover intronic/pre-mRNA and antisense-in-intron
signal; use `gene`/`transcript`/`all` if you consider anything inside a gene body
to be "already annotated."

**Stranded vs positional masking.** For stranded libraries the mask is applied
**per strand** by default: a `+` feature masks only the `+` strand, so
**antisense transcription lying directly over an annotated feature stays
discoverable** (and is labelled `antisense_to_gene`). Region metrics are then
computed on the masked coverage, so the host gene's sense signal over that
feature does not count toward the antisense candidate. This only applies to
stranded (forward/reverse) libraries; for unstranded input the tool
automatically falls back to **positional (strand-agnostic)** masking (a base
annotated on either strand is masked on both), and `--no-stranded-masking`
forces positional masking for any library. Under positional masking, antisense
over an annotated feature is not recovered — only antisense that does not overlap
a masked feature (e.g. intronic-antisense in `exon` mode, or intergenic).

> **Chasing intronic antisense transcription (e.g. an eRNA on an intronic
> enhancer)? Use `--annotation-mode gene`, not the default `exon`.**
> This is counter-intuitive, because `exon` mode masks *less*. But in `exon` mode
> the host gene's intronic pre-mRNA coverage stays in the data, and as soon as it
> clears `--min-depth` the **whole intron becomes one `+`-dominant candidate** that
> merely *contains* the antisense locus — so the antisense feature is no longer
> resolvable on its own, and its `context_label` comes back `intronic` rather than
> `antisense_to_gene`. In `gene` mode with stranded masking the host span is
> masked on its own strand only, which removes the competing sense signal and
> leaves the antisense locus as a clean single-strand candidate
> (dominant-strand fraction 1.0) at exactly its own coordinates, whatever the
> host's intronic depth. Both behaviours are pinned in
> `tests/test_discovery.py`. Note this only holds with stranded masking:
> `gene` mode plus `--no-stranded-masking` masks the intron on both strands and
> the locus is lost entirely.

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

# external annotation overlay (see below)
--annotate [LABEL=]FILE ...   # GFF3/GFF/GTF sources to annotate consensus coordinates with
--annotate-labels reg ...      # explicit column prefixes matching --annotate order
--annotate-feature-types enhancer promoter ...   # keep only these column-3 types
--annotate-nearest-window 10000                  # report nearest feature when nothing overlaps
--annotate-stranded           # require feature strand to match (default: strand-agnostic)
--no-annotate-names           # keep bare consensus_transcript_N IDs (default: add the type suffix)
```

### Outputs (`--out-prefix cohort`)
| File | Contents |
|------|----------|
| `cohort.consensus_regions.tsv` | Reproducible clusters: consensus class, n_samples, per-sample classes, union coordinates, mean metrics, provenance, **+ external-annotation columns** |
| `cohort.consensus_transcripts.gtf` | Reproducible novel loci as `consensus_transcript_N` — **`consensus_transcript_N-<feature_type>` when annotated** (union span; carries `n_samples`, `samples`, `member_region_ids`, `<label>_types`) |
| `cohort.reference_plus_consensus.gtf` | **(with `--reference-gtf`)** reference + consensus — the analysis-ready GTF for featureCounts on the original STAR BAMs |
| `cohort.consensus_summary.json` | Parameters and counts per consensus class, **+ the annotation cross-tab** |

Consensus coordinates use the **union span** of the clustered members.

> **Caveat:** reproducibility filtering will discard genuinely sample-specific
> biology (e.g. a transcript induced in a single condition). It is the right
> tool for finding *robust* novel transcripts, not a complete catalogue, and it
> assumes consistent chromosome naming across samples.

---

## Annotating consensus coordinates with external data

The consensus step produces a list of reproducible coordinates. Those
coordinates often fall on features that are already known — regulatory elements,
repeats, CAGE TSS peaks — and knowing that **upfront, in the consensus table**,
is far more useful than discovering it later in a browser.

`--annotate` takes one or more GFF3 / GFF / GTF files (optionally gzipped) and
attaches overlap information to every consensus region:

```bash
# Ensembl Regulatory Build (release 116, GRCh38):
#   https://ftp.ensembl.org/pub/release-116/regulation/homo_sapiens/GRCh38/annotation/
python merge_candidates.py \
  --tsv A.candidate_regions.tsv B.candidate_regions.tsv C.candidate_regions.tsv \
  --reference-gtf reference.gtf \
  --out-prefix cohort --min-samples 2 \
  --annotate regulatory=Homo_sapiens.GRCh38.regulatory_features.v116.gff3.gz \
             emar=Homo_sapiens.GRCh38.EMARs.v116.gff.gz
```

Each source contributes eight columns, prefixed with its label:

| Column | Contents |
|---|---|
| `<label>_n` | number of overlapping features |
| `<label>_types` | overlapping feature types with counts, e.g. `CTCF_binding_site:1,enhancer:1` |
| `<label>_ids` | feature IDs, e.g. `ENSR1_B33F;ENSR1_538P5` (capped at 10, then `+N_more`) |
| `<label>_overlap_bp` | bases of the region covered by any feature (union, not sum) |
| `<label>_overlap_frac` | that as a fraction of region length |
| `<label>_genes` | linked gene names — Ensembl **promoters** carry `gene_id`/`gene_name` |
| `<label>_nearest` | nearest feature when nothing overlaps, `type:id` |
| `<label>_nearest_distance` | gap in bp; `0` when overlapping, `-1` when nothing is in the window |

The label defaults to a short name derived from the filename (the Ensembl
`...regulatory_features.v116.gff3.gz` becomes `regulatory`); use `label=path` or
`--annotate-labels` to set it explicitly. Overlapping types are also written into
the consensus GTF as a `<label>_types` attribute, so the analysis-ready
annotation is self-documenting.

### Transcript names carry the annotation

When a consensus locus overlaps a feature, its transcript name gains that
feature type as a **suffix**, so the GTF — and every downstream featureCounts /
DE table built from it — says what the locus sits on without a lookup:

```
consensus_transcript_2-enhancer
consensus_transcript_9-promoter
consensus_transcript_1-CTCF_binding_site
consensus_transcript_17                     <- overlaps nothing, bare name kept
```

`gene_id` follows (`consensus_transcript_2-enhancer_gene`). Design points:

- **Suffix, not prefix**, so the `consensus_transcript_N` stem survives: IDs stay
  grouped by number and anything matching `^consensus_transcript_` keeps working.
  Base names are already unique, so appending can never collide.
- **The token is the first type in the `<label>_types` column** — i.e. the
  alphabetically-first overlapping type of the first `--annotate` source that
  hits. There is deliberately **no built-in list of feature types and no priority
  table**, so a new annotation source needs no code change. Source precedence is
  the order you pass `--annotate` in.
- Consequence worth knowing: alphabetical order means `CTCF_binding_site` beats
  `enhancer`, `open_chromatin_region` and `promoter` on loci that overlap several.
  If you would rather the names reflect the transcription-relevant elements, drop
  CTCF from consideration with
  `--annotate-feature-types enhancer promoter open_chromatin_region`.
- `--no-annotate-names` keeps bare IDs — use it when you need identifiers
  comparable with an earlier un-annotated run, since renaming changes the row
  labels of any count matrix built from the GTF.

### It annotates; it does not classify

The overlay runs **after** clustering and the class vote, and it can only add
columns — it never changes a call. That is deliberate:

- The classifier stays coverage-only and auditable. Every call remains explained
  by `reason_for_classification` alone.
- It keeps the annotation **independent of the calls**. "Our novel transcripts
  are enriched at enhancers" is only a finding if enhancer annotation was not an
  input to the calls; if it were, the enrichment would be guaranteed by
  construction and would mean nothing.

### Reading the numbers honestly

`consensus_summary.json` gains an `annotation` block with a **cross-tab of
consensus class against feature type**, plus each source's genomic footprint:

```json
"crosstab_passing": {
  "regulatory": {
    "by_consensus_class": {
      "reproducible_novel": { "n_regions": 24, "n_with_overlap": 24,
                              "enhancer": 12, "promoter": 7,
                              "CTCF_binding_site": 15, "open_chromatin_region": 7 },
      "recurrent_gDNA":     { "n_regions": 10 }
    }
  }
}
```

**A raw overlap count is not interpretable on its own.** The Ensembl Regulatory
Build is large: 237k enhancers covering ~121 Mb, i.e. roughly 4% of the genome
before you even account for candidate length. A meaningful fraction of *any*
interval set will touch an enhancer by chance.

Two things make it interpretable:

- **The cross-class contrast** (the built-in control). Every consensus class went
  through identical masking, discovery and clustering, so the difference between
  the `reproducible_novel` row and the `recurrent_gDNA` /
  `recurrent_multimapper_artifact` rows is informative in a way that either row
  alone is not.
- **`covered_bp_by_type`** in the summary, which is each feature type's union
  footprint — the input you need to build a proper null.

Note that a genome-wide uniform null would be *wrong* here: candidates are not
drawn uniformly from the genome. Annotated RNA was masked out upstream and
candidates only exist where there was continuous coverage, so the eligible
territory is a specific, non-random subset. A length- and context-matched shuffle
within the unmasked, covered genome is the right null if you need a p-value;
the cross-class contrast is the honest quick read.

### Notes

- **Strand-agnostic by default.** Ensembl regulatory features are mostly
  strandless (`.`), and where a strand exists (`CTCF_binding_site`) it denotes
  motif orientation, not transcription — so requiring a strand match would drop
  real hits. Use `--annotate-stranded` for sources where strand *does* mean
  transcription (e.g. CAGE TSS peaks).
- **Chromosome naming is normalised** (`chr1`↔`1`, `chrM`↔`MT`), so an
  Ensembl regulation file works against a UCSC-named cohort. If **no** name could
  ever match, the run **fails with an error** rather than writing a table of `NA`
  — a silent zero-overlap result looks exactly like a genuine "nothing is
  regulatory" answer, which is the worst possible outcome.
- **The overlay is not restricted to regulatory data.** Any GFF3/GTF works —
  repeat annotation, CAGE peaks, or a lncRNA catalogue that is *not* in the GTF
  you masked with (a candidate reproducing a known lncRNA is not novel, it is
  annotated in a database you did not use).
- Pure standard library (`bisect` + `gzip`) — no pysam, no numpy, so it runs
  natively anywhere the merge step runs. The full 365k-feature Ensembl file loads
  in a few seconds.

---

## Extracting novel-transcript sequences (FASTA)

`extract_novel_fasta.py` pulls the sequences of rescued novel transcripts out of
a genome FASTA, ready for BLAST / ORF / homology / annotation work:

```bash
# from a per-sample run:
python extract_novel_fasta.py --genome genome.fa \
  --gtf sample_analysis.unknown_transcripts.gtf --out sample_analysis.novel.fa

# from the cross-sample consensus:
python extract_novel_fasta.py --genome genome.fa \
  --gtf cohort.consensus_transcripts.gtf --out cohort.novel.fa
```

- Exons of each transcript are concatenated in coordinate order (single-exon
  novel models are just their span), and the sequence is **reverse-complemented
  for `-` strand** transcripts, so the FASTA is in 5'→3' transcript orientation.
- FASTA headers carry provenance: `>transcript_id gene_id=… strand=… loc=chr:start-end(strand) class=… length=…`.
- If you point it at a merged GTF (reference + rescued), use `--source
  gdna_rescue` / `--source gdna_rescue_consensus` to extract only the rescued
  features.
- A genome `.fai` index is built automatically if absent. Uses pysam when
  available, and a dependency-free indexed reader otherwise (so it runs on
  native Windows too).

---

## Recommended end-to-end pipeline

```
FastQC / fastp (trim)
  -> STAR  (2-pass; KEEP multimappers; SortedByCoordinate + samtools index)
  -> extRNA  detect_gdna_vs_novel.py   (per sample)
  -> extRNA  merge_candidates.py --reference-gtf reference.gtf   (cohort consensus)
       => cohort.reference_plus_consensus.gtf
       (optional: --annotate regulatory=Homo_sapiens.GRCh38.regulatory_features.v116.gff3.gz
        => enhancer / promoter / CTCF overlap per consensus locus, in the table)
       (optional: extract_novel_fasta.py --gtf cohort.consensus_transcripts.gtf
        => cohort.novel.fa for BLAST/ORF/annotation of the rescued transcripts)
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
  With `--annotate`, those IDs also carry the overlapping feature type
  (`consensus_transcript_7-enhancer`), so the count matrix is self-describing —
  but decide on annotation *before* quantifying, since the IDs are the row labels.
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
- **Unstranded libraries are out of scope.** The method compares transcription
  strands, which is only informative when the library preserves strand. In an
  unstranded library even genuine single-strand RNA maps to both strands ~50/50,
  so gDNA cannot be told from RNA by strand symmetry. The tool runs and warns on
  unstranded input, but its classifications are not meaningful there — this tool
  is designed for stranded (forward / reverse) RNA-seq only.
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
- Per-epigenome regulatory activity in the overlay: Ensembl ships
  `regulatory_activity.v116.tsv.gz` (ACTIVE/INACTIVE per cell type, keyed by the
  same `ENSR…` IDs the overlay already reports). "Enhancer active in a tissue
  matching the sample" is a much stronger eRNA story than bare enhancer overlap.
- Length- and context-matched interval shuffling within the unmasked, covered
  genome, to turn the annotation cross-tab into a real enrichment statistic.

---

## Package layout

```
detect_gdna_vs_novel.py     # thin CLI entry point (per-sample analysis)
merge_candidates.py         # thin CLI entry point (cross-sample consensus)
extract_novel_fasta.py      # thin CLI entry point (novel-transcript FASTA)
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
  overlay.py       # external annotation overlay for consensus regions (stdlib only)
  fasta.py         # genome-FASTA reader + novel-transcript sequence extraction
  cli.py           # argument parsing
tests/
  generate_test_data.py   # synthetic archetypes + synthetic BAM/GTF
  test_classify.py        # classifier unit tests (incl. multimapper artifact)
  test_discovery.py       # discovery/merging unit tests
  test_strandedness.py    # read -> strand mapping unit tests
  test_crosssample.py     # cross-sample consensus unit tests (polars)
  test_overlay.py         # external annotation overlay unit tests (stdlib)
  test_qc.py              # gDNA QC + read-assignment + MultiQC writer tests
  test_fasta.py           # FASTA reader + novel-transcript extraction tests
  test_pipeline.py        # end-to-end integration (needs pysam)
```
