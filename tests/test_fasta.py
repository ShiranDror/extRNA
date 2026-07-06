"""Unit tests for novel-transcript FASTA extraction (pysam-free path)."""

from gdna_rescue.fasta import (
    IndexedFasta,
    reverse_complement,
    read_transcripts_from_gtf,
    extract_transcript_fasta,
)


def _write_genome(tmp_path):
    # 3 chromosomes, wrapped at width 10 to exercise line-wrapped indexing.
    seq1 = "ACGTACGTAC" "GGGGCCCCTT" "AAAATTTTGG"   # chr1, 30 bp
    seq2 = "TTTTTTTTTT" "GGGGGGGGGG"                 # chr2, 20 bp
    fa = tmp_path / "genome.fa"
    fa.write_text(
        ">chr1\n" + "\n".join(seq1[i:i+10] for i in range(0, len(seq1), 10)) + "\n"
        ">chr2\n" + "\n".join(seq2[i:i+10] for i in range(0, len(seq2), 10)) + "\n"
    )
    return str(fa), seq1, seq2


def test_reverse_complement():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AAAC") == "GTTT"
    assert reverse_complement("acgtN") == "Nacgt"


def test_indexed_fasta_builds_index_and_fetches(tmp_path):
    fa_path, seq1, seq2 = _write_genome(tmp_path)
    fa = IndexedFasta(fa_path)
    # Whole chrom and sub-ranges (crossing line-wrap boundaries).
    assert fa.fetch("chr1", 0, 30) == seq1
    assert fa.fetch("chr1", 5, 15) == seq1[5:15]
    assert fa.fetch("chr2", 0, 20) == seq2
    assert fa.fetch("chr1", 25, 100) == seq1[25:]   # clamped to length
    assert fa.fetch("missing", 0, 10) == ""
    fa.close()


def test_extract_transcript_fasta_plus_and_minus(tmp_path):
    fa_path, seq1, seq2 = _write_genome(tmp_path)
    gtf = tmp_path / "novel.gtf"
    # + transcript on chr1 [1..10], - transcript on chr1 [11..20].
    gtf.write_text(
        'chr1\tgdna_rescue\ttranscript\t1\t10\t.\t+\t.\tgene_id "u1_gene"; transcript_id "u1"; classification "likely_novel_transcript";\n'
        'chr1\tgdna_rescue\texon\t1\t10\t.\t+\t.\tgene_id "u1_gene"; transcript_id "u1";\n'
        'chr1\tgdna_rescue\ttranscript\t11\t20\t.\t-\t.\tgene_id "u2_gene"; transcript_id "u2"; classification "likely_novel_transcript";\n'
        'chr1\tgdna_rescue\texon\t11\t20\t.\t-\t.\tgene_id "u2_gene"; transcript_id "u2";\n'
    )
    out = tmp_path / "novel.fa"
    n = extract_transcript_fasta(str(gtf), fa_path, str(out))
    assert n == 2

    records = {}
    name = None
    for line in open(out):
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = ""
        else:
            records[name] += line.strip()
    # + strand: genomic sequence as-is.
    assert records["u1"] == seq1[0:10]
    # - strand: reverse complement of the genomic span.
    assert records["u2"] == reverse_complement(seq1[10:20])


def test_source_filter(tmp_path):
    fa_path, seq1, seq2 = _write_genome(tmp_path)
    gtf = tmp_path / "merged.gtf"
    gtf.write_text(
        'chr1\tref\ttranscript\t1\t10\t.\t+\t.\tgene_id "g"; transcript_id "ref1";\n'
        'chr1\tref\texon\t1\t10\t.\t+\t.\tgene_id "g"; transcript_id "ref1";\n'
        'chr1\tgdna_rescue\ttranscript\t11\t20\t.\t+\t.\tgene_id "u_gene"; transcript_id "u1";\n'
        'chr1\tgdna_rescue\texon\t11\t20\t.\t+\t.\tgene_id "u_gene"; transcript_id "u1";\n'
    )
    txs = read_transcripts_from_gtf(str(gtf), source_filter="gdna_rescue")
    assert [tid for tid, _ in txs] == ["u1"]

    out = tmp_path / "novel.fa"
    n = extract_transcript_fasta(str(gtf), fa_path, str(out), source_filter="gdna_rescue")
    assert n == 1
    assert ">u1 " in open(out).read()
    assert "ref1" not in open(out).read()
