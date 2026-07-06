"""Extract transcript sequences from a genome FASTA given a GTF.

Used to pull the sequences of rescued novel transcripts (from
``*.unknown_transcripts.gtf`` or ``*.consensus_transcripts.gtf``) into a FASTA
for downstream BLAST / homology / ORF / annotation work.

Sequence handling:
  * exon features of a transcript are concatenated in coordinate order (a
    single-exon novel transcript is just its span);
  * the sequence is reverse-complemented for transcripts on the '-' strand, so
    the FASTA is in the transcript's 5'->3' orientation.

Genome access uses pysam.FastaFile when available (Linux/macOS/WSL), and falls
back to a dependency-free ``.fai``-indexed reader otherwise (so this module and
its tests run on native Windows too). Either backend needs, or builds, a ``.fai``
index next to the FASTA.
"""

from __future__ import annotations

import os
import re
from typing import Dict, List, Optional, Tuple

_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')
_COMPLEMENT = str.maketrans("ACGTNacgtnRYSWKMryswkm", "TGCANtgcanYRSWMKyrswmk")


def reverse_complement(seq: str) -> str:
    """Reverse-complement a nucleotide string (IUPAC-aware, case-preserving)."""
    return seq.translate(_COMPLEMENT)[::-1]


# --------------------------------------------------------------------------- #
# Genome access
# --------------------------------------------------------------------------- #

class IndexedFasta:
    """Dependency-free random-access reader for a (possibly multi-record) FASTA.

    Uses a samtools-style ``.fai`` index (built if absent). Assumes each record
    has a uniform line width, which standard genome FASTAs satisfy.
    """

    def __init__(self, path: str):
        self.path = path
        self.fai_path = path + ".fai"
        if not os.path.exists(self.fai_path):
            self._build_fai()
        self.index: Dict[str, Tuple[int, int, int, int]] = self._load_fai()
        self._fh = open(path, "rb")

    def _build_fai(self) -> None:
        with open(self.path, "rb") as fh, open(self.fai_path, "w") as out:
            offset = 0
            cur = None  # [name, length, seq_offset, linebases, linewidth]
            for raw in fh:
                if raw.startswith(b">"):
                    if cur is not None:
                        out.write("\t".join(map(str, cur)) + "\n")
                    name = raw[1:].split()[0].decode()
                    cur = [name, 0, offset + len(raw), 0, 0]
                elif cur is not None:
                    seq = raw.rstrip(b"\r\n")
                    if cur[3] == 0:
                        cur[3] = len(seq)       # linebases
                        cur[4] = len(raw)       # linewidth (incl newline bytes)
                    cur[1] += len(seq)          # length
                offset += len(raw)
            if cur is not None:
                out.write("\t".join(map(str, cur)) + "\n")

    def _load_fai(self) -> Dict[str, Tuple[int, int, int, int]]:
        index = {}
        with open(self.fai_path) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 5:
                    index[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))
        return index

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Return the sequence for 0-based half-open [start, end)."""
        if chrom not in self.index:
            return ""
        length, offset, linebases, linewidth = self.index[chrom]
        if linebases <= 0:
            return ""
        start = max(0, start)
        end = min(length, end)
        if end <= start:
            return ""

        def bytepos(pos):
            return offset + (pos // linebases) * linewidth + (pos % linebases)

        self._fh.seek(bytepos(start))
        raw = self._fh.read(bytepos(end) - bytepos(start))
        return raw.replace(b"\n", b"").replace(b"\r", b"").decode().upper()

    def close(self) -> None:
        self._fh.close()


class _PysamFasta:
    def __init__(self, ff):
        self._ff = ff

    def fetch(self, chrom: str, start: int, end: int) -> str:
        try:
            return self._ff.fetch(reference=chrom, start=start, end=end).upper()
        except (KeyError, ValueError):
            return ""

    def close(self) -> None:
        self._ff.close()


def open_fasta(path: str):
    """Open a genome FASTA with pysam if available, else the pure reader."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"genome FASTA not found: {path!r}")
    try:
        import pysam
    except ImportError:
        return IndexedFasta(path)
    if not os.path.exists(path + ".fai"):
        pysam.faidx(path)   # build the index if it is missing
    return _PysamFasta(pysam.FastaFile(path))


# --------------------------------------------------------------------------- #
# GTF -> transcripts
# --------------------------------------------------------------------------- #

def read_transcripts_from_gtf(
    path: str, source_filter: Optional[str] = None
) -> List[Tuple[str, dict]]:
    """Return [(transcript_id, record)] in first-seen order.

    record = {chrom, strand, source, exons:[(start0,end)], span, attrs}.
    ``source_filter`` (GTF column 2) restricts to transcripts from that source,
    e.g. 'gdna_rescue' / 'gdna_rescue_consensus' when reading a merged GTF.
    """
    tx: Dict[str, dict] = {}
    order: List[str] = []
    with open(path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            chrom, source, feature, s, e, _score, strand, _frame, attrs = p[:9]
            if feature not in ("exon", "transcript"):
                continue
            if source_filter is not None and source != source_filter:
                continue
            m = dict(_ATTR_RE.findall(attrs))
            tid = m.get("transcript_id")
            if not tid:
                continue
            try:
                start = int(s) - 1
                end = int(e)
            except ValueError:
                continue
            rec = tx.get(tid)
            if rec is None:
                rec = {"chrom": chrom, "strand": strand, "source": source,
                       "exons": [], "span": None, "attrs": m}
                tx[tid] = rec
                order.append(tid)
            if feature == "exon":
                rec["exons"].append((start, end))
            elif feature == "transcript":
                rec["span"] = (start, end)
    return [(tid, tx[tid]) for tid in order]


def extract_transcript_fasta(
    gtf_path: str,
    genome_path: str,
    out_path: str,
    source_filter: Optional[str] = None,
    line_width: int = 60,
) -> int:
    """Write a FASTA of transcript sequences from ``gtf_path``. Returns count."""
    fa = open_fasta(genome_path)
    transcripts = read_transcripts_from_gtf(gtf_path, source_filter=source_filter)
    written = 0
    try:
        with open(out_path, "w") as out:
            for tid, rec in transcripts:
                exons = sorted(rec["exons"])
                if not exons and rec["span"] is not None:
                    exons = [rec["span"]]
                if not exons:
                    continue
                seq = "".join(fa.fetch(rec["chrom"], s, e) for s, e in exons)
                if not seq:
                    continue
                if rec["strand"] == "-":
                    seq = reverse_complement(seq)
                attrs = rec["attrs"]
                cls = attrs.get("classification") or attrs.get("consensus_class", "")
                loc = f"{rec['chrom']}:{exons[0][0] + 1}-{exons[-1][1]}"
                header = (
                    f">{tid} gene_id={attrs.get('gene_id', '')} "
                    f"strand={rec['strand']} loc={loc}({rec['strand']}) "
                    f"class={cls} length={len(seq)}"
                )
                out.write(header + "\n")
                for i in range(0, len(seq), line_width):
                    out.write(seq[i:i + line_width] + "\n")
                written += 1
    finally:
        fa.close()
    return written
