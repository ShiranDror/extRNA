"""Unit tests for region discovery / merging (pysam-free)."""

import numpy as np
import pytest

from gdna_rescue.config import Config
from gdna_rescue.discovery import discover_intervals
from gdna_rescue.utils import find_true_runs, merge_runs_with_gap


@pytest.fixture
def cfg():
    return Config(min_depth=3, max_gap=50, min_region_length=100,
                  min_covered_fraction=0.5, min_covered_bases=0)


def test_find_true_runs():
    mask = np.array([0, 1, 1, 0, 0, 1, 0], dtype=bool)
    assert find_true_runs(mask) == [(1, 3), (5, 6)]


def test_merge_respects_gap():
    runs = [(0, 100), (120, 200), (400, 500)]
    merged = merge_runs_with_gap(runs, max_gap=50)
    assert merged == [(0, 200), (400, 500)]


def test_merge_does_not_bridge_annotation():
    runs = [(0, 100), (120, 200)]
    blocked = np.zeros(200, dtype=bool)
    blocked[105:115] = True  # an annotated feature sits inside the gap
    merged = merge_runs_with_gap(runs, max_gap=50, blocked_mask=blocked)
    assert merged == [(0, 100), (120, 200)]


def test_discover_basic_region(cfg):
    n = 2000
    plus = np.zeros(n, dtype=np.int32)
    minus = np.zeros(n, dtype=np.int32)
    plus[500:1000] = 10          # a clean 500 bp covered region
    annotated = np.zeros(n, dtype=bool)
    intervals = discover_intervals(plus, minus, annotated, cfg)
    assert (500, 1000) in intervals


def test_annotated_region_is_excluded(cfg):
    n = 2000
    plus = np.zeros(n, dtype=np.int32)
    minus = np.zeros(n, dtype=np.int32)
    plus[500:1000] = 10
    annotated = np.zeros(n, dtype=bool)
    annotated[500:1000] = True   # fully annotated -> nothing discovered
    intervals = discover_intervals(plus, minus, annotated, cfg)
    assert intervals == []


def test_short_regions_filtered(cfg):
    n = 2000
    plus = np.zeros(n, dtype=np.int32)
    minus = np.zeros(n, dtype=np.int32)
    plus[500:550] = 10           # 50 bp < min_region_length (100)
    intervals = discover_intervals(plus, minus, np.zeros(n, dtype=bool), cfg)
    assert intervals == []


def test_gap_merging_in_discovery(cfg):
    n = 2000
    plus = np.zeros(n, dtype=np.int32)
    minus = np.zeros(n, dtype=np.int32)
    plus[500:700] = 10
    plus[720:1000] = 10          # 20 bp gap < max_gap (50) -> merge
    intervals = discover_intervals(plus, minus, np.zeros(n, dtype=bool), cfg)
    assert (500, 1000) in intervals


def test_stranded_masking_discovers_antisense_over_exon(tmp_path):
    from gdna_rescue.gtf_io import parse_gtf
    from gdna_rescue.discovery import build_candidates_for_chrom, CTX_ANTISENSE
    from gdna_rescue.classify import LIKELY_NOVEL

    gtf = tmp_path / "a.gtf"
    gtf.write_text(
        'chr1\tsim\tgene\t1001\t2000\t.\t+\t.\tgene_id "g"; gene_name "g";\n'
        'chr1\tsim\texon\t1001\t2000\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'
    )
    ann = parse_gtf(str(gtf), "exon")
    cfg = Config()  # min_region_length 200, min_depth 10, min_covered_fraction 0.7

    n = 3000
    plus = np.zeros(n, dtype=np.int32)
    minus = np.zeros(n, dtype=np.int32)
    multi = np.zeros(n, dtype=np.int32)
    plus[1000:2000] = 25    # host mRNA on + (annotated exon)
    minus[1000:2000] = 15   # antisense transcription on - (novel)

    # Positional masking: exon masks BOTH strands -> nothing discovered there.
    pos = build_candidates_for_chrom(
        "chr1", plus.copy(), minus.copy(), multi.copy(), ann, cfg,
        stranded_masking=False,
    )
    assert pos == [] or all(not (c.start < 2000 and c.end > 1000) for c in pos)

    # Stranded masking: only + is masked over the + exon, so the - antisense
    # region IS discovered, labelled antisense, and rescued as novel.
    st = build_candidates_for_chrom(
        "chr1", plus.copy(), minus.copy(), multi.copy(), ann, cfg,
        stranded_masking=True,
    )
    over = [c for c in st if c.start < 2000 and c.end > 1000]
    assert len(over) == 1
    c = over[0]
    assert c.metrics.dominant_strand == "-"
    assert c.context_label == CTX_ANTISENSE
    assert c.label == LIKELY_NOVEL


def test_low_covered_fraction_filtered():
    cfg = Config(min_depth=3, max_gap=0, min_region_length=100,
                 min_covered_fraction=0.9, min_covered_bases=0)
    n = 2000
    plus = np.zeros(n, dtype=np.int32)
    # Sparse punctate coverage: only 30% of a would-be region is covered.
    plus[500:1000:3] = 10
    intervals = discover_intervals(plus, np.zeros(n, dtype=np.int32),
                                   np.zeros(n, dtype=bool), cfg)
    # With max_gap=0 punctate spikes never merge into a long region.
    assert all((e - s) < 100 or (plus[s:e] >= 3).mean() >= 0.9
               for s, e in intervals)
