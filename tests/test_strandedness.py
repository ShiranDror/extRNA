"""Unit tests for the read -> transcription-strand mapping (pysam-free).

``transcription_strand`` lives in bam_io but only touches read attributes, so we
exercise it with lightweight mock reads and never import pysam.
"""

from dataclasses import dataclass

from gdna_rescue.bam_io import transcription_strand


@dataclass
class MockRead:
    is_reverse: bool = False
    is_read1: bool = True
    is_paired: bool = False


def test_single_end_forward():
    # Forward library, SE read on + strand -> transcript +.
    assert transcription_strand(MockRead(is_reverse=False), "forward") == "+"
    assert transcription_strand(MockRead(is_reverse=True), "forward") == "-"


def test_single_end_reverse():
    # Reverse (dUTP) library: SE read strand is flipped.
    assert transcription_strand(MockRead(is_reverse=False), "reverse") == "-"
    assert transcription_strand(MockRead(is_reverse=True), "reverse") == "+"


def test_unstranded_uses_alignment_strand():
    assert transcription_strand(MockRead(is_reverse=False), "unstranded") == "+"
    assert transcription_strand(MockRead(is_reverse=True), "unstranded") == "-"


def test_paired_forward_read1_vs_read2():
    # forward: read1 same strand as transcript, read2 opposite.
    r1 = MockRead(is_reverse=False, is_read1=True, is_paired=True)
    r2 = MockRead(is_reverse=False, is_read1=False, is_paired=True)
    assert transcription_strand(r1, "forward") == "+"
    assert transcription_strand(r2, "forward") == "-"


def test_paired_reverse_read1_vs_read2():
    # reverse/dUTP: read1 opposite strand of transcript, read2 same.
    r1 = MockRead(is_reverse=False, is_read1=True, is_paired=True)
    r2 = MockRead(is_reverse=False, is_read1=False, is_paired=True)
    assert transcription_strand(r1, "reverse") == "-"
    assert transcription_strand(r2, "reverse") == "+"


def test_paired_reverse_read1_on_minus():
    r1 = MockRead(is_reverse=True, is_read1=True, is_paired=True)
    assert transcription_strand(r1, "reverse") == "+"
