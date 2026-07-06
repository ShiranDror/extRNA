"""gdna_rescue: detect likely gDNA contamination and rescue novel transcripts.

The package is layered so that the biological logic can be imported and tested
without pysam (which has no Windows wheels):

    config       - all tunable thresholds (Config dataclass)
    utils        - dependency-light helpers (logging, interval maths)
    gtf_io       - GTF parsing / interval indexing / masking
    classify     - per-region metrics + transparent rule-based classifier (numpy)
    discovery    - unannotated region discovery + context labelling (numpy)
    bam_io       - the ONLY pysam-dependent module (coverage extraction)
    strandedness - automatic library-strandedness inference
    writers      - TSV / GTF / JSON / BED / bedGraph output
    pipeline     - end-to-end orchestration
    cli          - argument parsing / entry point
"""

__version__ = "1.0.0"
