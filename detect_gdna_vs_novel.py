#!/usr/bin/env python3
"""Thin entry point so the tool can be run as `python detect_gdna_vs_novel.py`.

All logic lives in the ``gdna_rescue`` package.
"""

import sys

from gdna_rescue.cli import main

if __name__ == "__main__":
    sys.exit(main())
