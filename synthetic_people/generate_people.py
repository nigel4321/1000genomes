#!/usr/bin/env python3
"""Generate synthetic single-person VCFs (CLI shim).

Implementation lives in the ``syntheticgen`` package next to this file.
See README.md for usage and IMPLEMENTATION_PLAN.md for the roadmap.
"""

from __future__ import annotations

import sys
from pathlib import Path

# Make the sibling package importable when invoked as a plain script.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from syntheticgen.cli import main  # noqa: E402


if __name__ == "__main__":
    sys.exit(main())
