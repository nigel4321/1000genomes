"""Process-shared on-disk cache for tests that invoke ``cli.main`` end-to-end.

Why this exists
---------------

``cli.main`` unconditionally calls ``fetch_clinvar(args.cache_dir, args.build)``
near the top of every run (see ``syntheticgen/cli.py``). On a cache miss
``fetch_clinvar`` downloads the ~70 MB ClinVar VCF + tabix index from
NCBI. Tests historically passed each one a per-test
``tempfile.TemporaryDirectory()`` as ``--cache-dir``, so every test saw a
fresh empty cache and re-downloaded ClinVar.

A single CI run showed 15+ ``clinvar.vcf.gz`` downloads (+ 15 ``.tbi``)
across the unit-test job — ~1 GB of redundant transfer per workflow,
~7 s wall-time per fetch. The stderr message "Fetching ClinVar
(cached across runs)..." is misleading: the cache is only "cached
across runs" if the cache directory is reused.

How to use
----------

Tests that invoke ``cli.main`` end-to-end (and don't specifically need
to assert cache-miss behaviour) should pass
:data:`SHARED_TEST_CACHE_DIR` as ``--cache-dir`` instead of a per-test
tmpdir. The first test to fetch ClinVar populates it; every subsequent
test in the same process gets a cache hit and skips the download.

The path lives under ``tempfile.gettempdir()`` so it persists across
tests within a pytest process and gets cleared by the OS / CI runner
between sessions. It is intentionally NOT under the repo or under the
user's ``~/.cache`` (which is where the cli's runtime auto-fetch
writes), so a test run cannot collide with a real cli invocation.

When NOT to use this
--------------------

Tests that specifically exercise the cache-population code path —
``FetchReferenceFastaTest`` in ``test_reference.py``, the
``fetch_clinvar`` direct tests — must keep their per-test cache_dir
plus their existing ``urllib.request.urlopen`` mocks, because they're
asserting download-on-empty-cache behaviour.

Layout (matches ``clinvar.fetch_clinvar`` +
``reference.fetch_reference_fasta``)::

    <SHARED>/clinvar_<BUILD>.vcf.gz[.tbi]      — ClinVar VCF
    <SHARED>/reference/<BUILD>.fa[.fai]        — REF FASTA (only when
                                                  a test opts in to
                                                  fetching it)
"""

from __future__ import annotations

import tempfile
from pathlib import Path

SHARED_TEST_CACHE_DIR: Path = (
    Path(tempfile.gettempdir()) / "synthetic_people_test_cache"
)
