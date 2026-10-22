"""Machine-independent paths for the snrna/integration Python scripts.

Everything is resolved from this file's own location, so the scripts run
unchanged on any checkout of the repo (local workstation or cluster) instead of
assuming /private/groups/colquittlab/song-system-grn.

Scripts in a subdirectory put this directory on sys.path and import from here:

    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # -> snrna/integration
    from integration_paths import GG_HYBRID

Importing this module also puts finch-integration-toolkit (reciprocal_score,
zaremba_prep, composite_score) on sys.path. It is looked for, in order, at
$FINCH_INTEGRATION_TOOLKIT, a checkout next to this repo, and the cluster path.
"""
import os
import sys
from pathlib import Path

INTEGRATION_DIR = Path(__file__).resolve().parent
REPO_ROOT = INTEGRATION_DIR.parents[1]
DATASETS = INTEGRATION_DIR / "datasets"
GG_HYBRID = INTEGRATION_DIR / "gg_glutonly_hybrid"

_TOOLKIT_CANDIDATES = [
    os.environ.get("FINCH_INTEGRATION_TOOLKIT"),
    REPO_ROOT.parent / "finch-integration-toolkit",
    "/private/groups/colquittlab/finch-integration-toolkit",
]
TOOLKIT_DIR = next((Path(p) for p in _TOOLKIT_CANDIDATES if p and Path(p).is_dir()), None)
if TOOLKIT_DIR is not None and str(TOOLKIT_DIR) not in sys.path:
    sys.path.insert(1, str(TOOLKIT_DIR))
