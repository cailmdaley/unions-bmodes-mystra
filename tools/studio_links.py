"""Point Studio's results directories at lc's.

lc keeps an output at results/<universe>/<a>/<b>/<file>; Studio reads it at
analyses/<a>/analyses/<b>/results/<universe>/<file>. For each universe this links
the second to the first, one symlink per analysis. Root outputs need none.

    python tools/studio_links.py [<universe> ...]     # default: every universes/*.yaml
"""

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def link(universe):
    for spec in ROOT.glob("analyses/**/astra.yaml"):
        owner = spec.parent
        inner = [part for part in owner.relative_to(ROOT).parts if part != "analyses"]
        target = ROOT / "results" / universe / Path(*inner)
        place = owner / "results" / universe
        place.parent.mkdir(exist_ok=True)
        relative = os.path.relpath(target, place.parent)
        if place.is_symlink():
            place.unlink()
        place.symlink_to(relative)


if __name__ == "__main__":
    for universe in sys.argv[1:] or [p.stem for p in ROOT.glob("universes/*.yaml")]:
        link(universe)
