#!/usr/bin/env python
"""Build the self-contained docs/live site (wheel + static previews).

Idempotent, CI-friendly: it (re)builds the treeio wheel, regenerates the gallery
preview PNGs, copies both into docs/live (which is what gets published), and
verifies every example in examples.js has a matching preview.

Run from the repo root:
    python scripts/build_live.py
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LIVE = ROOT / "docs" / "live"
PREVIEWS = LIVE / "static" / "previews"
PACKAGES = LIVE / "packages"


def sh(cmd: list[str]) -> None:
    print("$ " + " ".join(cmd))
    subprocess.run(cmd, cwd=ROOT, check=True)


def build_wheel() -> Path:
    """Build a wheel into dist/ and move it into docs/live/packages."""
    # Prefer `uv build` (fast, hermetic); fall back to `python -m build`.
    try:
        sh(["uv", "build", "--wheel"])
    except FileNotFoundError:
        sh([sys.executable, "-m", "build", "--wheel"])
    dist = ROOT / "dist"
    whl = sorted(dist.glob("treeio-*-py3-none-any.whl"))[-1]
    PACKAGES.mkdir(parents=True, exist_ok=True)
    # clean stale wheels so the published folder only has the current version
    for old in PACKAGES.glob("treeio-*.whl"):
        old.unlink()
    dst = PACKAGES / whl.name
    shutil.copy2(whl, dst)
    print(f"  -> {dst.relative_to(ROOT)}")
    return dst


def generate_previews() -> None:
    """Regenerate gallery preview PNGs into docs/live/static/previews."""
    import matplotlib

    matplotlib.use("Agg")
    import sys as _sys
    sys.path.insert(0, str(ROOT))  # let `examples.gallery` import treeio
    sys.path.insert(0, str(ROOT / "examples"))
    import gallery  # noqa: E402  (imports treeio from ROOT)

    PREVIEWS.mkdir(parents=True, exist_ok=True)
    # map the gallery's OUT to our previews dir by monkeypatching before calling
    gallery.OUT = PREVIEWS
    gallery.random.seed(20)
    gallery.make_layout_gallery()
    gallery.make_trait()
    gallery.make_highlight()
    gallery.make_rings()
    gallery.make_heatmap()
    gallery.make_radial()
    gallery.make_boxed()
    gallery.make_symbols()
    print(f"  -> previews in {PREVIEWS.relative_to(ROOT)}")


def verify_examples() -> None:
    """Every example in examples.js must reference an existing preview PNG."""
    text = (LIVE / "examples.js").read_text()
    previews = sorted(p.name for p in PREVIEWS.glob("*.png"))
    missing = []
    for m in re.finditer(r"preview:\s*'([^']+)'", text):
        ref = Path(m.group(1))
        if not (LIVE / ref).exists():
            missing.append(str(ref))
    if missing:
        raise SystemExit("examples.js references missing previews: " + ", ".join(missing))
    print("  -> all example previews present")


def main() -> int:
    print("== building live site ==")
    build_wheel()
    generate_previews()
    verify_examples()
    print("ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
