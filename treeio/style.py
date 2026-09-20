#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-03-29 18:32

"""Colour scales, palettes and figure sizing shared by all backends."""

from __future__ import annotations

from typing import Dict, Sequence

from .tree import Tree

__all__ = [
    "color_by_value",
    "color_map",
    "named_palette",
    "suggest_figsize",
]

# deterministic colour palettes (used when no mapping is supplied)
_PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
]
_VIRIDIS = [
    "#440154", "#46327e", "#365c8d", "#277f8e", "#1fa187",
    "#4ac16d", "#a0da39", "#fde725",
]


def color_by_value(values: Dict[str, object]) -> Dict[str, str]:
    """Build a deterministic ``{name: hex}`` palette from a value mapping.

    Numeric values map to a continuous viridis gradient; categorical values
    cycle through a discrete palette.
    """
    items = [(k, v) for k, v in values.items() if v is not None]
    if not items:
        return {}
    numeric = [v for _, v in items if isinstance(v, (int, float))]
    if len(numeric) == len(items):
        lo, hi = min(numeric), max(numeric)
        span = (hi - lo) or 1.0
        return {name: _gradient((v - lo) / span) for name, v in items}
    cats = {}
    out = {}
    for name, v in items:
        key = str(v)
        if key not in cats:
            cats[key] = _PALETTE[len(cats) % len(_PALETTE)]
        out[name] = cats[key]
    return out


def color_map(values: Dict[str, object]) -> Dict[str, str]:
    """Alias of :func:`color_by_value` (trait -> colour)."""
    return color_by_value(values)


def named_palette(names: Sequence[str]) -> Dict[str, str]:
    """Assign a distinct colour to each of ``names``."""
    return {name: _PALETTE[i % len(_PALETTE)] for i, name in enumerate(names)}


def suggest_figsize(tree: Tree, layout: str = "rectangular",
                    base_width: int = 800, tip_slot: int = 14):
    """A sensible canvas size ``(width, height)`` for a tree.

    Height scales with the number of tips so points / labels do not overlap;
    polar layouts get a square canvas sized ~ sqrt(tips).
    """
    n = tree.nleaves or 1
    if layout in ("circular", "fan", "radial", "unrooted"):
        import math as _m
        side = int(min(max(400, _m.sqrt(n) * base_width / 5.0), 6000))
        return side, side
    height = int(min(max(300, n * tip_slot + 80), 20000))
    return int(base_width), height


def _gradient(t: float) -> str:
    """A viridis-like continuous palette."""
    t = min(max(t, 0.0), 1.0)
    pos = t * (len(_VIRIDIS) - 1)
    i = int(pos)
    frac = pos - i
    a, b = _VIRIDIS[i], _VIRIDIS[min(i + 1, len(_VIRIDIS) - 1)]
    ar, ag, ab = int(a[1:3], 16), int(a[3:5], 16), int(a[5:7], 16)
    br, bg, bb = int(b[1:3], 16), int(b[3:5], 16), int(b[5:7], 16)
    r = int(ar + (br - ar) * frac)
    g = int(ag + (bg - ag) * frac)
    bl = int(ab + (bb - ab) * frac)
    return f"#{r:02x}{g:02x}{bl:02x}"


def _fmt_num(x) -> str:
    f = float(x)
    if abs(f) >= 1000 or (abs(f) < 0.01 and f != 0):
        return f"{f:.1e}"
    if f == int(f):
        return str(int(f))
    return f"{f:.2f}"
