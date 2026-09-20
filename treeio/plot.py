#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-03-29 18:32

"""
Phylogenetic tree plotting (matplotlib backend).

Trees are drawn onto a real :class:`matplotlib.axes.Axes` using ``ax.plot`` /
``ax.scatter`` / ``ax.text`` / ``LineCollection``, so every matplotlib feature
(subplots, legends, colormaps, limits, ``savefig``, animations) works.  The
renderer is object-oriented and matplotlib-compatible.

Layouts
-------
* ``rectangular`` / ``roundrect`` -- L-shaped / rounded phylogram
* ``circular`` / ``fan``          -- radial / partial-sector tree
* ``radial``                      -- equal-angle ring dendrogram
* ``slanted``                     -- diagonal phylogram
* ``unrooted``                    -- equal-angle circular layout
* ``time``                        -- time-scaled tree (tips left, root right)

Examples
--------
>>> import matplotlib.pyplot as plt
>>> from treeio import render, draw_tree, treeplot
>>> tree = read("data/animals.nwk")
>>> ax = render(tree, backend="mpl", layout="circular", tip_labels=True)
>>> plt.savefig("tree.png")
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

from .tree import Tree
from .style import (_PALETTE, _VIRIDIS, _gradient, _fmt_num, color_by_value, color_map,
                    named_palette, suggest_figsize)

__all__ = [
    "tree_coords",
    "edge_segments",
    "draw_tree",
    "treeplot",
    "TreePlotter",
    "render",
    "draw",
    "plot",
    "save",
    "tree_theme",
    "add_scalebar",
    "facet_grid",
    "grid_of_trees",
    "register_backend",
    "register_layout",
    "unregister_backend",
    "unregister_layout",
    "layouts",
]

# user-registered render backends / layouts
_BACKEND_REGISTRY: dict = {}
_LAYOUT_REGISTRY: dict = {}


_BUILTIN_BACKENDS = {"mpl", "matplotlib", "ascii", "text", "cli"}
_BUILTIN_LAYOUTS = {"rectangular", "roundrect", "circular", "fan", "radial", "slanted", "unrooted", "time"}


def register_backend(name: str, func):
    """Register a custom render backend, usable via ``render(..., backend=name)``.

    ``func(tree, **kwargs)`` should return whatever the backend produces (e.g.
    a matplotlib axes, a string, a figure).  Built-in backend names cannot be
    overwritten; use :func:`unregister_backend` first to replace one.
    """
    if name in _BUILTIN_BACKENDS:
        raise ValueError(f"cannot overwrite built-in backend {name!r}")
    _BACKEND_REGISTRY[name] = func
    return func


def unregister_backend(name: str):
    """Remove a previously registered custom backend."""
    return _BACKEND_REGISTRY.pop(name, None)


def register_layout(name: str, func):
    """Register a custom layout function.

    ``func(tree)`` returns ``{node: (x, y)}`` (final Cartesian coordinates); it
    becomes usable as ``layout=name`` in :func:`draw_tree` / :func:`tree_coords`.
    Built-in layout names cannot be overwritten.
    """
    if name in _BUILTIN_LAYOUTS:
        raise ValueError(f"cannot overwrite built-in layout {name!r}")
    _LAYOUT_REGISTRY[name] = func
    return func


def unregister_layout(name: str):
    """Remove a previously registered custom layout."""
    return _LAYOUT_REGISTRY.pop(name, None)


def layouts():
    """Return the list of available layout names (built-in + registered)."""
    return sorted(set(_BUILTIN_LAYOUTS) | set(_LAYOUT_REGISTRY))


# ----------------------------------------------------------------------- --
# geometry (data coordinates, independent of canvas)
# ----------------------------------------------------------------------- --
def tree_coords(tree: Tree, layout: str = "rectangular", reverse_x: bool = False) -> Dict[Tree, Tuple[float, float]]:
    """Return ``{node: (x, y)}`` in data coordinates for ``layout``."""
    if layout in _LAYOUT_REGISTRY:
        return dict(_LAYOUT_REGISTRY[layout](tree))
    use_depth = layout == "radial"
    xpos = _x_positions(tree, reverse_x=reverse_x, use_depth=use_depth)
    tips = tree.get_tips()
    tipy = {tip: float(i) for i, tip in enumerate(tips)}
    coords: Dict[Tree, Tuple[float, float]] = {}
    for node in tree.traverse("postorder"):
        if not node.children:
            y = tipy[node]
        else:
            y = sum(coords[c][1] for c in node.children) / len(node.children)
        coords[node] = (xpos[node], y)
    if layout in ("circular", "fan", "radial", "unrooted"):
        return _polarize(coords, tree, layout)
    return coords


def _x_positions(tree: Tree, reverse_x: bool, use_depth: bool) -> Dict[Tree, float]:
    if use_depth:
        def depth(node):
            d = 0
            cur = node
            while cur is not None and cur.parent is not None:
                d += 1
                cur = cur.parent
            return float(d)
        out = {n: depth(n) for n in tree.traverse("preorder")}
    else:
        has_len = any(t.branch_length for t in tree.traverse("preorder") if t.branch_length)
        if has_len:
            out = {}
            stack = [(tree, 0.0)]
            while stack:
                node, x = stack.pop()
                out[node] = x
                for c in node._children:
                    d = c.branch_length if c.branch_length is not None else 0.0
                    stack.append((c, x + d))
        else:
            def depth(node):
                d = 0
                cur = node
                while cur is not None and cur._parent is not None:
                    d += 1
                    cur = cur._parent
                return float(d)
            out = {n: depth(n) for n in tree.traverse("preorder")}
    if reverse_x:
        mx = max(out.values(), default=0.0) or 1.0
        out = {n: mx - x for n, x in out.items()}
    return out


def _polarize(coords: Dict[Tree, Tuple[float, float]], tree: Tree, layout: str):
    n = tree.nleaves or 1
    xmax = max((c[0] for c in coords.values()), default=0.0) or 1.0
    out = {}
    if layout == "fan":
        sweep = math.pi
        start = math.pi / 2 - sweep / 2
        for node, (x, y) in coords.items():
            r = x / xmax
            theta = start + (y / n) * sweep
            out[node] = (r * math.cos(theta), r * math.sin(theta))
        return out
    for node, (x, y) in coords.items():
        r = x / xmax
        theta = ((y / n) * 2 * math.pi - math.pi / 2) if layout != "unrooted" else ((y / n) * 2 * math.pi)
        out[node] = (r * math.cos(theta), r * math.sin(theta))
    return out


def edge_segments(tree: Tree, coords=None, layout: str = "rectangular") -> List[List[Tuple[float, float]]]:
    """Return a list of polylines (each a list of points) for every edge."""
    if coords is None:
        coords = tree_coords(tree, layout)
def edge_segments(tree: Tree, coords=None, layout: str = "rectangular") -> List[List[Tuple[float, float]]]:
    """Return a list of polylines (each a list of points) for every edge."""
    if coords is None:
        coords = tree_coords(tree, layout)
    straight = layout in ("circular", "fan", "radial", "unrooted", "slanted") or layout in _LAYOUT_REGISTRY
    segs = []
    for node in tree.traverse("preorder"):
        if node.parent is None:
            continue
        xp, yp = coords[node.parent]
        xc, yc = coords[node]
        if straight:
            segs.append([(xp, yp), (xc, yc)])
        else:
            segs.append([(xp, yp), (xp, yc), (xc, yc)])
    return segs


def _mpl():
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as e:  # pragma: no cover
        raise ImportError("matplotlib is required for treeio.plot; install it with `pip install matplotlib`") from e
    return plt


def _figsize(tree: Tree, layout: str):
    return suggest_figsize(tree, layout)


def _new_axes(tree: Tree, layout: str, reverse_x: bool = False, ax=None, **kw):
    plt = _mpl()
    polar = layout in ("circular", "fan", "radial", "unrooted")
    if ax is None:
        w, h = _figsize(tree, layout)
        figsize = kw.pop("figsize", (w / 80.0, h / 80.0))
        ax = plt.subplots(figsize=figsize, subplot_kw={"aspect": "equal"} if polar else None)[1]
    if polar:
        ax.set_axis_off()
    return ax


# ----------------------------------------------------------------------- --
# draw onto a matplotlib Axes
# ----------------------------------------------------------------------- --
def draw_tree(
    tree: Tree,
    ax=None,
    layout: str = "rectangular",
    tip_labels: bool = True,
    tip_points: bool = True,
    node_points: bool = False,
    node_support: bool = False,
    tip_colors: Optional[Dict[str, str]] = None,
    tip_color_field: str = None,
    show_tip_legend: bool = True,
    edge_color: str = "#333333",
    edge_width: float = 1.2,
    label_size: int = 10,
    label_radial_offset: float = 0.04,
    axis_label_size: int = 9,
    tick_size: int = 8,
    reverse_x: bool = False,
    branch_color_field: str = None,
    cmap: str = "viridis",
    theme: str = "clean",
    show_colorbar: bool = False,
    colorbar_label: str = "value",
    **kwargs,
):
    """Draw ``tree`` onto ``ax`` (creating one if ``ax`` is None) and return it.

    Uses real ``ax`` primitives (``LineCollection`` + ``scatter`` + ``text``)
    so the returned axes is fully matplotlib-compatible.  Edges are batched
    into a single collection and all tips into one ``scatter``, so very large
    trees render quickly.

    ``branch_color_field`` colours each branch by a per-node value and adds a
    continuous colourbar if ``show_colorbar`` is true; ``theme`` applies a
    :func:`tree_theme` style.
    """
    plt = _mpl()
    from matplotlib.collections import LineCollection
    polar = layout in ("circular", "fan", "radial", "unrooted")
    if ax is None:
        figsize = kwargs.pop("figsize", None)
        w, h = _figsize(tree, layout)
        ax = plt.subplots(figsize=figsize or (w / 80.0, h / 80.0),
                          subplot_kw={"aspect": "equal"} if polar else None)[1]
    coords = tree_coords(tree, layout=layout, reverse_x=reverse_x)
    segments = edge_segments(tree, coords, layout=layout)

    mappable = None
    if branch_color_field is not None:
        import numpy as np
        import matplotlib
        vals = [c.get_data(branch_color_field) for c in tree.traverse("preorder") if c.parent is not None]
        arr = np.asarray([np.nan if v is None else float(v) for v in vals], dtype=float)
        if np.all(np.isnan(arr)):
            lc = LineCollection(segments, colors=edge_color, linewidths=edge_width,
                                capstyle="round", antialiaseds=True, zorder=1)
        else:
            lo = float(np.nanmin(arr)); hi = float(np.nanmax(arr))
            if not np.isfinite(lo) or lo == hi:
                lo, hi = 0.0, 1.0
            norm = matplotlib.colors.Normalize(vmin=lo, vmax=hi)
            lc = LineCollection(segments, array=arr, cmap=matplotlib.colormaps.get_cmap(cmap),
                                norm=norm, linewidths=edge_width, capstyle="round", zorder=1)
            mappable = lc
    else:
        lc = LineCollection(segments, colors=edge_color, linewidths=edge_width,
                            capstyle="round", antialiaseds=True, zorder=1)
    ax.add_collection(lc)

    tips = tree.get_tips()
    tip_x = [coords[t][0] for t in tips]
    tip_y = [coords[t][1] for t in tips]
    tip_c, tip_mappable, tip_legend = _tip_colors_for(tree, tips, tip_colors, tip_color_field, cmap, edge_color)
    if tip_points:
        ax.scatter(tip_x, tip_y, s=18, c=tip_c, zorder=3, edgecolors="white", linewidths=0.4)
    if tip_labels:
        if polar:
            import math
            n = len(tips) or 1
            # wrap-around-safe angular crowding drives the radial offset so
            # dense clusters push their labels outward (avoiding overlap).
            threshold = (2 * math.pi / n) * 0.6
            angles = [math.atan2(tip_y[i], tip_x[i]) for i in range(n)]
            offset_map = {}
            for i in range(n):
                a = angles[i]
                crowd = sum(1 for j in range(n)
                            if abs((angles[j] - a + math.pi) % (2 * math.pi) - math.pi) < threshold)
                offset_map[i] = label_radial_offset * (1 + 0.3 * max(0, crowd - 1))
        else:
            offset_map = {}
        for i, (t, x, y, c) in enumerate(zip(tips, tip_x, tip_y, tip_c)):
            off = offset_map.get(i, label_radial_offset)
            _write_tip_label(ax, x, y, t.name, polar, c, label_size, offset=off)

    nodes = tree.get_internal_nodes()
    if node_points:
        ax.scatter([coords[n][0] for n in nodes], [coords[n][1] for n in nodes], s=12, c=["#999999"], zorder=3)
    if node_support:
        for n in nodes:
            if n.support is not None:
                x, y = coords[n]
                ax.text(x, y + 0.05, str(n.support), fontsize=label_size - 2, color="#666666",
                        ha="center", va="bottom")

    if polar and theme == "void":
        ax.set_axis_off()
    elif not polar:
        tree_theme(ax, style=theme)
        ax.set_xlabel("branch length" if any(t.branch_length for t in tree.traverse("preorder") if t.branch_length) else "edges",
                      fontsize=axis_label_size)
        ax.tick_params(labelsize=tick_size, direction="out")
        if len(tips) <= 200:
            ax.set_yticks(list(range(len(tips))))
            ax.set_yticklabels([])
    ax.autoscale()
    # colourbar / legend for the tip field
    if tip_mappable is not None and show_tip_legend:
        cb = ax.figure.colorbar(tip_mappable, ax=ax)
        cb.set_label(tip_color_field, fontsize=9)
    elif tip_legend and show_tip_legend:
        import matplotlib.patches as mpatches
        handles = [mpatches.Patch(color=c, label=str(k)) for k, c in tip_legend]
        ax.legend(handles=handles, title=tip_color_field, fontsize=8, title_fontsize=9,
                  loc="upper left", bbox_to_anchor=(1.02, 1.0))
    if show_colorbar and mappable is not None:
        cb = ax.figure.colorbar(mappable, ax=ax)
        if colorbar_label:
            cb.set_label(colorbar_label)
        ax._treeio_mappable = mappable
    return ax


def _tip_colors_for(tree, tips, tip_colors, tip_color_field, cmap, default):
    """Return ``(colors, mappable, legend)`` for tip markers.

    If ``tip_color_field`` is given the colour is derived from each tip's value
    (continuous -> colourbar mappable; categorical -> discrete legend).
    """
    import numpy as np
    import matplotlib
    if tip_color_field is not None:
        values = {t.name: t.get_data(tip_color_field) for t in tips}
        known = {k: v for k, v in values.items() if v is not None}
        if known and all(isinstance(v, (int, float)) for v in known.values()):
            arr = np.array([float(values[t.name]) if values[t.name] is not None else np.nan for t in tips])
            lo = float(np.nanmin(arr)); hi = float(np.nanmax(arr))
            if not np.isfinite(lo) or lo == hi:
                lo, hi = 0.0, 1.0
            norm = matplotlib.colors.Normalize(vmin=lo, vmax=hi)
            cobj = matplotlib.colormaps.get_cmap(cmap)
            colors = cobj(norm(arr))
            mappable = matplotlib.cm.ScalarMappable(norm=norm, cmap=cobj)
            return colors, mappable, None
        cats = sorted(set(known.values()), key=str)
        pal = named_palette([str(c) for c in cats])
        colors = [pal.get(str(values[t.name]), default) for t in tips]
        return colors, None, [(c, pal[str(c)]) for c in cats]
    return [(tip_colors or {}).get(t.name, default) for t in tips], None, None


def _write_tip_label(ax, x, y, text, polar, color, size, offset=0.04):
    """Place a tip label; radial layouts rotate it so labels do not overlap."""
    if polar:
        import math
        theta = math.degrees(math.atan2(y, x))
        lx = x + math.cos(math.radians(theta)) * offset
        ly = y + math.sin(math.radians(theta)) * offset
        rot = theta
        if 90 < theta < 270:
            rot = theta - 180
            ha = "right"
        else:
            ha = "left"
        return ax.text(lx, ly, text, rotation=rot, rotation_mode="anchor", ha=ha, va="center",
                       fontsize=size, color=color, zorder=4)
    return ax.text(x, y, text, fontsize=size, color=color, va="center", ha="left", zorder=4)


# ----------------------------------------------------------------------- --
# object-oriented wrapper (matplotlib-style: ax.plot / ax.scatter)
# ----------------------------------------------------------------------- --
class TreePlotter:
    """An OO, matplotlib-compatible tree plotter.

    Holds a real :class:`matplotlib.axes.Axes` and exposes marker methods that
    delegate to the matching ``ax`` method (``.plot_edges()``,
    ``.scatter_tips()``, ``.scatter_nodes()``, ``.plot_tip_labels()``,
    ``.plot_node_labels()``, ``.legend()``, ``.savefig()``)."""

    def __init__(self, tree: Tree, ax=None, layout: str = "rectangular",
                 edge_color: str = "#333333", edge_width: float = 1.2, **kw):
        self.tree = tree
        self.layout = layout
        self.edge_color = edge_color
        self.edge_width = edge_width
        self.reverse_x = kw.pop("reverse_x", False)
        self.coords = tree_coords(tree, layout=layout, reverse_x=self.reverse_x)
        self._ax = _new_axes(tree, layout, self.reverse_x, ax=ax, **kw)

    @property
    def ax(self):
        """The underlying matplotlib axes (use any ``ax`` method freely)."""
        return self._ax

    def scatter_tips(self, size: float = 18, color=None, aes: str = None, cmap: str = "viridis", **kw):
        """Mark tips (``ax.scatter``).  ``aes`` colours by a per-tip field."""
        tips = self.tree.get_tips()
        c, _, _ = _tip_colors_for(self.tree, tips, None, aes, cmap, self.edge_color)
        return self._ax.scatter([self.coords[t][0] for t in tips],
                                [self.coords[t][1] for t in tips],
                                s=size, c=color or c, **kw)

    def scatter_nodes(self, size: float = 12, color="#999999", **kw):
        """Mark internal nodes (``ax.scatter``)."""
        nodes = self.tree.get_internal_nodes()
        return self._ax.scatter([self.coords[n][0] for n in nodes],
                                [self.coords[n][1] for n in nodes],
                                s=size, c=color, **kw)

    def plot_tip_labels(self, size: int = 10, aes: str = None, color=None,
                        radial_offset: float = 0.04, **kw):
        """Write tip labels; radial layouts rotate them (collision avoidance).

        ``aes`` colours the labels by a per-tip field."""
        polar = self.layout in ("circular", "fan", "radial", "unrooted")
        c, _, _ = _tip_colors_for(self.tree, self.tree.get_tips(), None, aes, "viridis", self.edge_color)
        out = []
        for i, tip in enumerate(self.tree.get_tips()):
            x, y = self.coords[tip]
            col = color or (c[i] if isinstance(c[i], str) else self.edge_color)
            out.append(_write_tip_label(self._ax, x, y, tip.name, polar, col, size, offset=radial_offset))
        return out

    def plot_node_labels(self, size: int = 8, color="#666666", show_support=True, **kw):
        """Write internal-node labels / support values (``ax.text``)."""
        out = []
        for node in self.tree.get_internal_nodes():
            label = node.name
            if (label in (None, "", "unknown")) and show_support and node.support is not None:
                label = str(node.support)
            if label in (None, "", "unknown"):
                continue
            x, y = self.coords[node]
            out.append(self._ax.text(x, y, label, fontsize=size, color=color, ha="center", va="center", **kw))
        return out

    def plot_edges(self, aes: str = None, color=None, width=None, cmap: str = "viridis", norm=None):
        """Draw the branches.

        If ``aes`` names a per-node field, each edge is coloured by the child
        node's value (a continuous colour scale); call :meth:`add_colorbar` to
        add a colourbar.  Otherwise a single ``LineCollection`` is used.
        """
        from matplotlib.collections import LineCollection
        segments = edge_segments(self.tree, self.coords, self.layout)
        if aes is not None:
            import numpy as np
            import matplotlib
            vals = [
                c.get_data(aes)
                for c in self.tree.traverse("preorder")
                if c.parent is not None
            ]
            arr = np.asarray([np.nan if v is None else float(v) for v in vals], dtype=float)
            lo = float(np.nanmin(arr))
            hi = float(np.nanmax(arr))
            if not np.isfinite(lo) or lo == hi:
                lo, hi = 0.0, 1.0
            if norm is None:
                norm = matplotlib.colors.Normalize(vmin=lo, vmax=hi)
            cmap_obj = matplotlib.colormaps.get_cmap(cmap)
            lc = LineCollection(segments, array=arr, cmap=cmap_obj, norm=norm,
                                linewidths=width or self.edge_width, capstyle="round")
            self._edge_mappable = lc
        else:
            self._edge_mappable = None
            lc = LineCollection(segments, colors=color or self.edge_color,
                                linewidths=width or self.edge_width, capstyle="round")
        self._ax.add_collection(lc)
        self._ax.autoscale()
        return lc

    def add_colorbar(self, label: str = "value", **kwargs):
        """Add a continuous colourbar for the ``aes``-coloured branches."""
        if getattr(self, "_edge_mappable", None) is None:
            raise RuntimeError("plot_edges(aes=...) must be called first")
        cb = self._ax.figure.colorbar(self._edge_mappable, ax=self._ax, **kwargs)
        if label:
            cb.set_label(label)
        return cb

    def legend_for(self, field: str, cmap: str = "viridis"):
        """Add a continuous colourbar or categorical legend for ``field``.

        Uses the tip values of ``field``; numeric -> colourbar, categorical -> a
        legend of categories.
        """
        import matplotlib
        tips = self.tree.get_tips()
        values = {t.name: t.get_data(field) for t in tips}
        known = {k: v for k, v in values.items() if v is not None}
        if known and all(isinstance(v, (int, float)) for v in known.values()):
            import numpy as np
            arr = np.array([float(values[t.name]) if values[t.name] is not None else np.nan for t in tips])
            lo = float(np.nanmin(arr)); hi = float(np.nanmax(arr))
            if not np.isfinite(lo) or lo == hi:
                lo, hi = 0.0, 1.0
            norm = matplotlib.colors.Normalize(vmin=lo, vmax=hi)
            mappable = matplotlib.cm.ScalarMappable(norm=norm, cmap=matplotlib.colormaps.get_cmap(cmap))
            cb = self._ax.figure.colorbar(mappable, ax=self._ax)
            cb.set_label(field, fontsize=9)
            return cb
        cats = sorted(set(known.values()), key=str)
        pal = named_palette([str(c) for c in cats])
        import matplotlib.patches as mpatches
        handles = [mpatches.Patch(color=pal[str(c)], label=str(c)) for c in cats]
        return self._ax.legend(handles=handles, title=field, fontsize=8, title_fontsize=9,
                               loc="upper left", bbox_to_anchor=(1.02, 1.0))

    def theme(self, style: str = "clean", **kwargs):
        """Apply a :func:`tree_theme` style to the axes."""
        tree_theme(self._ax, style=style, **kwargs)
        return self

    def scale_bar(self, unit: float = None, label: str = "branch length",
                  loc: Tuple[float, float] = None):
        """Draw a branch-unit scale bar on the axes."""
        add_scalebar(self._ax, self.coords, unit=unit, label=label, loc=loc)
        return self

    def legend(self, *args, **kwargs):
        return self._ax.legend(*args, **kwargs)

    def savefig(self, path, **kwargs):
        return self._ax.figure.savefig(path, **kwargs)


def treeplot(tree: Tree, ax=None, layout: str = "rectangular", **kwargs) -> TreePlotter:
    """Return a :class:`TreePlotter` (matplotlib-OO) for ``tree``."""
    return TreePlotter(tree, ax=ax, layout=layout, **kwargs)


# ----------------------------------------------------------------------- --
# unified entry points
# ----------------------------------------------------------------------- --
def render(tree: Tree, backend: str = "mpl", layout: str = "rectangular",
           path=None, savefig_kwargs=None, dpi: int = 300, format: str = None,
           tight: bool = True, **kwargs):
    """One entry point, several backends.

    ``backend`` is ``"mpl"`` (default; returns / saves a matplotlib axes) or
    ``"ascii"`` (returns an ASCII-art string).  If ``path`` is given the figure
    is written there with ``savefig`` (``dpi`` / ``format`` / ``tight`` control
    resolution, file format and bounding box).
    """
    backend = backend.lower().strip().lstrip("-")
    if backend in ("mpl", "matplotlib"):
        ax = draw_tree(tree, layout=layout, **kwargs)
        if path is not None:
            skw = dict(savefig_kwargs or {})
            skw.setdefault("bbox_inches", "tight" if tight else None)
            skw.setdefault("dpi", dpi)
            if format:
                skw.setdefault("format", format)
            ax.figure.savefig(path, **skw)
            return path
        return ax
    if backend in ("ascii", "text", "cli"):
        from .show import tree_to_ascii

        is_compact = kwargs.pop("is_compact", False)
        is_pruned = kwargs.pop("is_pruned", True)
        return tree_to_ascii(tree, is_compact, is_pruned)
    custom = _BACKEND_REGISTRY.get(backend)
    if custom is not None:
        return custom(tree, layout=layout, **kwargs)
    raise ValueError(f"unknown backend {backend!r}; choose from 'mpl', 'ascii', or a registered backend")


def draw(*args, **kwargs):
    """Alias of :func:`render` (matplotlib backend)."""
    return render(*args, **kwargs)


def plot(*args, **kwargs):
    """Alias of :func:`render` (matplotlib backend)."""
    return render(*args, **kwargs)


def save(tree: Tree, path: str, dpi: int = 300, format: str = None, tight: bool = True, **kwargs):
    """Draw ``tree`` and ``figure.savefig`` to ``path`` (publication-ready)."""
    from pathlib import Path
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    ax = render(tree, **kwargs)
    ax.figure.savefig(path, dpi=dpi, format=format, bbox_inches="tight" if tight else None)
    return path


# ----------------------------------------------------------------------- --
# publication-grade helpers
# ----------------------------------------------------------------------- --
def tree_theme(ax, style: str = "clean", grid: bool = False, spine_width: float = 0.8):
    """Apply a minimal publication theme to a tree axes.

    ``style`` may be ``"clean"`` (no top/right spines), ``"void"``
    (axis off), ``"minimal"`` (only a bottom spine) or ``"plain"`` (full box).
    """
    if style == "void":
        ax.set_axis_off()
        return ax
    if style == "minimal":
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.spines["bottom"].set_linewidth(spine_width)
        ax.set_yticks([])
    elif style == "plain":
        for s in ("top", "right", "bottom", "left"):
            ax.spines[s].set_visible(True)
            ax.spines[s].set_linewidth(spine_width)
    else:  # clean
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_linewidth(spine_width)
    ax.grid(grid)
    return ax


def _nice_unit(raw: float) -> float:
    """A round (1/2/5 * 10^k) value near ``raw``."""
    if raw <= 0:
        return 1.0
    import math as _m
    exp = _m.floor(_m.log10(raw))
    base = 10.0 ** exp
    for mult in (1, 2, 5, 10):
        if base * mult >= raw:
            return base * mult
    return base * 10


def add_scalebar(ax, coords, unit: float = None, label: str = "branch length",
                 loc: Tuple[float, float] = None):
    """Draw a branch-unit scale bar in data coordinates.

    ``coords`` is the node -> ``(x, y)`` mapping from :func:`tree_coords`.
    """
    import numpy as np
    xs = [c[0] for c in coords.values()]
    ys = [c[1] for c in coords.values()]
    lo_x, hi_x = min(xs), max(xs)
    lo_y, hi_y = min(ys), max(ys)
    span = (hi_x - lo_x) or 1.0
    if unit is None:
        unit = _nice_unit(span / 5.0)
    if loc is None:
        x0 = lo_x + 0.02 * span
        y0 = hi_y + 0.04 * (hi_y - lo_y + 1.0)
    else:
        x0, y0 = loc
    ax.plot([x0, x0 + unit], [y0, y0], color="#333333", lw=1.2)
    h = 0.02 * (hi_y - lo_y + 1.0)
    ax.plot([x0, x0], [y0 - h, y0 + h], color="#333333", lw=1.2)
    ax.plot([x0 + unit, x0 + unit], [y0 - h, y0 + h], color="#333333", lw=1.2)
    ax.text(x0 + unit / 2, y0 + 0.03 * (hi_y - lo_y + 1.0), _fmt_num(unit),
            ha="center", va="bottom", fontsize=8, color="#333333")
    if label:
        ax.text(x0, y0 - 1.2 * h, label, ha="left", va="top", fontsize=8, color="#333333")
    return ax


def facet_grid(tree: Tree, *fields: str, layout: str = "rectangular",
               figsize=None, cmap: str = "viridis", tip_labels: bool = True,
               share_y: bool = True, shared_legend: bool = True):
    """Draw a tree on the left and one per-tip data panel per ``field`` on the right.

    Returns the matplotlib ``Figure``.  Tips are aligned across panels.  When
    ``shared_legend`` is true a single colourbar / categorical legend is drawn
    at the figure edge.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.collections import LineCollection, PolyCollection
    import matplotlib
    fields = list(fields)
    ncols = 1 + len(fields)
    fig, axes = plt.subplots(1, ncols, figsize=figsize or ((4 + 2 * len(fields)), 6),
                             sharey=share_y)
    axes = list(axes) if ncols > 1 else [axes]
    coords = tree_coords(tree, layout=layout)
    tips = tree.get_tips()
    tip_y = [coords[t][1] for t in tips]

    # left panel: the tree
    tree_ax = axes[0]
    segs = edge_segments(tree, coords, layout=layout)
    tree_ax.add_collection(LineCollection(segs, colors="#333333", linewidths=1.2))
    tree_ax.scatter([coords[t][0] for t in tips], tip_y, s=14, c="#333333", zorder=3)
    if tip_labels:
        for t, x, y in zip(tips, [coords[t][0] for t in tips], tip_y):
            tree_ax.text(x, y, t.name, fontsize=8, va="center")
    tree_ax.set_ylabel("tree" if layout in ("circular", "unrooted") else "")
    tree_ax.autoscale()

    widths = 0.6
    legend = []
    for i, field in enumerate(fields):
        pax = axes[i + 1]
        values = {t.name: t.get_data(field) for t in tips}
        known = {k: v for k, v in values.items() if v is not None}
        numeric = bool(known) and all(isinstance(v, (int, float)) for v in known.values())
        if numeric:
            arr = np.array([float(values[t.name]) if values[t.name] is not None else np.nan for t in tips])
            lo = float(np.nanmin(arr)); hi = float(np.nanmax(arr))
            if not np.isfinite(lo) or lo == hi:
                lo, hi = 0.0, 1.0
            norm = matplotlib.colors.Normalize(vmin=lo, vmax=hi)
            colors = matplotlib.colormaps.get_cmap(cmap)(norm(arr))
        else:
            uniq = {v: i for i, v in enumerate(sorted(set(known.values()), key=str))}
            pal = named_palette([str(c) for c in uniq])
            colors = [pal.get(str(values[t.name]), "#cccccc") for t in tips]
            if shared_legend:
                legend.extend((c, pal[str(c)]) for c in uniq)
        # draw one rectangle per tip
        for j, (t, y) in enumerate(zip(tips, tip_y)):
            pax.add_patch(plt.Rectangle((0, y - widths / 2), widths, widths,
                                        facecolor=colors[j], edgecolor="none"))
        pax.set_xlim(0, widths)
        pax.set_ylim(min(tip_y) - 0.6, max(tip_y) + 0.6)
        pax.set_title(field, fontsize=9)
        pax.spines["top"].set_visible(False)
        pax.spines["right"].set_visible(False)
        pax.set_xticks([])
        pax.set_yticks([])

    # one shared colourbar at the right edge (first continuous field)
    if shared_legend and fields:
        first = fields[0]
        values = {t.name: t.get_data(first) for t in tips}
        known = {k: v for k, v in values.items() if v is not None}
        if known and all(isinstance(v, (int, float)) for v in known.values()):
            arr = np.array([float(values[t.name]) if values[t.name] is not None else np.nan for t in tips])
            lo = float(np.nanmin(arr)); hi = float(np.nanmax(arr))
            if not np.isfinite(lo) or lo == hi:
                lo, hi = 0.0, 1.0
            norm = matplotlib.colors.Normalize(vmin=lo, vmax=hi)
            mappable = matplotlib.cm.ScalarMappable(norm=norm, cmap=matplotlib.colormaps.get_cmap(cmap))
            cb = fig.colorbar(mappable, ax=axes[-1], fraction=0.03, pad=0.04)
            cb.ax.tick_params(labelsize=7)
            cb.set_label(first, fontsize=8)
        elif legend:
            import matplotlib.patches as mpatches
            seen = set()
            handles = [mpatches.Patch(color=c, label=str(k)) for k, c in legend if not (k in seen or seen.add(k))]
            if handles:
                fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8)

    fig.tight_layout()
    return fig


def grid_of_trees(trees, labels=None, layout: str = "rectangular", figsize=None,
                  share_y: bool = True, tip_labels: bool = True, ncols: int = None,
                  layouts=None, **kwargs):
    """Plot several trees in a grid sharing the y (tip) axis.

    Trees are arranged left-to-right (or wrapped to ``ncols`` columns); when
    ``share_y`` is true the tip y-axis is shared so tips align across panels.
    ``layouts`` may be a list giving the layout of each tree.  Returns the
    matplotlib ``Figure``.
    """
    import matplotlib.pyplot as plt
    trees = list(trees)
    if not trees:
        raise ValueError("no trees given")
    n = len(trees)
    if ncols is None:
        ncols = n
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize or (5 * ncols, 6 * nrows),
                             sharey=share_y, squeeze=False)
    max_tips = max(t.nleaves for t in trees)
    flat = axes.flatten()
    for i, t in enumerate(trees):
        ax = flat[i]
        lay = layouts[i] if layouts else layout
        draw_tree(t, ax=ax, layout=lay, tip_labels=tip_labels, tick_size=7, **kwargs)
        if labels is not None:
            ax.set_title(labels[i] if i < len(labels) else f"tree {i + 1}", fontsize=10)
        if share_y:
            ax.set_ylim(-0.5, max_tips - 0.5)
    for j in range(n, len(flat)):
        flat[j].set_axis_off()
    if share_y:
        flat[0].set_yticks([])
        flat[0].set_ylabel("")
    fig.tight_layout()
    return fig
