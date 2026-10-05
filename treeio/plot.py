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

from .style import (
    _fmt_num,
    color_by_value,
    color_map,
    named_palette,
    suggest_figsize,
)
from .tree import Tree

__all__ = [
    "TreePlotter",
    "add_boxed_labels",
    "add_ring",
    "add_rings",
    "add_scalebar",
    "box_label",
    "color_by_value",
    "color_map",
    "draw",
    "draw_tree",
    "edge_segments",
    "facet_grid",
    "gheatmap",
    "grid_of_trees",
    "highlight_clade",
    "layouts",
    "plot",
    "register_backend",
    "register_layout",
    "render",
    "save",
    "tree_coords",
    "tree_theme",
    "treeplot",
    "unregister_backend",
    "unregister_layout",
]

# user-registered render backends / layouts
_BACKEND_REGISTRY: dict = {}
_LAYOUT_REGISTRY: dict = {}


_BUILTIN_BACKENDS = {"mpl", "matplotlib", "ascii", "text", "cli"}
_BUILTIN_LAYOUTS = {
    "rectangular",
    "roundrect",
    "slanted",
    "ellipse",
    "circular",
    "fan",
    "radial",
    "unrooted",
    "equal_angle",
    "daylight",
    "time",
}
_POLAR_LAYOUTS = frozenset(
    {"circular", "fan", "radial", "unrooted", "equal_angle", "daylight", "ellipse"}
)


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
def tree_coords(
    tree: Tree,
    layout: str = "rectangular",
    reverse_x: bool = False,
    use_branch_length: bool | None = None,
) -> dict[Tree, tuple[float, float]]:
    """Return ``{node: (x, y)}`` in data coordinates for ``layout``.

    Supported layouts: ``rectangular``, ``roundrect``, ``slanted``, ``ellipse``,
    ``circular``, ``fan``, ``radial``, ``unrooted``, ``equal_angle``,
    ``daylight`` and ``time``.  Unrecognised layouts raise ``ValueError`` (they
    are never silently drawn as a rectangular tree).
    """
    if layout in _LAYOUT_REGISTRY:
        return dict(_LAYOUT_REGISTRY[layout](tree))
    if layout == "equal_angle":
        return _equal_angle_coords(tree, reverse_x=reverse_x)
    if layout == "daylight":
        return _daylight_coords(tree, reverse_x=reverse_x)
    if layout == "ellipse":
        use_depth = False
        coords = _rect_coords(
            tree,
            reverse_x=reverse_x,
            use_depth=use_depth,
            use_branch_length=use_branch_length,
        )
        return _polarize(coords, tree, "ellipse")
    if layout not in _BUILTIN_LAYOUTS:
        raise ValueError(f"unknown tree layout {layout!r}")
    use_depth = layout == "radial"
    coords = _rect_coords(
        tree,
        reverse_x=reverse_x,
        use_depth=use_depth,
        use_branch_length=use_branch_length,
    )
    if layout in ("circular", "fan", "radial", "unrooted"):
        return _polarize(coords, tree, layout)
    return coords


def _rect_coords(
    tree: Tree, reverse_x: bool, use_depth: bool, use_branch_length: bool | None = None
) -> dict[Tree, tuple[float, float]]:
    """Rectangular (x = branch length / depth, y = tip order) node coordinates."""
    xpos = _x_positions(
        tree,
        reverse_x=reverse_x,
        use_depth=use_depth,
        use_branch_length=use_branch_length,
    )
    tips = tree.get_tips()
    tipy = {tip: float(i) for i, tip in enumerate(tips)}
    coords: dict[Tree, tuple[float, float]] = {}
    for node in tree.traverse("postorder"):
        if not node.children:
            y = tipy[node]
        else:
            y = sum(coords[c][1] for c in node.children) / len(node.children)
        coords[node] = (xpos[node], y)
    return coords


def _x_positions(
    tree: Tree, reverse_x: bool, use_depth: bool, use_branch_length: bool | None = None
) -> dict[Tree, float]:
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
        # ``use_branch_length``: None -> auto (has_len), True -> force, False -> depth
        if use_branch_length is None:
            use_branch_length = any(
                t.branch_length for t in tree.traverse("preorder") if t.branch_length
            )
        has_len = use_branch_length
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


def _polarize(coords: dict[Tree, tuple[float, float]], tree: Tree, layout: str):
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
    if layout == "ellipse":
        # Elliptical phylogram (ggtree ``ellipse``): tips are swept over a
        # partial arc so the crown is left open, and the radius is stretched on
        # the x-axis so the outline is a true ellipse rather than the unit
        # circle of ``circular``.  This is genuinely distinct from ``circular``:
        # the same ``(y/n)*2pi - pi/2`` mapping is NOT reused.
        sweep = 1.5 * math.pi  # 270 deg arc -> 90 deg open crown
        start = math.pi + math.pi / 4  # anchor so the opening faces left
        stretch = 1.2  # x-axis elongation -> elliptical outline
        for node, (x, y) in coords.items():
            r = x / xmax
            theta = start + (y / n) * sweep
            out[node] = (r * stretch * math.cos(theta), r * math.sin(theta))
        return out
    for node, (x, y) in coords.items():
        r = x / xmax
        theta = (
            ((y / n) * 2 * math.pi - math.pi / 2)
            if layout != "unrooted"
            else ((y / n) * 2 * math.pi)
        )
        out[node] = (r * math.cos(theta), r * math.sin(theta))
    return out


def _equal_angle_coords(
    tree: Tree, reverse_x: bool = False
) -> dict[Tree, tuple[float, float]]:
    """Equal-angle unrooted layout (Meacham / PHYLIP / PAUP*).

    Each subtree is allocated an angular arc proportional to its number of tips;
    the root is at the centre (angle 0, radius 0) and the tips spread around a
    full circle.  Radius is the (unit-normalised) cumulative branch length, or the
    node depth when no branch lengths are present (otherwise every node collapses
    to the origin, radius 0).  ``reverse_x`` inverts the radius so the root moves
    to the outer edge and the tips collapse toward the centre (mirroring the
    rectangular ``reverse_x``).
    """
    has_bl = _has_branch_lengths(tree)
    xmax = (
        max((_x_for(tree, n, has_bl) for n in tree.traverse("preorder")), default=1.0)
        or 1.0
    )
    ang: dict[Tree, float] = {}

    def assign(node: Tree, theta0: float, theta1: float) -> None:
        ang[node] = (theta0 + theta1) / 2.0
        if node.is_leaf():
            return
        kids = node.children
        width = theta1 - theta0
        # arc width proportional to each child's tip count
        total = sum(len(c.get_tips()) for c in kids) or 1
        cursor = theta0
        for c in kids:
            frac = len(c.get_tips()) / total
            assign(c, cursor, cursor + width * frac)
            cursor += width * frac

    assign(tree, 0.0, 2 * math.pi)
    out = {}
    for node in tree.traverse("preorder"):
        r = _x_for(tree, node, has_bl) / xmax
        if reverse_x:
            # invert the radius so the root sits on the outer edge and the tips
            # collapse toward the centre (mirrors the rectangular reverse_x).
            r = 1.0 - r
        th = ang[node]
        out[node] = (r * math.cos(th), r * math.sin(th))
    return out


def _has_branch_lengths(tree: Tree) -> bool:
    """True if any node in ``tree`` carries a (truthy) branch length."""
    return any(t.branch_length for t in tree.traverse("preorder") if t.branch_length)


def _x_for(tree: Tree, node: Tree, has_bl: bool | None = None) -> float:
    """Cumulative branch length from the root to ``node``.

    When no branch lengths are present (``has_bl`` False, or auto-detected), fall
    back to the node's depth so the equal-angle / daylight layouts still fan out
    instead of collapsing every node onto the origin (radius 0).
    """
    if has_bl is None:
        has_bl = _has_branch_lengths(tree)
    d = 0.0
    cur = node
    if not has_bl:
        while cur is not None and cur.parent is not None:
            d += 1.0
            cur = cur.parent
        return d
    while cur is not None and cur.parent is not None:
        d += cur.branch_length if cur.branch_length is not None else 0.0
        cur = cur.parent
    return d


def _daylight_coords(
    tree: Tree, reverse_x: bool = False
) -> dict[Tree, tuple[float, float]]:
    """Daylight unrooted layout (PAUP*).

    Start from the equal-angle layout and iteratively re-balance every interior
    node: each node's child subtrees are swung as rigid bodies around the node
    so that the tips are distributed more evenly (the classic daylight
    improvement).  Edges stay connected because a whole subtree is rotated by
    the same angle about their common ancestor.
    """
    coords = _equal_angle_coords(tree, reverse_x=reverse_x)

    def rotate_about(parent, child, delta):
        """Rotate the subtree rooted at ``child`` rigidly around ``parent``.

        Angles are measured about the *parent* node, not the global origin, so a
        child subtree is swung as a rigid body about its connection point.  Edges
        keep their length and direction relative to the parent, and symmetric
        subtrees stay distinct instead of swinging onto one another.
        """
        px, py = coords[parent]
        c, s = math.cos(delta), math.sin(delta)
        for x in child.traverse("preorder"):
            dx = coords[x][0] - px
            dy = coords[x][1] - py
            coords[x] = (px + dx * c - dy * s, py + dx * s + dy * c)

    for _ in range(3):
        for node in tree.traverse("preorder"):
            kids = list(node.children)
            if len(kids) < 2:
                continue
            px, py = coords[node]
            # local (around-the-parent) angle of each child
            local = {k: math.atan2(coords[k][1] - py, coords[k][0] - px) for k in kids}
            ordered = sorted(kids, key=lambda k: local[k] % (2 * math.pi))
            angs = [local[k] % (2 * math.pi) for k in ordered]
            # find the biggest gap between consecutive children
            gaps = [
                (angs[(i + 1) % len(angs)] - angs[i]) % (2 * math.pi)
                for i in range(len(angs))
            ]
            gi = gaps.index(max(gaps))
            start = angs[(gi + 1) % len(angs)]
            even = (2 * math.pi) / len(kids)
            # place each child (and its subtree) evenly, in angular order
            for j, k in enumerate(ordered):
                target = (start + j * even) % (2 * math.pi)
                cur = angs[j]
                delta = (target - cur) % (2 * math.pi)
                if delta > math.pi:
                    delta -= 2 * math.pi
                rotate_about(node, k, delta)
    return coords


def edge_segments(
    tree: Tree, coords=None, layout: str = "rectangular"
) -> list[list[tuple[float, float]]]:
    """Return a list of polylines (each a list of points) for every edge.

    ``roundrect`` renders L-shaped phylogram edges but with rounded corners at
    each elbow, matching the ``ggtree`` layout of the same name.  ``circular`` /
    ``fan`` / ``radial`` / ``unrooted`` / ``slanted`` draw straight edges.
    """
    if coords is None:
        coords = tree_coords(tree, layout)
    straight = (
        layout in _POLAR_LAYOUTS or layout == "slanted" or layout in _LAYOUT_REGISTRY
    )
    segs = []
    for node in tree.traverse("preorder"):
        if node.parent is None:
            continue
        xp, yp = coords[node.parent]
        xc, yc = coords[node]
        if straight:
            segs.append([(xp, yp), (xc, yc)])
        elif layout == "roundrect":
            segs.append(_rounded_lpath(xp, yp, xc, yc))
        else:
            segs.append([(xp, yp), (xp, yc), (xc, yc)])
    return segs


def _rounded_lpath(
    xp: float, yp: float, xc: float, yc: float, radius: float | None = None, n: int = 6
) -> list[tuple[float, float]]:
    """An L-shaped edge ``(xp,yp)->(xp,yc)->(xc,yc)`` with a rounded elbow.

    The corner at ``(xp, yc)`` is replaced by a quarter arc of the given
    ``radius`` (default: 40% of the shorter leg, clamped), drawn with ``n``
    points.  Used by the ``roundrect`` layout.
    """
    dx = xc - xp
    dy = yc - yp
    # default radius: fraction of the shorter leg (so the arc never flips)
    if radius is None:
        radius = 0.4 * min(abs(dx), abs(dy))
    r = min(radius, abs(dx) * 0.5, abs(dy) * 0.5)
    if r <= 1e-12:
        # degenerate (zero-length leg): plain L-shape
        return [(xp, yp), (xp, yc), (xc, yc)]
    sx = -1.0 if dx < 0 else 1.0
    sy = -1.0 if dy < 0 else 1.0
    pts = [(xp, yp)]
    # approach point on the vertical leg, ending r short of the corner
    pts.append((xp, yc - sy * r))
    # quarter arc centred at (xp + sx*r, yc - sy*r)
    import math

    cx, cy = xp + sx * r, yc - sy * r
    for i in range(1, n + 1):
        t = (i / n) * (math.pi / 2)
        # param so that t=0 -> (xp, yc - sy*r) and t=pi/2 -> (xp + sx*r, yc)
        px = cx - sx * r * math.cos(t)
        py = cy + sy * r * math.sin(t)
        pts.append((px, py))
    pts.append((xp + sx * r, yc))
    pts.append((xc, yc))
    return pts


def _mpl():
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as e:  # pragma: no cover
        raise ImportError(
            "matplotlib is required for treeio.plot; install it with `pip install matplotlib`"
        ) from e
    return plt


def _figsize(tree: Tree, layout: str):
    return suggest_figsize(tree, layout)


def _new_axes(tree: Tree, layout: str, reverse_x: bool = False, ax=None, **kw):
    plt = _mpl()
    polar = layout in _POLAR_LAYOUTS
    if ax is None:
        w, h = _figsize(tree, layout)
        figsize = kw.pop("figsize", (w / 80.0, h / 80.0))
        ax = plt.subplots(
            figsize=figsize, subplot_kw={"aspect": "equal"} if polar else None
        )[1]
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
    tip_colors: dict[str, str] | None = None,
    tip_color_field: str | None = None,
    show_tip_legend: bool = True,
    edge_color: str = "#333333",
    edge_width: float = 1.2,
    label_size: int = 10,
    label_radial_offset: float = 0.04,
    axis_label_size: int = 9,
    tick_size: int = 8,
    reverse_x: bool = False,
    branch_color_field: str | None = None,
    branch_palette: dict[str, str] | None = None,
    cmap: str = "viridis",
    theme: str = "clean",
    show_colorbar: bool = False,
    colorbar_label: str = "value",
    use_branch_length: bool | None = None,
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

    polar = layout in _POLAR_LAYOUTS
    if ax is None:
        figsize = kwargs.pop("figsize", None)
        w, h = _figsize(tree, layout)
        ax = plt.subplots(
            figsize=figsize or (w / 80.0, h / 80.0),
            subplot_kw={"aspect": "equal"} if polar else None,
        )[1]
    # ``use_branch_length=False`` draws a cladogram (branch.length='none'), which
    # also fixes the "collapsed core" circular artifact seen on skewed trees by
    # spacing tips evenly instead of scaling radius by a single longest path.
    coords = tree_coords(
        tree, layout=layout, reverse_x=reverse_x, use_branch_length=use_branch_length
    )
    segments = edge_segments(tree, coords, layout=layout)

    mappable = None
    branch_legend = None
    if branch_color_field is not None:
        import matplotlib
        import numpy as np

        childs = [c for c in tree.traverse("preorder") if c.parent is not None]
        vals = [c.get_data(branch_color_field) for c in childs]
        numeric = vals and all(
            v is not None and isinstance(v, (int, float, np.number)) for v in vals
        )
        if vals and numeric:
            arr = np.asarray(
                [np.nan if v is None else float(v) for v in vals], dtype=float
            )
            if np.all(np.isnan(arr)):
                lc = LineCollection(
                    segments,
                    colors=edge_color,
                    linewidths=edge_width,
                    capstyle="round",
                    antialiaseds=True,
                    zorder=1,
                )
            else:
                lo = float(np.nanmin(arr))
                hi = float(np.nanmax(arr))
                if not np.isfinite(lo) or lo == hi:
                    lo, hi = 0.0, 1.0
                norm = matplotlib.colors.Normalize(vmin=lo, vmax=hi)
                lc = LineCollection(
                    segments,
                    array=arr,
                    cmap=matplotlib.colormaps.get_cmap(cmap),
                    norm=norm,
                    linewidths=edge_width,
                    capstyle="round",
                    zorder=1,
                )
                mappable = lc
        else:
            # categorical branch colouring (discrete classes drawn as solid segs)
            alls = [str(v) for v in vals if v is not None]
            cats = sorted(set(alls), key=str)
            pal = branch_palette or named_palette(cats)
            colors = [
                pal.get(str(v), edge_color) if v is not None else edge_color
                for v in vals
            ]
            lc = LineCollection(
                segments,
                colors=colors,
                linewidths=edge_width,
                capstyle="round",
                antialiaseds=True,
                zorder=1,
            )
            import matplotlib.patches as mpatches

            branch_legend = [(c, pal[str(c)]) for c in cats]
    else:
        lc = LineCollection(
            segments,
            colors=edge_color,
            linewidths=edge_width,
            capstyle="round",
            antialiaseds=True,
            zorder=1,
        )
    ax.add_collection(lc)
    if branch_legend and show_tip_legend:
        import matplotlib.patches as mpatches

        handles = [mpatches.Patch(color=c, label=str(k)) for k, c in branch_legend]
        ax.legend(
            handles=handles,
            title=branch_color_field,
            fontsize=8,
            title_fontsize=9,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
        )

    tips = tree.get_tips()
    tip_x = [coords[t][0] for t in tips]
    tip_y = [coords[t][1] for t in tips]
    tip_c, tip_mappable, tip_legend = _tip_colors_for(
        tree, tips, tip_colors, tip_color_field, cmap, edge_color
    )
    if tip_points:
        # a single PathCollection is batched by matplotlib (already the compact
        # representation); nothing more is needed here.
        ax.scatter(
            tip_x, tip_y, s=18, c=tip_c, zorder=3, edgecolors="none", linewidths=0
        )
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
                crowd = sum(
                    1
                    for j in range(n)
                    if abs((angles[j] - a + math.pi) % (2 * math.pi) - math.pi)
                    < threshold
                )
                offset_map[i] = label_radial_offset * (1 + 0.3 * max(0, crowd - 1))
        else:
            offset_map = {}
        for i, (t, x, y, c) in enumerate(zip(tips, tip_x, tip_y, tip_c)):
            off = offset_map.get(i, label_radial_offset)
            _write_tip_label(ax, x, y, t.name, polar, c, label_size, offset=off)

    nodes = tree.get_internal_nodes()
    if node_points:
        ax.scatter(
            [coords[n][0] for n in nodes],
            [coords[n][1] for n in nodes],
            s=12,
            c=["#999999"],
            zorder=3,
        )
    if node_support:
        for n in nodes:
            if n.support is not None:
                x, y = coords[n]
                ax.text(
                    x,
                    y + 0.05,
                    str(n.support),
                    fontsize=label_size - 2,
                    color="#666666",
                    ha="center",
                    va="bottom",
                )

    if polar and theme == "void":
        ax.set_axis_off()
    elif not polar:
        tree_theme(ax, style=theme)
        ax.set_xlabel(
            "branch length"
            if any(
                t.branch_length for t in tree.traverse("preorder") if t.branch_length
            )
            else "edges",
            fontsize=axis_label_size,
        )
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
        ax.legend(
            handles=handles,
            title=tip_color_field,
            fontsize=8,
            title_fontsize=9,
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
        )
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
    import matplotlib
    import numpy as np

    if tip_color_field is not None:
        values = {t.name: t.get_data(tip_color_field) for t in tips}
        known = {k: v for k, v in values.items() if v is not None}
        if known and all(isinstance(v, (int, float)) for v in known.values()):
            arr = np.array(
                [
                    float(values[t.name]) if values[t.name] is not None else np.nan
                    for t in tips
                ]
            )
            lo = float(np.nanmin(arr))
            hi = float(np.nanmax(arr))
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

        # ``atan2`` returns (-180, 180]; normalise to [0, 360) so the left half
        # (90..270) includes both the upper-left AND lower-left arcs.  Without
        # the normalisation the lower-left labels never got flipped and rendered
        # mirrored.
        theta = math.degrees(math.atan2(y, x)) % 360.0
        lx = x + math.cos(math.radians(theta)) * offset
        ly = y + math.sin(math.radians(theta)) * offset
        rot = theta
        if 90 < theta < 270:
            rot = theta - 180
            ha = "right"
        else:
            ha = "left"
        return ax.text(
            lx,
            ly,
            text,
            rotation=rot,
            rotation_mode="anchor",
            ha=ha,
            va="center",
            fontsize=size,
            color=color,
            zorder=4,
        )
    return ax.text(
        x, y, text, fontsize=size, color=color, va="center", ha="left", zorder=4
    )


def box_label(
    ax,
    x,
    y,
    text,
    fill="#a7d7a7",
    textcolor="#1a1a1a",
    fontsize=9,
    boxstyle="round,pad=0.35",
    edgecolor="#333333",
    edge_width=0.8,
    round=True,
    ha="center",
    va="center",
    rotation=0,
    zorder=5,
):
    """Draw text inside a coloured, rounded box at ``(x, y)`` (data coords).

    This is the ``ggtree`` collapsed-clade box label.  The box is rendered by
    matplotlib's native text ``bbox`` so it is exactly sized to the text (and
    scales correctly for any axis): ``round=True`` draws a rounded rectangle,
    ``round=False`` a square one.  Returns the text object.
    """
    style = boxstyle if round else "square,pad=0.35"
    return ax.text(
        x,
        y,
        text,
        fontsize=fontsize,
        color=textcolor,
        ha=ha,
        va=va,
        rotation=rotation,
        zorder=zorder,
        bbox={
            "boxstyle": style,
            "facecolor": fill,
            "edgecolor": edgecolor,
            "linewidth": edge_width,
        },
    )


def add_boxed_labels(
    tree,
    labels,
    ax=None,
    layout="rectangular",
    fill="#a7d7a7",
    fontsize=9,
    round=True,
    coords=None,
    **kwargs,
):
    """Draw boxed labels at named nodes/clades (``ggtree`` clade-box labels).

    ``labels`` is a ``{node_or_name: text}`` mapping (or ``{node_or_name:
    (text, fill)}``).  Each target node is located by name (tip or internal) or
    by a :class:`Tree` node and a coloured rounded box with ``text`` is drawn at
    its coordinates.  Returns the axes.
    """
    if ax is None:
        import matplotlib.pyplot as plt

        ax = plt.gca()
    if coords is None:
        coords = tree_coords(tree, layout=layout)
    for target, spec in labels.items():
        node = (
            target
            if isinstance(target, Tree)
            else (
                tree.get_node_by_label(str(target))
                or tree.get_tip_by_label(str(target))
            )
        )
        if node is None:
            continue
        if isinstance(spec, (tuple, list)):
            text, f = spec
        else:
            text, f = spec, fill
        x, y = coords[node]
        box_label(ax, x, y, str(text), fill=f, fontsize=fontsize, round=round, **kwargs)
    return ax


# ----------------------------------------------------------------------- --
# clade highlight boxes (ggtree geom_hilight analogue)
# ----------------------------------------------------------------------- --
def _clade_bounds(tree, coords, clade, layout):
    """Data bounds covering ``clade``'s tips.

    Rectangular layouts return ``(x0, x1, y0, y1)``; the polar layouts return
    ``(r0, r1, theta0, theta1)`` in degrees (the shortest arc enclosing the
    clade's tips, handling wrap-around).
    """
    import math

    tips = clade.get_tips()
    if not tips:
        return None
    if layout in _POLAR_LAYOUTS:
        r0 = min(math.hypot(*coords[t]) for t in tips)
        r1 = max(math.hypot(*coords[t]) for t in tips)
        angles = sorted(
            math.degrees(math.atan2(coords[t][1], coords[t][0])) % 360.0 for t in tips
        )
        gaps = [angles[i + 1] - angles[i] for i in range(len(angles) - 1)]
        gaps.append(360.0 - (angles[-1] - angles[0]))
        gi = gaps.index(max(gaps))
        a0 = angles[(gi + 1) % len(angles)]
        a1 = angles[gi]
        if a1 < a0:
            a1 += 360.0
        return (r0, r1, a0, a1)
    # rectangular: the box starts at the clade's split node (x of the clade
    # itself), not at the root, matching ggtree geom_hilight.
    x0 = coords[clade][0]
    x1 = max(coords[t][0] for t in tips)
    y0 = min(coords[t][1] for t in tips)
    y1 = max(coords[t][1] for t in tips)
    return (x0, x1, y0, y1)


def highlight_clade(
    tree,
    clade,
    ax=None,
    layout="rectangular",
    fill="#2e8b57",
    alpha=0.35,
    edgecolor="none",
    lw=0,
    zorder=0.6,
    pad=0.45,
    coords=None,
    extend_center=False,
    **kwargs,
):
    """Shade the region behind a clade (``geom_hilight`` analogue).

    ``clade`` may be a :class:`Tree` node, a node name, or an iterable of tip
    names (their MRCA).  Rectangular layouts draw a filled rectangle covering
    the clade tips; polar layouts draw an annular wedge.  Returns the patch.
    """
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    if ax is None:
        ax = plt.gca()
    if isinstance(clade, str):
        clade = tree.get_node_by_label(clade) or tree.get_tip_by_label(clade)
    elif isinstance(clade, (list, tuple, set)):
        clade = tree.get_mrca(*list(clade))
    if clade is None:
        raise ValueError("clade not found in tree")
    if coords is None:
        coords = tree_coords(tree, layout=layout)
    b = _clade_bounds(tree, coords, clade, layout)
    if b is None:
        return None
    if layout in _POLAR_LAYOUTS:
        r0, r1, a0, a1 = b
        if r1 <= r0:
            r1 = r0 + 1e-3
        if extend_center:
            # full-sector wedge from the origin out to the clade's tip radius
            patch = mpatches.Wedge(
                (0, 0),
                r1,
                a0,
                a1,
                facecolor=fill,
                alpha=alpha,
                edgecolor=edgecolor,
                linewidth=lw,
                zorder=zorder,
                **kwargs,
            )
        else:
            patch = mpatches.Wedge(
                (0, 0),
                r1,
                a0,
                a1,
                width=r1 - r0,
                facecolor=fill,
                alpha=alpha,
                edgecolor=edgecolor,
                linewidth=lw,
                zorder=zorder,
                **kwargs,
            )
    else:
        x0, x1, y0, y1 = b
        # box starts at the clade's split node (x0), extends past its tips
        left = x0
        right = x1 + 0.02 * max(abs(x1), 1.0)
        patch = mpatches.Rectangle(
            (left, y0 - pad),
            right - left,
            (y1 - y0) + 2 * pad,
            facecolor=fill,
            alpha=alpha,
            edgecolor=edgecolor,
            linewidth=lw,
            zorder=zorder,
            **kwargs,
        )
    ax.add_patch(patch)
    return patch


# ----------------------------------------------------------------------- --
# concentric tip-data rings around a tree (ggtree ring / annoRing analogue)
# ----------------------------------------------------------------------- --
def _tip_angles(tree, i, n, layout, coords=None):
    """Angular extent (degrees) of tip index ``i`` in layout ``layout``.

    The wedge for a tip is pivoted onto the tip's *actual* polar angle
    (``atan2`` of its polarised coordinates) rather than the uniform
    ``(i/n)*360 - 90`` grid.  That grid is only correct for the layouts whose
    tips are evenly spaced on the circle starting at -90 (``circular`` /
    ``radial``).  For ``unrooted`` the base angle is 0, and for
    ``equal_angle`` / ``daylight`` the tips are not uniformly spaced at all, so
    the hard-coded grid puts the ring/heatmap wedges where no tip lives.

    Tips are sorted by their real angle and each tip is given the arc that runs
    from just before itself to just before its angular neighbour (cyclic), so
    the sectors tile the whole circle once -- no gaps, no overlaps -- and each
    wedge sits on top of the tip it colours.  Uniform layouts reduce to the old
    grid (each tip at the leading edge of a ``360/n`` sector), so ``circular``
    and ``radial`` behaviour is unchanged.

    The tiny ``1e-6`` back-shift is invisible and removes the floating-point
    knife-edge: it guarantees the wedge midpoint lands a hair *before* the
    ``tip + 180/n`` ideal, so callers that measure the sector/tip offset by
    wrapping to ``[0, 360)`` (see ``test_unrooted_ring_aligns_with_tips``) get a
    stable result instead of ``+/-1e-14`` float noise flipping the sign.
    """
    if layout == "fan":
        sweep = 180.0
        start = 90.0 - sweep / 2
        return start + (i / n) * sweep, start + ((i + 1) / n) * sweep
    if coords is None:
        coords = tree_coords(tree, layout)
    tips = tree.get_tips()
    m = len(tips) or 1
    angles = [
        math.degrees(math.atan2(coords[t][1], coords[t][0])) % 360.0 for t in tips
    ]
    # order tips by their actual angle (cyclic).  Each tip owns the arc from just
    # before itself to just before its next angular neighbour, so the whole ring
    # is tiled once with each tip sitting at its own sector's leading edge.
    order = sorted(range(m), key=lambda k: angles[k])
    pos = order.index(i)
    eps = 1e-6
    a0 = (angles[i] - eps) % 360.0
    nxt = order[(pos + 1) % m]
    width = (angles[nxt] - a0) % 360.0
    if width < 1e-9:  # degenerate: two tips on the same ray -- give an even slice
        width = 360.0 / m
    return a0, a0 + width


def add_ring(
    tree,
    field,
    ax=None,
    layout="circular",
    cmap="viridis",
    discrete=False,
    ring_width=0.06,
    start=0.12,
    pad=0.02,
    colorbar=False,
    colorbar_label=None,
    palette=None,
    norm=None,
    zorder=1.0,
    coords=None,
    **kwargs,
):
    """Draw one concentric ring of per-tip ``field`` values around a tree.

    The ring is an annulus centred on the origin whose segments are coloured by
    each tip's ``field`` value.  Continuous values map through ``cmap`` (add a
    colourbar with ``colorbar=True``); categorical values use a discrete palette.

    Rings are stacked outward: the first ring starts at the tree's outer tip
    radius (``start`` gap beyond it), and each subsequent ring is drawn farther
    out.  Returns ``(mappable_or_None, legend_handles)``.
    """
    import matplotlib
    import numpy as np
    from matplotlib.patches import Wedge

    if ax is None:
        import matplotlib.pyplot as plt

        ax = plt.gca()
    if coords is None:
        coords = tree_coords(tree, layout=layout)
    tips = tree.get_tips()
    n = len(tips) or 1
    values = {t.name: t.get_data(field) for t in tips}
    known = {k: v for k, v in values.items() if v is not None}
    if not known:
        return None, []
    base_outer = max(math.hypot(*coords[t]) for t in tips) or 1.0
    outer = getattr(ax, "_treeio_ring_r", base_outer + start)
    w = ring_width
    numeric = all(isinstance(v, (int, float)) for v in known.values())
    if numeric and not discrete:
        arr = np.array(
            [
                float(values[t.name]) if values[t.name] is not None else np.nan
                for t in tips
            ]
        )
        lo = float(np.nanmin(arr))
        hi = float(np.nanmax(arr))
        if not np.isfinite(lo) or lo == hi:
            lo, hi = 0.0, 1.0
        if norm is None:
            norm = matplotlib.colors.Normalize(vmin=lo, vmax=hi)
        cmap_obj = matplotlib.colormaps.get_cmap(cmap)
        colors = cmap_obj(norm(arr))
        mappable = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_obj)
        for i, tip in enumerate(tips):
            if values[tip.name] is None:
                continue
            a0, a1 = _tip_angles(tree, i, n, layout, coords)
            ax.add_patch(
                Wedge(
                    (0, 0),
                    outer + w,
                    a0,
                    a1,
                    width=w,
                    facecolor=colors[i],
                    edgecolor="none",
                    zorder=zorder,
                    **kwargs,
                )
            )
        ax._treeio_ring_r = outer + w + pad
        if colorbar:
            cb = ax.figure.colorbar(mappable, ax=ax, fraction=0.03, pad=0.04)
            cb.set_label(colorbar_label or field, fontsize=8)
        return mappable, None
    # discrete / categorical
    cats = sorted(set(known.values()), key=str)
    pal = palette or named_palette([str(c) for c in cats])
    for i, tip in enumerate(tips):
        if values[tip.name] is None:
            continue
        a0, a1 = _tip_angles(tree, i, n, layout, coords)
        ax.add_patch(
            Wedge(
                (0, 0),
                outer + w,
                a0,
                a1,
                width=w,
                facecolor=pal[str(values[tip.name])],
                edgecolor="none",
                zorder=zorder,
                **kwargs,
            )
        )
    ax._treeio_ring_r = outer + w + pad
    import matplotlib.patches as mpatches

    handles = [mpatches.Patch(color=c, label=str(k)) for k, c in pal.items()]
    return None, handles


def add_rings(
    tree,
    fields,
    ax=None,
    layout="circular",
    cmap="viridis",
    ring_width=0.06,
    start=0.12,
    pad=0.02,
    **kwargs,
):
    """Draw several concentric ``fields`` rings around a tree.

    ``fields`` is a list of field names; the first is drawn innermost.
    ``cmap`` may be a single name or a list matching ``fields``.  Returns ``ax``.
    """
    if isinstance(cmap, (list, tuple)):
        assert len(cmap) == len(fields), "cmap list must match fields"
        cmaps = list(cmap)
    else:
        cmaps = [cmap] * len(fields)
    # reset the stacking anchor for each call
    if hasattr(ax, "_treeio_ring_r"):
        del ax._treeio_ring_r
    for f, cm in zip(fields, cmaps):
        add_ring(
            tree,
            f,
            ax=ax,
            layout=layout,
            cmap=cm,
            ring_width=ring_width,
            start=start,
            pad=pad,
            **kwargs,
        )
    return ax


def gheatmap(
    tree,
    columns,
    ax=None,
    layout="rectangular",
    cmap="viridis",
    discrete=False,
    cell_width=0.8,
    cell_gap=0.0,
    group_gap=0.4,
    group_by=None,
    tip_labels=False,
    label_size=9,
    colorbar=False,
    colorbar_label="value",
    palette=None,
    zorder=2.0,
    coords=None,
    **kwargs,
):
    """Draw a per-tip heatmap matrix beside a tree (``gheatmap`` analogue).

    ``columns`` is a list of field names attached to tips; each is one column
    (rectangular layout) or one concentric ring (circular layout).  In the
    rectangular layout the tree is drawn first and the coloured matrix is placed
    immediately to its right, one row per tip, with optional gaps between
    ``group_by`` column groups.  Returns ``(ax, column_specs)`` where
    ``column_specs`` holds each column's left edge (rectangular) so the caller
    can draw column-axis labels.

    Continuous values share a single ``cmap``/norm across all columns;
    categorical values use a discrete palette (``discrete=True``).
    """
    import matplotlib
    import numpy as np
    from matplotlib.patches import Rectangle, Wedge

    if ax is None:
        import matplotlib.pyplot as plt

        ax = plt.gca()
    if coords is None:
        coords = tree_coords(tree, layout=layout)
    tips = tree.get_tips()
    n = len(tips) or 1
    # gather per-column values (aligned to tips; missing -> None)
    allvals = []
    numeric = True
    cats = set()
    for col in columns:
        vals = {t.name: t.get_data(col) for t in tips}
        allvals.append(vals)
        known = [v for v in vals.values() if v is not None]
        if not known:
            continue
        if not all(isinstance(v, (int, float, np.number)) for v in known):
            numeric = False
        cats.update(str(v) for v in known)
    if numeric and not discrete:
        # shared continuous colour scale
        known = [float(v) for vals in allvals for v in vals.values() if v is not None]
        lo = min(known) if known else 0.0
        hi = max(known) if known else 1.0
        if not np.isfinite(lo) or lo == hi:
            lo, hi = 0.0, 1.0
        norm = matplotlib.colors.Normalize(vmin=lo, vmax=hi)
        cmap_obj = matplotlib.colormaps.get_cmap(cmap)
        mappable = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap_obj)
        colfunc = lambda v: cmap_obj(norm(float(v)))
    else:
        cats = sorted(cats, key=str)
        pal = palette or named_palette(cats)
        mappable = None
        colfunc = lambda v: pal[str(v)]

    if layout in ("circular", "fan", "radial", "unrooted"):
        # concentric rings -- each column is drawn at a progressively larger radius
        base_outer = max(math.hypot(*coords[t]) for t in tips) or 1.0
        w = cell_width / 8.0
        for vals in allvals:
            outer = getattr(ax, "_treeio_ring_r", base_outer + 0.12)
            for i, tip in enumerate(tips):
                v = vals.get(tip.name)
                if v is None:
                    continue
                a0, a1 = _tip_angles(tree, i, n, layout, coords)
                ax.add_patch(
                    Wedge(
                        (0, 0),
                        outer + w,
                        a0,
                        a1,
                        width=w,
                        facecolor=colfunc(v),
                        edgecolor="none",
                        zorder=zorder,
                    )
                )
            ax._treeio_ring_r = outer + w + 0.02
        if colorbar and mappable is not None:
            cb = ax.figure.colorbar(mappable, ax=ax, fraction=0.03, pad=0.04)
            cb.set_label(colorbar_label, fontsize=8)
        return ax, None

    # rectangular: draw the matrix to the right of the tips
    tip_y = [coords[t][1] for t in tips]
    x0 = max(coords[t][0] for t in tips)
    left = x0 + 0.05 * max(abs(x0), 1.0)
    specs = []
    len(columns)
    for ci, vals in enumerate(allvals):
        # group separator
        if group_by and ci > 0 and group_by[ci - 1] is not None:
            left += group_gap
        w = cell_width
        specs.append({"left": left, "width": w, "field": columns[ci]})
        for i, tip in enumerate(tips):
            v = vals.get(tip.name)
            if v is None:
                continue
            ax.add_patch(
                Rectangle(
                    (left, tip_y[i] - 0.4),
                    w,
                    0.8,
                    facecolor=colfunc(v),
                    edgecolor="none",
                    zorder=zorder,
                )
            )
        left += w + cell_gap
    if tip_labels:
        for i, tip in enumerate(tips):
            ax.text(
                left + 0.6,
                tip_y[i],
                tip.name,
                fontsize=label_size,
                va="center",
                ha="left",
            )
    if colorbar and mappable is not None:
        cb = ax.figure.colorbar(mappable, ax=ax, fraction=0.03, pad=0.04)
        cb.set_label(colorbar_label, fontsize=8)
    ax.set_xlim((-0.05 * max(abs(x0), 1.0), left + (cell_width if tip_labels else 0)))
    return ax, specs


# ----------------------------------------------------------------------- --
# object-oriented wrapper (matplotlib-style: ax.plot / ax.scatter)
# ----------------------------------------------------------------------- --
class TreePlotter:
    """An OO, matplotlib-compatible tree plotter.

    Holds a real :class:`matplotlib.axes.Axes` and exposes marker methods that
    delegate to the matching ``ax`` method (``.plot_edges()``,
    ``.scatter_tips()``, ``.scatter_nodes()``, ``.plot_tip_labels()``,
    ``.plot_node_labels()``, ``.legend()``, ``.savefig()``)."""

    def __init__(
        self,
        tree: Tree,
        ax=None,
        layout: str = "rectangular",
        edge_color: str = "#333333",
        edge_width: float = 1.2,
        **kw,
    ):
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

    def scatter_tips(
        self,
        size: float = 18,
        color=None,
        aes: str | None = None,
        cmap: str = "viridis",
        **kw,
    ):
        """Mark tips (``ax.scatter``).  ``aes`` colours by a per-tip field."""
        tips = self.tree.get_tips()
        c, _, _ = _tip_colors_for(self.tree, tips, None, aes, cmap, self.edge_color)
        return self._ax.scatter(
            [self.coords[t][0] for t in tips],
            [self.coords[t][1] for t in tips],
            s=size,
            c=color or c,
            **kw,
        )

    def scatter_nodes(self, size: float = 12, color="#999999", **kw):
        """Mark internal nodes (``ax.scatter``)."""
        nodes = self.tree.get_internal_nodes()
        return self._ax.scatter(
            [self.coords[n][0] for n in nodes],
            [self.coords[n][1] for n in nodes],
            s=size,
            c=color,
            **kw,
        )

    def plot_tip_labels(
        self,
        size: int = 10,
        aes: str | None = None,
        color=None,
        radial_offset: float = 0.04,
        **kw,
    ):
        """Write tip labels; radial layouts rotate them (collision avoidance).

        ``aes`` colours the labels by a per-tip field."""
        polar = self.layout in _POLAR_LAYOUTS
        c, _, _ = _tip_colors_for(
            self.tree, self.tree.get_tips(), None, aes, "viridis", self.edge_color
        )
        out = []
        for i, tip in enumerate(self.tree.get_tips()):
            x, y = self.coords[tip]
            col = color or (c[i] if isinstance(c[i], str) else self.edge_color)
            out.append(
                _write_tip_label(
                    self._ax, x, y, tip.name, polar, col, size, offset=radial_offset
                )
            )
        return out

    def plot_node_labels(self, size: int = 8, color="#666666", show_support=True, **kw):
        """Write internal-node labels / support values (``ax.text``)."""
        out = []
        for node in self.tree.get_internal_nodes():
            label = node.name
            if (
                (label in (None, "", "unknown"))
                and show_support
                and node.support is not None
            ):
                label = str(node.support)
            if label in (None, "", "unknown"):
                continue
            x, y = self.coords[node]
            out.append(
                self._ax.text(
                    x,
                    y,
                    label,
                    fontsize=size,
                    color=color,
                    ha="center",
                    va="center",
                    **kw,
                )
            )
        return out

    def plot_edges(
        self,
        aes: str | None = None,
        color=None,
        width=None,
        cmap: str = "viridis",
        norm=None,
    ):
        """Draw the branches.

        If ``aes`` names a per-node field, each edge is coloured by the child
        node's value (a continuous colour scale); call :meth:`add_colorbar` to
        add a colourbar.  Otherwise a single ``LineCollection`` is used.
        """
        from matplotlib.collections import LineCollection

        segments = edge_segments(self.tree, self.coords, self.layout)
        if aes is not None:
            import matplotlib
            import numpy as np

            vals = [
                c.get_data(aes)
                for c in self.tree.traverse("preorder")
                if c.parent is not None
            ]
            arr = np.asarray(
                [np.nan if v is None else float(v) for v in vals], dtype=float
            )
            lo = float(np.nanmin(arr))
            hi = float(np.nanmax(arr))
            if not np.isfinite(lo) or lo == hi:
                lo, hi = 0.0, 1.0
            if norm is None:
                norm = matplotlib.colors.Normalize(vmin=lo, vmax=hi)
            cmap_obj = matplotlib.colormaps.get_cmap(cmap)
            lc = LineCollection(
                segments,
                array=arr,
                cmap=cmap_obj,
                norm=norm,
                linewidths=width or self.edge_width,
                capstyle="round",
            )
            self._edge_mappable = lc
        else:
            self._edge_mappable = None
            lc = LineCollection(
                segments,
                colors=color or self.edge_color,
                linewidths=width or self.edge_width,
                capstyle="round",
            )
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

            arr = np.array(
                [
                    float(values[t.name]) if values[t.name] is not None else np.nan
                    for t in tips
                ]
            )
            lo = float(np.nanmin(arr))
            hi = float(np.nanmax(arr))
            if not np.isfinite(lo) or lo == hi:
                lo, hi = 0.0, 1.0
            norm = matplotlib.colors.Normalize(vmin=lo, vmax=hi)
            mappable = matplotlib.cm.ScalarMappable(
                norm=norm, cmap=matplotlib.colormaps.get_cmap(cmap)
            )
            cb = self._ax.figure.colorbar(mappable, ax=self._ax)
            cb.set_label(field, fontsize=9)
            return cb
        cats = sorted(set(known.values()), key=str)
        pal = named_palette([str(c) for c in cats])
        import matplotlib.patches as mpatches

        handles = [mpatches.Patch(color=pal[str(c)], label=str(c)) for c in cats]
        return self._ax.legend(
            handles=handles,
            title=field,
            fontsize=8,
            title_fontsize=9,
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
        )

    def theme(self, style: str = "clean", **kwargs):
        """Apply a :func:`tree_theme` style to the axes."""
        tree_theme(self._ax, style=style, **kwargs)
        return self

    def scale_bar(
        self,
        unit: float | None = None,
        label: str = "branch length",
        loc: tuple[float, float] | None = None,
    ):
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
def _apply_fast(ax, threshold: int = 20000):
    """Rasterise any collection with a very large number of path vertices.

    A single dense layer (e.g. tens of thousands of coloured rings / wedges /
    per-node patches) bloats a vector PDF with millions of path operators and
    colour-state switches.  Rasterising just those huge collections turns them
    into one crisp image (rendered at the current savefig DPI) while everything
    else -- branches, labels -- stays editable vector.
    """
    n = 0
    for coll in ax.collections:
        total = 0
        if not hasattr(coll, "get_paths"):
            continue
        for p in coll.get_paths():
            total += len(p.vertices)
        if total > threshold:
            coll.set_rasterized(True)
            n += 1
    return n


def render(
    tree: Tree,
    backend: str = "mpl",
    layout: str = "rectangular",
    path=None,
    savefig_kwargs=None,
    dpi: int = 300,
    format: str | None = None,
    tight: bool = True,
    **kwargs,
):
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
            _apply_fast(ax)
            from pathlib import Path

            Path(path).parent.mkdir(parents=True, exist_ok=True)
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
    raise ValueError(
        f"unknown backend {backend!r}; choose from 'mpl', 'ascii', or a registered backend"
    )


def draw(*args, **kwargs):
    """Alias of :func:`render` (matplotlib backend)."""
    return render(*args, **kwargs)


def plot(*args, **kwargs):
    """Alias of :func:`render` (matplotlib backend)."""
    return render(*args, **kwargs)


def save(
    tree: Tree,
    path: str,
    dpi: int = 300,
    format: str | None = None,
    tight: bool = True,
    **kwargs,
):
    """Draw ``tree`` and ``figure.savefig`` to ``path`` (publication-ready).

    Very dense layers are automatically drawn as a crisp image so large trees
    open quickly while branches and labels stay editable vector.
    """
    from pathlib import Path

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    return render(
        tree,
        backend="mpl",
        path=path,
        dpi=dpi,
        format=format,
        tight=tight,
        **kwargs,
    )


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
    base = 10.0**exp
    for mult in (1, 2, 5, 10):
        if base * mult >= raw:
            return base * mult
    return base * 10


def add_scalebar(
    ax,
    coords,
    unit: float | None = None,
    label: str = "branch length",
    loc: tuple[float, float] | None = None,
):
    """Draw a branch-unit scale bar in data coordinates.

    ``coords`` is the node -> ``(x, y)`` mapping from :func:`tree_coords`.
    """

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
    ax.text(
        x0 + unit / 2,
        y0 + 0.03 * (hi_y - lo_y + 1.0),
        _fmt_num(unit),
        ha="center",
        va="bottom",
        fontsize=8,
        color="#333333",
    )
    if label:
        ax.text(
            x0, y0 - 1.2 * h, label, ha="left", va="top", fontsize=8, color="#333333"
        )
    return ax


def facet_grid(
    tree: Tree,
    *fields: str,
    layout: str = "rectangular",
    figsize=None,
    cmap: str = "viridis",
    tip_labels: bool = True,
    share_y: bool = True,
    shared_legend: bool = True,
):
    """Draw a tree on the left and one per-tip data panel per ``field`` on the right.

    Returns the matplotlib ``Figure``.  Tips are aligned across panels.  When
    ``shared_legend`` is true a single colourbar / categorical legend is drawn
    at the figure edge.
    """
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.collections import LineCollection

    fields = list(fields)
    ncols = 1 + len(fields)
    fig, axes = plt.subplots(
        1, ncols, figsize=figsize or ((4 + 2 * len(fields)), 6), sharey=share_y
    )
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
        numeric = bool(known) and all(
            isinstance(v, (int, float)) for v in known.values()
        )
        if numeric:
            arr = np.array(
                [
                    float(values[t.name]) if values[t.name] is not None else np.nan
                    for t in tips
                ]
            )
            lo = float(np.nanmin(arr))
            hi = float(np.nanmax(arr))
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
            pax.add_patch(
                plt.Rectangle(
                    (0, y - widths / 2),
                    widths,
                    widths,
                    facecolor=colors[j],
                    edgecolor="none",
                )
            )
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
            arr = np.array(
                [
                    float(values[t.name]) if values[t.name] is not None else np.nan
                    for t in tips
                ]
            )
            lo = float(np.nanmin(arr))
            hi = float(np.nanmax(arr))
            if not np.isfinite(lo) or lo == hi:
                lo, hi = 0.0, 1.0
            norm = matplotlib.colors.Normalize(vmin=lo, vmax=hi)
            mappable = matplotlib.cm.ScalarMappable(
                norm=norm, cmap=matplotlib.colormaps.get_cmap(cmap)
            )
            cb = fig.colorbar(mappable, ax=axes[-1], fraction=0.03, pad=0.04)
            cb.ax.tick_params(labelsize=7)
            cb.set_label(first, fontsize=8)
        elif legend:
            import matplotlib.patches as mpatches

            seen = set()
            handles = [
                mpatches.Patch(color=c, label=str(k))
                for k, c in legend
                if not (k in seen or seen.add(k))
            ]
            if handles:
                fig.legend(
                    handles=handles,
                    loc="center left",
                    bbox_to_anchor=(1.0, 0.5),
                    fontsize=8,
                )

    fig.tight_layout()
    return fig


def grid_of_trees(
    trees,
    labels=None,
    layout: str = "rectangular",
    figsize=None,
    share_y: bool = True,
    tip_labels: bool = True,
    ncols: int | None = None,
    layouts=None,
    **kwargs,
):
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
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=figsize or (5 * ncols, 6 * nrows),
        sharey=share_y,
        squeeze=False,
    )
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
