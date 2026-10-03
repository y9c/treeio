"""Reproduce the ggtree-style gallery figures with treeio.

Covers:
  1. Layout gallery (rectangular / roundrect / slanted / circular / fan ...)
  2. Circular trait-coloured trees + colorbar
  3. Clade-highlight boxes (rectangular, circular, unrooted)
  4. Concentric data rings (continuous + discrete)
  5. Tree + heatmap matrix (gheatmap)
  6. Radial tree: categorical branch colours + full-sector highlighted clades
     + node support labels
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math
import random
from pathlib import Path

from treeio import (
    Tree, attach, draw_tree, highlight_clade, add_ring, add_rings, gheatmap,
    tree_coords, named_palette,
)

OUT = Path(__file__).parent / "gallery"
OUT.mkdir(exist_ok=True)


def build(n, skew=0.7):
    """A random binary tree with sensible branch lengths."""
    pool = [Tree(f"t{i}", branch_length=random.uniform(0.5, 1.5)) for i in range(n)]
    while len(pool) > 2:
        a = pool.pop(random.randrange(len(pool)))
        b = pool.pop(random.randrange(len(pool)))
        c = Tree("internal", branch_length=random.uniform(0.3, 1.0) * skew)
        c.extend_children([a, b])
        pool.append(c)
    if len(pool) == 2:
        r = Tree("root")
        r.extend_children(pool)
        return r
    return pool[0]


def tidy(ax, polar=False):
    ax.set_axis_off()
    if polar:
        ax.axis("equal")
    return ax


def make_layout_gallery():
    fig, axes = plt.subplots(3, 3, figsize=(11, 11))
    grid = [
        [("rectangular", "A"), ("roundrect", "B"), ("slanted", "C")],
        [("slanted", "D"), ("roundrect", "E"), ("rectangular", "F")],
        [("circular", "G"), ("fan", "H"), ("roundrect", "I")],
    ]
    for r in range(3):
        for c in range(3):
            lay, letter = grid[r][c]
            ax = axes[r][c]
            draw_tree(build(36), ax=ax, layout=lay, tip_labels=False, tip_points=False)
            tidy(ax, polar=lay in ("circular", "fan", "radial", "unrooted"))
            ax.set_title(letter, loc="left", fontsize=12)
    fig.savefig(OUT / "01_layouts.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def make_trait():
    t = build(46)
    attach(t, {tip.name: round(random.uniform(3.5, 5.0), 2) for tip in t.get_tips()},
           key="trait")
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    for ax in axes:
        draw_tree(t, ax=ax, layout="fan", tip_labels=True, tip_points=True,
                  tip_color_field="trait", cmap="turbo", show_tip_legend=False)
        tidy(ax, polar=True)
    import matplotlib
    mappable = matplotlib.cm.ScalarMappable(
        norm=matplotlib.colors.Normalize(3.5, 5.0), cmap=matplotlib.colormaps["turbo"])
    fig.colorbar(mappable, ax=axes[0], fraction=0.03, pad=0.04).set_label("trait", fontsize=9)
    for i, ax in enumerate(axes):
        ax.set_title(chr(65 + i), loc="left", fontsize=13)
    fig.savefig(OUT / "02_trait.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def make_highlight():
    t = build(44)
    tips = t.get_tips()
    ca = t.get_mrca(tips[4], tips[10])
    cb = t.get_mrca(tips[22], tips[29])
    configs = [
        ("rectangular", "#2e8b57", "#4a80b5", False),
        ("circular", "#2e8b57", "#4a80b5", False),
        ("unrooted", "#4a80b5", "#2e8b57", False),
        ("rectangular", "#4a80b5", "#2e8b57", False),
        ("rectangular", "#e08080", "#40c0c0", False),
        ("rectangular", "#4a80b5", "#2e8b57", True),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    for idx, (lay, cA, cB, ext) in enumerate(configs):
        ax = axes[idx // 3][idx % 3]
        draw_tree(t, ax=ax, layout=lay, tip_labels=False, tip_points=False)
        highlight_clade(t, ca, ax=ax, layout=lay, fill=cA, alpha=0.4, extend_center=ext)
        highlight_clade(t, cb, ax=ax, layout=lay, fill=cB, alpha=0.4, extend_center=ext)
        tidy(ax, polar=lay in ("circular", "fan", "radial", "unrooted"))
        ax.set_title(chr(65 + idx), loc="left", fontsize=13)
    fig.savefig(OUT / "03_highlight.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def make_rings():
    t = build(72)
    attach(t, {tip.name: round(math.sin(i / 5.0) * 1.8, 2)
               for i, tip in enumerate(t.get_tips())}, key="continuous")
    letters = "abcdefgh"
    attach(t, {tip.name: letters[i % 8] for i, tip in enumerate(t.get_tips())}, key="d1")
    attach(t, {tip.name: letters[(i + 3) % 8] for i, tip in enumerate(t.get_tips())}, key="d2")
    fig, ax = plt.subplots(figsize=(10, 9))
    draw_tree(t, ax=ax, layout="circular", tip_labels=False, tip_points=False)
    add_ring(t, "continuous", ax=ax, cmap="magma", ring_width=0.06, pad=0.02)
    add_ring(t, "continuous", ax=ax, cmap="magma", ring_width=0.06)
    add_ring(t, "d1", ax=ax, discrete=True, palette=named_palette(list(letters)), ring_width=0.06)
    add_ring(t, "d2", ax=ax, discrete=True, palette=named_palette(list(letters)), ring_width=0.06)
    tidy(ax, polar=True)
    import matplotlib
    mappable = matplotlib.cm.ScalarMappable(
        norm=matplotlib.colors.Normalize(-1.8, 1.8), cmap=matplotlib.colormaps["magma"])
    fig.colorbar(mappable, ax=ax, fraction=0.03, pad=0.06).set_label("continuous\nvalue", fontsize=9)
    pal = named_palette(list(letters))
    handles = [mpatches.Patch(color=pal[k], label=k) for k in letters]
    ax.legend(handles=handles, title="discrete\nvalue", loc="upper left",
              bbox_to_anchor=(1.05, 0.98), fontsize=7, title_fontsize=8)
    fig.savefig(OUT / "04_rings.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def make_heatmap():
    t = build(26)
    cols = [f"c{i}" for i in range(20)]
    for col in cols:
        attach(t, {tip.name: round(random.uniform(0, 10), 1) for tip in t.get_tips()},
               key=col)
    fig, ax = plt.subplots(figsize=(14, 8))
    draw_tree(t, ax=ax, layout="rectangular", tip_labels=True, tip_points=False)
    gheatmap(t, cols, ax=ax, layout="rectangular", cmap="magma", tip_labels=False)
    ax.set_axis_off()
    fig.savefig(OUT / "05_heatmap.png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    # circular ring variant
    fig, ax = plt.subplots(figsize=(8, 8))
    draw_tree(t, ax=ax, layout="circular", tip_labels=False, tip_points=False)
    gheatmap(t, cols, ax=ax, layout="circular", cmap="magma")
    ax.set_axis_off()
    fig.savefig(OUT / "05_heatmap_ring.png", dpi=120, bbox_inches="tight")
    plt.close(fig)


def make_radial():
    t = build(24)
    tips = t.get_tips()
    for nd in t.get_internal_nodes():
        nd.support = random.choice([98, 100, 97, 90, 95])
    # Partition the tree into three disjoint clades for categorical branch colour.
    ca = t.get_mrca(tips[3], tips[8])
    cb = t.get_mrca(tips[12], tips[17])
    cc = t.get_mrca(tips[20], tips[23])
    for nd in t.traverse("preorder"):
        grp = None
        if ca and ca.is_ancestor_of(nd):
            grp = "clade A"
        elif cb and cb.is_ancestor_of(nd):
            grp = "clade B"
        elif cc and cc.is_ancestor_of(nd):
            grp = "clade C"
        nd.set_data("grp", grp)
    fig, ax = plt.subplots(figsize=(9, 9))
    draw_tree(t, ax=ax, layout="circular", tip_labels=True, label_size=7,
              tip_points=False, branch_color_field="grp", edge_width=2.0,
              node_support=True, label_radial_offset=0.05)
    highlight_clade(t, ca, ax=ax, layout="circular", fill="#bfe6e6", alpha=0.5, extend_center=True)
    highlight_clade(t, cb, ax=ax, layout="circular", fill="#f5c6c6", alpha=0.5, extend_center=True)
    highlight_clade(t, cc, ax=ax, layout="circular", fill="#dae9f2", alpha=0.5, extend_center=True)
    ax.set_axis_off()
    fig.savefig(OUT / "06_radial.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


def make_symbols():
    """Categorical per-tip symbol tree + legend (CG/GB/TC/TG/Tr/Tw)."""
    t = build(60)
    cats = ["CG", "GB", "TC", "TG", "Tr", "Tw"]
    attach(t, {tip.name: random.choice(cats) for tip in t.get_tips()}, key="stat")
    fig, ax = plt.subplots(figsize=(12, 8))
    draw_tree(t, ax=ax, layout="rectangular", tip_labels=True, tip_points=True,
              tip_color_field="stat", show_tip_legend=True)
    ax.set_axis_off()
    fig.savefig(OUT / "08_symbols.png", dpi=120, bbox_inches="tight")
    plt.close(fig)


def make_boxed():
    """Collapsed-clade boxed tree (green clade boxes + blue N + rotated Y)."""
    from treeio import box_label, add_boxed_labels

    yeast = [Tree(f"ADH{i}", branch_length=1.0) for i in range(1, 5)]
    fungi = Tree("Fungi", branch_length=1.0)
    fungi.extend_children(yeast)
    human = [Tree(f"hADH{i}", branch_length=1.0) for i in (1, 2)]
    ynode = Tree("Y", branch_length=0.6, support=100)
    ynode.extend_children(human)
    primates = Tree("primates", branch_length=1.0)
    primates.extend_children([ynode])
    insect = Tree("ADHX", branch_length=1.0)
    nematode = Tree("ADHY", branch_length=1.0)
    metazoa = Tree("metazoa", branch_length=1.0)
    metazoa.extend_children([primates, insect, nematode])
    t = Tree("N")
    t.extend_children([fungi, metazoa])

    coords = tree_coords(t, layout="rectangular")
    fig, ax = plt.subplots(figsize=(9, 6))
    draw_tree(t, ax=ax, layout="rectangular", tip_labels=False, tip_points=False, edge_width=1.4)
    add_boxed_labels(t, {
        "Fungi": "Fungi",
        "ADH1": "yeast", "ADH2": "yeast", "ADH3": "yeast", "ADH4": "yeast",
        "primates": "primates",
        "hADH1": "human", "hADH2": "human",
        "ADHX": "insect", "ADHY": "nematode",
    }, ax=ax, layout="rectangular", fill="#a7d7a7", fontsize=8)
    for name in ("N", "metazoa"):
        x, y = coords[t.get_node_by_label(name)]
        box_label(ax, x, y, "N", fill="#5b8bd0", textcolor="white", fontsize=8)
    x, y = coords[t.get_node_by_label("Y")]
    box_label(ax, x, y, "Y", fill="#5b8bd0", textcolor="white", fontsize=8, rotation=90, round=False)
    ax.text(x + 0.4, y, "100", fontsize=8, color="#333333", ha="left", va="center")
    for tip in t.get_tips():
        x, y = coords[tip]
        ax.text(x + 1.3, y, tip.name, fontsize=8, color="#333333", ha="left", va="center")
    ax.set_axis_off()
    fig.savefig(OUT / "07_boxed.png", dpi=120, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    random.seed(20)
    make_layout_gallery()
    make_trait()
    make_highlight()
    make_rings()
    make_heatmap()
    make_radial()
    make_boxed()
    make_symbols()
    print(f"wrote gallery to {OUT}")