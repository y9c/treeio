/*
 * treeio live gallery — every example is a self-contained, browser-runnable
 * snippet.  Each entry has: id, title, desc, preview (static PNG for instant
 * first paint) and code (executed by runner.js when live WASM is ready).
 */
window.TREEIO_EXAMPLES = [
  {
    id: 'layouts',
    title: 'Layout gallery (A–I)',
    desc: 'The full set of tree layouts: rectangular, roundrect, slanted, circular and fan.',
    preview: 'static/previews/01_layouts.png',
    code: `
from treeio import Tree, draw_tree
import matplotlib.pyplot as plt
import random
random.seed(20)

def build(n):
    pool = [Tree(f"t{i}", branch_length=random.uniform(0.5, 1.5)) for i in range(n)]
    while len(pool) > 2:
        a = pool.pop(random.randrange(len(pool)))
        b = pool.pop(random.randrange(len(pool)))
        c = Tree("internal", branch_length=random.uniform(0.3, 1.0))
        c.extend_children([a, b])
        pool.append(c)
    r = Tree("root"); r.extend_children(pool); return r

fig, axes = plt.subplots(3, 3, figsize=(11, 11))
grid = [["rectangular", "roundrect", "slanted"],
        ["slanted", "roundrect", "rectangular"],
        ["circular", "fan", "roundrect"]]
for r in range(3):
    for c in range(3):
        ax = axes[r][c]
        lay = grid[r][c]
        draw_tree(build(36), ax=ax, layout=lay, tip_labels=False, tip_points=False)
        ax.set_axis_off()
        if lay in ("circular", "fan"):
            ax.axis("equal")
        ax.set_title(chr(65 + r * 3 + c), loc="left", fontsize=12)
plt.tight_layout()
`
  },
  {
    id: 'trait',
    title: 'Trait-coloured circular tree',
    desc: 'Tips and branches coloured by a continuous trait with a colorbar (ggtree geom_tippoint + scale_color_viridis).',
    preview: 'static/previews/02_trait.png',
    code: `
from treeio import Tree, attach, draw_tree
import matplotlib.pyplot as plt
import random
random.seed(1)

def build(n):
    pool = [Tree(f"t{i}", branch_length=random.uniform(0.4, 1.2)) for i in range(n)]
    while len(pool) > 2:
        a = pool.pop(random.randrange(len(pool)))
        b = pool.pop(random.randrange(len(pool)))
        c = Tree("internal", branch_length=random.uniform(0.2, 0.7))
        c.extend_children([a, b])
        pool.append(c)
    r = Tree("root"); r.extend_children(pool); return r

t = build(46)
attach(t, {tip.name: round(random.uniform(3.5, 5.0), 2) for tip in t.get_tips()}, key="trait")
fig, ax = plt.subplots(figsize=(7, 7))
draw_tree(t, ax=ax, layout="fan", tip_labels=True, tip_points=True,
          tip_color_field="trait", cmap="turbo", show_tip_legend=False)
ax.set_axis_off()
import matplotlib
m = matplotlib.cm.ScalarMappable(
    norm=matplotlib.colors.Normalize(3.5, 5.0), cmap=matplotlib.colormaps["turbo"])
fig.colorbar(m, ax=ax, fraction=0.03, pad=0.04).set_label("trait", fontsize=9)
`
  },
  {
    id: 'highlight',
    title: 'Clade highlight boxes',
    desc: 'Shade clades with highlight_clade() (geom_hilight): rectangular, circular and unrooted, including full sectors.',
    preview: 'static/previews/03_highlight.png',
    code: `
from treeio import Tree, draw_tree, highlight_clade
import matplotlib.pyplot as plt
import random
random.seed(9)

def build(n):
    pool = [Tree(f"t{i}", branch_length=random.uniform(0.4, 1.2)) for i in range(n)]
    while len(pool) > 2:
        a = pool.pop(random.randrange(len(pool)))
        b = pool.pop(random.randrange(len(pool)))
        c = Tree("internal", branch_length=random.uniform(0.2, 0.7))
        c.extend_children([a, b])
        pool.append(c)
    r = Tree("root"); r.extend_children(pool); return r

t = build(44)
tips = t.get_tips()
ca = t.get_mrca(tips[4], tips[10])
cb = t.get_mrca(tips[22], tips[29])
fig, axes = plt.subplots(1, 2, figsize=(13, 6))
draw_tree(t, ax=axes[0], layout="rectangular", tip_labels=False, tip_points=False)
highlight_clade(t, ca, ax=axes[0], layout="rectangular", fill="#2e8b57", alpha=0.4)
highlight_clade(t, cb, ax=axes[0], layout="rectangular", fill="#4a80b5", alpha=0.4)
axes[0].set_axis_off()
draw_tree(t, ax=axes[1], layout="circular", tip_labels=False, tip_points=False)
highlight_clade(t, ca, ax=axes[1], layout="circular", fill="#2e8b57", alpha=0.4)
highlight_clade(t, cb, ax=axes[1], layout="circular", fill="#4a80b5", alpha=0.4)
axes[1].set_axis_off(); axes[1].axis("equal")
plt.tight_layout()
`
  },
  {
    id: 'rings',
    title: 'Concentric data rings',
    desc: 'Stack continuous and discrete trait rings around a circular tree (ring / geom_annoRing).',
    preview: 'static/previews/04_rings.png',
    code: `
from treeio import Tree, attach, draw_tree, add_ring, named_palette
import matplotlib.pyplot as plt
import math

def build(n):
    pool = [Tree(f"t{i}", branch_length=0.5 + (i % 5) * 0.2) for i in range(n)]
    while len(pool) > 2:
        a = pool.pop(random.randrange(len(pool)))
        b = pool.pop(random.randrange(len(pool)))
        c = Tree("internal", branch_length=0.4)
        c.extend_children([a, b])
        pool.append(c)
    r = Tree("root"); r.extend_children(pool); return r

import random
random.seed(3)
t = build(72)
attach(t, {tip.name: round(math.sin(i / 5.0) * 1.8, 2) for i, tip in enumerate(t.get_tips())}, key="continuous")
letters = "abcdefgh"
attach(t, {tip.name: letters[i % 8] for i, tip in enumerate(t.get_tips())}, key="d1")
attach(t, {tip.name: letters[(i + 3) % 8] for i, tip in enumerate(t.get_tips())}, key="d2")
fig, ax = plt.subplots(figsize=(10, 9))
draw_tree(t, ax=ax, layout="circular", tip_labels=False, tip_points=False)
add_ring(t, "continuous", ax=ax, cmap="magma", ring_width=0.06, pad=0.02)
add_ring(t, "continuous", ax=ax, cmap="magma", ring_width=0.06)
add_ring(t, "d1", ax=ax, discrete=True, palette=named_palette(list(letters)), ring_width=0.06)
add_ring(t, "d2", ax=ax, discrete=True, palette=named_palette(list(letters)), ring_width=0.06)
ax.set_axis_off(); ax.axis("equal")
import matplotlib
a = matplotlib.cm.ScalarMappable(
    norm=matplotlib.colors.Normalize(-1.8, 1.8), cmap=matplotlib.colormaps["magma"])
fig.colorbar(a, ax=ax, fraction=0.03, pad=0.06).set_label("continuous\nvalue", fontsize=9)
`
  },
  {
    id: 'heatmap',
    title: 'Tree + heatmap matrix',
    desc: 'gheatmap() draws a per-tip data matrix beside (rectangular) or around (circular) the tree.',
    preview: 'static/previews/05_heatmap.png',
    code: `
from treeio import Tree, attach, draw_tree, gheatmap
import matplotlib.pyplot as plt
import random
random.seed(5)

def build(n):
    pool = [Tree(f"t{i}", branch_length=random.uniform(0.4, 1.4)) for i in range(n)]
    while len(pool) > 2:
        a = pool.pop(random.randrange(len(pool)))
        b = pool.pop(random.randrange(len(pool)))
        c = Tree("internal", branch_length=random.uniform(0.2, 0.8))
        c.extend_children([a, b])
        pool.append(c)
    r = Tree("root"); r.extend_children(pool); return r

t = build(26)
cols = [f"c{i}" for i in range(20)]
for col in cols:
    attach(t, {tip.name: round(random.uniform(0, 10), 1) for tip in t.get_tips()}, key=col)
fig, ax = plt.subplots(figsize=(14, 8))
draw_tree(t, ax=ax, layout="rectangular", tip_labels=True, tip_points=False)
gheatmap(t, cols, ax=ax, layout="rectangular", cmap="magma", tip_labels=False)
ax.set_axis_off()
`
  },
  {
    id: 'heatmap_ring',
    title: 'Tree + heatmap (circular rings)',
    desc: 'gheatmap() also draws the per-tip matrix as concentric rings around a circular tree.',
    preview: 'static/previews/05_heatmap_ring.png',
    code: `
from treeio import Tree, attach, draw_tree, gheatmap
import matplotlib.pyplot as plt
import random
random.seed(6)

def build(n):
    pool = [Tree(f"t{i}", branch_length=random.uniform(0.4, 1.4)) for i in range(n)]
    while len(pool) > 2:
        a = pool.pop(random.randrange(len(pool)))
        b = pool.pop(random.randrange(len(pool)))
        c = Tree("internal", branch_length=random.uniform(0.2, 0.8))
        c.extend_children([a, b])
        pool.append(c)
    r = Tree("root"); r.extend_children(pool); return r

t = build(26)
cols = [f"c{i}" for i in range(20)]
for col in cols:
    attach(t, {tip.name: round(random.uniform(0, 10), 1) for tip in t.get_tips()}, key=col)
fig, ax = plt.subplots(figsize=(8, 8))
draw_tree(t, ax=ax, layout="circular", tip_labels=False, tip_points=False)
gheatmap(t, cols, ax=ax, layout="circular", cmap="magma")
ax.set_axis_off(); ax.axis("equal")
`
  },
  {
    id: 'radial',
    title: 'Radial tree: category branches + support',
    desc: 'Circular tree with categorical branch colours, full-width clade sectors and node support labels.',
    preview: 'static/previews/06_radial.png',
    code: `
from treeio import Tree, draw_tree, highlight_clade
import matplotlib.pyplot as plt
import random
random.seed(2)

def build(n):
    pool = [Tree(f"t{i}", branch_length=random.uniform(0.4, 1.2)) for i in range(n)]
    while len(pool) > 2:
        a = pool.pop(random.randrange(len(pool)))
        b = pool.pop(random.randrange(len(pool)))
        c = Tree("internal", branch_length=random.uniform(0.2, 0.7))
        c.extend_children([a, b])
        pool.append(c)
    r = Tree("root"); r.extend_children(pool); return r

t = build(24)
tips = t.get_tips()
for nd in t.get_internal_nodes():
    nd.support = random.choice([98, 100, 97, 90, 95])
ca = t.get_mrca(tips[3], tips[8])
cb = t.get_mrca(tips[12], tips[17])
cc = t.get_mrca(tips[20], tips[23])
for nd in t.traverse("preorder"):
    grp = None
    if ca and ca.is_ancestor_of(nd): grp = "clade A"
    elif cb and cb.is_ancestor_of(nd): grp = "clade B"
    elif cc and cc.is_ancestor_of(nd): grp = "clade C"
    nd.set_data("grp", grp)
fig, ax = plt.subplots(figsize=(9, 9))
draw_tree(t, ax=ax, layout="circular", tip_labels=True, label_size=7,
          tip_points=False, branch_color_field="grp", edge_width=2.0, node_support=True)
highlight_clade(t, ca, ax=ax, layout="circular", fill="#bfe6e6", alpha=0.5, extend_center=True)
highlight_clade(t, cb, ax=ax, layout="circular", fill="#f5c6c6", alpha=0.5, extend_center=True)
highlight_clade(t, cc, ax=ax, layout="circular", fill="#dae9f2", alpha=0.5, extend_center=True)
ax.set_axis_off(); ax.axis("equal")
`
  },
  {
    id: 'boxed',
    title: 'Collapsed-clade boxed tree',
    desc: 'Collapsed clades shown as coloured boxed labels (geom_label) with a rotated branch callout.',
    preview: 'static/previews/07_boxed.png',
    code: `
from treeio import Tree, draw_tree, box_label, add_boxed_labels, tree_coords
import matplotlib.pyplot as plt

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
  "Fungi": "Fungi", "ADH1": "yeast", "ADH2": "yeast",
  "ADH3": "yeast", "ADH4": "yeast", "primates": "primates",
  "hADH1": "human", "hADH2": "human", "ADHX": "insect", "ADHY": "nematode",
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
`
  },
  {
    id: 'symbols',
    title: 'Categorical symbol tree',
    desc: 'Per-tip points coloured by a discrete stat with an auto legend (geom_tippoint + discrete scale).',
    preview: 'static/previews/08_symbols.png',
    code: `
from treeio import Tree, attach, draw_tree
import matplotlib.pyplot as plt
import random
random.seed(3)

def build(n):
    pool = [Tree(f"t{i}", branch_length=random.uniform(0.5, 1.5)) for i in range(n)]
    while len(pool) > 2:
        a = pool.pop(random.randrange(len(pool)))
        b = pool.pop(random.randrange(len(pool)))
        c = Tree("internal", branch_length=random.uniform(0.3, 0.9))
        c.extend_children([a, b])
        pool.append(c)
    r = Tree("root"); r.extend_children(pool); return r

t = build(60)
attach(t, {tip.name: random.choice(["CG", "GB", "TC", "TG", "Tr", "Tw"]) for tip in t.get_tips()}, key="stat")
fig, ax = plt.subplots(figsize=(12, 8))
draw_tree(t, ax=ax, layout="rectangular", tip_labels=True, tip_points=True,
          tip_color_field="stat", show_tip_legend=True)
ax.set_axis_off()
`
  },
];
