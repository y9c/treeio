#!/usr/bin/env python
"""Minimal working example: read the bundled data and render a figure."""

import matplotlib
matplotlib.use("Agg")

from treeio import read, draw_tree, add_scalebar, tree_coords, color_by_value, attach

# 1) read the bundled sample tree (auto-detected format)
tree = read("data/animals.nwk")
print("tips:", tree.tip_names)
print("nodes:", tree.nnodes, " leaves:", tree.nleaves)

# 2) attach a trait, and colour tip points by it (precomputed -> no extra legend)
attach(tree, {"raccoon": "mammal", "bear": "mammal", "monkey": "primate",
              "dog": "primate", "weasel": "primate", "cat": "primate",
              "sea_lion": "mammal", "seal": "mammal"}, key="class")
cmap = color_by_value(tree.get_tipdata("class"))

# 3) draw a clean publication figure: branches coloured by length + one colorbar
ax = draw_tree(tree, layout="rectangular", tip_labels=True, tip_points=True,
               tip_colors=cmap, branch_color_field="branch_length",
               show_colorbar=True, colorbar_label="branch length", theme="clean")
add_scalebar(ax, tree_coords(tree), label="substitutions")
ax.figure.savefig("examples/animals.png", dpi=200, bbox_inches="tight")

print("wrote examples/animals.png")
