"""Reproduce the collapsed-clade boxed tree (ggtree style).

Rectangular tree where collapsed clades are shown as coloured rounded-box
labels (``Fungi``, ``metazoa``, ``primates``, ``yeast``, ...), a blue boxed
``N`` at the root, a rotated ``Y`` branch label with support ``100``, and tip
gene names (``ADH1``...).
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from treeio import Tree, add_boxed_labels, box_label, draw_tree, tree_coords


def gene_tree():
    """Build the ADH gene-family tree with collapsed clades."""
    # Fungi clade: four yeast tips, each collapsed under its gene name.
    yeast = [Tree(f"ADH{i}", branch_length=1.0) for i in range(1, 5)]
    fungi = Tree("Fungi", branch_length=1.0)
    fungi.extend_children(yeast)

    # primates -> human clade (two tips) on a branch labelled "Y", support 100.
    human = [Tree(f"hADH{i}", branch_length=1.0) for i in (1, 2)]
    ynode = Tree("Y", branch_length=0.6, support=100)
    ynode.extend_children(human)
    primates = Tree("primates", branch_length=1.0)
    primates.extend_children([ynode])

    insect = Tree("ADHX", branch_length=1.0)
    nematode = Tree("ADHY", branch_length=1.0)
    metazoa = Tree("metazoa", branch_length=1.0)
    metazoa.extend_children([primates, insect, nematode])

    root = Tree("N")
    root.extend_children([fungi, metazoa])
    return root


def main():
    t = gene_tree()
    coords = tree_coords(t, layout="rectangular")

    fig, ax = plt.subplots(figsize=(9, 6))
    draw_tree(
        t,
        ax=ax,
        layout="rectangular",
        tip_labels=False,
        tip_points=False,
        edge_width=1.4,
    )

    GREEN = "#a7d7a7"
    # Boxed collapsed-clade labels drawn at each named node (green).
    boxed = {
        "Fungi": "Fungi",
        "ADH1": "yeast",
        "ADH2": "yeast",
        "ADH3": "yeast",
        "ADH4": "yeast",
        "primates": "primates",
        "hADH1": "human",
        "hADH2": "human",
        "ADHX": "insect",
        "ADHY": "nematode",
    }
    add_boxed_labels(t, boxed, ax=ax, layout="rectangular", fill=GREEN, fontsize=8)

    # Blue boxed "N" at the root and at the metazoa node (metazoa is blue, not green).
    for name in ("N", "metazoa"):
        x, y = coords[t.get_node_by_label(name)]
        box_label(ax, x, y, "N", fill="#5b8bd0", textcolor="white", fontsize=8)

    # Rotated "Y" branch label + support value at the human clade root.
    x, y = coords[t.get_node_by_label("Y")]
    box_label(
        ax,
        x,
        y,
        "Y",
        fill="#5b8bd0",
        textcolor="white",
        fontsize=8,
        rotation=90,
        round=False,
    )
    ax.text(x + 0.4, y, "100", fontsize=8, color="#333333", ha="left", va="center")

    # Gene names as plain tip text (right of each leaf box, matching the ref).
    for tip in t.get_tips():
        x, y = coords[tip]
        ax.text(
            x + 1.3, y, tip.name, fontsize=8, color="#333333", ha="left", va="center"
        )

    ax.set_axis_off()
    fig.savefig(
        "/home/yec/Coding/treeio/examples/gallery/07_boxed.png",
        dpi=120,
        bbox_inches="tight",
    )
    plt.close(fig)
    print("wrote examples/gallery/07_boxed.png")


if __name__ == "__main__":
    main()
