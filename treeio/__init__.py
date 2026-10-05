"""
Top-level package for treeio.

A pure-Python library for reading, writing, manipulating and plotting
phylogenetic trees.  The core engine works on the standard library alone;
plotting uses matplotlib (imported lazily on first draw).
"""

__author__ = "Ye Chang"
__email__ = "yech1990@gmail.com"
__version__ = "0.0.0.dev12"

from .anno import (
    attach,
    attach_to_nodes,
    attach_to_tips,
    get_nodedata,
    get_tipdata,
)
from .io import (
    detect_format,
    read,
    read_many,
    register_format,
    unregister_format,
    write,
)
from .jplace import read_jplace
from .json_io import read_json, write_json
from .newick import read_newick, read_newicks, write_newick
from .nexml import read_nexml, read_nexmls, write_nexml
from .nexus import read_nexus, read_nexuses, write_nexus
from .phylip import read_phylip, read_phylips, write_phylip
from .phylogeny import build_tree, distance_matrix, neighbor_joining
from .phyloxml import read_phyloxml, read_phyloxmls, write_phyloxml
from .plot import (
    TreePlotter,
    add_boxed_labels,
    add_ring,
    add_rings,
    add_scalebar,
    box_label,
    color_by_value,
    color_map,
    draw,
    draw_tree,
    edge_segments,
    facet_grid,
    gheatmap,
    grid_of_trees,
    highlight_clade,
    layouts,
    named_palette,
    plot,
    register_backend,
    register_layout,
    render,
    save,
    suggest_figsize,
    tree_coords,
    tree_theme,
    treeplot,
    unregister_backend,
    unregister_layout,
)
from .tree import Tree, fortify
from .treeio import convert_format, convert_string
from .vendor import read_iqtree, read_mrbayes, read_mrbayeses, read_raxml

__all__ = [
    "Tree",
    "TreePlotter",
    "add_boxed_labels",
    "add_ring",
    "add_rings",
    "add_scalebar",
    "attach",
    "attach_to_nodes",
    "attach_to_tips",
    "box_label",
    "build_tree",
    "color_by_value",
    "color_map",
    "convert_format",
    "convert_string",
    "detect_format",
    "distance_matrix",
    "draw",
    "draw_tree",
    "edge_segments",
    "facet_grid",
    "fortify",
    "get_nodedata",
    "get_tipdata",
    "gheatmap",
    "grid_of_trees",
    "highlight_clade",
    "layouts",
    "named_palette",
    "neighbor_joining",
    "plot",
    "read",
    "read_iqtree",
    "read_jplace",
    "read_json",
    "read_many",
    "read_mrbayes",
    "read_mrbayeses",
    "read_newick",
    "read_newicks",
    "read_nexml",
    "read_nexmls",
    "read_nexus",
    "read_nexuses",
    "read_phylip",
    "read_phylips",
    "read_phyloxml",
    "read_phyloxmls",
    "read_raxml",
    "register_backend",
    "register_format",
    "register_layout",
    "render",
    "save",
    "suggest_figsize",
    "tree_coords",
    "tree_theme",
    "treeplot",
    "unregister_backend",
    "unregister_format",
    "unregister_layout",
    "write",
    "write_json",
    "write_newick",
    "write_nexml",
    "write_nexus",
    "write_phylip",
    "write_phyloxml",
]
