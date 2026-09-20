# -*- coding: utf-8 -*-

"""
Top-level package for treeio.

A pure-Python library for reading, writing, manipulating and plotting
phylogenetic trees.  The core engine works on the standard library alone;
plotting uses matplotlib (imported lazily on first draw).
"""

__author__ = "Ye Chang"
__email__ = "yech1990@gmail.com"
__version__ = "0.0.0.dev6"

from .tree import Tree, fortify
from .newick import read_newick, read_newicks, write_newick
from .nexus import read_nexus, read_nexuses, write_nexus
from .json_io import read_json, write_json
from .phyloxml import read_phyloxml, read_phyloxmls, write_phyloxml
from .phylip import read_phylip, read_phylips, write_phylip
from .jplace import read_jplace
from .nexml import read_nexml, read_nexmls, write_nexml
from .vendor import read_mrbayes, read_mrbayeses, read_iqtree, read_raxml
from .io import read, read_many, write, detect_format, register_format, unregister_format
from .treeio import convert_format, convert_string
from .plot import (
    render, draw, plot, save, draw_tree, treeplot, TreePlotter, tree_coords, edge_segments,
    tree_theme, add_scalebar, facet_grid, grid_of_trees,
    register_backend, register_layout, unregister_backend, unregister_layout, layouts,
    color_by_value, color_map, named_palette, suggest_figsize,
)
from .anno import (
    attach, attach_to_tips, attach_to_nodes, get_tipdata, get_nodedata,
)


__all__ = [
    "Tree",
    "fortify",
    "read",
    "read_many",
    "write",
    "detect_format",
    "convert_format",
    "convert_string",
    "read_newick",
    "read_newicks",
    "write_newick",
    "read_nexus",
    "read_nexuses",
    "write_nexus",
    "read_json",
    "write_json",
    "read_phyloxml",
    "read_phyloxmls",
    "write_phyloxml",
    "read_phylip",
    "read_phylips",
    "write_phylip",
    "read_jplace",
    "read_nexml",
    "read_nexmls",
    "write_nexml",
    "read_mrbayes",
    "read_mrbayeses",
    "read_iqtree",
    "read_raxml",
    "draw",
    "plot",
    "save",
    "render",
    "color_by_value",
    "color_map",
    "named_palette",
    "suggest_figsize",
    "attach",
    "attach_to_tips",
    "attach_to_nodes",
    "get_tipdata",
    "get_nodedata",
    "draw_tree",
    "treeplot",
    "TreePlotter",
    "tree_coords",
    "edge_segments",
    "tree_theme",
    "add_scalebar",
    "facet_grid",
    "grid_of_trees",
    "register_format",
    "unregister_format",
    "register_backend",
    "register_layout",
    "unregister_backend",
    "unregister_layout",
    "layouts",
]   