#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-03-27 15:58

"""Read and write PHYLIP tree format.

PHYLIP trees are newick bodies, optionally preceded by a header line of two
integers ``<ntax> <ntrees>`` (the number of taxa and trees).  ``read_phylip``
accepts trees with or without the header; ``write_phylip`` emits it by default.

Examples
--------
>>> from treeio import read_phylip, write_phylip
>>> tree = read_phylip("3 1\n((A:0.1,B:0.2):0.3,C:0.4);")
>>> tree.tip_names
['A', 'B', 'C']
"""

from __future__ import annotations

import re
from typing import List

from .tree import Tree
from .newick import read_newicks, write_newick

# ``<ntax> <ntrees>`` optional header
_HEADER_RE = re.compile(r"^\s*\d+\s+\d+\s*$")

__all__ = ["read_phylip", "read_phylips", "write_phylip"]


def _strip_header(text: str) -> str:
    """Remove the optional ``<ntax> <ntrees>`` leading header line."""
    lines = text.split("\n", 1)
    head = lines[0].strip()
    if len(lines) == 2 and _HEADER_RE.match(head):
        return lines[1]
    return text


def read_phylip(phylip_string: str) -> Tree:
    """Read the first tree from a PHYLIP string."""
    trees = read_phylips(phylip_string)
    if not trees:
        raise ValueError("no tree found in PHYLIP input")
    return trees[0]


def read_phylips(phylip_string: str) -> List[Tree]:
    """Read every tree in a PHYLIP string."""
    body = _strip_header(phylip_string)
    return read_newicks(body)


def write_phylip(tree, header: bool = True) -> str:
    """Write a tree as PHYLIP text (optionally with a ``<ntax> <ntrees>`` header)."""
    if isinstance(tree, Tree):
        trees = [tree]
    else:
        trees = list(tree)
    nwk = "\n".join(write_newick(t) for t in trees)
    ntax = tree.nleaves if isinstance(tree, Tree) else trees[0].nleaves
    if header:
        return f"{ntax} {len(trees)}\n{nwk}\n"
    return nwk + "\n"
