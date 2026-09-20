#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-03-27 15:58

"""Read the Jplace format (pplacer / EPA placement output).

Jplace is a JSON document holding a reference tree plus placements of query
sequences onto its edges.  ``read_jplace`` returns the reference :class:`Tree`
with the placements and fields exposed as attributes on the tree (``tree...``).

Examples
--------
>>> from treeio import read_jplace
>>> tree = read_jplace('{"tree": "(A,B);", "fields": ["like_weight_ratio"], "placements": []}')
>>> tree.tip_names
['A', 'B']
"""

from __future__ import annotations

import json

from .tree import Tree
from .newick import read_newick

__all__ = ["read_jplace"]


def read_jplace(jplace_string: str) -> Tree:
    """Parse a Jplace document and return the reference Tree.

    The parsed placements (a ``list`` keyed by placement) are stored on the
    tree as ``tree.placements`` (a dict with ``fields``, ``placements``,
    ``metadata`` and ``version``), so downstream code can inspect or plot
    them without losing the raw data.
    """
    data = json.loads(jplace_string)
    if not isinstance(data, dict) or "tree" not in data:
        raise ValueError("not a valid Jplace document (missing 'tree' key)")
    tree = read_newick(data["tree"], annotations=True)
    tree.placements = {
        "fields": data.get("fields", []),
        "placements": data.get("placements", []),
        "metadata": data.get("metadata", {}),
        "version": data.get("version", 1),
    }
    return tree
