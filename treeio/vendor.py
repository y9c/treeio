#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-03-27 15:58

"""Readers for popular phylogenetic-software output formats.

These are thin adapters over the core parsers, tuned for how each tool writes
its trees:

* ``read_mrbayes`` -- MrBayes NEXUS (a ``begin trees`` block, plus the tool's
  ``begin mrbayes`` block is captured as tree metadata).
* ``read_iqtree``  -- iQ-TREE ``.treefile`` (newick with node supports).
* ``read_raxml``   -- RAxML output trees (newick with node supports).
"""

from __future__ import annotations

import re

from .tree import Tree
from .nexus import read_nexuses
from .newick import read_newick

__all__ = ["read_mrbayes", "read_mrbayeses", "read_iqtree", "read_raxml"]


def read_mrbayes(text: str) -> Tree:
    """Read the first tree from a MrBayes NEXUS file."""
    trees = read_mrbayeses(text)
    if not trees:
        raise ValueError("no tree found in MrBayes output")
    return trees[0]


def read_mrbayeses(text: str) -> list:
    """Read every tree from a MrBayes NEXUS file."""
    trees = read_nexuses(text)
    meta = _mrbayes_block(text)
    for t in trees:
        if meta:
            t.set_data("mrbayes", meta)
    return trees


def _mrbayes_block(text: str) -> dict:
    """Capture key MrBayes settings from the ``begin mrbayes`` block."""
    m = re.search(r"BEGIN\s+MRBAYES;(.*?)END;", text, re.S | re.I)
    if not m:
        return {}
    body = m.group(1)
    out = {}
    for key in ("lset", "prset", "mcmc", "unlink", "outgroup"):
        km = re.search(r"\b" + key + r"\s+([^;]+);", body, re.I)
        if km:
            out[key] = km.group(1).strip()
    return out


def read_iqtree(treefile_text: str) -> Tree:
    """Read the tree from an iQ-TREE ``.treefile`` (newick)."""
    return read_newick(treefile_text, annotations=True)


def read_raxml(tree_text: str) -> Tree:
    """Read a RAxML output tree (newick with node supports)."""
    return read_newick(tree_text)
