#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-03-27 15:58

"""
Read and write newick format.
"""

from .tree import Tree


def read_newick(nwk_string: str) -> Tree:
    """Read newick file into Tree object."""
    return Tree()


def write_newick(tree: Tree) -> str:
    """Write Tree object into string in newick format."""
    nwk_string = ""
    return nwk_string
