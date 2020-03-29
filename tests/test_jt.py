#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `treeio` package."""


from treeio import jt

with open("../data/animals.json") as f:
    tree = jt.read_json(f.read())
    print(tree)
