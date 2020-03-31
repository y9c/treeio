#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for tree sub-module in `treeio` package."""

import unittest

from treeio import show
from treeio import Tree


class TestTreeClass(unittest.TestCase):
    def test_show(self):
        node_a = Tree("Alpha")
        node_b = Tree("Beta", parent=node_a)
        node_c = Tree("Gamma", parent=node_a)
        node_d = Tree("Delta", parent=node_b)
        node_e = Tree("Theta", parent=node_b)

        tree = node_a
        tree_shown = show.tree2ascii(tree, False, False)
        print(tree_shown)
        self.assertEqual(
            tree_shown,
            "               ┌ Delta \n"
            "       ┌─ Beta ┤\n"
            " Alpha ┤       └ Theta \n"
            "       │\n"
            "       └ Gamma ",
        )


if __name__ == "__main__":
    unittest.main()
