#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `treeio` package."""

import unittest

from treeio import show
from treeio import Tree
from treeio import jt


class TestTreeShow(unittest.TestCase):
    def test_show(self):
        node_a = Tree("Alpha")
        node_b = Tree("Beta")
        node_c = Tree("Gamma")
        node_d = Tree("Delta")
        node_e = Tree("Theta")

        tree = node_a
        node_a.children = [node_b, node_c]
        node_b.children = [node_d, node_e]

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
