"""Tests for tree sub-module in `treeio` package."""

import unittest

from treeio import Tree, show


class TestTreeClass(unittest.TestCase):
    def test_show(self):
        node_a = Tree("Alpha")
        node_b = Tree("Beta", parent=node_a)
        Tree("Gamma", parent=node_a)
        Tree("Delta", parent=node_b)
        Tree("Theta", parent=node_b)

        tree = node_a
        tree_shown = show.tree_to_ascii(tree, False, False)
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
