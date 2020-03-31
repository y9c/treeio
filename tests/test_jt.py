#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `treeio` package."""

import unittest

from treeio import Tree
from treeio import jt


class TestJsonIO(unittest.TestCase):
    def test_read_json(self):
        with open("./data/animals.json") as file_json:
            trees = jt.read_json(file_json.read())
            for tree in trees:
                self.assertIsInstance(tree, Tree)

    def test_write_json(self):
        node_a = Tree("Alpha")
        node_b = Tree("Beta")
        node_c = Tree("Gamma")
        node_d = Tree("Delta")
        node_e = Tree("Theta")
        node_a.children = [node_b, node_c]
        node_b.children = [node_d, node_e]
        tree = node_a
        json_string = jt.write_json([tree])
        self.assertEqual(
            json_string,
            '{"name": "Alpha", "branch_length": null, "support": null, "children": [{"name": "Beta", "branch_length": null, "support": null, "children": [{"name": "Delta", "branch_length": null, "support": null}, {"name": "Theta", "branch_length": null, "support": null}]}, {"name": "Gamma", "branch_length": null, "support": null}]}',
        )


if __name__ == "__main__":
    unittest.main()
