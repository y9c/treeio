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


if __name__ == "__main__":
    unittest.main()
