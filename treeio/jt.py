#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-03-27 15:58

"""
Read and write json format.
"""

import json
from typing import List

from .tree import Tree


def read_json(
    json_string: str,
    name_key="name",
    child_key="children",
    dist_key="branch_length",
    supp_key="support",
) -> List[Tree]:
    """Return a json object in the format desribed below

    ```json
    {
      "name": "root",
      "branch_length": 0,
      "children": [
        {
          "name": "unknown",
          "branch_length": 0.846,
          "children": [
            { "name": "raccoon", "branch_length": 19.19959 },
            { "name": "bear", "branch_length": 6.80041 }
          ]
        },
        {
          "name": "unknown",
          "branch_length": 3.87382,
          "children": [

          ...

          ]
        },
        { "name": "dog", "branch_length": 25.46154 }
      ]
    }
    ```
    """

    data = json.loads(json_string)
    # assert type(data) is dict, "Only single tree is supported."

    tree_scratch = Tree("scratch")

    def _parse_node(obj, tree_cur):
        """Recursively search for values of key in JSON tree."""
        if isinstance(obj, dict):
            node = Tree(
                name=obj.get(name_key),
                dist=obj.get(dist_key),
                supp=obj.get(supp_key),
            )
            tree_cur.append_child(node)
            if child_key in obj:
                tree_cur = node
                _parse_node(obj[child_key], tree_cur)
        elif isinstance(obj, list):
            for item in obj:
                _parse_node(item, tree_cur)

    _parse_node(data, tree_scratch)

    return [t.isolated() for t in tree_scratch.children]


def write_json(
    trees: List[Tree],
    name_key="name",
    child_key="children",
    dist_key="branch_length",
    supp_key="support",
) -> str:
    """Return a json object in the format desribed below
    """

    def _record_node(node):
        attr_key = ["name", "dist", "supp"]
        attr_values = [name_key, dist_key, supp_key]
        data = {v: getattr(node, k) for k, v in zip(attr_key, attr_values)}
        children = [_record_node(child) for child in node.children]
        if children:
            data[child_key] = children
        return data

    data = _record_node(trees[0])

    json_string = json.dumps(data)
    return json_string


if __name__ == "__main__":
    with open("../data/animals.json") as f:
        TREE = read_json(f.read())
        print(TREE[0])
