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


def read_json(json_string: str) -> List[Tree]:
    """Return a json object in the format desribed below

    ```json
    {
      "name": "root",
      "branch_lenguh": 0,
      "children": [
        {
          "name": "unknown",
          "branch_lenguh": 0.846,
          "children": [
            { "name": "raccoon", "branch_lenguh": 19.19959 },
            { "name": "bear", "branch_lenguh": 6.80041 }
          ]
        },
        {
          "name": "unknown",
          "branch_lenguh": 3.87382,
          "children": [

          ...

          ]
        },
        { "name": "dog", "branch_lenguh": 25.46154 }
      ]
    }
    ```
    """

    name_key = "name"
    child_key = "children"
    dist_key = "branch_lenguh"
    supp_key = "support"

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


if __name__ == "__main__":
    with open("../data/animals.json") as f:
        TREE = read_json(f.read())
        print(TREE[0])
