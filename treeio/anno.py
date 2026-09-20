#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-03-27 15:58

"""
Annotation & data attachment (the ``treeio`` signature feature).

``treeio``'s strength is associating per-tip and per-node data (bootstrap
values, substitution rates, dates, trait values, ...) with a tree.  This
module provides a small data-table attachment API on top of the per-node
attributes that :class:`treeio.Tree` already supports::

    from treeio import read, attach
    tree = read("tree.nwk")
    attach(tree, {"Homo": "human", "Pan": "chimpanzee"}, key="species")
    tree.get_tipdata("species")

It works with a plain ``dict`` (name -> value), a list of ``dict`` rows (as
from a CSV / DNAML), or a two-sequence mapping of names to values.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional

from .tree import Tree


def attach(
    tree: Tree,
    data,
    key: Optional[str] = None,
    match: str = "name",
) -> Tree:
    """Attach external data to the nodes of ``tree``.

    ``data`` may be one of:

    * a ``dict`` ``{label: value}``
    * a list of ``dict`` rows where each row has a label column and value(s)
    * a ``dict`` ``{label: {attr: value}}`` for multi-attribute rows

    ``key`` names the attribute to store a scalar value under (ignored for
    dict rows that supply their own attribute names).  ``match`` selects the
    node attribute to match on (``'name'`` for tips/internal labels).

    Returns ``tree`` for chaining.
    """
    index = tree.label_index()
    if isinstance(data, dict):
        for label, value in data.items():
            node = index.get(str(label))
            if node is None:
                continue
            if isinstance(value, dict):
                for k, v in value.items():
                    node.set_data(k, v)
            elif key is not None:
                node.set_data(key, value)
            else:
                # store under its own name as attribute
                node.set_data(str(label), value)
    elif isinstance(data, (list, tuple)):
        for row in data:
            if not isinstance(row, dict):
                continue
            label = row.get(match) or row.get("name") or row.get("id")
            if label is None:
                continue
            node = index.get(str(label))
            if node is None:
                continue
            for k, v in row.items():
                if k != match and k not in ("name", "id"):
                    node.set_data(k, v)
    else:
        raise TypeError(f"unsupported data type {type(data).__name__}")
    return tree


def attach_to_tips(tree: Tree, data, key: str = None) -> Tree:
    """Attach data only to tip nodes (see :func:`attach`)."""
    return _attach_filtered(tree, data, key, terminal=True)


def attach_to_nodes(tree: Tree, data, key: str = None) -> Tree:
    """Attach data only to internal nodes (see :func:`attach`)."""
    return _attach_filtered(tree, data, key, terminal=False)


def _attach_filtered(tree: Tree, data, key, terminal):
    mapping = {}
    if isinstance(data, dict):
        index = tree.label_index()
        for label, value in data.items():
            node = index.get(str(label))
            if node is not None and node.is_leaf() == terminal:
                mapping[node] = value
    # apply via attach-like semantics on the filtered subset
    from .tree import Tree

    for node, value in mapping.items():
        if isinstance(value, dict):
            for k, v in value.items():
                node.set_data(k, v)
        elif key is not None:
            node.set_data(key, value)
        else:
            node.set_data(node.name, value)
    return tree


def get_tipdata(tree: Tree, key: str, default=None) -> Dict[str, object]:
    """Return ``{tip_name: value}`` for a tip-level annotation."""
    return tree.get_tipdata(key, default)


def get_nodedata(tree: Tree, key: str, default=None) -> Dict[str, object]:
    """Return ``{node_name: value}`` for an internal-node annotation."""
    return tree.get_nodedata(key, default)


__all__ = [
    "attach",
    "attach_to_tips",
    "attach_to_nodes",
    "get_tipdata",
    "get_nodedata",
]
