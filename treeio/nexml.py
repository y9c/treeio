#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-03-27 15:58

"""Read and write NeXML.

NeXML is the community XML interchange standard for phylogenetic data.  This
module reads the ``otus``/``trees`` (``nex:Tree``) block into :class:`Tree`
objects, mapping ``otu`` labels to tip names, ``edge`` lengths to branch
lengths and ``meta`` annotations to per-node attributes, and writes a
compatible document back out.

Examples
--------
>>> from treeio import read_nexml
>>> tree = read_nexml('<nexml xmlns="http://www.nexml.org/2009">'
...                   '<otus id="o"><otu id="t1" label="A"/><otu id="t2" label="B"/></otus>'
...                   '<trees id="ts"><tree id="tr"><node id="n1" otu="t1"/>'
...                   '<node id="n2" otu="t2"/><node id="n0"/>'
...                   '<root id="n0"/><edge id="e1" source="n0" target="n1"/>'
...                   '</tree></trees></nexml>')
>>> tree.tip_names
['A', 'B']
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from typing import Dict, List

from .tree import Tree

__all__ = ["read_nexml", "read_nexmls", "write_nexml"]

_NS = "http://www.nexml.org/2009"


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _to_scalar(text):
    """Best-effort scalar cast for a meta value."""
    if text is None:
        return None
    t = str(text).strip()
    try:
        return float(t)
    except (ValueError, TypeError):
        pass
    if t.lower() == "true":
        return True
    if t.lower() == "false":
        return False
    return text


def _fmt(x):
    f = float(x)
    return str(int(f)) if f == int(f) else repr(f)


def _label(tag: str) -> str:
    return tag if tag.startswith("{") else "{" + _NS + "}" + tag


def _find_all(elem, name: str):
    return elem.iter(_label(name))


def read_nexml(xml_string: str) -> Tree:
    """Read the first NeXML tree into a :class:`Tree`."""
    trees = read_nexmls(xml_string)
    if not trees:
        raise ValueError("no tree found in NeXML data")
    return trees[0]


def read_nexmls(xml_string: str) -> List[Tree]:
    """Read every NeXML ``tree`` into a list of :class:`Tree`."""
    root = ET.fromstring(xml_string)
    otus: Dict[str, str] = {}
    for otu in _find_all(root, "otu"):
        otus[otu.get("id")] = otu.get("label") or otu.get("id")
    trees = []
    for tree in _find_all(root, "tree"):
        t = _parse_tree(tree, otus)
        if t is not None:
            trees.append(t)
    return trees


def _parse_tree(tree_elem, otus):
    # map node id -> Tree, and edges (source -> [(target, length)])
    nodes: Dict[str, Tree] = {}
    edges = []
    root_id = None
    for ed in _find_all(tree_elem, "edge"):
        edges.append((ed.get("source"), ed.get("target"), ed.get("length")))
    for nd in _find_all(tree_elem, "node"):
        nid = nd.get("id")
        otu = nd.get("otu")
        node = Tree(name=(otus.get(otu, "unknown") if otu else "unknown"))
        for meta in nd.iter(_label("meta")):
            key = meta.get("property")
            if key:
                setattr(node, key, _to_scalar(meta.get("content", "")))
        if _find_all(nd, "branch_length") or nd.get("branch_length"):
            bl = list(_find_all(nd, "branch_length"))
            if bl:
                node.branch_length = float(bl[0].text)
        nodes[nid] = node
    # attach children
    for source, target, length in edges:
        if source in nodes and target in nodes:
            child = nodes[target]
            if length is not None:
                child.branch_length = float(length)
            nodes[source].append_child(child)
    # find the root (node referenced by <root/>)
    for rt in _find_all(tree_elem, "root"):
        root_id = rt.get("id")
    if root_id is None:
        # fall back to the single node that is not a target of any edge
        targets = {t for _, t, _ in edges}
        for nid, node in nodes.items():
            if nid not in targets:
                root_id = nid
                break
    if root_id is None:
        return None
    return nodes[root_id]


def write_nexml(tree: Tree, name: str = "Tree") -> str:
    """Write a :class:`Tree` as a NeXML document string."""
    root = ET.Element(_label("nexml"))
    root.set("xmlns", _NS)
    root.set("xmlns:nex", _NS)
    otus = ET.SubElement(root, _label("otus")); otus.set("id", "otus1")
    # assign otu ids to tips
    tips = tree.get_tips()
    otu_ids = {tip: f"otu{i + 1}" for i, tip in enumerate(tips)}
    for tip in tips:
        otu = ET.SubElement(otus, _label("otu"))
        otu.set("id", otu_ids[tip])
        otu.set("label", tip.name)
    trees = ET.SubElement(root, _label("trees")); trees.set("id", "trees1")
    tree_el = ET.SubElement(trees, _label("tree"))
    tree_el.set("id", name)
    nids = {node: f"n{index}" for index, node in enumerate(tree.traverse("preorder"))}
    # nodes: internal first then leaves / record order doesn't matter to readers
    for node in tree.traverse("preorder"):
        nd = ET.SubElement(tree_el, _label("node"))
        nd.set("id", nids[node])
        if node.is_leaf():
            nd.set("otu", otu_ids[node])
        if node.support is not None:
            meta = ET.SubElement(nd, _label("meta"))
            meta.set("property", "support")
            meta.set("content", _fmt(node.support))
    for node in tree.traverse("preorder"):
        if node.parent is None:
            continue
        edge = ET.SubElement(tree_el, _label("edge"))
        edge.set("id", "e" + nids[node])
        edge.set("source", nids[node.parent])
        edge.set("target", nids[node])
        if node.branch_length is not None:
            edge.set("length", repr(node.branch_length))
    rt = ET.SubElement(tree_el, _label("root"))
    rt.set("id", nids[tree])
    ET.indent(root)
    return ET.tostring(root, encoding="unicode") + "\n"
