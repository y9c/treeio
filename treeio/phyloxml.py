#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-03-27 15:58

"""
Read and write PhyloXML.

PhyloXML is an XML interchange format for phylogenetic trees.  This module
reads the ``clade`` hierarchy into a :class:`Tree`, preserving ``name``,
``branch_length``, ``confidence`` (support) and generic ``<property>``
annotations, and writes a compatible document back out.

Examples
--------
>>> from treeio import read_phyloxml, write_phyloxml
>>> xml = ('<?xml version=\'1.0\'?>'
...        '<phyloxml><phylogeny rooted="true"><clade>'
...        '<clade><name>A</name><branch_length>0.1</branch_length></clade>'
...        '<clade><name>B</name><branch_length>0.2</branch_length></clade>'
...        '</clade></phylogeny></phyloxml>')
>>> tree = read_phyloxml(xml)
>>> tree.tip_names
['A', 'B']
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from typing import List

from .tree import Tree


def read_phyloxml(xml_string: str) -> Tree:
    """Read the first phylogeny of a PhyloXML string into a Tree."""
    trees = read_phyloxmls(xml_string)
    if not trees:
        raise ValueError("no phylogeny found in PhyloXML data")
    return trees[0]


def read_phyloxmls(xml_string: str) -> List[Tree]:
    """Read every ``<phylogeny>`` element into a list of Trees."""
    root = ET.fromstring(xml_string)
    trees = []
    for phylo in root.iter():
        tag = _local(phylo.tag)
        if tag == "phylogeny":
            t = _parse_clade(phylo)
            t._is_rooted = phylo.get("rooted", "false").lower() == "true"
            trees.append(t)
    return trees


def _parse_clade(elem) -> Tree:
    node = Tree()
    for child in elem:
        tag = _local(child.tag)
        if tag == "name":
            node.name = child.text
        elif tag == "branch_length":
            node.branch_length = float(child.text)
        elif tag == "confidence":
            node.support = _to_float(child.text)
        elif tag == "clade":
            node.append_child(_parse_clade(child))
        elif tag == "property":
            key = child.get("ref") or child.get("name")
            if key:
                setattr(node, key, _cast_scalar(child.text))
    return node


def _cast_scalar(text):
    """Best-effort scalar cast for a property value (numbers become numbers)."""
    if text is None:
        return None
    t = text.strip()
    if t == "":
        return text
    try:
        return float(t)
    except (ValueError, TypeError):
        pass
    if t.lower() == "true":
        return True
    if t.lower() == "false":
        return False
    return text


def write_phyloxml(
    tree: Tree, rooted: bool = True, name: str = None, properties: List[str] = None
) -> str:
    """Write a Tree as a PhyloXML document string."""
    phylo = ET.Element("phylogeny")
    phylo.set("rooted", "true" if rooted else "false")
    if name:
        phylo.set("name", name)
    phylo.append(_clade_elem(tree, properties))

    root = ET.Element("phyloxml")
    root.append(phylo)
    ET.indent(root)
    return ET.tostring(root, encoding="unicode") + "\n"


def _clade_elem(node: Tree, properties: List[str]) -> ET.Element:
    clade = ET.Element("clade")
    if node.name:
        name = ET.SubElement(clade, "name")
        name.text = node.name
    if node.branch_length is not None:
        bl = ET.SubElement(clade, "branch_length")
        bl.text = _fmt(node.branch_length)
    if node.support is not None:
        conf = ET.SubElement(clade, "confidence")
        conf.set("type", "support")
        conf.text = _fmt(node.support)
    if properties:
        for key in properties:
            value = node.get_data(key)
            if value is not None:
                prop = ET.SubElement(clade, "property")
                prop.set("ref", key)
                prop.text = str(value)
    for c in node.children:
        clade.append(_clade_elem(c, properties))
    return clade


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _to_float(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return text


def _fmt(x) -> str:
    f = float(x)
    return str(int(f)) if f == int(f) else repr(f)


__all__ = ["read_phyloxml", "read_phyloxmls", "write_phyloxml"]
