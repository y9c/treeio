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

import math
import xml.etree.ElementTree as ET

from .tree import Tree


def read_phyloxml(xml_string: str) -> Tree:
    """Read the first phylogeny of a PhyloXML string into a Tree."""
    trees = read_phyloxmls(xml_string)
    if not trees:
        raise ValueError("no phylogeny found in PhyloXML data")
    return trees[0]


def read_phyloxmls(xml_string: str) -> list[Tree]:
    """Read every ``<phylogeny>`` element into a list of Trees."""
    root = ET.fromstring(xml_string)
    trees = []
    for phylo in root.iter():
        tag = _local(phylo.tag)
        if tag == "phylogeny":
            # The ``<phylogeny>`` element wraps exactly one real ``<clade>``;
            # parse that child directly rather than the wrapper element itself,
            # otherwise a spurious anonymous root node would be injected.
            clade = next((c for c in phylo if _local(c.tag) == "clade"), None)
            if clade is None:
                continue
            t = _parse_clade(clade)
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
            # child is freshly built (no parent, never a duplicate), so attach
            # directly to avoid ``append_child``'s per-child linear scan on wide nodes.
            child_node = _parse_clade(child)
            node._children.append(child_node)
            child_node._parent = node
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
    tree: Tree,
    rooted: bool = True,
    name: str | None = None,
    properties: list[str] | None = None,
) -> str:
    """Write a Tree as a PhyloXML document string."""
    phylo = ET.Element("phylogeny")
    phylo.set("rooted", "true" if rooted else "false")
    if name:
        phylo.set("name", name)
    phylo.append(_clade_elem(tree, properties))

    root = ET.Element("phyloxml")
    root.append(phylo)
    # ``ET.indent`` recurses in the stdlib, so it would overflow Python's
    # recursion limit on a very deep tree.  Use the same formatting produced by
    # ``ET.indent`` but computed iteratively.
    _indent(root)
    return ET.tostring(root, encoding="unicode") + "\n"


def _clade_elem(node: Tree, properties: list[str]) -> ET.Element:
    """Build a ``<clade>`` element for ``node`` iteratively.

    An explicit ``(node, parent_element)`` stack is used instead of recursion so
    a very deep tree does not exhaust Python's recursion limit.
    """

    def build(n: Tree) -> ET.Element:
        clade = ET.Element("clade")
        if n.name:
            name = ET.SubElement(clade, "name")
            name.text = n.name
        if n.branch_length is not None:
            bl = ET.SubElement(clade, "branch_length")
            bl.text = _fmt(n.branch_length)
        if n.support is not None:
            conf = ET.SubElement(clade, "confidence")
            conf.set("type", "support")
            conf.text = _fmt(n.support)
        if properties:
            for key in properties:
                value = n.get_data(key)
                if value is not None:
                    prop = ET.SubElement(clade, "property")
                    prop.set("ref", key)
                    prop.text = str(value)
        return clade

    build_root = build(node)
    stack: list[tuple[Tree, ET.Element]] = []
    # push reversed so a LIFO pop appends children left-to-right
    for c in reversed(node._children):
        stack.append((c, build_root))
    while stack:
        n, parent_elem = stack.pop()
        elem = build(n)
        parent_elem.append(elem)
        for c in reversed(n._children):
            stack.append((c, elem))
    return build_root


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _indent(elem: ET.Element, space: str = "  ") -> None:
    """Insert newlines/indentation like :func:`xml.etree.ElementTree.indent`.

    The stdlib ``indent`` recurses per level, so a very deep tree would raise
    ``RecursionError``.  This reproduces the same output iteratively.
    """
    if len(elem) == 0:
        return
    indentations = ["\n"]
    stack = [(elem, 0)]
    while stack:
        el, lvl = stack.pop()
        child_level = lvl + 1
        try:
            ci = indentations[child_level]
        except IndexError:
            ci = indentations[lvl] + space
            indentations.append(ci)
        if not el.text or not el.text.strip():
            el.text = ci
        children = list(el)
        for ch in children:
            if not ch.tail or not ch.tail.strip():
                ch.tail = ci
            if len(ch):
                stack.append((ch, child_level))
        # dedent the last child to the current level (mirrors stdlib behaviour)
        if children and (children[-1].tail is None or not children[-1].tail.strip()):
            children[-1].tail = indentations[lvl]


def _to_float(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return text


def _fmt(x) -> str:
    f = float(x)
    if not math.isfinite(f):
        return str(f)
    return str(int(f)) if f == int(f) else repr(f)


__all__ = ["read_phyloxml", "read_phyloxmls", "write_phyloxml"]
