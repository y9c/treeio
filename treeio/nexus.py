#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-03-27 15:58

"""
Read and write NEXUS format.

Supports the ``TREES`` block including a ``TRANSLATE`` table, one or more
``TREE`` definitions, and BEAST-style embodied annotations (``[&height=...,
rate=...]``) which are captured as per-node attributes.  Writing produces a
standard ``#NEXUS`` file readable by PAUP*, Mesquite, BEAST, RAxML and R
``ape``.

Examples
--------
>>> from treeio import read_nexus
>>> tree = read_nexus(
...     "#NEXUS\nBEGIN TREES;\n  TRANSLATE 1 A, 2 B, 3 C;\n  TREE t = (1,2);\nEND;\n")
>>> tree.tip_names
['A', 'B', 'C']
"""

from __future__ import annotations

import re
from typing import List, Optional

from .tree import Tree
from .newick import read_newick, write_newick

_BLOCK_RE = re.compile(r"BEGIN\s+(\w+);(.*?)END;", re.S | re.I)
_TREE_RE = re.compile(r"TREE\s*\*?\s*([A-Za-z0-9_.-]*)\s*=\s*(.*?);", re.S | re.I)
_TRANSLATE_RE = re.compile(r"TRANSLATE\b(.*?)(?:;|BEGIN|END)", re.S | re.I)
_ENTRY_RE = re.compile(r"(\w+)\s*([^,;]+?)(?=[,;]|$)")
_NEEDS_QUOTE_RE = re.compile(r"[\s,;()\[\]{}\\s:'\"]")
_ID_RE = re.compile(r"[\w.-]+")


def _parse_translate(block: str) -> dict:
    mapping = {}
    for m in _ENTRY_RE.finditer(block):
        num, label = m.group(1), m.group(2).strip()
        if len(label) >= 2 and label[0] == label[-1] and label[0] in "'\"":
            label = label[1:-1].replace(label[0] * 2, label[0])
        mapping[num] = label
    return mapping


def _substitute_translate(nwk: str, mapping: dict) -> str:
    """Replace integer translate ids with quoted taxon labels.

    Only whole-number tokens that are not part of a branch length (i.e. not
    preceded by ``:`` or ``.`` or a digit) are substituted, so lengths like
    ``19.19959`` are left untouched.
    """
    def repl(match):
        tok = match.group(0)
        label = mapping.get(tok)
        if label is not None:
            return "'" + label.replace("'", "''") + "'"
        return tok

    return re.sub(r"(?<![\w.\d:])\d+", repl, nwk)


def read_nexus(nexus_string: str, tree_name: Optional[str] = None) -> Tree:
    """Read the first tree from a NEXUS string.  See :func:`read_nexuses`."""
    trees = read_nexuses(nexus_string, tree_name=tree_name)
    if not trees:
        raise ValueError("No TREE found in NEXUS input")
    return trees[0]


def read_nexuses(nexus_string: str, tree_name: Optional[str] = None) -> List[Tree]:
    """Read every tree in a NEXUS string, honouring TRANSLATE and capturing
    embodied ``[&key=value]`` node annotations."""
    mapping: dict = {}
    trees: List[Tree] = []

    for block in _BLOCK_RE.finditer(nexus_string):
        if block.group(1).upper() != "TREES":
            continue
        body = block.group(2)
        tm = _TRANSLATE_RE.search(body)
        if tm:
            mapping = _parse_translate(tm.group(1))
        for tree_match in _TREE_RE.finditer(body):
            name = tree_match.group(1).strip()
            if tree_name and name and name != tree_name:
                continue
            body_str = tree_match.group(2).strip()
            nwk = _substitute_translate(body_str, mapping)
            t = read_newick(nwk, annotations=True)
            if name:
                # the TREE name is a file-level title, not a node label
                t.set_data("tree_title", name)
            trees.append(t)
    return trees


def write_nexus(trees, title: str = "Tree") -> str:
    """Write one or more trees as a NEXUS file with a TAXA + TREES block.

    ``trees`` may be a single :class:`Tree` or a list of them.
    """
    if isinstance(trees, Tree):
        trees = [trees]
    tip_labels = _union_tip_labels(trees)

    lines = ["#NEXUS"]
    lines.append("BEGIN TAXA;")
    lines.append(f"  DIMENSIONS NTAX={len(tip_labels)};")
    lines.append("  TAXLABELS " + " ".join(_quote_label(n) for n in tip_labels) + ";")
    lines.append("END;")
    lines.append("BEGIN TREES;")
    trans = [(str(i + 1), n) for i, n in enumerate(tip_labels)]
    if trans:
        lines.append(
            "  TRANSLATE " + ", ".join(f"{num} {_quote_label(n)}" for num, n in trans) + ";"
        )
    for i, t in enumerate(trees):
        label = title if len(trees) == 1 else f"{title}_{i + 1}"
        lines.append(f"  TREE {label} = [&R] " + _write_ref(t, trans) + ";")
    lines.append("END;")
    return "\n".join(lines)


def _union_tip_labels(trees: List[Tree]) -> List[str]:
    out, seen = [], set()
    for t in trees:
        for name in t.get_tip_names():
            if name not in seen:
                seen.add(name)
                out.append(name)
    return out


def _write_ref(node: Tree, trans: list) -> str:
    if node.is_leaf():
        s = _ref_label(node.name, trans)
    else:
        inner = ",".join(_write_ref(c, trans) for c in node.children)
        s = "(" + inner + ")"
        if node.name and node.name != "unknown":
            s += node.name
        elif node.support is not None:
            s += str(node.support)
    if node.branch_length is not None:
        s += ":" + repr(node.branch_length)
    s += _write_annotations(node)
    return s


def _write_annotations(node: Tree) -> str:
    """Emit BEAST-style ``[&key=value, ...]`` node annotations (nothing if none)."""
    from .newick import _format_anno

    ann = getattr(node, "_annotations", None)
    if not ann:
        return ""
    parts = [f"{k}={_format_anno(v)}" for k, v in ann.items()]
    if not parts:
        return ""
    return "[&" + ", ".join(parts) + "]"


def _ref_label(name: str, trans: list) -> str:
    for num, label in trans:
        if label == name:
            return num
    return _quote_label(name)


def _quote_label(name: str) -> str:
    if _NEEDS_QUOTE_RE.search(name) or name.lower() in {"end", "begin", "translate"}:
        return "'" + name.replace("'", "''") + "'"
    return name


__all__ = ["read_nexus", "read_nexuses", "write_nexus"]
