#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-03-27 15:58

"""
Read and write newick format.

A dependency-free, robust Newick parser/writer.  It understands:

* ``(A:0.1,B:0.2):0.3;``         branch lengths
* ``((A,B)95,C)90;``             support values (numeric internal labels)
* ``'A ' , "B:1"``               single- / double-quoted labels with escapes
* ``[comment]``                  bracketed comments (ignored)
* polytomies, unrooted trees and missing branch lengths (``:`` with no value)

Examples
--------
>>> from treeio import read_newick, write_newick
>>> t = read_newick("((raccoon:19.2,bear:6.8):0.85,dog:25.5);")
>>> t.tip_names
['raccoon', 'bear', 'dog']
>>> round(t.get_tree_length(), 2)
52.35
"""

from __future__ import annotations

import math
from typing import List, Optional

from .tree import Tree

# characters that terminate an unquoted token
_SPECIAL = set(":,();[]'\"")


def _is_number(s: str) -> bool:
    try:
        float(s)
        return True
    except (ValueError, TypeError):
        return False


class _NewickParser:
    """Recursive-descent Newick parser over a string.

    When ``annotations`` is true, embodied ``[&key=value,...]`` comments are
    captured and stored as attributes on the node they immediately follow
    (the BEAST / NEXUS annotation idiom).
    """

    def __init__(self, text: str, annotations: bool = False):
        self.text = text
        self.i = 0
        self.n = len(text)
        self.annotations = annotations

    # -- lexer helpers ------------------------------------------------ --
    def _skip(self) -> None:
        """Skip whitespace and bracketed comments (but not ``[&...]`` if we
        are capturing annotations)."""
        t, n = self.text, self.n
        while self.i < n:
            c = t[self.i]
            if c in " \t\n\r\f\v":
                self.i += 1
            elif c == "[":
                if self.annotations and self.i + 1 < n and t[self.i + 1] == "&":
                    break
                # skip a (possibly nested) bracketed comment
                depth = 1
                self.i += 1
                while self.i < n and depth:
                    if t[self.i] == "[":
                        depth += 1
                    elif t[self.i] == "]":
                        depth -= 1
                    self.i += 1
            else:
                break

    def _peek(self) -> str:
        return self.text[self.i] if self.i < self.n else ""

    def _quoted(self) -> str:
        """Read a single/double quoted label ('' or \"\" escapes an inner quote)."""
        q = self.text[self.i]
        self.i += 1
        out = []
        t, n = self.text, self.n
        while self.i < n:
            c = t[self.i]
            if c == q:
                if self.i + 1 < n and t[self.i + 1] == q:
                    out.append(q)
                    self.i += 2
                    continue
                self.i += 1
                break
            out.append(c)
            self.i += 1
        return "".join(out)

    def _token(self) -> str:
        """Read an unquoted token (label or number), stopping at specials."""
        start = self.i
        t, n = self.text, self.n
        while self.i < n and t[self.i] not in _SPECIAL and t[self.i] not in " \t\n\r\f\v":
            self.i += 1
        return t[start:self.i]

    # -- grammar ------------------------------------------------------- --
    def parse(self) -> Tree:
        self._skip()
        # BEAST-style tree-level markers ([&R], [&U], ...) precede the root
        # node and belong to the tree, not to a node.
        leading = self._consume_leading_annotations()
        tree = self._subtree()
        for ann in leading:
            for key, value in ann.items():
                setattr(tree, key, value)
        self._skip()
        if self._peek() == ";":
            self.i += 1
        if tree.parent is not None:
            tree.isolated()
        return tree

    def _consume_leading_annotations(self) -> list:
        if not self.annotations:
            return []
        out = []
        while True:
            self._skip()
            if self._peek() == "[" and self.i + 1 < self.n and self.text[self.i + 1] == "&":
                out.append(self._read_annotation())
            else:
                break
        return out

    def _subtree(self) -> Tree:
        self._skip()
        if self._peek() == "(":
            self.i += 1
            children = []
            while True:
                self._skip()
                children.append(self._subtree())
                self._skip()
                if self._peek() == ",":
                    self.i += 1
                    continue
                break
            self._skip()
            if self._peek() == ")":
                self.i += 1
            node = Tree()
            for c in children:
                node.append_child(c)
            self._label_branch(node)
            return node

        node = Tree()
        self._label_branch(node)
        return node

    def _label_branch(self, node: Tree) -> None:
        """Parse an optional node label and optional branch length."""
        self._skip()
        label = None
        if self._peek() in ("'", '"'):
            label = self._quoted()
        else:
            tok = self._token()
            if tok:
                label = tok
        if label is not None:
            node.name = label
            if _is_number(label):
                try:
                    node.support = float(label)
                except (ValueError, TypeError):
                    pass

        self._skip()
        if self._peek() == ":":
            self.i += 1
            num = self._token()
            if num:
                try:
                    node.branch_length = float(num)
                except (ValueError, TypeError):
                    raise ValueError(f"invalid branch length {num!r} in newick")

        # embodied BEAST-style annotations, e.g. [&rate=0.5, height=2.0]
        if self.annotations:
            self._capture_annotations(node)

    def _capture_annotations(self, node: Tree) -> None:
        while True:
            self._skip()
            if self._peek() == "[" and self.i + 1 < self.n and self.text[self.i + 1] == "&":
                for key, value in self._read_annotation().items():
                    setattr(node, key, value)
            else:
                return

    def _read_annotation(self) -> dict:
        """Read one annotation block into a dict.

        Supports BEAST-style ``[&key=value,...]`` and NHX-style
        ``[&&NHX:key=value,...]``.
        """
        self.i += 2  # skip [&
        if self.text.startswith("&NHX:", self.i):
            self.i += 5  # skip &NHX:  (we already consumed the first '&')
        body = []
        while self.i < self.n and self.text[self.i] != "]":
            body.append(self.text[self.i])
            self.i += 1
        if self.i < self.n:
            self.i += 1  # skip ]
        out = {}
        for item in self._split_kv("".join(body)):
            key, _, value = item.partition("=")
            if not key:
                continue
            key = key.strip().strip("'")
            out[key] = _cast(value.strip()) if value else True
        return out

    @staticmethod
    def _split_kv(s: str):
        """Split a comma-separated key=value string, respecting quotes & braces."""
        parts, buf, depth, in_q = [], [], 0, False
        for ch in s:
            if ch in "'\"":
                in_q = not in_q
                buf.append(ch)
            elif ch in "({[":
                depth += 1
                buf.append(ch)
            elif ch in ")}]":
                depth -= 1
                buf.append(ch)
            elif ch == "," and depth == 0 and not in_q:
                parts.append("".join(buf))
                buf = []
            else:
                buf.append(ch)
        if buf:
            parts.append("".join(buf))
        return parts


def _cast(value: str):
    """Best-effort scalar cast for an annotation value."""
    if value.startswith("'") and value.endswith("'"):
        return value[1:-1].replace("''", "'")
    if value.startswith('"') and value.endswith('"'):
        return value[1:-1].replace('""', '"')
    try:
        return float(value)
    except (ValueError, TypeError):
        pass
    if value.lower() == "true":
        return True
    if value.lower() == "false":
        return False
    return value


def _format_number(x) -> str:
    """Format a number for newick output, keeping it short but lossless."""
    if isinstance(x, bool):
        x = int(x)
    if isinstance(x, int):
        return str(x)
    f = float(x)
    if not math.isfinite(f):
        raise ValueError(f"cannot format non-finite branch length {x!r}")
    if f == int(f) and abs(f) < 1e16:
        return str(int(f))
    return repr(f)


def _needs_quote(label: str) -> bool:
    return any(c in _SPECIAL or c in " \t\n\r\f\v" for c in label)


def write_newick(
    tree: Tree, include_dist: bool = True, include_support: bool = True,
    annotations: Optional[list] = None,
) -> str:
    """Write Tree object into a string in newick format.

    Parameters
    ----------
    tree : Tree
        The tree (rooted) to serialize.
    include_dist : bool
        Whether to emit branch lengths.
    include_support : bool
        Whether to emit support values as internal-node labels when a node
        has no explicit name.
    annotations : optional list
        Names of node attributes to emit as BEAST-style ``[&key=value]``
        embodied data on each node.
    """
    return _write_node(tree, include_dist, include_support, annotations) + ";"


def _write_node(
    node: Tree, include_dist: bool, include_support: bool, annotations: Optional[list]
) -> str:
    label = _write_label(node, include_support)
    if node.is_leaf():
        s = label
    else:
        child = ",".join(
            _write_node(c, include_dist, include_support, annotations)
            for c in node.children
        )
        s = "(" + child + ")" + label
    if include_dist and node.branch_length is not None:
        s += ":" + _format_number(node.branch_length)
    s += _write_annotations(node, annotations)
    return s


def _write_annotations(node: Tree, annotations: Optional[list]) -> str:
    if not annotations:
        return ""
    parts = []
    for key in annotations:
        value = node.get_data(key)
        if value is None:
            continue
        parts.append(f"{key}={_format_anno(value)}")
    if not parts:
        return ""
    return "[&" + ", ".join(parts) + "]"


def _format_anno(value) -> str:
    if isinstance(value, str):
        return "'" + value.replace("'", "''") + "'"
    if isinstance(value, bool):
        return "true" if value else "false"
    return _format_number(value)


def _write_label(node: Tree, include_support: bool) -> str:
    label = node.name
    if label is None or label in ("", "unknown"):
        label = str(node.support) if include_support and node.support is not None else ""
    if label == "":
        return label
    return f"'{label}'" if _needs_quote(label) else label


def read_newick(nwk_string: str, annotations: bool = False) -> Tree:
    """Read a newick string into a Tree object.

    If ``annotations`` is true, embodied ``[&key=value,...]`` data are parsed
    into per-node attributes (BEAST / NEXUS style).  Only the first tree is
    returned; use :func:`read_newicks` for a list.
    """
    parser = _NewickParser(nwk_string, annotations=annotations)
    return parser.parse()


def read_newicks(nwk_string: str) -> List[Tree]:
    """Read one or more newick trees (separated by ``;``) into a list."""
    trees: List[Tree] = []
    rest = nwk_string
    # pull tree-by-tree using a parser that stops after one top-level element
    idx = 0
    length = len(nwk_string)
    while True:
        parser = _NewickParser(nwk_string)
        parser.i = idx
        parser._skip()
        if parser.i >= length:
            break
        # only start a new tree on a parenthesised group, a lone node label, or
        # a bracketed root marker; anything else is trailing garbage -> stop.
        head = nwk_string[parser.i:parser.i + 1]
        if head not in ("(", "[") and not _is_number(head):
            break
        trees.append(parser.parse())
        idx = parser.i
    return trees


if __name__ == "__main__":
    nwk = "((raccoon:19.19959,bear:6.80041):0.84600," \
          "((sea_lion:11.99700,seal:12.00300):7.52973," \
          "((monkey:100.85930,cat:47.14069):20.59201,weasel:18.87953):2.09460):3.87382,dog:25.46154);"
    t = read_newick(nwk)
    print(t.tip_names)
    print(t.nleaves, t.nnodes, t.is_binary())
    print(f"tree length = {t.get_tree_length():.5f}")
    s = write_newick(t)
    print(s)
    t2 = read_newick(s)
    print("round trip", t2 == t)
    print(t)


__all__ = ["read_newick", "read_newicks", "write_newick"]