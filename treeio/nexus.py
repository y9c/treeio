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

from .newick import read_newick
from .tree import Tree

_BLOCK_RE = re.compile(r"BEGIN\s+(\w+);(.*?)END;", re.DOTALL | re.IGNORECASE)
_TREE_RE = re.compile(
    r"TREE\s*\*?\s*([A-Za-z0-9_.-]*)\s*=\s*(.*?);", re.DOTALL | re.IGNORECASE
)
_TRANSLATE_RE = re.compile(
    r"TRANSLATE\b(.*?)(?:;|BEGIN|END)", re.DOTALL | re.IGNORECASE
)
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


def read_nexus(nexus_string: str, tree_name: str | None = None) -> Tree:
    """Read the first tree from a NEXUS string.  See :func:`read_nexuses`."""
    trees = read_nexuses(nexus_string, tree_name=tree_name)
    if not trees:
        raise ValueError("No TREE found in NEXUS input")
    return trees[0]


def read_nexuses(nexus_string: str, tree_name: str | None = None) -> list[Tree]:
    """Read every tree in a NEXUS string, honouring TRANSLATE and capturing
    embodied ``[&key=value]`` node annotations."""
    mapping: dict = {}
    trees: list[Tree] = []

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
    # label -> translate id; a dict makes the per-tip ``_ref_label`` lookup O(1)
    # instead of a linear scan through the translate table on every node.
    trans = {n: str(i + 1) for i, n in enumerate(tip_labels)}
    if trans:
        lines.append(
            "  TRANSLATE "
            + ", ".join(f"{num} {_quote_label(n)}" for n, num in trans.items())
            + ";"
        )
    for i, t in enumerate(trees):
        label = title if len(trees) == 1 else f"{title}_{i + 1}"
        lines.append(f"  TREE {label} = [&R] " + _write_ref(t, trans) + ";")
    lines.append("END;")
    return "\n".join(lines)


def _union_tip_labels(trees: list[Tree]) -> list[str]:
    out, seen = [], set()
    for t in trees:
        for name in t.get_tip_names():
            if name not in seen:
                seen.add(name)
                out.append(name)
    return out


def _write_ref(node: Tree, trans: dict) -> str:
    """Serialize ``node`` as a (translated) newick string, iteratively.

    A shared output list plus an explicit task stack is used instead of the
    previous recursion, so a very deep tree does not exhaust Python's recursion
    limit and the per-level ``"(" + inner + ")"`` string copying that made deep
    trees quadratic is avoided.
    """
    out: list = []
    stack = [("node", node)]
    while stack:
        kind, data = stack.pop()
        if kind == "node":
            n = data
            if not n._children:
                out.append(_ref_label(n.name, trans))
                _write_ref_suffix(out, n)
            else:
                out.append("(")
                children = n._children
                tasks = []
                for i, c in enumerate(children):
                    tasks.append(("node", c))
                    if i != len(children) - 1:
                        tasks.append(("comma", None))
                tasks.append(("post", n))
                stack.extend(reversed(tasks))
        elif kind == "comma":
            out.append(",")
        else:  # "post"
            n = data
            out.append(")")
            if n.name and n.name != "unknown":
                out.append(n.name)
            elif n.support is not None:
                out.append(str(n.support))
            _write_ref_suffix(out, n)
    return "".join(out)


def _write_ref_suffix(out: list, node: Tree) -> None:
    """Append ``[:length]`` and any ``[&key=value]`` block for ``node``."""
    if node.branch_length is not None:
        out.append(":" + repr(node.branch_length))
    s = _write_annotations(node)
    if s:
        out.append(s)


def _write_annotations(node: Tree) -> str:
    """Emit BEAST-style ``[&key=value, ...]`` node annotations (nothing if none)."""
    from .newick import _format_anno

    ann = getattr(node, "_annotations", None)
    if not ann:
        return ""
    parts = []
    for k, v in ann.items():
        if v is None:
            continue
        parts.append(f"{k}={_format_anno(v)}")
    if not parts:
        return ""
    return "[&" + ", ".join(parts) + "]"


def _ref_label(name: str, trans: dict) -> str:
    num = trans.get(name)
    if num is not None:
        return num
    return _quote_label(name)


def _quote_label(name: str) -> str:
    if _NEEDS_QUOTE_RE.search(name) or name.lower() in {"end", "begin", "translate"}:
        return "'" + name.replace("'", "''") + "'"
    return name


__all__ = ["read_nexus", "read_nexuses", "write_nexus"]
