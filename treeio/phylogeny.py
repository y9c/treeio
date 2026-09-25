#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-03-27 15:58

"""
Distance-based tree building (aligned sequences -> distance matrix -> NJ tree).

This closes the one capability gap vs. Biopython: building a tree from an
*alignment*, rather than only parsing serialised trees.  It is pure-Python
(stdlib only), so Biopython is **not** required -- but sequences may be given as
raw strings, ``(name, sequence)`` pairs, a ``{name: sequence}`` dict, or
Biopython ``SeqRecord`` objects (detected via a ``.seq`` attribute).

Workflow
--------
1. ``mafft`` (or any aligner) produces an alignment -- an external tool.
2. ``distance_matrix(sequences, model=...)`` ---- pairwise distances.
3. ``neighbor_joining(labels, matrix)`` ---- the NJ tree.

Examples
--------
>>> from treeio import build_tree
>>> seqs = {"A": "ACGT", "B": "ACGA", "C": "AGGT"}
>>> tree = build_tree(seqs, model="p")
>>> tree.tip_names
['A', 'B', 'C']
"""

from __future__ import annotations

import math
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

from .tree import Tree

__all__ = ["distance_matrix", "neighbor_joining", "build_tree"]


# ---------------------------------------------------------------------- --
# sequence handling
# ---------------------------------------------------------------------- --
def _seq_text(seq) -> str:
    """Best-effort sequence text from a raw str, dict value or SeqRecord."""
    if hasattr(seq, "seq"):
        return str(seq.seq)
    return str(seq)


def _seq_label(seq, index: int) -> str:
    if hasattr(seq, "id"):
        return str(seq.id)
    if hasattr(seq, "name"):
        return str(seq.name)
    return str(index)


def _normalise(sequences) -> Tuple[List[str], List[str]]:
    """Return ``(labels, texts)`` for many input forms."""
    if isinstance(sequences, dict):
        items = list(sequences.items())
        return [str(k) for k, _ in items], [_seq_text(v) for _, v in items]
    items = list(sequences)
    if items and isinstance(items[0], (tuple, list)) and len(items[0]) == 2 and isinstance(items[0][1], str):
        return [str(k) for k, _ in items], [_seq_text(v) for _, v in items]
    return [str(_seq_label(s, i)) for i, s in enumerate(items)], [_seq_text(s) for s in items]


# ---------------------------------------------------------------------- --
# distance models
# ---------------------------------------------------------------------- --
def _p_distance(a: str, b: str) -> float:
    """Raw p-distance: fraction of comparable, differing sites.

    Gaps / unknown '?' are ignored; length is taken from the shorter sequence.
    """
    total = diff = 0
    for x, y in zip(a, b):
        if x in "-?" or y in "-?":
            continue
        total += 1
        if x != y:
            diff += 1
    return (diff / total) if total else 0.0


def _jukes_cantor(p: float) -> float:
    """Jukes-Cantor corrected distance from a p-distance."""
    if p <= 0.0:
        return 0.0
    if p >= 0.75:
        return float("inf")  # saturated
    return -0.75 * math.log(1.0 - 4.0 * p / 3.0)


_MODELS = {
    "p": _p_distance,
    "p-distance": _p_distance,
    "raw": _p_distance,
    "jc": lambda a, b: _jukes_cantor(_p_distance(a, b)),
    "jc69": lambda a, b: _jukes_cantor(_p_distance(a, b)),
    "jukes-cantor": lambda a, b: _jukes_cantor(_p_distance(a, b)),
}


def distance_matrix(sequences, model: str = "p") -> Tuple[List[str], List[List[float]]]:
    """Return ``(labels, matrix)`` of pairwise distances.

    ``model`` may be ``"p"`` (raw p-distance, default) or ``"jc"`` / ``"jc69"``
    (Jukes-Cantor corrected).  All sequences must be aligned (equal length).
    """
    if model not in _MODELS:
        raise ValueError(f"unknown distance model {model!r}; choose from {sorted(set(_MODELS))}")
    labels, texts = _normalise(sequences)
    n = len(texts)
    func = _MODELS[model]
    matrix = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = func(texts[i], texts[j])
            matrix[i][j] = matrix[j][i] = d
    return labels, matrix


# ---------------------------------------------------------------------- --
# neighbour joining (pure Python)
# ---------------------------------------------------------------------- --
def neighbor_joining(labels: Sequence[str], matrix: Sequence[Sequence[float]]) -> Tree:
    """Build a Tree from a pairwise distance matrix using neighbour joining.

    ``labels`` are the taxon names, ``matrix`` is a symmetric distance matrix.
    Returns a (rooted) :class:`Tree` whose internal nodes are unnamed and whose
    branch lengths come from the NJ limb-length estimates.
    """
    n = len(labels)
    if len(matrix) != n or any(len(row) != n for row in matrix):
        raise ValueError("distance matrix must be n x n for the given labels")
    if n == 1:
        return Tree(name=labels[0])

    active = list(range(n))
    nodes = [Tree(name=labels[i]) for i in range(n)]
    dist = [[float(matrix[i][j]) for j in range(n)] for i in range(n)]

    while len(active) > 1:
        m = len(active)
        if m == 2:
            i, j = active
            half = dist[i][j] / 2.0
            nodes[i].branch_length = half
            nodes[j].branch_length = half
            root = Tree()
            root.append_child(nodes[i])
            root.append_child(nodes[j])
            return root

        # row sums over the active set
        row_sum = {a: sum(dist[a][b] for b in active) for a in active}
        # Q matrix, pick min pair
        best = None
        for x in range(m):
            for y in range(x + 1, m):
                a, b = active[x], active[y]
                q = (m - 2) * dist[a][b] - row_sum[a] - row_sum[b]
                if best is None or q < best[0]:
                    best = (q, a, b)
        _, i, j = best

        limb_i = 0.5 * dist[i][j] + (row_sum[i] - row_sum[j]) / (2.0 * (m - 2))
        limb_j = dist[i][j] - limb_i
        limb_i = max(limb_i, 0.0)
        limb_j = max(limb_j, 0.0)

        u = Tree()
        nodes[i].branch_length = limb_i
        nodes[j].branch_length = limb_j
        u.append_child(nodes[i])
        u.append_child(nodes[j])

        # new active clusters (drop i, j; add u)
        new_active = [a for a in active if a not in (i, j)]
        keys = new_active + ["_u"]
        k2i = {a: x for x, a in enumerate(new_active)}
        new_dist = [[0.0] * len(keys) for _ in keys]
        for x, a in enumerate(new_active):
            for y, b in enumerate(new_active):
                if a != b:
                    new_dist[x][y] = new_dist[y][x] = dist[a][b]
            # new cluster <-> u
            d = 0.5 * (dist[i][a] + dist[j][a] - dist[i][j])
            new_dist[x][-1] = new_dist[-1][x] = max(d, 0.0)

        active = list(range(len(keys)))
        dist = new_dist
        nodes = [nodes[a] for a in new_active] + [u]

    # unreachable for n>=2 (loop always returns at m==2)
    raise RuntimeError("unexpected neighbour-joining state")


# ---------------------------------------------------------------------- --
# convenience: sequences -> tree
# ---------------------------------------------------------------------- --
def build_tree(
    sequences,
    model: str = "p",
    rooted: bool = True,
) -> Tree:
    """Build a Tree from an alignment / sequences using distance + neighbour joining.

    ``sequences`` may be a ``{name: sequence}`` dict, an iterable of
    ``(name, sequence)`` pairs, a list of raw sequence strings, or Biopython
    ``SeqRecord`` objects.  Sequences must be aligned (equal length).

    ``model`` selects the distance model (``"p"`` or ``"jc"`` / ``"jc69"``).
    Returns a :class:`Tree`.
    """
    labels, matrix = distance_matrix(sequences, model=model)
    tree = neighbor_joining(labels, matrix)
    if not rooted:
        tree.unroot()
    return tree
