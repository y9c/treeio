#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-03-29 18:32
# Tree style and algorithm inspired by the Haskell snippet at:
# https://doisinkidney.com/snippets/drawing-trees.html
# https://rosettacode.org/wiki/Visualize_a_tree#Simple_decorated-outline_tree

"""
Textually visualized tree, with vertically-centered parent nodes.
"""

from typing import Callable, Tuple, List
from functools import reduce
from itertools import chain, takewhile

from .tree import Tree
from .jt import read_json


def tree2ascii(tree: Tree, is_compact: bool, is_pruned: bool) -> str:
    """
    Monospaced UTF8 left-to-right text tree.

    In a compact or expanded format,
    with any lines containing no nodes optionally pruned out.

    Bool -> Bool -> Tree a -> String
    Demo output.

                ┌ Epsilon
                ├─── Zeta
        ┌─ Beta ┼──── Eta
        │       │         ┌───── Mu
        │       └── Theta ┤
    Alpha ┤                 └───── Nu
        ├ Gamma ────── Xi ─ Omicron
        │       ┌─── Iota
        └ Delta ┼── Kappa
                └─ Lambda
    """

    def padding(s, n=1):
        return " " * n + s + " " * n

    def compose(g):
        """Right to left function composition."""
        return lambda f: lambda x: g(f(x))

    def intercalate(x: List[str]) -> Callable:
        """The concatenation of xs interspersed with copies of x."""
        return (
            lambda xs: list(
                chain.from_iterable(
                    reduce(lambda a, v: a + [x, v], xs[1:], [xs[0]])
                )
            )
            if xs
            else []
        )

    def lmrFromStrings(xs):
        """Lefts, Mid, Rights."""
        i = len(xs) // 2
        ls, rs = xs[0:i], xs[i:]
        return ls, rs[0], rs[1:]

    def stringsFromLMR(lmr):
        ls, m, rs = lmr
        return ls + [m] + rs

    def fghOverLMR(f, g, h):
        def go(lmr):
            ls, m, rs = lmr
            return ([f(x) for x in ls], g(m), [h(x) for x in rs])

        return lambda lmr: go(lmr)

    # leftPad :: Int -> String -> String
    def leftPad(n):
        return lambda s: (" " * n) + s

    def treeFix(l, m, r):
        def cfix(x):
            return lambda xs: x + xs

        return compose(stringsFromLMR)(fghOverLMR(cfix(l), cfix(m), cfix(r)))

    def levels(tree):
        def go(x):
            v = x
            while len(v) > 0:
                yield v
                v = [i for t in v if not t.is_leaf() in t for i in t.children]

        return [[padding(t.name) for t in x] for x in go([tree])]

    def foldr(tree_dict):
        """
        Right to left reduction of a list, using the binary operator f,
        and starting with an initial accumulator value.
        """

        level_widths = reduce(
            lambda a, xs: a + [max(len(x) for x in xs)], levels(tree_dict), [],
        )

        def g(f, w):
            def go(wsTree):
                x = padding(wsTree.name)
                nChars = len(x)
                _x = ("─" * (w - nChars)) + x
                xs = wsTree.children
                lng = len(xs)

                def linked(s):
                    c = s[0]
                    t = s[1:]
                    return (
                        _x + "┬" + t
                        if "┌" == c
                        else (
                            _x + "┤" + t
                            if "│" == c
                            else (_x + "┼" + t if "├" == c else (_x + "┴" + t))
                        )
                    )

                # LEAF ------------------------------------
                if 0 == lng:
                    return ([], _x, [])

                # SINGLE CHILD ----------------------------
                elif 1 == lng:

                    def lineLinked(z):
                        return _x + "─" + z

                    rightAligned = leftPad(1 + w)
                    return fghOverLMR(rightAligned, lineLinked, rightAligned)(
                        f(xs[0])
                    )

                # CHILDREN --------------------------------
                else:
                    rightAligned = leftPad(w)
                    lmrs = [f(x) for x in xs]
                    return fghOverLMR(rightAligned, linked, rightAligned)(
                        lmrFromStrings(
                            intercalate([] if is_compact else ["│"])(
                                [treeFix(" ", "┌", "│")(lmrs[0])]
                                + [
                                    treeFix("│", "├", "│")(x)
                                    for x in lmrs[1:-1]
                                ]
                                + [treeFix("│", "└", " ")(lmrs[-1])]
                            )
                        )
                    )

            return lambda wsTree: go(wsTree)

        return reduce(g, level_widths[::-1], None)(tree_dict)

    tree_lines = stringsFromLMR(foldr(tree))

    # Modify tree lines if is not compact and is pruned.
    if not is_compact and is_pruned:
        tree_lines = [s for s in tree_lines if any(c not in "│ " for c in s)]
    return "\n".join(tree_lines)


if __name__ == "__main__":
    with open("./data/animals.json") as f:
        json_string = f.read()
    treex = read_json(json_string)[0]
    print(tree2ascii(treex, True, True))
    # "Fully compacted (parents not all centered):"
    print(tree2ascii(treex, True, False))
    # "Expanded with vertically centered parents:"
    print(tree2ascii(treex, False, False))
    # "Centered parents with nodeless lines pruned out:"
    print(tree2ascii(treex, False, True))
