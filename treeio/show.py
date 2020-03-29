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

    # lmrFromStrings :: [String] -> ([String], String, [String])
    def lmrFromStrings(xs):
        """Lefts, Mid, Rights."""
        i = len(xs) // 2
        ls, rs = xs[0:i], xs[i:]
        return ls, rs[0], rs[1:]

    # stringsFromLMR :: ([String], String, [String]) -> [String]
    def stringsFromLMR(lmr):
        ls, m, rs = lmr
        return ls + [m] + rs

    # fghOverLMR
    # :: (String -> String)
    # -> (String -> String)
    # -> (String -> String)
    # -> ([String], String, [String])
    # -> ([String], String, [String])
    def fghOverLMR(f, g, h):
        def go(lmr):
            ls, m, rs = lmr
            return ([f(x) for x in ls], g(m), [h(x) for x in rs])

        return lambda lmr: go(lmr)

    # leftPad :: Int -> String -> String
    def leftPad(n):
        return lambda s: (" " * n) + s

    # treeFix :: (Char, Char, Char) -> ([String], String, [String])
    #                               ->  [String]
    def treeFix(l, m, r):
        def cfix(x):
            return lambda xs: x + xs

        return compose(stringsFromLMR)(fghOverLMR(cfix(l), cfix(m), cfix(r)))

    def lmrBuild(w, f):
        def go(wsTree):
            nChars, x = wsTree["root"]
            _x = ("─" * (w - nChars)) + x
            xs = wsTree["nest"]
            lng = len(xs)

            # linked :: String -> String
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
                            + [treeFix("│", "├", "│")(x) for x in lmrs[1:-1]]
                            + [treeFix("│", "└", " ")(lmrs[-1])]
                        )
                    )
                )

        return lambda wsTree: go(wsTree)

    # foldr :: (a -> b -> b) -> b -> [a] -> b
    def foldr(f):
        """
        Right to left reduction of a list, using the binary operator f,
        and starting with an initial accumulator value.
        """

        def g(x, a):
            return f(a, x)

        return lambda acc: lambda xs: reduce(g, xs[::-1], acc)

    def measured(x) -> Tuple[int, str]:
        """Value of a tree node tupled with string length."""
        s = " " + str(x) + " "
        return len(s), s

    def go(x):
        return Node(measured(x["root"]))([go(v) for v in x["nest"]])

    measured_tree = go(tree)
    print(measured_tree)

    print(levels(measured_tree))

    level_widths = reduce(
        lambda a, xs: a + [max(x[0] for x in xs)], levels(measured_tree), [],
    )
    print(level_widths)

    tree_lines = stringsFromLMR(
        foldr(lmrBuild)(None)(level_widths)(measured_tree)
    )
    # Modify tree lines if is not compact and is pruned.
    if not is_compact and is_pruned:
        tree_lines = [s for s in tree_lines if any(c not in "│ " for c in s)]
    return "\n".join(tree_lines)


# GENERIC -------------------------------------------------

# Node :: a -> [Tree a] -> Tree a
def Node(v):
    """Contructor for a Tree node which connects a
        value of some kind to a list of zero or
        more child trees.
    """
    return lambda xs: {"root": v, "nest": xs}


# compose (<<<) :: (b -> c) -> (a -> b) -> a -> c
def compose(g):
    """Right to left function composition."""
    return lambda f: lambda x: g(f(x))


# concatMap :: (a -> [b]) -> [a] -> [b]
def concatMap(f):
    """A concatenated list over which a function has been mapped.
        The list monad can be derived by using a function f which
        wraps its output in a list,
        (using an empty list to represent computational failure).
    """
    return lambda xs: list(chain.from_iterable(map(f, xs)))


# fmapTree :: (a -> b) -> Tree a -> Tree b
def fmapTree(f: Callable):
    """A new tree holding the results of
        applying f to each root in
        the existing tree.
    """

    def go(x):
        return Node(f(x["root"]))([go(v) for v in x["nest"]])

    return lambda tree: go(tree)


# intercalate :: [a] -> [[a]] -> [a]
# intercalate :: String -> [String] -> String
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


# iterate :: (a -> a) -> a -> Gen [a]
def iterate(f):
    """An infinite list of repeated
        applications of f to x.
    """

    def go(x):
        v = x
        while True:
            yield v
            v = f(v)

    return lambda x: go(x)


# levels :: Tree a -> [[a]]
def levels(tree):
    """A list of the nodes at each level of the tree."""
    return list(
        map_(map_(root))(takewhile(bool, iterate(concatMap(nest))([tree])))
    )


# map :: (a -> b) -> [a] -> [b]
def map_(f):
    """The list obtained by applying f
        to each element of xs.
    """
    return lambda xs: list(map(f, xs))


# nest :: Tree a -> [Tree a]
def nest(t):
    """Accessor function for children of tree node."""
    return t["nest"] if "nest" in t else None


# root :: Tree a -> a
def root(t):
    """Accessor function for data of tree node."""
    return t["root"] if "root" in t else None


# main :: IO ()
# tree1 :: Tree Int
tree1 = Node(1)(
    [
        Node(2)([Node(4)([Node(7)([])]), Node(5)([])]),
        Node(3)([Node(6)([Node(8)([]), Node(9)([])])]),
    ]
)

# tree :: Tree String
tree2 = Node("Alpha")(
    [
        Node("Beta")(
            [
                Node("Epsilon")([]),
                Node("Zeta")([]),
                Node("Eta")([]),
                Node("Theta")([Node("Mu")([]), Node("Nu")([])]),
            ]
        ),
        Node("Gamma")(
            [Node("Xiiiiiiiiiiiiiiiiiiiiii")([Node("Omicron")([])])]
        ),
        Node("Delta")(
            [Node("Iota")([]), Node("Kappa")([]), Node("Lambda")([])]
        ),
    ]
)


def main():
    print(tree2ascii(tree2, True, True))
    # "Fully compacted (parents not all centered):"
    print(tree2ascii(tree2, True, False))
    # "Expanded with vertically centered parents:"
    print(tree2ascii(tree2, False, False))
    # "Centered parents with nodeless lines pruned out:"
    print(tree2ascii(tree2, False, True))
