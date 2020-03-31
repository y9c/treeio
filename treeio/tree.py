#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-03-27 00:28

"""Tree class."""


from __future__ import annotations

from typing import Optional, Iterable, List
from itertools import chain

from .utils import dedup


class Tree:
    """
    Tree class is used to store a tree object.

                  ┌ Theta
          ┌─ Beta ┤
    Alpha ┤       └ Delta
          │
          └ Gamma
    """

    def __init__(self, name="unknown", dist=None, supp=None, **kwargs):
        """Init."""
        self.name: str = name
        self.dist: Optional[float] = dist
        self.supp: Optional[float] = supp
        self._parent: Optional[Tree] = None
        self._children: List[Tree] = []
        # support any value
        # Is this a good feature or not?
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __iter__(self):
        for v in chain(*map(iter, self.children)):
            yield v
        yield self

    def __repr__(self):
        """Print."""
        return f"<Tree: {self.name}>"

    def __str__(self):
        """ Print tree in console by ascii art."""
        # return self.as_ascii()
        from .show import tree2ascii

        return tree2ascii(self, False, True)

    @property
    def parent(self) -> Optional[Tree]:
        """Get the parent of tree node."""
        return self._parent

    @parent.setter
    def parent(self, value: Tree) -> None:
        if isinstance(value, type(self)) or value is None:
            self._parent = value
            # should not use `.children`, or will set twice
            value._children.append(self)
        else:
            raise ValueError("xxx")

    @parent.deleter
    def parent(self):
        del self._parent

    @property
    def children(self) -> List[Tree]:
        """Get the children of tree node."""
        return self._children

    @children.setter
    def children(self, value: Iterable[Tree]) -> None:
        if hasattr(value, "__iter__") and all(
            isinstance(n, type(self)) for n in value
        ):
            self._children = dedup(value)
            for node in value:
                # should not use `.parent`, or will set twice
                node._parent = self
        else:
            raise ValueError("xxx")

    @children.deleter
    def children(self):
        del self._children

    def append_child(self, tree: Tree):
        if tree not in self.children:
            self.children += [tree]
        return self

    def extend_children(self, tree: List[Tree]):
        self.children += tree
        return self

    def remove_child(self, tree):
        assert self.parent is self, "The input node is not a child node."
        tree._parent = None
        self._children = [c for c in self.children if c is not tree]
        return self

    def isolated(self):
        """Isolate tree, turn into root node."""
        if self.is_root():
            return self
        p = self._parent
        p._children = [c for c in p.children if c is not self]
        self._parent = None
        return self

    def is_leaf(self):
        """Chech node is a leaf(terminal node) or not."""
        return len(self.children) == 0

    def is_root(self):
        """Chech node is a root(starting node) or not."""
        return self.parent is None


if __name__ == "__main__":
    tree = Tree("A")
    clade = Tree("B")
    clade.children = [Tree("X"), Tree("Y")]
    tree.children = [clade, Tree("C"), Tree("D")]
    print(tree)
