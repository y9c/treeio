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


class Tree:
    """
    Tree class is used to store a tree object.
    """

    def __init__(self, name="tree"):
        """Init."""
        self.name = name
        self._parent: Optional[Tree] = None
        self._children: List[Tree] = []

    def __repr__(self):
        """Print."""
        return f"<Tree: {self.name}>"

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
        if hasattr(value, "__iter__") and all(isinstance(n, type(self)) for n in value):
            self._children = list(value)
            for node in value:
                # should not use `.parent`, or will set twice
                node._parent = self
        else:
            raise ValueError("xxx")

    @children.deleter
    def children(self):
        del self._children


# Alias
TreeNode = Tree
