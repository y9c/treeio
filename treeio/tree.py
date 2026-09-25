#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-03-27 00:28

"""Tree class.

The :class:`Tree` object stores a rooted / unrooted phylogenetic tree as a
collection of nodes linked by ``parent`` / ``children`` pointers.  The design
is deliberately inspired by the R ``treeio``/``ape`` data model and the Python
``toytree``/``ete3`` API while staying dependency free.

A node holds a handful of slot attributes (``name``, ``dist``, ``supp``) but
also accepts arbitrary keyword attributes -- the same capability that lets
``treeio`` attach bootstrap values, substitution rates, dates, or any custom
annotation to a node::

    node = Tree("Homo", branch_length=0.1, rate=0.05, origin="Africa")

Attributes are stored *per node*.  Convenience helpers (:meth:`get_data`,
:meth:`set_data`, :meth:`get_tipdata`, ...) make it easy to work with
node- or tip-level annotation tables, rather like ``treeio``'s tipdata /
nodedata.

Examples
--------
Build a tree by hand::

    root = Tree("root")
    a, b, c = Tree("A"), Tree("B", branch_length=0.5), Tree("C", branch_length=0.7)
    a.children = [b, c]
    root.children = [a, Tree("D", branch_length=1.0)]

Or parse one directly from a newick string::

    from treeio import read_newick
    tree = read_newick("(A:0.1,B:0.2):0.3,C:0.4;")
"""

from __future__ import annotations

from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Union


class Tree:
    """
    Tree class is used to store a tree object.

                  ┌ Theta
          ┌─ Beta ┤
    Alpha ┤       └ Delta
          │
          └ Gamma
    """

    # first-class slots (``__slots__`` keeps per-node memory small); arbitrary
    # annotation attributes are consolidated into a single lazy ``_annotations``
    # dict instead of scattering extra keys into each node's ``__dict__``.
    __slots__ = ("name", "branch_length", "support", "_parent", "_children",
                 "_annotations", "__weakref__")

    # slots that are treated as first-class attributes.
    _SLOTS = ("name", "branch_length", "support", "_parent", "_children")

    def __init__(self, name="unknown", branch_length=None, support=None, **kwargs):
        """Create a node.

        Each argument is stored as an attribute on the node.  ``branch_length``
        is the branch length leading to this node (``None`` if unknown) and
        ``support`` the support / bootstrap value (``None`` if unknown).

        Keyword arguments become arbitrary per-node annotations (stored in a
        single dict).  If ``parent`` is supplied the node is attached
        underneath it.
        """
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "branch_length", branch_length)
        object.__setattr__(self, "support", support)
        object.__setattr__(self, "_parent", None)
        object.__setattr__(self, "_children", [])
        object.__setattr__(self, "_annotations", None)

        # support any value, attached as annotation.
        for key, value in kwargs.items():
            if key == "parent":
                self.parent = value
            else:
                setattr(self, key, value)

    def __setattr__(self, name: str, value) -> None:
        """Store slots / properties directly; everything else in ``_annotations``."""
        try:
            object.__setattr__(self, name, value)
        except AttributeError:
            ann = object.__getattribute__(self, "_annotations")
            if ann is None:
                ann = {}
                object.__setattr__(self, "_annotations", ann)
            ann[name] = value

    def __getattr__(self, name: str):
        """Look up annotation attributes stored in ``_annotations``."""
        ann = object.__getattribute__(self, "_annotations")
        if ann is not None and name in ann:
            return ann[name]
        raise AttributeError(name)


    # ---------------------------------------------------------------- ---
    # iteration / string
    # ---------------------------------------------------------------- ---
    def __iter__(self) -> Iterator["Tree"]:
        """Iterate over all descendant nodes (post-order), self last."""
        stack = [(self, False)]
        while stack:
            node, visited = stack.pop()
            if visited:
                yield node
            else:
                stack.append((node, True))
                # push children in reverse so they are yielded left-to-right
                for ch in reversed(node._children):
                    stack.append((ch, False))

    def __len__(self) -> int:
        """Number of nodes that are descendants of (and including) self."""
        count = 0
        stack = [self]
        while stack:
            node = stack.pop()
            count += 1
            stack.extend(node._children)
        return count

    def __contains__(self, node) -> bool:
        """``node in tree`` is true when ``node`` is a descendant of ``self``."""
        return self.is_ancestor_of(node)

    def __bool__(self) -> bool:
        return True

    def __repr__(self) -> str:
        return f"<Tree: {self.name}>"

    def __str__(self) -> str:
        """Print tree in console by ascii art."""
        from .show import tree_to_ascii

        return tree_to_ascii(self, False, True)

    def __eq__(self, other) -> bool:
        if not isinstance(other, Tree):
            return NotImplemented
        stack = [(self, other)]
        while stack:
            a, b = stack.pop()
            if a.name != b.name or a.branch_length != b.branch_length or a.support != b.support:
                return False
            if len(a._children) != len(b._children):
                return False
            stack.extend(zip(a._children, b._children))
        return True

    def __hash__(self) -> int:  # used by set/dict if hashing needed
        return id(self)

    # ---------------------------------------------------------------- ---
    # parent / children
    # ---------------------------------------------------------------- ---
    @property
    def parent(self) -> Optional["Tree"]:
        """Get the parent of tree node."""
        return self._parent

    @parent.setter
    def parent(self, value: Optional["Tree"]) -> None:
        if value is None:
            self._detach_from_parent()
            return
        if isinstance(value, type(self)):
            if self._parent is not None and self._parent is not value:
                self._detach_from_parent()
            self._parent = value
            if not any(c is self for c in value._children):
                value._children.append(self)
        else:
            raise ValueError("parent must be a Tree or None")

    @parent.deleter
    def parent(self) -> None:
        self._detach_from_parent()

    def _detach_from_parent(self) -> None:
        """Remove self from the current parent's children (keeps graph consistent)."""
        p = self._parent
        self._parent = None
        if p is not None:
            p._children = [c for c in p._children if c is not self]

    @property
    def children(self) -> "_ChildrenView":
        """Read-only view of the children (use :meth:`append_child` to mutate)."""
        return _ChildrenView(self._children)

    @children.setter
    def children(self, value: Iterable["Tree"]) -> None:
        if hasattr(value, "__iter__") and all(
            isinstance(n, type(self)) for n in value
        ):
            self._children = []
            for node in value:
                self.append_child(node)
        else:
            raise ValueError("children must be an iterable of Tree nodes")

    @children.deleter
    def children(self) -> None:
        for c in list(self._children):
            c._parent = None
        self._children = []

    def append_child(self, tree: "Tree"):
        """Add a child node (detaching it from any previous parent)."""
        if tree is self:
            raise ValueError("cannot attach a node to itself")
        if any(c is tree for c in self._children):
            return self
        if tree._parent is not None:
            tree._detach_from_parent()
        tree._parent = self
        self._children.append(tree)
        return self

    def extend_children(self, tree: Iterable["Tree"]):
        """Add several children and return self."""
        for t in tree:
            self.append_child(t)
        return self

    def remove_child(self, tree: "Tree"):
        """Detach ``tree`` from ``self``.

        Raises ``ValueError`` if ``tree`` is not a direct child of ``self``.
        """
        if not any(c is tree for c in self._children):
            raise ValueError("The input node is not a child node.")
        tree._parent = None
        self._children = [c for c in self._children if c is not tree]
        return self

    def isolated(self) -> "Tree":
        """Isolate tree, turn into root node."""
        self._detach_from_parent()
        return self

    # ---------------------------------------------------------------- ---
    # topology queries
    # ---------------------------------------------------------------- ---
    def is_leaf(self) -> bool:
        """Check node is a leaf (terminal node) or not."""
        return len(self.children) == 0

    def is_root(self) -> bool:
        """Check node is a root (starting node) or not."""
        return self.parent is None

    def is_internal(self) -> bool:
        """Check node is an internal (non-terminal, non-root) node."""
        return not self.is_leaf()

    def is_binary(self) -> bool:
        """Check the (sub)tree is fully bifurcating (no multifurcation)."""
        if self.is_leaf():
            return True
        if len(self.children) != 2:
            return False
        return all(c.is_binary() for c in self.children)

    # alias used by toytree / ete3
    is_bifurcating = is_binary

    def is_rooted(self) -> bool:
        """Return ``True`` (a Tree is always rooted in this model).

        Unrooted trees are represented with a degree-3 root node; helpers
        like :meth:`unroot` produce that representation.
        """
        return True

    def get_tips(self) -> List["Tree"]:
        """Return leaf nodes (terminal taxa) in depth-first order."""
        return [n for n in self if n.is_leaf()]

    get_leaves = get_tips

    def get_internal_nodes(self) -> List["Tree"]:
        """Return all non-terminal nodes in depth-first order."""
        return [n for n in self if n.is_internal()]

    def get_nodes(self) -> List["Tree"]:
        """Return all nodes in depth-first (post-order) order."""
        return list(self)

    def get_names(self) -> List[str]:
        """Return the names of every node in post-order."""
        return [n.name for n in self]

    def get_tip_names(self) -> List[str]:
        """Return tip labels in the order they appear (depth-first)."""
        return [n.name for n in self.get_tips()]

    @property
    def tips(self) -> List["Tree"]:
        return self.get_tips()

    @property
    def tip_name2node(self) -> Dict[str, "Tree"]:
        return {n.name: n for n in self.get_tips()}

    @property
    def tip_names(self) -> List[str]:
        """Tip labels in depth-first order."""
        return self.get_tip_names()

    @property
    def node_names(self) -> List[str]:
        """All node labels in post-order."""
        return [n.name for n in self]

    @property
    def node_index(self) -> Dict[str, int]:
        """Map node name -> 0-based post-order index."""
        return {n.name: i for i, n in enumerate(self)}

    @property
    def nleaves(self) -> int:
        return len(self.get_tips())

    @property
    def ntips(self) -> int:
        return self.nleaves

    @property
    def nnodes(self) -> int:
        return len(self)

    @property
    def total_node_count(self) -> int:
        return self.nnodes

    @property
    def nternal_node_count(self) -> int:
        return len(self.get_internal_nodes())

    def is_monophyletic(self, *names) -> bool:
        """Return ``True`` if *all* ``names`` form a single clade."""
        tips = {n.name for n in self.get_tips()}
        target = set(names)
        if not target <= tips:
            return False
        # find node containing exactly target
        for node in self.get_nodes():
            below = {n.name for n in node.get_tips()}
            if below == target:
                return True
        return False

    def is_outgroup(self, *names) -> bool:
        """True when the clade of ``names`` is sibling to everything else."""
        tips = {n.name for n in self.get_tips()}
        target = set(names)
        if not target <= tips:
            return False
        # an outgroup splits the root into exactly two children
        if len(self.children) != 2:
            return False
        for c in self.children:
            if {n.name for n in c.get_tips()} == target:
                return True
        return False

    # ---------------------------------------------------------------- ---
    # traversal (iterators)
    # ---------------------------------------------------------------- ---
    def traverse(self, order: str = "preorder") -> Iterator["Tree"]:
        """Iterate over descendant nodes.

        order: one of ``'preorder'``, ``'postorder'``, ``'levelorder'``.
        The node ``self`` is included.  (Iterative, so it copes with very
        deep trees without hitting the recursion limit.)
        """
        if order == "postorder":
            stack = [(self, False)]
            while stack:
                node, visited = stack.pop()
                if visited:
                    yield node
                else:
                    stack.append((node, True))
                    for ch in reversed(node._children):
                        stack.append((ch, False))
        elif order == "levelorder":
            frontier = [self]
            while frontier:
                nxt = []
                for node in frontier:
                    yield node
                    nxt.extend(node._children)
                frontier = nxt
        elif order == "preorder":
            stack = [self]
            while stack:
                node = stack.pop()
                yield node
                # push children in reverse so they are visited left-to-right
                for ch in reversed(node._children):
                    stack.append(ch)
        else:
            raise ValueError(f"unknown order {order!r}")

    def is_ancestor_of(self, node: "Tree") -> bool:
        """Return ``True`` if ``node`` is ``self`` or a descendant of ``self``."""
        cur = node
        while cur is not None:
            if cur is self:
                return True
            cur = cur._parent
        return False

    def walk(self, order: str = "preorder") -> Iterator["Tree"]:
        """Alias of :meth:`traverse`."""
        return self.traverse(order)

    # ---------------------------------------------------------------- ---
    # search / ancestry
    # ---------------------------------------------------------------- ---
    def label_index(self) -> Dict[str, "Tree"]:
        """Build ``{node_name: node}`` (first occurrence wins) in one pass.

        Use this for bulk label lookups (e.g. data attachment) instead of
        calling :meth:`get_node_by_label` repeatedly.
        """
        index: Dict[str, "Tree"] = {}
        for node in self.traverse("preorder"):
            name = node.name
            if name not in index:
                index[name] = node
        return index

    def get_node_by_label(self, name: str) -> Optional["Tree"]:
        """Return the (first) node whose ``name`` equals ``name``."""
        for n in self:
            if n.name == name:
                return n
        return None

    def get_tip_by_label(self, name: str) -> Optional["Tree"]:
        """Return the tip whose label equals ``name``."""
        for n in self.get_tips():
            if n.name == name:
                return n
        return None

    def get_ancestors(self, node: "Tree") -> List["Tree"]:
        """Return the chain of ancestors from ``node`` up to the root."""
        out, cur = [], node
        while cur is not None and cur is not self:
            out.append(cur)
            cur = cur.parent
        if cur is not None:
            out.append(self)
        return out

    def get_ancestors_of(self, node: "Tree") -> List["Tree"]:
        """Alias of :meth:`get_ancestors` (own ancestors, not including node)."""
        out, cur = [], node.parent
        while cur is not None:
            out.append(cur)
            cur = cur.parent
        return out

    def get_descendants(self, node: "Tree") -> List["Tree"]:
        """Return all strict descendants of ``node`` (excluding ``node``)."""
        return [n for n in node if n is not node]

    def get_nodes_below(self, node: "Tree") -> List["Tree"]:
        """Return the clade rooted at ``node`` (including it)."""
        return list(node)

    def get_common_ancestor(self, node1: "Tree", node2: "Tree") -> Optional["Tree"]:
        """Return the most recent common ancestor of two nodes."""
        anc1 = set(self._ancestor_set(node1))
        cur = node2
        while cur is not None:
            if cur in anc1:
                return cur
            cur = cur.parent
        return None

    def get_mrca(self, *nodes) -> Optional["Tree"]:
        """Return the most recent common ancestor of a set of nodes.

        Accepts ``Tree`` nodes or names (strings).
        """
        resolved = [
            n if isinstance(n, Tree) else self.get_node_by_label(n)
            for n in nodes
        ]
        if not resolved or any(n is None for n in resolved):
            return None
        out = resolved[0]
        for n in resolved[1:]:
            out = self.get_common_ancestor(out, n)
        return out

    # "mrca" as a method alias, following toytree/ete3 naming
    mrca = get_mrca

    def _ancestor_set(self, node: "Tree") -> List["Tree"]:
        out = []
        cur = node
        while cur is not None:
            out.append(cur)
            cur = cur.parent
        return out

    # ---------------------------------------------------------------- ---
    # distances
    # ---------------------------------------------------------------- ---
    def get_distance(self, node1: "Tree", node2: "Tree") -> float:
        """Return the patristic distance between two nodes.

        Branch lengths are summed along the path; missing (``None``) branch
        lengths are treated as 0.
        """
        if not isinstance(node1, Tree) or not isinstance(node2, Tree):
            raise ValueError("get_distance expects Tree nodes")

        anc1 = set(self._ancestor_set(node1))
        cur = node2
        lca = None
        while cur is not None:
            if cur in anc1:
                lca = cur
                break
            cur = cur.parent
        if lca is None:
            raise ValueError("nodes are not in the same tree")

        d = 0.0
        cur = node1
        while cur is not lca:
            d += self._dist(cur)
            cur = cur.parent
        cur = node2
        while cur is not lca:
            d += self._dist(cur)
            cur = cur.parent
        return d

    def get_cophenetic_distance(self) -> List[List[float]]:
        """Return the pairwise patristic distance matrix for all tips."""
        tips = self.get_tips()
        n = len(tips)
        mat = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                d = self.get_distance(tips[i], tips[j])
                mat[i][j] = mat[j][i] = d
        return mat

    def get_depth_to_root(self, node: "Tree") -> float:
        """Distance from ``node`` up to the root."""
        return self.get_distance(node, self)

    def get_tree_length(self) -> float:
        """Sum of all branch lengths in the (sub)tree."""
        return sum(self._dist(n) for n in self if n.parent is not None)

    def get_edge_lengths(self) -> List[float]:
        """Return the branch length of every non-root node (post-order)."""
        return [self._dist(n) for n in self if n.parent is not None]

    @staticmethod
    def _dist(node: "Tree") -> float:
        d = node.branch_length
        return 0.0 if d is None else float(d)

    # ---------------------------------------------------------------- ---
    # data / annotation attachment
    # ---------------------------------------------------------------- ---
    def get_data(self, key: str, default=None):
        """Get an annotation value.

        Annotations are looked up first (so a key that coincides with a method
        name, e.g. ``root``, still returns the stored annotation), then the
        genuine data slots (``name`` / ``branch_length`` / ``support``) and
        on-the-fly coordinate attributes.  Class attributes such as methods are
        never returned as data.
        """
        ann = self._annotations
        if ann is not None and key in ann:
            return ann[key]
        if key in ("name", "branch_length", "support"):
            return object.__getattribute__(self, key)
        if key.startswith("_"):
            # e.g. _x / _y set during layout, stored as annotations (never a method)
            return getattr(self, key, default)
        return default

    def set_data(self, key: str, value) -> None:
        """Set an annotation value."""
        setattr(self, key, value)

    def get_tipdata(self, key: str, default=None) -> Dict[str, object]:
        """Return ``{tip_name: value}`` for a tip-level annotation."""
        return {n.name: n.get_data(key, default) for n in self.get_tips()}

    def get_nodedata(self, key: str, default=None) -> Dict[str, object]:
        """Return ``{node_name: value}`` for an internal-node annotation."""
        return {n.name: n.get_data(key, default) for n in self.get_internal_nodes()}

    def set_nodedata(self, mapping: Dict[str, object]) -> None:
        """Set a per-node annotation from ``{node_name: value}``."""
        for name, value in mapping.items():
            node = self.get_node_by_label(name)
            if node is not None:
                node.set_data(name, value)

    # ---------------------------------------------------------------- ---
    # topology manipulation
    # ---------------------------------------------------------------- ---
    def ladderize(self, reverse: bool = False) -> "Tree":
        """Reorder children so the tree is 'laddered' (sorted by clade size).

        Mutates the tree in place and returns ``self``.
        """
        self._ladderize(reverse)
        return self

    def _ladderize(self, reverse: bool) -> None:
        for c in self.children:
            c._ladderize(reverse)
        self.children = sorted(
            self.children,
            key=lambda c: (c.nleaves, c.name),
            reverse=reverse,
        )

    def rotate(self, node: "Tree" = None) -> "Tree":
        """Rotate (reverse the order of) a node's children in place.

        If ``node`` is ``None`` the root's children are rotated.  Returns self.
        """
        if node is None:
            node = self
        node.children = list(reversed(node.children))
        return self

    def root(self, outgroup: Union["Tree", str, Sequence] = None) -> "Tree":
        """Re-root the tree so ``outgroup`` is a basal (direct child) clade.

        ``outgroup`` may be a :class:`Tree` node, a single tip name, or a
        sequence of names forming a clade.  The tree is modified in place and
        ``self`` remains the root.  Returns self.
        """
        if outgroup is None:
            return self
        out = self._resolve_clade(outgroup)
        if out is None:
            raise ValueError("outgroup not found in tree")

        # walk the path out -> root, collecting every sibling clade encountered
        pieces = []
        node = out
        while node.parent is not None:
            p = node.parent
            p._children = [c for c in p._children if c is not node]
            pieces.extend(c for c in p.children if c is not node)
            node = p
            if node is self:
                break

        # reconnect the remaining clades under a single "rest" subtree
        if len(pieces) == 1:
            rest = pieces[0]
        elif len(pieces) > 1:
            rest = Tree("internal", branch_length=0.0)
            rest.children = pieces
        else:
            rest = None

        if rest is None:
            self.children = [out]
        else:
            self.children = [out, rest]
        return self

    def _resolve_clade(self, outgroup) -> Optional["Tree"]:
        """Resolve an outgroup spec into a node.

        Accepts a :class:`Tree`, a tip/name string, or a sequence of names.
        """
        if isinstance(outgroup, Tree):
            return outgroup
        if isinstance(outgroup, str):
            return self.get_tip_by_label(outgroup) or self.get_node_by_label(
                outgroup
            )
        if isinstance(outgroup, (list, tuple, set)):
            return self.get_mrca(*list(outgroup))
        return None

    def unroot(self) -> "Tree":
        """Collapse the root into a single (possibly polytomous) network node.

        Branch lengths on the two root edges are averaged so overall tree
        length is preserved.  Returns self.
        """
        if len(self.children) == 2:
            c1, c2 = self.children
            new_dist = (self._dist(c1) + self._dist(c2)) / 2.0
            c1.branch_length = new_dist
            c2.branch_length = new_dist
        return self

    def drop_tips(self, names: Iterable[str]) -> "Tree":
        """Remove tips whose name is in ``names`` and prune emptied clades.

        Mutates and returns ``self``.
        """
        banned = set(names)
        kept = [t for t in self.get_tips() if t.name not in banned]
        keptset = set(kept)
        if not keptset:
            self._children = []
            self.name = "unknown"
            self.branch_length = None
            return self
        self._prune_to(keptset)
        return self

    def _prune_to(self, keptset: set) -> None:
        """Keep only descendants that lead to a tip in ``keptset``."""
        if self.is_leaf():
            return
        self._children = [c for c in self._children if c._has_kept(keptset)]
        for c in self._children:
            c._parent = self
            c._prune_to(keptset)
        self._splice_single_child()

    def _has_kept(self, keptset: set) -> bool:
        """True if any tip of this (sub)tree is in ``keptset``."""
        return any(t in keptset for t in self.get_tips())

    def _splice_single_child(self) -> None:
        # Collapse any node (including the root) that has exactly one child by
        # promoting the child upward (its incoming branch length is added to the
        # child's, so total tree length is preserved).
        if len(self._children) == 1:
            child = self._children[0]
            if child.branch_length is not None or self.branch_length is not None:
                child.branch_length = (self._dist(self) + self._dist(child)) or None
            p = self._parent
            if p is not None:
                for i, c in enumerate(p._children):
                    if c is self:
                        p._children[i] = child
                        break
                child._parent = p
            else:
                # self is the root: adopt the child's identity and children
                self.name = child.name
                self.support = child.support
                for k, v in _iter_attrs(child):
                    setattr(self, k, v)
                self._children = list(child._children)
                for gc in self._children:
                    gc._parent = self
        for c in self._children:
            c._splice_single_child()

    def prune(self, names: Iterable[str]) -> "Tree":
        """Keep only the tips named in ``names``; everything else is dropped.

        Mutates and returns ``self``.
        """
        return self.drop_tips(
            set(self.get_tip_names()) - set(names)
        )

    def collapse(self, node: "Tree") -> "Tree":
        """Collapse the clade under ``node`` into a single node.

        ``node`` becomes a tip while keeping its name/support.  Returns self.
        """
        node._children = []
        return self

    def resolve_polytomies(self) -> "Tree":
        """Insert (dummy, length-0) internal nodes to make the tree
        bifurcating.  Mutates and returns ``self``.
        """
        for c in self._children:
            c.resolve_polytomies()
        if len(self._children) > 2:
            kids = self._children
            first = kids[0]
            tail = kids[1]
            for k in kids[2:]:
                dummy = Tree("internal", branch_length=0.0)
                dummy._parent = self
                dummy._children = [tail, k]
                tail._parent = dummy
                k._parent = dummy
                tail = dummy
            self._children = [first, tail]
            first._parent = self
            tail._parent = self
        return self

    def get_dummy_polytomy(self) -> "Tree":
        """Return a copy of the tree's root polytomy resolved to a dummy node."""
        copy = self.copy()
        root = copy
        if len(root._children) > 2:
            kids = root._children
            dummy = Tree("internal", branch_length=0.0)
            dummy._children = list(kids)
            for k in kids:
                k._parent = dummy
            root._children = [dummy]
            dummy._parent = root
        return copy

    # ---------------------------------------------------------------- ---
    # copy / export
    # ---------------------------------------------------------------- ---
    def copy(self, deep: bool = True) -> "Tree":
        """Independent deep copy of the subtree rooted at self."""
        new = Tree()
        # iterative clone (copes with very deep trees without recursion)
        new.name = self.name
        new.branch_length = self.branch_length
        new.support = self.support
        _copy_attrs(self, new)
        stack = [(self, new)]
        while stack:
            src, dst = stack.pop()
            for child in src.children:
                nc = Tree()
                nc.name = child.name
                nc.branch_length = child.branch_length
                nc.support = child.support
                _copy_attrs(child, nc)
                nc._parent = dst
                dst._children.append(nc)
                stack.append((child, nc))
        return new

    def clone(self) -> "Tree":
        """Alias of :meth:`copy`."""
        return self.copy()

    def get_subtree(self, node: "Tree") -> "Tree":
        """Return a copy of the clade rooted at ``node``."""
        return node.copy()

    # ---------------------------------------------------------------- ---
    # coordinates for plotting
    # ---------------------------------------------------------------- ---
    def get_y_and_x(self, x_as_branch: bool = True, root_dist: float = 0.0):
        """Assign ``x`` / ``y`` coordinates for a rectangular layout.

        ``x`` is the cumulative branch length from the root (or the node
        index if ``x_as_branch`` is False); ``y`` is the tip index position
        (internal nodes get the midpoint of their children).  Mutates the
        nodes in place.  Returns ``self`` so it can be used fluently.
        """
        self._assign_layout(x_as_branch, root_dist)
        return self

    def _assign_layout(self, x_as_branch, root_dist) -> None:
        leaves = self.get_tips()
        for i, leaf in enumerate(leaves):
            leaf._y = i

        if x_as_branch:
            # accumulate branch length from the root downward (iterative)
            stack = [(self, root_dist)]
            while stack:
                node, x = stack.pop()
                node._x = x
                for c in node.children:
                    stack.append((c, x + self._dist(c)))
        else:
            # x is simply the post-order node index
            for i, node in enumerate(self.traverse("postorder")):
                node._x = float(i)

        # y for internal nodes is the midpoint of their children.  Iterate in
        # post-order (children first) so each child's _y is already computed.
        for node in self.get_internal_nodes():
            node._y = (node.children[0]._y + node.children[-1]._y) / 2.0
        return self

    # ---------------------------------------------------------------- ---
    # misc
    # ---------------------------------------------------------------- ---
    def get_node2index(self) -> Dict["Tree", int]:
        """Map every node to a 0-based index in post-order."""
        return {n: i for i, n in enumerate(self)}

    def get_index2node(self) -> List["Tree"]:
        return list(self)

    def edge_vector(self) -> List[tuple]:
        """Return list of ``(parent_index, child_index)`` in post-order."""
        idx = self.get_node2index()
        return [(idx[n._parent], idx[n]) for n in self if n._parent is not None]

    def to_table(self, layout: str = "rectangular") -> List[dict]:
        """Return the tree as a list of node records (a data-frame style table).

        Each row describes one node with ``id``, ``parent``, ``label``,
        ``is_tip``, ``branch_length``, ``support``, ``x`` and ``y`` coordinates
        plus every per-node annotation attribute.  This is the analogue of
        the classic ``fortify`` data-frame view for phylogenetic trees.
        """
        return fortify(self, layout=layout)

    def as_newick(self, *args, **kwargs) -> str:
        """Serialize to newick string."""
        from .newick import write_newick

        return write_newick(self, *args, **kwargs)

    def save(self, path: str, layout: str = "rectangular", dpi: int = 300,
             format: str = None, tight: bool = True, **kwargs):
        """Draw this tree and save a publication-ready figure to ``path``."""
        from .plot import save as _save

        return _save(self, path, layout=layout, dpi=dpi, format=format, tight=tight, **kwargs)

    def render(self, backend: str = "mpl", layout: str = "rectangular", **kwargs):
        """Render this tree with the chosen backend.

        ``backend`` may be ``'mpl'`` (returns a matplotlib axes) or
        ``'ascii'`` (returns ASCII art).  See :func:`treeio.render`.
        """
        from .plot import render as _render

        return _render(self, backend=backend, layout=layout, **kwargs)

    @classmethod
    def from_alignment(cls, sequences, model: str = "p", rooted: bool = True) -> "Tree":
        """Build a Tree from an alignment via distance + neighbour joining.

        ``sequences`` may be a ``{name: sequence}`` dict, iterable of
        ``(name, sequence)`` pairs, raw sequence strings, or Biopython
        ``SeqRecord`` objects (aligned, equal length).  ``model`` selects the
        distance model (``"p"`` or ``"jc"`` / ``"jc69"``).

        Equivalent to :func:`treeio.build_tree`.
        """
        from .phylogeny import build_tree

        return build_tree(sequences, model=model, rooted=rooted)

    def as_dict(self, **fields) -> dict:
        """Serialize to a nested python dict (deep-tree safe, iterative).

        By default ``name``, ``branch_length`` and ``support`` are emitted for
        every node, and ``children`` follows the subtree.  Custom annotation
        fields can be requested via keyword ``name->attr``, e.g.
        ``as_dict(rate="rate")``.
        """
        field_map = {"name": "name", "branch_length": "branch_length", "support": "support"}
        field_map.update(fields)

        def make(node):
            return {key: node.get_data(attr) for key, attr in field_map.items()}

        root = make(self)
        # iterative post-order assembly (no recursion)
        out = {id(self): root}
        stack = [(self, root)]
        while stack:
            node, dic = stack.pop()
            if not node._children:
                continue
            kids = []
            for c in node._children:
                cd = make(c)
                kids.append(cd)
                out[id(c)] = cd
                stack.append((c, cd))
            dic["children"] = kids
        return root


def fortify(tree: Tree, layout: str = "rectangular") -> List[dict]:
    """Return the tree as a list of node records (a data-frame style table).

    Columns
    -------
    id, parent, label, is_tip, branch_length, support, x, y
    plus every per-node annotation attribute.

    ``layout`` accepts ``rectangular`` / ``circular`` to choose how ``x`` and
    ``y`` are computed.  This mirrors the standard
    ``fortify``.
    """
    # compute x/y coordinates (rectangular layout by default)
    tree.get_y_and_x(x_as_branch=True)
    index = tree.get_node2index()
    columns = ["id", "parent", "label", "is_tip", "branch_length", "support", "x", "y"]
    rows: List[dict] = []
    for node in tree.traverse("postorder"):
        row = {
            "id": index[node],
            "parent": index.get(node.parent),
            "label": node.name,
            "is_tip": node.is_leaf(),
            "branch_length": node.branch_length,
            "support": node.support,
            "x": float(getattr(node, "_x", 0.0)),
            "y": float(getattr(node, "_y", 0.0)),
        }
        # append per-node annotation attributes
        for key, value in _iter_attrs(node):
            if key not in row:
                row[key] = value
        rows.append(row)
    return rows


def _iter_attrs(node):
    ann = node._annotations
    if ann:
        for key, value in ann.items():
            yield key, value


def _copy_attrs(src, dst):
    for key, value in _iter_attrs(src):
        setattr(dst, key, value)


class _ChildrenView:
    """A read-only sequence view of a node's children.

    Prevents accidental mutation (``tree.children.append(x)``) that would
    corrupt the parent/child invariants; use :meth:`Tree.append_child`.
    """
    __slots__ = ("_lst",)

    def __init__(self, lst):
        object.__setattr__(self, "_lst", lst)

    def __len__(self):
        return len(self._lst)

    def __iter__(self):
        return iter(self._lst)

    def __getitem__(self, i):
        return self._lst[i]

    def __contains__(self, x):
        return x in self._lst

    def __reversed__(self):
        return reversed(self._lst)

    def __eq__(self, other):
        if hasattr(other, "__iter__") and not isinstance(other, str):
            return list(self._lst) == list(other)
        return self._lst == other

    def __repr__(self):
        return repr(self._lst)


if __name__ == "__main__":
    tree = Tree("A")
    clade = Tree("B")
    clade.children = [Tree("X"), Tree("Y")]
    tree.children = [clade, Tree("C"), Tree("D")]
    print(tree)