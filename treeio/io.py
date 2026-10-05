#
# Copyright © 2020 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2020-03-27 15:58

"""
Unified tree I/O with automatic format detection.

This is the Python analogue of the R ``treeio`` ``read.tree`` / ``write.tree``
dispatcher (and ``toytree``'s ``toytree.tree`` auto-detection).  It reads and
writes the common tree file formats from a single, guess-free-ish interface:

* newick  (``.nwk``, ``.newick``, ``.tree``)
* nexus   (``.nex``, ``.nexus``)
* json    (``.json``)
* phyloxml (``.xml``)

Formats can be forced with the ``format`` argument or inferred from the file
extension / first character of the payload.

Examples
--------
>>> from treeio import read, write
>>> tree = read("data/animals.nwk")
>>> write(tree, "out/out.nexus")
"""

from __future__ import annotations

import re as _re
from pathlib import Path

from .jplace import read_jplace
from .json_io import read_json, write_json
from .newick import read_newicks, write_newick
from .nexml import read_nexmls, write_nexml
from .nexus import read_nexuses, write_nexus
from .phylip import read_phylips, write_phylip
from .phyloxml import read_phyloxmls, write_phyloxml
from .tree import Tree
from .vendor import read_iqtree, read_mrbayeses, read_raxml

_PHYLIP_HEADER_RE = _re.compile(r"^\s*\d+\s+\d+\s*$")

_FORMATS = {
    "newick": ("nwk", "newick", "tree", "nhx"),
    "nexus": ("nex", "nexus", "nxs"),
    "json": ("json",),
    "phyloxml": ("xml", "phyloxml"),
    "phylip": ("phy", "phylip"),
    "jplace": ("jplace", "jpl"),
    "nexml": ("nexml",),
    "mrbayes": ("mb", "mrbayes", "mrb"),
    "iqtree": ("treefile", "iqtree"),
    "raxml": ("raxml", "besttree", "bipartitions"),
}

# user-registered formats: name -> {"read", "write", "extensions"}
_FORMAT_REGISTRY: dict = {}


def register_format(name: str, read=None, write=None, extensions=()):
    """Register a custom tree format so ``read`` / ``write`` / path detection work.

    ``read`` is ``callable(text) -> list[Tree]``, ``write`` is
    ``callable(tree) -> str``, and ``extensions`` are file suffixes to route to
    this format.  Built-in format names cannot be overwritten.
    """
    if name in _FORMATS:
        raise ValueError(f"cannot overwrite built-in format {name!r}")
    _FORMAT_REGISTRY[name] = {
        "read": read,
        "write": write,
        "extensions": [e.lstrip(".") for e in extensions],
    }
    return name


def unregister_format(name: str):
    """Remove a previously registered custom format."""
    return _FORMAT_REGISTRY.pop(name, None)


def _extension_map() -> dict:
    """Merged extension -> format map (built-in + registered)."""
    mapping = {ext: fmt for fmt, exts in _FORMATS.items() for ext in exts}
    for name, info in _FORMAT_REGISTRY.items():
        for ext in info.get("extensions", []):
            mapping[ext] = name
    return mapping


def detect_format(source: str) -> str:
    """Best-effort detection of the tree format of ``source``."""
    s = source.lstrip()
    if s.startswith("#NEXUS"):
        # only call it MrBayes when there is an actual ``begin mrbayes`` block
        if _re.search(r"BEGIN\s+MRBAYES\s*;", s, _re.IGNORECASE):
            return "mrbayes"
        return "nexus"
    # Determine the root element's *local* name across an optional ``<?xml?>``
    # declaration, leading markup (comments / DOCTYPE) and a namespace prefix
    # (e.g. ``<ns0:nexml>``).  The old ``startswith("<nexml")`` check missed a
    # declaration or a prefixed root, letting NeXML fall through to ``phyloxml``.
    decl = _re.match(r"\s*<\?xml[^>]*\?>", s)
    if decl:
        s = s[decl.end() :].lstrip()
    tag_m = _re.match(r"<([A-Za-z_][\w.-]*:)?([A-Za-z_][\w.-]*)", s)
    root_local = tag_m.group(2) if tag_m else None
    if root_local == "nexml":
        return "nexml"
    if root_local == "phyloxml":
        return "phyloxml"
    if s.startswith(("<?xml", "<")):
        return "phyloxml"
    if s.startswith("{"):
        # JSON may be a plain tree or a Jplace placement document
        return "jplace" if '"tree"' in s else "json"
    if _PHYLIP_HEADER_RE.match(s.split("\n", 1)[0]):
        return "phylip"
    if s.startswith(("(", "[")):
        # [&R] rooted marker -> newick-ish
        return "newick"
    return "newick"


def _format_from_path(path: str) -> str | None:
    ext = Path(path).suffix.lstrip(".").lower()
    return _extension_map().get(ext)


def read(source: str | Path | Tree, format: str | None = None) -> Tree:
    """Read a tree from a file path or raw string.

    If ``source`` is an existing file path it is read from disc; otherwise the
    string is treated as the tree payload itself.  Returns the first tree.
    Passing an already-parsed :class:`Tree` is a no-op (returns it unchanged).
    """
    if isinstance(source, Tree):
        return source
    payload = _read_payload(source)
    fmt = format or _format_from_path(str(source)) or detect_format(payload)
    trees = _read_trees(payload, fmt)
    if not trees:
        raise ValueError(f"no tree found in {source!r} ({fmt})")
    return trees[0]


def read_many(source: str | Path | Tree, format: str | None = None) -> list[Tree]:
    """Read every tree found in ``source`` (file path or raw string).

    Passing a single :class:`Tree`` returns ``[tree]``.
    """
    if isinstance(source, Tree):
        return [source]
    payload = _read_payload(source)
    fmt = format or _format_from_path(str(source)) or detect_format(payload)
    return _read_trees(payload, fmt)


def _read_payload(source: str | Path) -> str:
    s = str(source)
    # Structured data (JSON / newick / ...) is never a path: return it verbatim
    # before any filesystem stat, which would raise on over-long strings.
    stripped = s.lstrip()
    if stripped[:1] in ("{", "[", "(", '"'):
        return s
    if s.endswith((";", "};")):
        return s
    p = Path(source)
    try:
        if s and p.exists() and p.is_file():
            return p.read_text()
    except OSError:
        return s
    return s


def _read_trees(payload: str, fmt: str) -> list[Tree]:
    if fmt == "newick":
        return read_newicks(payload)
    if fmt == "nexus":
        return read_nexuses(payload)
    if fmt == "json":
        return read_json(payload)
    if fmt == "phyloxml":
        return read_phyloxmls(payload)
    if fmt == "phylip":
        return read_phylips(payload)
    if fmt == "jplace":
        return [read_jplace(payload)]
    if fmt == "nexml":
        return read_nexmls(payload)
    if fmt == "mrbayes":
        return read_mrbayeses(payload)
    if fmt == "iqtree":
        return [read_iqtree(payload)]
    if fmt == "raxml":
        return [read_raxml(payload)]
    info = _FORMAT_REGISTRY.get(fmt)
    if info is not None and info.get("read") is not None:
        return list(info["read"](payload))
    raise ValueError(f"unknown format {fmt!r}")


def write(
    tree: Tree | list[Tree],
    target: str | Path | None = None,
    format: str | None = None,
    **kwargs,
) -> str | None:
    """Write a tree to a file or return its serialised string.

    If ``target`` is ``None`` the serialised string is returned.  The format
    is taken from ``format`` or inferred from ``target``'s extension.
    """
    fmt = format or (_format_from_path(str(target)) if target else None) or "newick"
    payload = _serialize(tree, fmt, **kwargs)
    if target is not None:
        Path(target).parent.mkdir(parents=True, exist_ok=True)
        Path(target).write_text(payload)
        return None
    return payload


def _serialize(tree, fmt: str, **kwargs) -> str:
    if fmt == "newick":
        return write_newick(tree, **kwargs)
    if fmt == "nexus":
        return write_nexus(tree)
    if fmt == "json":
        return write_json(tree)
    if fmt == "phyloxml":
        return write_phyloxml(tree)
    if fmt == "phylip":
        return write_phylip(tree)
    if fmt == "nexml":
        return write_nexml(tree)
    info = _FORMAT_REGISTRY.get(fmt)
    if info is not None and info.get("write") is not None:
        return info["write"](tree)
    raise ValueError(f"unknown format {fmt!r}")


__all__ = [
    "detect_format",
    "read",
    "read_many",
    "register_format",
    "unregister_format",
    "write",
]
