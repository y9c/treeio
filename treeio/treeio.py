# -*- coding: utf-8 -*-

"""Main module -- high level format conversion."""

from typing import Optional

from .io import read, write


def convert_format(input_path, output_path, format_in: Optional[str] = None,
                   format_out: Optional[str] = None):
    """Convert a tree file from one format to another.

    Formats are inferred from the file extensions unless explicitly given.
    """
    tree = read(input_path, format=format_in)
    write(tree, output_path, format=format_out)
    return output_path


def convert_string(tree_string: str, out_format: str = "newick",
                   in_format: Optional[str] = None) -> str:
    """Convert a tree from a string to another format (string)."""
    tree = read(tree_string, format=in_format)
    return write(tree, format=out_format)
