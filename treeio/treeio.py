# -*- coding: utf-8 -*-

"""Main module."""

from pathlib import Path
from .newick import read_newick, write_newick


def convert_format(input_path, output_path):
    """Convert format."""
    input_string = Path(input_path).read_text()
    tree = read_newick(input_string)
    output_string = write_newick(tree)
    Path(output_path).write_text(output_string)
