r"""
Console script for treeio.

Examples
--------
Convert between formats::

    treeio convert -i tree.nwk -o tree.nex

Print a tree to the terminal::

    treeio show -i tree.nwk

Render a tree to SVG::

    treeio plot -i tree.nwk -o tree.svg --layout rectangular
"""

from __future__ import annotations

from pathlib import Path

import click

from .io import read, read_many, write
from .plot import render
from .show import tree_to_ascii


@click.group(invoke_without_command=True)
@click.pass_context
def cli(ctx: click.Context) -> None:
    """Parse, transform, display and plot phylogenetic trees."""
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


@cli.command()
@click.option("--input", "-i", "input_path", required=True, help="Path of input file.")
@click.option(
    "--output", "-o", "output_path", required=True, help="Path of output file."
)
@click.option("--format-in", "fmt_in", default=None, help="Force input format.")
@click.option("--format-out", "fmt_out", default=None, help="Force output format.")
def convert(
    input_path: str,
    output_path: str,
    fmt_in: str | None = None,
    fmt_out: str | None = None,
) -> None:
    """Convert tree formats (newick, nexus, json, phyloxml)."""
    try:
        tree = read(input_path, format=fmt_in)
        write(tree, output_path, format=fmt_out)
    except Exception as err:
        raise click.ClickException(str(err)) from err
    click.echo(f"Wrote {output_path}")


@cli.command()
@click.option("--input", "-i", "input_path", required=True, help="Path of input file.")
@click.option("--index", default=0, help="Tree index (0-based) for multi-tree files.")
def show(input_path: str, index: int) -> None:
    """Render a tree as ASCII art in the terminal."""
    try:
        tree = read_many(input_path)[index]
    except Exception as err:
        raise click.ClickException(str(err)) from err
    click.echo(tree_to_ascii(tree, False, True))


@cli.command("plot")
@click.option("--input", "-i", "input_path", required=True, help="Path of input file.")
@click.option(
    "--output",
    "-o",
    "output_path",
    default=None,
    help="Path of output file (stdout if omitted and ascii).",
)
@click.option(
    "--layout",
    "layout",
    default="rectangular",
    type=click.Choice(
        [
            "rectangular",
            "roundrect",
            "circular",
            "fan",
            "radial",
            "slanted",
            "unrooted",
            "time",
        ]
    ),
)
@click.option(
    "--backend", "backend", default="mpl", type=click.Choice(["mpl", "ascii"])
)
@click.option("--no-tip-labels", is_flag=True, default=False)
@click.option("--no-tip-points", is_flag=True, default=False)
@click.option("--support", is_flag=True, default=False, help="Show support values.")
@click.option("--width", default=800, type=int)
@click.option("--height", default=600, type=int)
def plot_tree(
    input_path: str,
    output_path: str | None,
    layout: str,
    backend: str,
    no_tip_labels: bool,
    no_tip_points: bool,
    support: bool,
    width: int,
    height: int,
) -> None:
    """Render a tree (mpl figure / ascii text)."""
    try:
        tree = read(input_path)
        opts: dict[str, object] = {
            "layout": layout,
            "width": width,
            "height": height,
            "tip_labels": not no_tip_labels,
            "tip_points": not no_tip_points,
            "node_support": support,
        }
        if backend == "ascii":
            out = tree_to_ascii(tree, False, True)
            if output_path:
                Path(output_path).write_text(out)
                click.echo(f"Wrote {output_path}")
            else:
                click.echo(out)
            return
        if not output_path:
            output_path = "tree.png"
        render(tree, backend="mpl", path=output_path, **opts)
        click.echo(f"Wrote {output_path}")
    except Exception as err:
        raise click.ClickException(str(err)) from err


if __name__ == "__main__":
    cli()
