# -*- coding: utf-8 -*-

r"""
Console script for treeio.

test latex math

$$\pi \times x = \frac{x}{y+z}$$
"""

import click
from .treeio import convert_format


@click.group()
@click.option("--debug/--no-debug", default=False)
def cli(debug=False):
    """Entry for command group."""
    click.echo("Debug mode is %s" % ("on" if debug else "off"))


@cli.command()
@click.option(
    "--input", "-i", "input_path", required=True, help="Path of input file."
)
@click.option(
    "--output",
    "-o",
    "output_path",
    required=True,
    help="Path of output file.",
)
def convert(input_path, output_path):
    """Convert tree formats."""
    click.echo(f"input path {input_path}!")
    click.echo(f"input path {output_path}!")
    convert_format(input_path, output_path)


@cli.command()
@click.option("--count", default=1, help="Number of greetings.")
@click.option(
    "--name", prompt="What is your name?\n", help="The person to greet."
)
def hello(count, name):
    """Program that greets NAME for a total of COUNT times."""
    for _ in range(count):
        click.echo("Hello %s!" % name)


if __name__ == "__main__":
    cli()
