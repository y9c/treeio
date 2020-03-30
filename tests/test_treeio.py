#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `treeio` package."""


import unittest
from click.testing import CliRunner

from treeio import cli


class TestCommandLine(unittest.TestCase):
    def test_command_line_interface(self):
        """Test the CLI."""
        runner = CliRunner()
        result = runner.invoke(cli.cli)
        self.assertEqual(result.exit_code, 0)
        help_result = runner.invoke(cli.cli, ["--help"])
        self.assertEqual(help_result.exit_code, 0)


if __name__ == "__main__":
    unittest.main()
