"""Tests for the starter command and shared logging."""

import subprocess
import sys
import unittest


class CommandTests(unittest.TestCase):
    def test_module_entry_point_and_logging(self):
        result = subprocess.run(
            [sys.executable, "-m", "python_template", "--name", "Ada"],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "Hello, Ada!\n")
        self.assertIn("INFO python_template.cli: Generating greeting", result.stderr)

    def test_help(self):
        result = subprocess.run(
            [sys.executable, "-m", "python_template", "--help"],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--name", result.stdout)
        self.assertEqual(result.stderr, "")
