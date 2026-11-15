"""Example command-line entry point."""

import argparse

from molecular_utilities.logging import configure_logging, get_logger

NAMESPACE = "python_template"
logger = get_logger("cli", namespace=NAMESPACE)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--name", default="world", help="name to greet")
    args = parser.parse_args(argv)
    configure_logging(namespace=NAMESPACE)
    logger.info("Generating greeting")
    print(f"Hello, {args.name}!")
    return 0
