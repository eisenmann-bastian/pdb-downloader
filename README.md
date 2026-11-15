# Python Project Template

A minimal Python starter with a command-line entry point, tests, examples, and
shared logging through [molecular-utilities](https://github.com/eisenmann-bastian/molecular-utilities).
Replace this introduction with your project's purpose after creating a repository.

## Contents

- [About](#about)
- [Getting started](#getting-started)
- [Usage](#usage)
- [Development](#development)
- [Roadmap](#roadmap)
- [Contributing](#contributing)
- [Contact and support](#contact-and-support)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## About

The starter provides an installable `python_template` package, a small CLI,
and an example using the same logging utilities as the molecular projects.
Its example command prints a greeting; replace it with your own functionality.

## Getting started

Requires Python 3.10 or newer, Git, and network access to install the GitHub
dependency.

```sh
git clone https://github.com/eisenmann-bastian/python-template.git
cd python-template
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -e .
```

For a new project created from this template:

1. Rename the project, Python package and console script in `pyproject.toml`.
2. Update the package imports, module entry point, tests, example scripts and
   the logging namespace consistently.
3. Replace the README description and links with your project information and
   review the copyright notice in `LICENSE`.
4. Run the tests and example below.

## Usage

```sh
python-template --name "world"
python -m python_template --name "world"
```

The command prints the greeting to standard output. Diagnostic messages use
`molecular-utilities` and are sent to standard error. The application entry
point configures logging once; library code obtains a namespaced child logger.

## Development

Run the test suite with an installed environment:

```sh
python -m unittest discover -s tests -v
```

The scripts in `examples/` follow the molecular projects' workflow:

```sh
./examples/setup-venv.sh          # Create and install into .venv
./examples/run-example.sh         # Create .venv if absent, then run the CLI
./examples/run-tests.sh           # Uses PYTHON or python3
```

Set `PYTHON` to choose the interpreter used by the setup and test scripts.
The setup script refuses to overwrite an existing environment; the example
runner uses an existing `.venv` without modifying it.

## Roadmap

- [ ] Describe the goals and milestones for your project.

## Contributing

Open an issue to discuss a change, or fork the repository and submit a pull
request. Include tests and update the documentation when behavior changes.

## Contact and support

- GitHub: [@eisenmann-bastian](https://github.com/eisenmann-bastian)
- LinkedIn: [Bastian Eisenmann](https://www.linkedin.com/in/bastianeisenmann/)
- [Buy me a coffee](https://buymeacoffee.com/bastian.eisenmann)

## License

Distributed under the [BSD 3-Clause License](LICENSE).

## Acknowledgments

- README structure inspired by [Best-README-Template](https://github.com/othneildrew/Best-README-Template).
- Shared logging provided by [molecular-utilities](https://github.com/eisenmann-bastian/molecular-utilities).
