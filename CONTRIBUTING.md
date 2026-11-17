# Contributing to TLS Finder

Thank you for your interest in TLS Finder. Bug reports, questions, feature
requests and code contributions are all welcome.

## Reporting bugs and asking questions

Open an issue at <https://github.com/AAKoksoy/TLS-Finder/issues>. For bugs,
please include:

- the TLS Finder version (`tls-finder --version`) and Python version,
- the command or code you ran and the parameters used,
- the full error message, and
- if possible, a small example dataset (synthetic data from
  `tls_finder.make_synthetic_tissue()` is ideal) that reproduces the problem.

Please do not upload identifiable patient data.

## Contributing code

1. Fork the repository and create a branch from `main`.
2. Install in development mode: `pip install -e ".[test]"`.
3. Make your change and add tests for it in `tests/`.
4. Run `pytest` and make sure all tests pass.
5. Add a line to the "Unreleased" section of `CHANGELOG.md`.
6. Open a pull request describing what you changed and why.

Changes to the detection logic must keep the regression tests against the
original implementation (`tests/reference.py`) passing, or clearly explain
in the pull request why results should change.

## Support and maintenance

TLS Finder is maintained by Ayse A. Koksoy. Please use GitHub issues rather
than email so that answers are visible to other users.

## Code of conduct

Please be respectful and constructive in all interactions.
