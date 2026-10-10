# Contributing to ocvfit

Thanks for your interest. Bug reports, questions and pull requests are welcome on
[GitHub Issues](https://github.com/catnap042/ocvfit/issues).

## Reporting a problem

Include the ocvfit version (`ocvfit --version`), Python version and OS, the exact command
or code, and the `status`, `errors` and `warnings` from the JSON output. If you can, attach
a small synthetic CSV that reproduces the problem. Do not post measured data you are not
allowed to share.

## Development

1. Fork and create a feature branch.
2. `pip install -e ".[test,mcp]" "ruff>=0.16,<0.17"`
3. Add tests for new behaviour. Fitting changes need a parameter-recovery test on
   synthetic data with known truth.
4. Run `ruff check src tests examples skills`, `ruff format src tests examples skills`
   and `pytest`.
5. Do not add measured data unless it is openly licensed and attributed. Prefer
   synthetic data from published OCP functions.
6. The CLI/MCP settings in `ocvfit.workflow.SETTINGS` are frozen on purpose: changing
   them changes published results, so it needs a changelog entry and a version bump.

Add a line under "Unreleased" in [CHANGELOG.md](CHANGELOG.md) for user-visible changes.
By contributing you agree that your contribution is licensed under the MIT license.
