# Bundling real or proprietary cell data

Requests to add measured curves from real cells to the examples, the tests or a
download helper are out of scope.

## Why this is out of scope

Measured data usually belongs to someone else, such as a company, a customer or a
lab, and its licence is rarely clear. The project also does not need it. Synthetic
curves with a known answer test the method more strictly than measured curves, where
the truth is unknown. Links to public datasets with an explicit licence are fine and
are listed in `docs/python-api.md`.

If a bug only shows up on your own data, reduce it to a synthetic or public example
before you file it.
