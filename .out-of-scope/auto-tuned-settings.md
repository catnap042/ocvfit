# Auto-tuning settings per dataset

Options that search smoothing windows, weights or loss functions per dataset to get
the lowest RMSE are out of scope.

## Why this is out of scope

The same input must give the same output in every release, so that two people and
an AI agent all get the same answer. Per-dataset tuning makes results depend on the
search and invites "nice" numbers. Frozen settings change only through a versioned
CHANGELOG entry.
