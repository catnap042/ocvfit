# Physics-based (P2D/SPM) simulation or rate-dependent fitting

Adding a pseudo-two-dimensional or single-particle model, or fitting higher-rate
curves with transport physics, is out of scope.

## Why this is out of scope

ocvfit does electrode alignment of quasi-equilibrium curves. Its strength is that
it has few parameters and that each of them can be checked for identifiability.
Transport models bring many more parameters that OCV data cannot constrain, and
mature open-source simulators already cover them. Use one of those and pass it the
electrode state that ocvfit reports.
