# Hosted service or usage telemetry

A web upload service, a hosted API, or anonymous usage statistics are out of scope.

## Why this is out of scope

Battery test data is often confidential. A tool that never opens a network
connection can be used in places where an upload is not allowed, and it is easy to
audit. Telemetry would break that promise even if it were opt-in.
