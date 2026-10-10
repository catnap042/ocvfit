# Security policy

## Supported versions

Only the latest release (and the `main` branch) receives fixes.

## Scope

ocvfit is an offline numerical library with a command-line tool and an optional local MCP
server (stdio). It makes no network connections, runs no shell commands and does not use
`pickle` or `eval`. Relevant reports include, for example: reading or writing files outside
the paths given to a command or MCP tool, bypassing `OCVFIT_MCP_ROOT`, crashes or runaway
run times from crafted input files, or a vulnerable dependency pin.

## Reporting a vulnerability

Please report vulnerabilities privately through GitHub: open the repository's
**Security** tab and choose **Report a vulnerability**. Do not open a public issue for a
security problem. If private reporting is not available, open an issue that only asks for
a private contact channel, without technical details.

You can expect an acknowledgement within a week. Fixes are released as a new version and
noted in the [changelog](CHANGELOG.md).
