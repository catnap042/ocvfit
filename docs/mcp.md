[README](../README.md) · [Skill](skill.md) · [MCP](mcp.md) · [CLI](cli.md) · [Python API](python-api.md) · [Examples](examples.md) · [Method](method.md)

# MCP server

An optional [Model Context Protocol](https://modelcontextprotocol.io) server exposes the
same commands as tools over stdio. It uses the official `mcp` Python SDK (1.x `FastMCP`
or 2.x `MCPServer`).

```bash
pip install "ocvfit[mcp] @ git+https://github.com/catnap042/ocvfit"
ocvfit-mcp   # stdio server; normally started by the client
```

| tool | CLI equivalent |
|---|---|
| `check_data` | `ocvfit check` |
| `fit_cell` | `ocvfit fit` |
| `degradation_modes` | `ocvfit modes` |
| `differential_analysis` | `ocvfit ica` |
| `full_report` | `ocvfit report` |
| `literature_references` | `ocvfit references` |

Tools take file paths (use absolute paths; relative ones resolve against the server's
working directory) and an `out_dir`, and return the same JSON document as the CLI plus
`absolute_files` for the figures. The server reads only the files it is given, writes
only into `out_dir`, makes no network connections and runs no shell commands. To confine
it to one directory tree, set `OCVFIT_MCP_ROOT`; input paths and `out_dir` outside that
directory are then refused (in the config below: `"env": {"OCVFIT_MCP_ROOT": "/path/to/data"}`).

Claude Code:

```bash
claude mcp add --transport stdio --scope user ocvfit -- ocvfit-mcp
```

Claude Desktop (`claude_desktop_config.json`) and Cursor (`~/.cursor/mcp.json` or the
project's `.cursor/mcp.json`) use the same `mcpServers` block:

```json
{
  "mcpServers": {
    "ocvfit": {
      "command": "/path/to/venv/bin/ocvfit-mcp",
      "args": []
    }
  }
}
```

Use the absolute path of `ocvfit-mcp` (`which ocvfit-mcp`) when the client does not
inherit your shell's `PATH`.
