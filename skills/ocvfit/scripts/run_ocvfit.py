#!/usr/bin/env python3
"""Fixed entry point for the ocvfit skill.

Forwards all arguments to the installed ``ocvfit`` command-line tool, e.g.::

    python run_ocvfit.py report --pe lfp --ne graphite --bol bol.csv --aged aged.csv --out results

It contains no analysis code of its own: all settings are frozen inside ocvfit.
"""

import sys

try:
    from ocvfit.cli import main
except ImportError:
    sys.stderr.write(
        'ocvfit is not installed. Install it with:\n  pip install "ocvfit[plot] @ '
        'git+https://github.com/catnap042/ocvfit.git"\n'
    )
    sys.exit(1)

sys.exit(main())
