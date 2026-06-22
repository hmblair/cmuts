#!/usr/bin/env python3
"""
Inject the top-level `cmuts -h` output into a markdown file, so the docs stay
in sync with the cmuts dispatcher's print-usage (the single source of truth).

The help text is read from stdin and written between:
    <!-- BEGIN AUTO-GENERATED USAGE -->
    <!-- END AUTO-GENERATED USAGE -->

Usage:
    cmuts -h | python generate_usage.py <output.md>
"""

import re
import sys
from pathlib import Path

BEGIN_MARKER = "<!-- BEGIN AUTO-GENERATED USAGE -->"
END_MARKER = "<!-- END AUTO-GENERATED USAGE -->"


def main():
    if len(sys.argv) != 2:
        print(f"Usage: cmuts -h | {sys.argv[0]} <output.md>", file=sys.stderr)
        sys.exit(1)

    help_text = sys.stdin.read().rstrip("\n")
    block = f"```text\n{help_text}\n```"

    filepath = Path(sys.argv[1])
    content = filepath.read_text()

    pattern = re.compile(
        re.escape(BEGIN_MARKER) + r".*?" + re.escape(END_MARKER),
        re.DOTALL,
    )
    if not pattern.search(content):
        print(f"Error: USAGE markers not found in {filepath}", file=sys.stderr)
        print(f"Expected: {BEGIN_MARKER} ... {END_MARKER}", file=sys.stderr)
        sys.exit(1)

    new_section = f"{BEGIN_MARKER}\n{block}\n{END_MARKER}"
    filepath.write_text(pattern.sub(new_section, content))
    print(f"Updated {filepath}")


if __name__ == "__main__":
    main()
