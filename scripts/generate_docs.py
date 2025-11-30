#!/usr/bin/env python3
"""
Generate markdown documentation for command line options from JSON config.

Usage:
    python generate_docs.py <input.json>              # Output to stdout
    python generate_docs.py <input.json> <output.md>  # Update file in-place

When updating a file in-place, the script replaces content between:
    <!-- BEGIN AUTO-GENERATED CLI OPTIONS -->
    <!-- END AUTO-GENERATED CLI OPTIONS -->
"""

import json
import sys
import re
from pathlib import Path

BEGIN_MARKER = "<!-- BEGIN AUTO-GENERATED CLI OPTIONS -->"
END_MARKER = "<!-- END AUTO-GENERATED CLI OPTIONS -->"


def format_arg_doc(arg: dict) -> str:
    """Format a single argument as markdown."""
    short_name = arg.get("short", "")
    long_name = arg["long"]
    help_text = arg["help"]
    arg_type = arg["type"]

    # Build the option string
    if short_name:
        option_str = f"`{short_name}, {long_name}`"
    else:
        option_str = f"`{long_name}`"

    # Add default value info if present
    if "default" in arg:
        default = arg["default"]
        if arg_type == "vector<string>":
            if not default:
                default_str = ""
            else:
                default_str = f" (default: {', '.join(default)})"
        elif arg_type == "string":
            if default == "":
                default_str = ""
            else:
                default_str = f" (default: {default})"
        elif arg_type == "int" and default == 2147483647:
            default_str = ""  # Don't show INT_MAX
        else:
            default_str = f" (default: {default})"
    else:
        if arg_type == "bool":
            default_str = ""
        else:
            default_str = " (required)"

    return f"**{option_str}** : {help_text}{default_str}\n"


def generate_docs(config: dict) -> str:
    """Generate markdown documentation from config."""
    lines = ["## Command Line Options\n"]

    for group in config.get("groups", []):
        group_name = group.get("name", "Options")
        args = group.get("args", [])

        if not args:
            continue

        lines.append(f"### {group_name}\n")

        for arg in args:
            lines.append(format_arg_doc(arg))

        lines.append("")  # Blank line between groups

    return "\n".join(lines).rstrip() + "\n"


def update_file(filepath: Path, new_content: str) -> None:
    """Update file by replacing content between markers."""
    with open(filepath, "r") as f:
        content = f.read()

    # Find and replace content between markers
    pattern = re.compile(
        re.escape(BEGIN_MARKER) + r".*?" + re.escape(END_MARKER),
        re.DOTALL
    )

    if not pattern.search(content):
        print(f"Error: Markers not found in {filepath}", file=sys.stderr)
        print(f"Expected: {BEGIN_MARKER} ... {END_MARKER}", file=sys.stderr)
        sys.exit(1)

    new_section = f"{BEGIN_MARKER}\n{new_content}{END_MARKER}"
    updated_content = pattern.sub(new_section, content)

    with open(filepath, "w") as f:
        f.write(updated_content)

    print(f"Updated {filepath}")


def main():
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print(f"Usage: {sys.argv[0]} <input.json> [output.md]", file=sys.stderr)
        sys.exit(1)

    input_path = Path(sys.argv[1])

    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    with open(input_path, "r") as f:
        config = json.load(f)

    docs_content = generate_docs(config)

    if len(sys.argv) == 3:
        output_path = Path(sys.argv[2])
        update_file(output_path, docs_content)
    else:
        print(docs_content)


if __name__ == "__main__":
    main()
