#!/usr/bin/env python3
"""
Code generator for argument parsing headers.

Usage:
    python generate_args.py <input.json> <output.hpp> <prefix>

This script reads argument definitions from a JSON file and generates
C++ header macros for use with the argparse-based argument system.
"""

import json
import sys
from pathlib import Path
from typing import Any


TYPE_MAP = {
    "bool": "bool",
    "int": "int",
    "float": "float",
    "string": "std::string",
    "vector<string>": "std::vector<std::string>",
}


def cpp_type(json_type: str) -> str:
    """Convert JSON type string to C++ type."""
    if json_type not in TYPE_MAP:
        raise ValueError(f"Unknown type: {json_type}")
    return TYPE_MAP[json_type]


def cpp_default_value(json_type: str, default: Any) -> str:
    """Convert JSON default value to C++ literal."""
    if json_type == "bool":
        return "true" if default else "false"
    elif json_type == "int":
        return str(default)
    elif json_type == "float":
        return f"{default}f"
    elif json_type == "string":
        return f'"{default}"'
    elif json_type == "vector<string>":
        if not default:
            return "std::vector<std::string>{}"
        items = ", ".join(f'"{item}"' for item in default)
        return f"std::vector<std::string>{{{items}}}"
    else:
        raise ValueError(f"Unknown type: {json_type}")


def generate_member_declaration(arg: dict) -> str:
    """Generate a single Arg<T> member declaration."""
    ctype = cpp_type(arg["type"])
    return f"Arg<{ctype}> {arg['name']}"


def generate_initializer(arg: dict, is_first_in_group: bool, group_name: str) -> str:
    """Generate a single initializer list entry."""
    name = arg["name"]
    short_name = arg.get("short", "")
    long_name = arg["long"]
    help_text = arg["help"].replace('"', '\\"')

    has_default = "default" in arg
    group = group_name if is_first_in_group else ""

    if has_default:
        # Always pass group param to avoid ambiguity with string defaults
        default_val = cpp_default_value(arg["type"], arg["default"])
        return f'{name}(_parser, "{short_name}", "{long_name}", "{help_text}", {default_val}, "{group}")'
    else:
        # Required argument (no default) - use GroupTag{} to disambiguate
        if group:
            return f'{name}(_parser, "{short_name}", "{long_name}", "{help_text}", GroupTag{{}}, "{group}")'
        else:
            return f'{name}(_parser, "{short_name}", "{long_name}", "{help_text}")'


def generate_header(config: dict, prefix: str) -> str:
    """Generate the complete header file content."""
    upper_prefix = prefix.upper()

    lines = [
        "// AUTO-GENERATED FILE - DO NOT EDIT",
        f"// Generated from JSON configuration for {prefix}",
        "",
        f"#ifndef {upper_prefix}_ARGS_GENERATED_HPP",
        f"#define {upper_prefix}_ARGS_GENERATED_HPP",
        "",
    ]

    # Generate program name macros
    program = config.get("program", {})
    program_name = program.get("name", prefix)

    # Handle conditional suffixes for program name
    if program.get("mpi_suffix") or program.get("debug_suffix"):
        lines.append("// Program name with conditional suffixes")
        lines.append(f"#ifdef MPI_BUILD")
        lines.append(f"#ifdef DEBUG")
        lines.append(f'#define {upper_prefix}_PROGRAM_NAME "{program_name} MPI (DEBUG)"')
        lines.append(f"#else")
        lines.append(f'#define {upper_prefix}_PROGRAM_NAME "{program_name} MPI"')
        lines.append(f"#endif")
        lines.append(f"#else")
        lines.append(f"#ifdef DEBUG")
        lines.append(f'#define {upper_prefix}_PROGRAM_NAME "{program_name} (DEBUG)"')
        lines.append(f"#else")
        lines.append(f'#define {upper_prefix}_PROGRAM_NAME "{program_name}"')
        lines.append(f"#endif")
        lines.append(f"#endif")
    else:
        lines.append(f'#define {upper_prefix}_PROGRAM_NAME "{program_name}"')

    lines.append("")

    # Program version
    if program.get("version_suffix"):
        lines.append(f'#define {upper_prefix}_PROGRAM_VERSION ({upper_prefix}_PROGRAM_NAME " " VERSION)')
    else:
        lines.append(f'#define {upper_prefix}_PROGRAM_VERSION VERSION')

    lines.append("")

    # Collect all arguments
    all_args = []
    for group in config.get("groups", []):
        group_name = group.get("name", "")
        args = group.get("args", [])
        for i, arg in enumerate(args):
            all_args.append({
                "arg": arg,
                "is_first_in_group": (i == 0),
                "group_name": group_name if i == 0 else ""
            })

    # Generate member declarations macro
    lines.append("// Member declarations - include in class body")
    lines.append(f"#define {upper_prefix}_ARG_MEMBERS \\")
    for i, entry in enumerate(all_args):
        decl = generate_member_declaration(entry["arg"])
        if i < len(all_args) - 1:
            lines.append(f"    {decl}; \\")
        else:
            lines.append(f"    {decl};")

    lines.append("")

    # Generate initializer list macro
    lines.append("// Constructor initializer list - use after Program base initializer")
    lines.append(f"#define {upper_prefix}_ARG_INIT \\")
    for i, entry in enumerate(all_args):
        init = generate_initializer(
            entry["arg"],
            entry["is_first_in_group"],
            entry["group_name"]
        )
        if i < len(all_args) - 1:
            lines.append(f"    {init}, \\")
        else:
            lines.append(f"    {init}")

    lines.append("")
    lines.append(f"#endif // {upper_prefix}_ARGS_GENERATED_HPP")
    lines.append("")

    return "\n".join(lines)


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <input.json> <output.hpp> <prefix>", file=sys.stderr)
        sys.exit(1)

    input_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])
    prefix = sys.argv[3]

    if not input_path.exists():
        print(f"Error: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    with open(input_path, "r") as f:
        config = json.load(f)

    header_content = generate_header(config, prefix)

    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        f.write(header_content)

    print(f"Generated {output_path}")


if __name__ == "__main__":
    main()
