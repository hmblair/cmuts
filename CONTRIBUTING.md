# Contributing to cmuts

Thank you for your interest in contributing to cmuts! This document provides guidelines for contributing code to the project.

## Code Style Guide

### C++ Conventions

#### Naming

| Element | Convention | Example |
|---------|------------|---------|
| Classes | PascalCase | `DataStream`, `BamFile`, `CramIterator` |
| Functions | snake_case | `read_bgzf()`, `get_file()` |
| Internal functions | Leading underscore | `_exists()`, `_throw_if_not_exists()` |
| Private members | Leading underscore | `_data`, `_pos`, `_remaining` |
| Type aliases | `_t` suffix | `base_t`, `seq_t`, `qual_t` |
| Constants | UPPER_SNAKE_CASE | `BASES`, `MAX_MAPQ`, `BGZF_BUFFER` |
| Template parameters | Single uppercase or descriptive | `T`, `dtype` |

#### Code Organization

- Use section headers to organize code:
  ```cpp
  //
  // Section Name
  //
  ```

- Group includes by category:
  1. Project headers (`#include "common.hpp"`)
  2. Standard library (`#include <vector>`)
  3. External dependencies (`#include <htslib/sam.h>`)

#### Formatting

- Use `.clang-format` for consistent formatting
- Run `clang-format -i src/cpp/*.cpp src/cpp/*.hpp` before committing
- 4-space indentation (no tabs)
- 100-character line limit

#### Error Handling

- Use `CMUTS_THROW(msg)` macro for errors with location info
- Use `CMUTS_THROW_IF(cond, msg)` for conditional errors
- Provide descriptive error messages with context

```cpp
// Good
CMUTS_THROW("Failed to read " + std::to_string(size) + " bytes from " + filename);

// Avoid
throw std::runtime_error("read failed");
```

### Python Conventions

#### Naming

| Element | Convention | Example |
|---------|------------|---------|
| Classes | PascalCase | `ProbingData`, `NormScheme` |
| Functions | snake_case | `normalize()`, `get_reactivity()` |
| Constants | UPPER_SNAKE_CASE | `IX_DEL`, `BOLD` |
| Private | Leading underscore | `_internal_method()` |

#### Type Hints

- Use type hints for all public functions
- Use `from __future__ import annotations` for forward references
- Prefer `list[T]` over `List[T]` (Python 3.9+)

```python
def normalize(
    file: str,
    fasta: str,
    opts: Opts
) -> ProbingData:
    ...
```

#### Documentation

- Use Google-style docstrings for public functions
- Include Args, Returns, and Raises sections

```python
def normalize(file: str, fasta: str, opts: Opts) -> ProbingData:
    """Normalize chemical probing data.

    Args:
        file: Path to HDF5 input file.
        fasta: Path to reference FASTA file.
        opts: Normalization options.

    Returns:
        Normalized probing data.

    Raises:
        ValueError: If file format is invalid.
    """
```

## Development Setup

### Building from Source

```bash
# Clone with submodules
git clone --recursive https://github.com/hmblair/cmuts.git
cd cmuts

# Build C++ components
./configure
cmake --build build --parallel

# Install Python package in development mode
pip install -e ".[dev]"
```

### Running Tests

```bash
# C++ tests (after building with BUILD_TESTING=ON)
ctest --test-dir build --output-on-failure

# Python tests
pytest tests/python -v
```

### Code Formatting

```bash
# C++ formatting
clang-format -i src/cpp/*.cpp src/cpp/*.hpp

# Python formatting and linting
ruff check src/python/cmuts
ruff format src/python/cmuts

# Type checking
mypy src/python/cmuts
```

## Pull Request Process

1. Create a feature branch from `dev`
2. Make your changes following the style guide
3. Add tests for new functionality
4. Ensure all tests pass
5. Update documentation if needed
6. Submit a pull request to `dev`

## License

By contributing to cmuts, you agree that your contributions will be licensed under the MIT License.
