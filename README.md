## Overview

`cmuts` is a program for counting mutations and computing reactivity profiles in MaP-seq experiments. It features
* Fast, compiled C++ code with native multithreading support
* Streamed IO and direct output to compressed HDF5 files
* Handling of arbitrary-length ambiguous deletions, including mutation-informed deletion spreading

In a picture:

![cmuts-overview](./docs/figures/overview.png)

Detailed documentation can be found at [hmblair.github.io/cmuts](https://hmblair.github.io/cmuts).
