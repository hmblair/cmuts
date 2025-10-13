
## Purpose

`cmuts visualize` may be used to overlay a reactivity profile generated via `cmuts` onto a 3D RNA structure.

## Usage

The general syntax is
```bash
cmuts visualize \
  --file "$PROFILES" \
  --dataset "$NAME" \
  --fasta "$FASTA" \
  --cif "$CIF"
```
The FASTA file is necessary to perform alignment of the sequence which was probed against the sequence of the 3D structure, which may be different e.g. due to flanking sequences.
