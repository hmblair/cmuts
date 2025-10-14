
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
and an example can be found under `./examples/visualize`.

The FASTA file is necessary to perform alignment of the sequence which was probed against the sequence of the 3D structure, which may be different e.g. due to flanking sequences.

## Command Line Options

### Core Options

**`--bin`** : Path to the ChimeraX executable (default: ChimeraX)

### Selection Options

**`--index`** : Index of the sequence in the FASTA and the reactivity in the HDF5 file (default: 0)

**`--chain`** : Color this chain of the structure (default: ALL)

**`--trim-5p`** : Trim this many bases from the 5' end of the reactivity data (default: 0)

**`--trim-3p`** : Trim this many bases from the 3' end of the reactivity data (default: 0)

### Output Options

**`--color`** : The ChimeraX color to use for highly reactive bases (default: indianred)

## Movies

Once the structure is loaded into ChimeraX, you can record a movie of it rotating by running

```
movie record
turn y 0.5 720
movie stop
movie encode PATH-TO-MOVIE
```

