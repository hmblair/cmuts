## Purpose

`cmuts generate` creates synthetic test data for validating the `cmuts` pipeline. It generates random reference sequences and simulated aligned reads with known modifications, allowing you to verify that `cmuts core` correctly identifies mutations, insertions, and deletions.

## Usage

The basic syntax is:

```bash
cmuts generate \
  --length 100 \
  --references 10 \
  --queries 50 \
  --out-fasta references.fa \
  --out-sam alignments.sam \
  --out-h5 expected.h5 \
  --min-mapq 10 \
  --min-phred 10 \
  --max-length 100 \
  --max-indel-length 10 \
  --quality-window 2
```

This generates:

1. **Reference sequences** (`--out-fasta`): Random nucleotide sequences of the specified length
2. **Aligned reads** (`--out-sam`): Simulated alignments with random modifications
3. **Expected counts** (`--out-h5`): Ground truth modification counts for validation

You can then run `cmuts core` on the generated data and compare its output to the expected counts to verify correctness.

## Example Workflow

```bash
# Generate test data
cmuts generate \
  --length 200 \
  --references 5 \
  --queries 100 \
  --out-fasta test_refs.fa \
  --out-sam test_reads.sam \
  --out-h5 expected_counts.h5 \
  --min-mapq 10 \
  --min-phred 10 \
  --max-length 200 \
  --max-indel-length 10 \
  --quality-window 2

# Run cmuts core on the generated data
cmuts core \
  -f test_refs.fa \
  -o actual_counts.h5 \
  test_reads.sam

# Compare actual_counts.h5 with expected_counts.h5
```

<!-- BEGIN AUTO-GENERATED CLI OPTIONS -->
## Command Line Options

### Generation

**`--length`** : The length of the reference sequences. (required)

**`--queries`** : The number of queries to generate per reference. (required)

**`--references`** : The number of references to generate. (required)

**`--seed`** : Random seed for reproducibility. If not provided, uses time-based seed. (default: -1)


### Output files

**`--out-fasta`** : The file to store the references in. (required)

**`--out-sam`** : The file to store the queries in. (required)

**`--out-h5`** : The file to store the expected modifications in. (required)


### Quality thresholds

**`--min-mapq`** : The minimum quality to consider a read. (required)

**`--min-phred`** : The minimum quality to consider a base. (required)

**`--min-length`** : The smallest query sequences to consider when counting modifications. (default: 2)

**`--max-length`** : The longest query sequences to consider when counting modifications. (required)

**`--max-indel-length`** : Skip indels longer than this. (required)

**`--collapse`** : The minimum number of bases between modifications to consider them consecutive. (default: 2)

**`--quality-window`** : The number of neighbouring bases to consider when calculating PHRED scores. (required)


### Mutation type filters

**`--no-mismatches`** : Do not count mismatches as modifications.

**`--no-insertions`** : Do not count insertions as modifications.

**`--no-deletions`** : Do not count deletions as modifications.
<!-- END AUTO-GENERATED CLI OPTIONS -->
