# `cmuts hmm`

## Purpose

Computing [mismatch](cmuts-hmm.md#mismatches), [insertion](cmuts-hmm.md#insertions), [deletion](cmuts-hmm.md#deletions), and [termination](cmuts-hmm.md#terminations) rates from alignment files.

## Requires

- One or more coordinate-sorted alignment files of single-end or merged reads. SAM, BAM, and CRAM formats are all supported
- The FASTA library

Alignments must be single-end to avoid double-counting; a read carrying the paired flag causes the program to exit early. Please merge reads upstream before passing data to `cmuts hmm`.

## Model Parameters

The pair HMM is configured by five internal parameters. `--dump-params` writes them in the form `--params` reads.

```sh
cmuts hmm --dump-params > params.txt
```

You may specify a subset of the parameters to modify only them.

## Bandwidth

The `--band` argument sets the bandwidth of the pair HMM, which determines the length of the longest deletion the HMM can relocate. The runtime grows linearly with the bandwidth, on top of a fixed cost. At zero, the HMM degenerates and ambiguous deletions are counted wherever the aligner places them.

DMS data is insensitive to the bandwidth, since methylations mostly cause mismatches. On the other hand, SHAPE data requires a bandwidth of at least the default of 2 to locate ambiguous deletions correctly. Larger bandwidths lead to rapidly diminishing returns.

## Sequencing Errors

Illumina reports a Phred quality per base, which is the probability that the sequencer introduced a substitution error at that base. Substitution errors serve as an alternative explanation for substitutions otherwise attributed to the reverse transcriptase.

Ultima reports a Phred quality per homopolymer run, which is the probability that the run was called one base too long or too short. Run-length errors serve as alternative explanations for insertions and deletions otherwise attributed to the reverse transcriptase. Since Ultima does not report substitution error probabilities, `--ultima-phred` sets one for every base.

A file may hold both kinds of read. `cmuts hmm` selects the error handling of each read from its tags.

The run-length tags are `tp`, a signed byte per base, and `t0`, a string of one score per base. They are read in the order of the stored sequence, which is the order Ultima's aligner writes with `--tp-output-order=reference`. A tool that copies the tags from the unaligned reads, as `minimap2 -y` does, leaves them in the order the bases were sequenced on a reverse-strand read, and `--ultima-tag-order sequenced` reads them in that order instead. `cmuts align` writes them in the order of the stored sequence, so its output needs no option. A warning is printed when the command line on an `@PG` line gives a `--tp-output-order` that disagrees with the order chosen, and when primary reads whose `RG` tag names a read group with a flow order carry no usable `tp` tag.

## Error Checking

`cmuts hmm` verifies the FASTA against the alignment header, by comparing each sequence's name and length, and its MD5 checksum when present. Any mismatch between the header and the FASTA ends the run early. This behavior is configurable via the `--verify` flag.

```{note}
The length check is required to avoid buffer overflows and cannot be disabled.
```

Both the DNA and the RNA forms of the FASTA sequence are accepted during the checksum comparison.

## Rejected Reads

Reads which do not pass the configured filters are rejected, meaning they are not processed by the HMM and do not contribute to the final rates. A secondary alignment, a read with MAPQ 255, a read missing a sequence, or a read with no CIGAR is automatically rejected.

```{warning}
`minimap2` and `bowtie2` assign MAPQ scores differently. The default threshold of 20 is tuned for `minimap2`-aligned data, but often rejects high-quality reads for `bowtie2`-aligned data. A lower threshold of 10 is recommended for `bowtie2`.
```

A read to which the HMM cannot assign a nonzero alignment probability is also rejected; this can only occur if you set one or more transition probabilities to zero.

Unmapped reads are rejected by nature of having no reference with which to compute each rate.

## Rate Computation

For each base $b$, the pair HMM computes five expectations under the posterior distribution over alignments:

| Symbol | Expectation |
| --- | --- |
| $C_b$ | coverage |
| $T_b$ | terminations |
| $M_b$ | mismatches |
| $I_b$ | insertion openings |
| $D_b$ | deletion openings |

Each of the four rates computed by `cmuts hmm` is the ratio of one expectation to the number of reads that could have produced the event.

### Mismatches

The mismatch rate is the probability that the wrong base was incorporated at base $b$, over every read pairing it.

$$
r_{\mathrm{M},b} = \frac{M_b}{C_b}
$$

### Insertions

The insertion opening rate is the probability that one or more bases were inserted directly 5' of base $b$. A read that terminates at base $b$ cannot open an insertion past it, so the denominator is the reads that pair base $b$ and continue 5'.

$$
r_{\mathrm{I},b} = \frac{I_b}{C_b - T_b}
$$

### Deletions

The deletion opening rate is the probability that base $b$ is the 3'-most base of a run of one or more deleted bases. A deletion at base $b$ is opened from base $b+1$, so the denominator is the reads that pair base $b+1$ and continue 5'.

$$
r_{\mathrm{D},b} = \frac{D_b}{C_{b+1} - T_{b+1}}
$$

### Terminations

The termination hazard rate is the probability that a read ends at base $b$, given that it reached it.

$$
r_{\mathrm{T},b} = \frac{T_b}{C_b}
$$

## Error Computation

Each rate is accompanied by a per-base estimate for its error, computed as the binomial standard error of the rate at the depth of its denominator. This is purely the statistical error introduced by finite read depths; it does not account for experimental or systemic errors.

## Output

<!-- BEGIN GENERATED cmuts-hmm DATASETS -->
| Dataset | Source |
| --- | --- |
| [`mismatches/rate`](format.md#mismatchesrate) | Estimated by the HMM. |
| [`mismatches/error`](format.md#mismatcheserror) | The binomial standard error of the rate at the position's depth. |
| [`insertions/rate`](format.md#insertionsrate) | Estimated by the HMM. |
| [`insertions/error`](format.md#insertionserror) | The binomial standard error of the rate at the position's depth. |
| [`deletions/rate`](format.md#deletionsrate) | Estimated by the HMM. |
| [`deletions/error`](format.md#deletionserror) | The binomial standard error of the rate at the position's depth. |
| [`terminations/rate`](format.md#terminationsrate) | Estimated by the HMM. |
| [`terminations/error`](format.md#terminationserror) | The binomial standard error of the rate at the position's depth. |
| [`coverage`](format.md#coverage) | Estimated by the HMM. |
| [`sequence`](format.md#sequence) | Tokenized from the FASTA. |
| [`reads/length-histogram`](format.md#readslength-histogram) | The read length reported in the alignment. |
| [`reads/quality-histogram`](format.md#readsquality-histogram) | The mapping quality reported in the alignment. |
| [`reads/primary/counted`](format.md#readsprimarycounted) | The number of primary reads the HMM successfully processed. |
| [`reads/primary/rejected`](format.md#readsprimaryrejected) | The number of primary reads rejected by a filter or by the HMM. |
| [`reads/supplementary/counted`](format.md#readssupplementarycounted) | The number of supplementary reads the HMM successfully processed. |
| [`reads/supplementary/rejected`](format.md#readssupplementaryrejected) | The number of supplementary reads rejected by a filter or by the HMM. |
| [`reads/unmapped`](format.md#readsunmapped) | The number of unmapped reads in the alignment. |
| [`pairwise/correlation`](format.md#pairwisecorrelation) | Estimated by the HMM. |
| [`pairwise/conditional`](format.md#pairwiseconditional) | Estimated by the HMM. |
| [`pairwise/coverage`](format.md#pairwisecoverage) | Estimated by the HMM. |
<!-- END GENERATED cmuts-hmm DATASETS -->

## CLI Options

<!-- BEGIN GENERATED cmuts-hmm OPTIONS -->
### Arguments

| Argument | Description |
| --- | --- |
| `BAM...` | coordinate-sorted alignments |

### Input and output

| Option | Description |
| --- | --- |
| `-f, --fasta FASTA` | reference sequences (required) |
| `-o, --output HDF5` | write results to this file (required) |
| `--overwrite` | replace the output file if it already exists |
| `--verify CHECKS` | identity checks to make against the FASTA (name, checksum, none; default name,checksum) |

### Filtering

| Option | Description |
| --- | --- |
| `-q, --min-mapq N` | reject alignments below this mapping quality (0 to 254; default 20) |
| `--min-length N` | reject reads shorter than this (default: no limit) |
| `--max-length N` | reject reads longer than this (default: no limit) |
| `-s, --strand STRANDS` | reject alignments not on these strands (forward, reverse; default forward,reverse) |
| `--alignment-type TYPES` | reject alignments not of these types (primary, supplementary; default primary) |

### Counting

| Option | Description |
| --- | --- |
| `--band N` | reference positions the pair HMM may look either side of the CIGAR (default 2) |
| `--min-phred Q` | assign the maximum sequencing error to bases below this quality (0 to 255; default 0) |
| `--ultima-phred Q` | the substitution Phred score of every base of an Ultima read, whose quality string describes run lengths only (0 to 255; default 40) |
| `--ultima-tag-order ORDER` | whether the tp and t0 values of a reverse-strand read follow its stored sequence, as Ultima's aligner and cmuts align write them, or the order its bases were sequenced, as minimap2 -y leaves them (stored, sequenced; default stored) |
| `--pairwise STATS` | write these statistics of how often two positions are modified together (correlation, conditional, none; default none) |
| `--min-depth D` | write NaN rates and errors for positions below this depth (default 1) |
| `--nan-5p N` | write NaN rates and errors for this many bases at the 5' end (default 0) |
| `--nan-3p N` | write NaN rates and errors for this many bases at the 3' end (default 0) |
| `--params FILE` | read the pair HMM's rates from this file (default: built in) |

### Performance

| Option | Description |
| --- | --- |
| `-j, --workers N` | threads running the processing step (default 1) |
| `--decode-threads N` | htslib threads for BGZF decompression (default 0) |

### Information

| Option | Description |
| --- | --- |
| `-h, --help` | show this help and exit |
| `-V, --version` | show the version and exit |
| `--dump-params` | write the rates in the form --params reads and exit |

### Advanced

Accepted, and left out of `--help`.

| Option | Description |
| --- | --- |
| `--queue-capacity N` | reads in transit at once (default 4096) |
| `--batch N` | reads transferred per queue operation (default 64) |
| `--live-refs N` | references in flight (default 64) |
| `--dump-options` | describe every argument as JSON and exit |
| `--dump-layout` | describe the output format as JSON and exit |
<!-- END GENERATED cmuts-hmm OPTIONS -->
