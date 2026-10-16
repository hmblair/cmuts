# Run length errors of Ultima reads in cmuts hmm

Notes from the sessions of 15 October 2026 on the branch `ultima-error-handling`.

## The problem

`cmuts hmm` reads each base's quality as the chance that the base was misread, and folds that chance into the match emission of the pair HMM. Quality is not used for gaps.

An Ultima sequencer measures the length of each run of one letter. Its errors are runs called one base too long or too short, and its quality string describes those errors, not misread bases. The `tp` tag gives, per base, the alternative run length that the base's quality is the chance of. The chance of an alternative length is the sum over the bases of the run that refer to it. The `t0` tag gives the chance that a base of another letter was missed between two runs. Sources: Almogy et al., bioRxiv 2022.05.29.493900, Methods, "Base calling"; GATK `FlowBasedRead.java`; and the data, where `tp` is palindromic in over 99 percent of runs.

So `cmuts hmm` counts the sequencer's length errors inside homopolymers as reverse transcription deletions. The benchmark team measured the result, in `max-analysis/`, as a loss of about 0.005 AUROC against the old tool on the Ultima data, with the deletion channel inside homopolymers near chance. The old tool avoided it by dropping bases below Q10, whose quality window reached across short runs and vetoed the deletion.

## The benchmark

Everything below is scored with `tools/mapseq/count.py` and `tools/mapseq/compare.py` from the paper repository, against the single-chain PDB structures, paired over PDB entries, at full depth. 2A3 counts mismatches and deletions; DMS counts mismatches only. Intervals are 95 percent. The data is under `/Users/hmblair/academic/papers/cmuts/data/raw-data/ribonanza-pdb/`, with the Ultima alignments in `alignments/ultima/` and the Illumina alignments in `alignments/illumina/`. The reference the programs take is the library with U replaced by T. Timings are one worker, user seconds, on the Illumina 2A3 BAM and on the first 60 references of the Ultima 2A3 BAM.

Main scores 0.6441 on Ultima 2A3, 0.7205 on Ultima DMS, 0.6302 on Illumina 2A3, and 0.6799 on Illumina DMS. Main takes 3.3 s on the Illumina BAM and 1.0 s on the Ultima slice.

## Three designs

### The five-state model

The HMM gains two states for a tagged read: a read base the sequencer added, and a reference base it dropped. Each is one base long, entered from a match and returning to a match. Their chances come from the tags, divided over the placements each error can take. The counts remain the posterior of the reverse transcription states. Every base of a tagged read takes a constant substitution quality, set by `--flow-phred` with a default of 40.

Untagged reads run main's loop unchanged, with output bitwise identical to main. On Ultima it gains 0.011 on 2A3, with interval 0.0004 to 0.0217, and 0.007 on DMS, with interval -0.0003 to 0.0144. The tagged path costs 1.7 times main per read, against a budget of 1.3. Not storing the sequencer states, passing the arms by pointer, precomputing per cell, and flushing denormals changed nothing. The notes of the first session say that most rows of a read carry no sequencer error chance, so a step that falls back to main's kernel on those rows is the lever not yet tried.

This implementation is in the stash at the top of the list, with the tag reading and decoding code reworded in plain terms. No tests, generator option, docs, or changelog entry exist for it.

### Reweighting the counts after the alignment

The alignment stays as main computes it. When the counts are added up, each deletion or insertion of exactly one base is divided between the reverse transcriptase and the sequencer in the ratio of the reverse transcription rate to the tag chance, and only the first part is counted. The one-base part of a deletion is read off the forward pass with one extra load, so the passes do not change.

On Ultima it gains 0.0074 on 2A3, with interval -0.0013 to 0.0164, and loses 0.0010 on DMS, with interval -0.0073 to 0.0051. It keeps 80 percent of the deletion count and 47 percent of the insertion count that main reports on the slice. Illumina output is identical to main. Untuned, it cost 1.6 times main on the slice, from reading the tags over the whole read, calling a lookup in another file four times per row, and reweighting every cell rather than only those with a nonzero chance.

The DMS result shows what this design cannot do. The DMS score counts mismatches only, and a mismatch at the end of a run can also be one run too short next to one run too long. The five-state model can move posterior mass from the mismatch to that pair of length errors, because the sequencer errors take part in the alignment. This design cannot, because the alignment is fixed before the tags are read. The change was discarded.

### The run model

The read is scored one run at a time. A cell pairs one read run with a prefix of the reference and holds three states: the run is complete after a pairing, a deletion is open after it, or an insertion is open after it. The step from one run to the next walks the bases of the run against the reference in a small matrix of its own, with the three states and the weights of the base model, so the per-base rates and the quality of each base enter as they do in main. A run gets two extra rows for the template bases the sequencer can fail to count, and the step is left at each possible template length with the chance from the tags. A base of another letter the sequencer missed between runs is a move over a one-base reference run with the `t0` chance. A deletion that opens after the last pairing of a run leaves the step in the open deletion state and extends only in the next step, so no path is counted twice.

The model is the generative model of the paper with the observation unit changed from a base to a run. On a per-base sequencer it reduces to main exactly: with main's band, the output on the Illumina 2A3 BAM is identical to main's in every channel to nine decimal places, and all 5051 tests pass.

With a band of two runs on either side of the CIGAR's run, it gains 0.0081 on Ultima 2A3, with interval -0.0018 to 0.0191, and 0.0088 on Ultima DMS, with interval 0.0004 to 0.0174. On Illumina it loses 0.0034 on 2A3 and 0.0017 on DMS. A band in runs is wide in bases wherever runs are long, and it admits alternative paths that main's band of two bases excludes, such as a deletion of two bases and an insertion of three in place of two mismatches and one insertion. Those paths blur the events. The benchmark team saw the same loss when they widened main's band.

With main's band of two bases around the CIGAR path, it gains 0.0106 on Ultima 2A3, with interval 0.0014 to 0.0204, and 0.0097 on Ultima DMS, with interval 0.0018 to 0.0181. Both intervals exclude zero. This is the best accuracy of the three designs.

Its cost is 3.2 times main on the Illumina BAM, 10.4 s against 3.3 s, and 4.9 times main on the Ultima slice, 5.1 s against 1.0 s. On Illumina the cell count equals main's, so the whole ratio is overhead per run: seeding the step from the row above, leaving it, clearing its backward matrix, and the loop over missed bases, each over the five cells of the band. A run is one or two bases on average, so that overhead is about the size of the work itself. On Ultima the step is left once per template length, five times rather than once, and the two extra rows per run are more than the run itself for most runs. Rewriting the step loops with precomputed per-base weights and plain pointers changed nothing, because the time was never in the address arithmetic. With the per-run overhead cut down, the floor is about twice main on Illumina and three times on Ultima. The extra rows and the five exits are the design.

On a Sherlock `biochem` node with 24 cores, with both trees built with gcc 12.4.0, the ratios are the same. Processor seconds over the whole run, the mean of three repetitions at 24 workers: main 12.7 against 40.0 on the Illumina BAM, and 4.3 against 21.3 on the Ultima slice, which is 3.1 and 5.0 times. Wall seconds at 24 workers: 0.8 against 2.0 on Illumina and 0.3 against 1.0 on the slice, where both inputs are too small to occupy 24 workers for long. One worker: 12.4 against 39.9 on Illumina and 4.1 against 21.6 on the slice. The job and its inputs are under `/scratch/users/hmblair/cmuts-timing/`.

## What was learned

- The sequencer errors must take part in the alignment. Reweighting the counts afterwards recovers the 2A3 gain only.
- The band must stay in bases. A band in runs costs accuracy on both platforms.
- The run model is the exact model, and it reduces to main where the sequencer is per base. Its cost is structural: a run is too short a unit on this data for per-run work to be cheap.
- The five-state model has the same gains within its intervals at a cost of 1.7 times on Ultima reads only. Its remaining lever is to run main's kernel on rows with no sequencer error chance.
- Design time must include a cell count against main before any code is written. Both prototypes of this session could have been costed from their design.

## Where things are

- The branch `ultima-error-handling` is at main's last commit. The working tree holds the run model, uncommitted: `include/runs.h` and `src/runs.c` split a sequence into runs; `include/flow.h` and `src/flow.c` decode the tags per run; `include/runhmm.h` and `src/runhmm.c` are the model; `src/bam.c` reads the tags and the read groups of the header; `src/tally.c`, `src/pipeline.c`, and `include/refctx.h` call the new model in place of `phmm_run`; `apps/hmm/options.c` adds `--flow-phred`. `src/phmm.c` is unchanged and unused by the tally.
- The stash at the top of the list holds the five-state model with the reworded tag code.
- The benchmark outputs of every design, the slice, the DNA reference, and a worktree of main with its binary are in the session's temporary directory, which does not persist. The worktree is registered in the repository's metadata until it is removed with `git worktree remove`.
- The data and tool paths are in the benchmark section above.
