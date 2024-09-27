#include "htsutils.h"


//
// Helper functions
//
//

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

static inline char* parseMD(const bam1_t *aln) {
    uint8_t* md_ptr = bam_aux_get(aln, "MD");

    if (md_ptr == NULL) {
        printf("No MD tag was found.");
        exit(EXIT_FAILURE);
    }

    return bam_aux2Z(md_ptr);
}

static inline int8_t baseToInt(char base) {
    switch (base) {
        case 'A':
            return IX_A;
        case 'C':
            return IX_C;
        case 'G':
            return IX_G;
        case 'T':
            return IX_T;
        case 'U':
            return IX_U;
        default:
            fprintf(stderr, "Error in basetoInt; '%c' was passed.\n", base);
            exit(EXIT_FAILURE);
    } 
}

static inline int getDecimalDigit(const char* str, int i) {

    int digit = 0;
    while(isdigit(str[i])) {
        digit = digit * 10 + (str[i] - '0');
        i++;
    }

    return digit;
}

// get the number of consecutive digits in the string
// str starting from position i.
static inline int getNumDigits(const char* str, int i) {

    int n = 0;
    while(isdigit(str[i+n])) {
        n++;
    }

    return n;
}

// get the the base at position ix in the query sequence
static inline int8_t getQuerySeqBase(uint8_t *querySequence, int ix) {

    // the query sequence value, in BAM format
    uint8_t bamIX = bam_seqi(querySequence, ix);
    // convert from BAM format (1,2,4,8) to standart
    // format (0,1,2,3)
    switch (bamIX) {
        case BAM_A:
            return IX_A;
        case BAM_C:
            return IX_C;
        case BAM_G:
            return IX_G;
        case BAM_T:
            return IX_T;
        default:
            return -1;
    }

}

static inline bool isSubsequence(
    char* s1,
    char* s2,
    int len1,
    int len2
) {

    // Traverse s2 and match characters with s1
    int i = 0;
    for (int j = 0; j < len2; j++) {
        if (s1[i] == s2[j]) {
            i++;
        }
    }

    // If we've matched all characters in s1, it's a subsequence
    return (i >= len1);

}



//
// For working with indexed FASTA files
//
//

IndexedFASTA openIndexedFASTA(const char *filename) {

    IndexedFASTA ixFASTA;

    ixFASTA.filename = filename;

    if (access(filename, R_OK) != 0) {
        fprintf(stderr, "There was an issue opening the file %s.", filename);
    }

    ixFASTA.fai = fai_load(filename);
    if (!ixFASTA.fai) {
        fprintf(stderr, "Failed to load indexed FASTA file %s.\n", filename);
        exit(EXIT_FAILURE);
    }

    ixFASTA.numSequences = faidx_nseq(ixFASTA.fai);

    return ixFASTA;

}

char *getSequence(IndexedFASTA ixFASTA, int ix) {

    if (ix < 0 || ix >= ixFASTA.numSequences) {
        fprintf(stderr, "The index (%d) is larger than the number of sequences in the FASTA file %s (%d).\n", ix, ixFASTA.filename, ixFASTA.numSequences);
        return NULL;
    }

    // Get the sequence name for the given sequence number
    const char *name = faidx_iseq(ixFASTA.fai, ix);
    if (!name) {
        fprintf(stderr, "Failed to get the name of sequence %s.\n", name);
        return NULL;
    }

    int len = 0;
    char *sequence = fai_fetch(ixFASTA.fai, name, &len);
    if (!sequence) {
        fprintf(stderr, "Failed to fetch sequence %s.\n", name);
        return NULL;
    }

    return sequence;


}

void closeIndexedFASTA(IndexedFASTA ixFASTA) {
    fai_destroy(ixFASTA.fai);
}


//
// For working with indexed SAM/BAM/CRAM files
//


IndexedBAM openIndexedBAM(const char* filename) {

    IndexedBAM ixBAM;
    ixBAM.filename = filename;

    // Open the ixBAM file
    ixBAM.file = sam_open(filename, "r");
    if (ixBAM.file == NULL) {
        fprintf(stderr, "Error: Failed to open the ixBAM file \"%s\".\n", filename);
        exit(EXIT_FAILURE);
    }

    // Open the ixBAM index file
    ixBAM.index = sam_index_load(ixBAM.file, filename);
    if (ixBAM.index == NULL) {
        fprintf(stderr, "Error: Failed to open ixBAM index associated with \"%s\".\n", filename);
        sam_close(ixBAM.file);
        exit(EXIT_FAILURE);
    }

    // Read the header of the ixBAM file
    ixBAM.header = sam_hdr_read(ixBAM.file);
    if (ixBAM.header == NULL) {
        fprintf(stderr, "Error: Failed to read header for \"%s\".\n", filename);
        sam_close(ixBAM.file);
        exit(EXIT_FAILURE);
    }

    // Get the number of reference sequences
    ixBAM.numReferences = ixBAM.header->n_targets;
    // Get the number of alignments and unmapped reads;
    // also keep track of the longest sequence
    ixBAM.numAlignments = 0;
    ixBAM.numUnmappedReads = 0;
    ixBAM.maxReferenceLength = 0;
    for (int i = 0; i < ixBAM.numReferences; i++) {
        uint64_t mapped = 0, unmapped = 0;
        if (hts_idx_get_stat(ixBAM.index, i, &mapped, &unmapped) < 0) {
            continue;
        }
        ixBAM.numAlignments += mapped;
        ixBAM.numUnmappedReads += unmapped;
        // If the current sequence is longer than the previous max,
        // then update the max length
        if (ixBAM.header->target_len[i] > ixBAM.maxReferenceLength) {
            ixBAM.maxReferenceLength = ixBAM.header->target_len[i];
        }
    }

    // Allocate space for one alignment
    ixBAM.alignment = bam_init1();

    return ixBAM;

}

void closeIndexedBAM(IndexedBAM ixBAM) {

    sam_close(ixBAM.file);
    bam_hdr_destroy(ixBAM.header);
    hts_idx_destroy(ixBAM.index);
    bam_destroy1(ixBAM.alignment);

}

CountingFlags defaultCountingFlags() {

    CountingFlags cFlags;

    cFlags.minMapQ = DEFAULT_MIN_MAPPING_QUALITY;
    cFlags.minBaseQ = DEFAULT_MIN_NUCLEOTIDE_QUALITY;
    cFlags.maxDelLength = DEFAULT_MAX_DELETION_LENGTH;
    cFlags.minQueryLength = DEFAULT_MIN_QUERY_LENGTH;
    cFlags.collapseMutations = DEFAULT_COLLAPSE_MUTATIONS;
    cFlags.spreadDeletions = DEFAULT_SPREAD_DELETIONS;
    cFlags.numNeighboursToCheck = DEFAULT_NEIGHBOURS_TO_CHECK;

    return cFlags;

}

bool verifyCountingFlags(CountingFlags cFlags) {

    if (cFlags.minMapQ < 0 || cFlags.minMapQ > 100) {
        fprintf(stderr, "The minimum mapping quailty must be between 0 and 100.\n");
        return false;
    }
    if (cFlags.minBaseQ < 0 || cFlags.minBaseQ > 100) {
        fprintf(stderr, "The minimum base quailty must be between 0 and 100.\n");
        return false;
    }
    if (cFlags.maxDelLength < 0) {
        fprintf(stderr, "The maximum deletion length must be non-negative.\n");
        return false;
    }
    if (cFlags.minQueryLength < 0) {
        fprintf(stderr, "The minimum query length must be non-negative.\n");
        return false;
    }
    if (cFlags.collapseMutations < 0) {
        fprintf(stderr, "The distance to collapose mutations over must be non-negative.\n");
        return false;
    }
    if (cFlags.numNeighboursToCheck< 0) {
        fprintf(stderr, "The number of neighbours to quality-check must be non-negative.\n");
        return false;
    }

    return true;

}



//
// Functions for determining the position of an ambiguous insertion
// or deletion. The "classic" versions are those which are used by
// rf-count, and are included for compatability reasons.
//



static int getAmbiguousStartClassic(
    int deletionStart,
    int deletionEnd,
    int deletionLength,
    const char *referenceSequence,
    int referenceLength
) {

    // If we are already at the sart of the sequence,
    // do nothing
    if (deletionStart == 0) {
        return deletionStart;
    }
    // Allocate enough memory to store the expanding deletion
    char *deletion = calloc(deletionEnd, sizeof(char));
    if (!deletion) {
        fprintf(stderr, "There was an error allocating memory for the deletion buffer.\n");
        exit(EXIT_FAILURE);
    }

    // Since we are trying to find the start, we are running 
    // backwards along the deletion.
    for (int i = 0; i < deletionLength; i++) {
        deletion[i] = referenceSequence[deletionEnd - i - 1];
    }

    deletion[deletionLength] = referenceSequence[deletionStart-1];

    int curr_ix = 0;
    while (deletion[deletionLength] == deletion[curr_ix]) {

        // If we have reached the start of the reference sequence,
        // then we should go no more
        if (deletionStart == 0) {
            free(deletion);
            return deletionStart;
        }

        deletionStart--;
        deletionLength++;
        deletion[deletionLength] = referenceSequence[deletionStart-1];
        curr_ix++;

        // Temporary to comply with rf-count
        if (curr_ix == 10) {
            free(deletion);
            return deletionStart;
        }

    }

    free(deletion);
    return deletionStart;

}

static int getAmbiguousEndClassic(
    int deletionStart,
    int deletionEnd,
    int deletionLength,
    const char *referenceSequence,
    int referenceLength
) {

    // If we are already at the end of the sequence,
    // do nothing
    if (deletionEnd == referenceLength) {
        return deletionEnd;
    }
    // Allocate enough memory to store the expanding deletion
    char *deletion = calloc(referenceLength - deletionStart, sizeof(char));
    if (!deletion) {
        fprintf(stderr, "There was an error allocating memory for the deletion buffer.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < deletionLength; i++) {
        deletion[i] = referenceSequence[i + deletionStart];
    }
    
    deletion[deletionLength] = referenceSequence[deletionEnd];

    int curr_ix = 0;
    while (deletion[deletionLength] == deletion[curr_ix]) {

        // If we have reached the end of the reference sequence,
        // then we should go no more
        if (deletionEnd == referenceLength) {
            free(deletion);
            return deletionEnd;
        }

        deletionLength++;
        deletionEnd++;
        deletion[deletionLength] = referenceSequence[deletionEnd];
        curr_ix++;

        // Temporary to comply with rf-count
        if (curr_ix == 10) {
            free(deletion);
            return deletionEnd;
        }

    }

    free(deletion);
    return deletionEnd ;

}

static inline int firstMatchInString(char c, const char* str, int M, int N) {
    // Loop and return the position after the first
    // match in the string if it is found
    for (int i = M; i < N; i++) {
        if (str[i] == c) {
            return i+1;
        }
    }
    // If there is no match, return -1.
    return -1;
}

static inline int lastMatchInString(char c, const char* str, int M, int N) {
    // Loop and return the position before the last
    // match in the string if it is found
    for (int i = N; i < M; i--) {
        if (str[i] == c) {
            return i-1;
        }
    }
    // If there is no match, return -1.
    return -1;
}

static int getAmbiguousEnd(
    int deletionStart,
    int deletionEnd,
    const char *referenceSequence,
    int referenceLength
) {

    // If we are already at the end of the sequence,
    // do nothing
    if (deletionEnd == referenceLength) {
        return deletionEnd;
    }

    // We wish to find the final base for which the deletion could
    // have occured at. Given a string S, with deletion starting at
    // position M and ending at position N, we wish to find the largest
    // K such that S[K] could be deleted to yield the same deletion.
    // If D is the string with the deletion, then we make use of
    // the fact that S[K] is a valid deletion if and only if 
    // D[M:K] is a subsequence of S[M:K-1].

    // We can find the largest such K in an efficient manner.
    // Check if D[N] is in S[M:N]. While this is so,
    //  - Advance M to the base after the first position in 
    //    S where D[N] occurs, and
    //  - Update N -> N + 1.
    // The final N is the furthest base in the 3' direction
    // that can be deleted to yield the same deletion. If N
    // ever reaches the end of the reference sequence, we 
    // must terminate.

    int M = deletionStart;
    int N = deletionEnd;

    char c = referenceSequence[N];
    M = firstMatchInString(c, referenceSequence, M, N);
    while (M != -1 && N != referenceLength) {
        N++;
        c = referenceSequence[N];
        M = firstMatchInString(c, referenceSequence, M, N);
    }

    return N;

}

static int getAmbiguousStart(
    int deletionStart,
    int deletionEnd,
    const char *referenceSequence,
    int referenceLength
) {

    // If we are already at the sart of the sequence,
    // do nothing
    if (deletionStart == 0) {
        return deletionStart;
    }

    // The algorithm for finding the first base that could
    // be feasilby deleted is dual to that for finding the
    // last such base.
    int M = deletionStart;
    int N = deletionEnd;

    char c = referenceSequence[M-1];
    N = lastMatchInString(c, referenceSequence, M, N);
    while (N != -1 && M != 0) {
        M--;
        c = referenceSequence[M-1];
        N = lastMatchInString(c, referenceSequence, M, N);
    }

    return M;

}


//
// The main mutation and insertion counting function
//


static int accumulateMutations(
    float *mutations,
    float *insertions,
    char *referenceSequence,
    IndexedBAM ixBAM,
    CountingFlags cFlags
) {

    // Guard based on the mapping quality
    uint8_t mappingQ = ixBAM.alignment->core.qual;
    if (mappingQ < cFlags.minMapQ) {
        return QUERY_SKIPPED;
    }

    // Guard based on the length of the aligned read
    int queryLength = ixBAM.alignment->core.l_qseq;
    if (queryLength < cFlags.minQueryLength) {
        return QUERY_SKIPPED;
    }

    // Get the current reference number
    int32_t refIX = ixBAM.alignment->core.tid;
    // Get the length of the reference sequence
    hts_pos_t referenceLength = ixBAM.header->target_len[refIX];
    // Initialise the position in the reference sequence based
    // on the offset of the alignment 
    hts_pos_t referencePosition = ixBAM.alignment->core.pos;
    // Initialise the position in the query sequence to zero
    hts_pos_t queryPosition = 0;

    // Get the CIGAR string, in the htslib format
    uint32_t *cigar = bam_get_cigar(ixBAM.alignment);
    // Get the query sequence, in the htslib format
    uint8_t *querySequence = bam_get_seq(ixBAM.alignment);
    // Get the MD tag, in human-readable format
    const char* mdTag = parseMD(ixBAM.alignment);
    // Get the per-nucleotide quality scores
    uint8_t *quality = bam_get_qual(ixBAM.alignment);
    // Get the PHRED scores of each nucleotide
    uint8_t *perNucMappingQ = bam_get_qual(ixBAM.alignment);

    // An index for tracking the current position in the MD tag
    int mdIX = 0;
    // Since the MD tag is in the inner loop, if the current
    // op of the MD tag is longer than the current op of the
    // CIGAR string, we need to keep track of the difference
    int mdOpPosition = 0;

    // For the purpose of collapsing nearby mutations,
    // we need to keep track of when the last mutation
    // that occured was. For now, we initialise it large
    // enough so that if the first mutation is on the
    // first base, there is not an error.
    int8_t lastMut = -cFlags.collapseMutations;
    int8_t obase_l = -1;
    int8_t mbase_l = -1;

    for (int cigIX = 0; cigIX < ixBAM.alignment->core.n_cigar; cigIX++) {

        int cigarOp = bam_cigar_op(cigar[cigIX]);
        hts_pos_t cigarOpLength = bam_cigar_oplen(cigar[cigIX]);

        // If the CIGAR string indicates that there is a mis/match,
        // we loop over the current portion of the MD tag to see
        // if any of the mis/matches are mutations.
        if (cigarOp == BAM_CMATCH || cigarOp == BAM_CDIFF) {

            while (cigarOpLength > 0) {

                // Get the current mdTag tag operation as an
                // alphanumeric value
                char mdChar = mdTag[mdIX];

                // If the operation is an integer n, then there are no
                // mutations for the next n bases, so we may skip forward
                // in the CIGAR string by this amount.
                if (isdigit(mdChar)) {

                    // Get the remaining length of the current MD tag
                    // operation
                    int mdOpLength = getDecimalDigit(mdTag, mdIX) - mdOpPosition;
                    // If the MD oplen is smaller than the CIGAR oplen,
                    // we can process this entire section of the
                    // mdTag tag
                    if (mdOpLength <= cigarOpLength) {
                        
                        // Advance the appropriate amount into the
                        // CIGAR string
                        cigarOpLength -= mdOpLength;

                        while (mdOpLength > 0) {

                            int8_t base = getQuerySeqBase(querySequence,queryPosition);
                            // If the base is ACGU, increment the
                            // corresponding coverage count
                            if (base != -1) {
                                (mutations[referencePosition * N_BASES * N_DELBASES + base * N_DELBASES + base])++;
                            }
                            // Update both the reference and query positions
                            referencePosition++;
                            queryPosition++;
                            mdOpLength--;

                        }

                        // Move along the MD tag until we reach the next
                        // operation
                        mdIX += getNumDigits(mdTag, mdIX);
                        // Zero out how far into the next MD operation
                        // we are
                        mdOpPosition = 0;

                    }

                    // Else, take note of how much of the MD tag
                    // was seen in this round, but do not move the
                    // MD tag index.
                    else {

                        mdOpPosition += cigarOpLength;

                        while (cigarOpLength > 0) {

                            // If the base is ACGU, increment the
                            // corresponding coverage count
                            int8_t base = getQuerySeqBase(querySequence, queryPosition);
                            if (base != -1) {
                                (mutations[referencePosition * N_BASES * N_DELBASES + base * N_DELBASES + base])++;
                            }
                            // Update both the reference and query positions
                            referencePosition++;
                            queryPosition++;
                            cigarOpLength--;

                        }


                    }

                // If the operation is a letter, then there is a mutation
                // from the base indicated by the letter to the base at
                // the current position in the query sequence.
                } else if (isalpha(mdChar)) {

                    // Get the original base
                    int8_t obase = baseToInt(mdChar);
                    // Get the mutated base
                    int8_t mbase = getQuerySeqBase(querySequence, queryPosition);

                    // Check the mapping quality of the base and
                    // of any neighbouring bases, if specified.
                    bool ofQuality = true;
                    int minNeighhourToCheck = MAX(queryPosition - cFlags.numNeighboursToCheck, 0);
                    int maxNeighhourToCheck = MIN(queryPosition + cFlags.numNeighboursToCheck + 1, queryLength);

                    for (int pos = minNeighhourToCheck; pos < maxNeighhourToCheck; pos++) {

                        if (perNucMappingQ[pos] < cFlags.minBaseQ) {
                            ofQuality = false;
                            break;
                        }

                    }

                    // If mbase is -1, then the mutated base is not one of
                    // ACGU. We do what RNAFramework does here, which is
                    // to ignore the mutation. Else, we increment the
                    // associated section of the mutation count.
                    if (mbase != -1 && ofQuality) { 

                        // If the last mutation is within cFlags.collapseMutations
                        // of the current mutation, we remove the previous
                        // mutation and only keep the current one.
                        if (referencePosition - lastMut < cFlags.collapseMutations) {
                            (mutations[lastMut * N_BASES * N_DELBASES + obase_l * N_DELBASES + mbase_l])--;
                        }
                        (mutations[referencePosition * N_BASES * N_DELBASES + obase * N_DELBASES + mbase])++;
                        lastMut = referencePosition;
                        obase_l = obase;
                        mbase_l = mbase;
                    } 
                    // In order for the position to count towards the 
                    // coverage, we act as if there is no mutation.
                    else {
                        (mutations[referencePosition * N_BASES * N_DELBASES + obase * N_DELBASES + obase])++;
                    }
                    // Move forward by one position in the reference and
                    // query sequences, as well as in the MD tag and the
                    // current CIGAR operation.
                    referencePosition++;
                    queryPosition++;
                    mdIX++;
                    cigarOpLength--;

                    // If the next operation in the MD tag is a 0, then
                    // the CIGAR string and MD tag will become out of
                    // sync, so we must check for this, and advance if
                    // necessary
                    if (mdTag[mdIX] == '0') {
                        mdIX++;
                    }

                // If the operation is a caret, then there is a deletion
                // in the MD tag which is not present in the CIGAR
                // string, so we raise an error.
                } else if (mdChar == '^') {

                    fprintf(stderr, "There is a deletion in the MD tag which is unaccounted for in the CIGAR string.");
                    exit(EXIT_FAILURE);

                }

            }

        // If the CIGAR string indicates that there is an insertion,
        // then we can handle them without reference to the MD tag
        // since they do not show up there.
        } else if (cigarOp == BAM_CINS) {

            // Get the inserted base
            int8_t base = getQuerySeqBase(querySequence, queryPosition);
            // If the inserted base is not ACGU, then we do
            // not perform any operations other than incrementing
            // the query position. We follow rf-count and place the
            // insertion at the reference base left of the insertion.
            if (base != -1) {
                insertions[(referencePosition - 1) * N_BASES + base]++;
            }
            // Incrememnt the query position by one for each
            // insertion. Do not incremement the reference sequence
            // as insertions do not affect it.
            queryPosition += cigarOpLength;
            cigarOpLength = 0;

        // If the CIGAR string indicates that there is a deletion,
        // then there should be a corresponding set of deletions 
        // in the MD tag and so we can loop over both strings 
        // simultaneously
        } else if (cigarOp == BAM_CDEL) {

            // Move forward by one in the MD tag to get 
            // past the caret
            if (mdTag[mdIX] != '^') {
                fprintf(stderr, "There is a deletion in the CIGAR string which is unaccounted for in the MD tag\n.");
                exit(EXIT_FAILURE);
            }
            mdIX++;

            if (cigarOpLength <= cFlags.maxDelLength) {

                // Get the start, end, and length of the deletion
                // as detected by the alignment
                int deletionLength = cigarOpLength;
                int deletionStart = referencePosition;
                int deletionEnd = referencePosition + deletionLength;
                // Traverse the sequence backwards to find the
                // furthest position where the deletion could be
                // in the 5' direction
                int trueDeletionStart = getAmbiguousStart(
                    deletionStart,
                    deletionEnd,
                    referenceSequence,
                    referenceLength
                );
                // Traverse the sequence forwards to find the
                // furthest position where the deletion could be
                // in the 3' direction
                int trueDeletionEnd = getAmbiguousEnd(
                    deletionStart,
                    deletionEnd,
                    referenceSequence,
                    referenceLength
                );
                // Get the length of the ambiguous region
                int trueDeletionLength = trueDeletionEnd - trueDeletionStart;

                if (cFlags.spreadDeletions) {

                    // Spread the deletion evenly among the bases determined
                    // to be part of the deletion
                    float delValue = (float) 1 / trueDeletionLength;
                    for (int pos = trueDeletionStart; pos < trueDeletionEnd; pos++) {
                        int base = baseToInt(referenceSequence[pos]);
                        mutations[pos * N_BASES * N_DELBASES + base * N_DELBASES + IX_DEL] += delValue;
                    }

                } else {

                    // Place the deletion at the 3'-most base
                    int base = baseToInt(referenceSequence[trueDeletionEnd - 1]);
                    mutations[(trueDeletionEnd - 1) * N_BASES * N_DELBASES + base * N_DELBASES + IX_DEL]++;

                }

            }

            // Move forward in the reference position and the
            // MD tag
            referencePosition += cigarOpLength;
            mdIX += cigarOpLength;
            cigarOpLength = 0;

            // If the next operation in the MD tag is a 0, then
            // the CIGAR string and MD tag will become out of
            // sync, so we must check for this, and advance if
            // necessary.
            if (mdTag[mdIX] == '0') {
                mdIX++;
            }

        } else if (cigarOp == BAM_CSOFT_CLIP) {

            queryPosition += cigarOpLength;
            cigarOpLength = 0;

        } else if (cigarOp == BAM_CHARD_CLIP) {

            // Do nothing.
            cigarOpLength = 0;

        } else if (cigarOp == BAM_CREF_SKIP) {

            referencePosition += cigarOpLength;
            cigarOpLength = 0;

        } else {

            fprintf(stderr, "An unknown CIGAR operation (%d) was detected.\n", cigarOp);
            exit(EXIT_FAILURE);

        }

    }

    return QUERY_PROCESSED;

}


AlignmentCount countMutations(
    float *mutations,
    float *insertions,
    char *referenceSequence,
    IndexedBAM ixBAM,
    int ix, 
    CountingFlags cFlags
) {

    // Create an object to store the number of alignments
    // processed and the number skipped
    AlignmentCount count;
    count.numAlignments = 0;
    count.numSkipped = 0;

    // Create an iterator to iterate all alignments associated with
    // the current reference sequence
    hts_itr_t *iter = sam_itr_queryi(ixBAM.index, ix, 0, HTS_POS_MAX);
    if (iter == NULL) {
        fprintf(stderr, "There was an error creating the iterator for reference sequence %d.\n", ix);
        exit(EXIT_FAILURE);
    }

    // Loop over the alignments and count how many were processed
    while (sam_itr_next(ixBAM.file, iter, ixBAM.alignment) >= 0) {
        count.numSkipped += accumulateMutations(
            mutations,
            insertions,
            referenceSequence,
            ixBAM,
            cFlags
        );
        count.numAlignments++;
    }

    return count;

}
