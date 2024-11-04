#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

char *randomSequence(
    int length
) {
    
    char *sequence = malloc(sizeof(char) * (length + 1));
    for (int i = 0; i < length; i++) {
        sequence[i] = "ACGU"[rand() % 4];
    }
    sequence[length] = '\0';

    return sequence;

}

char *randomMutation(
    char *sequence,
    int numMutations
) {

    int seqlen = strlen(sequence);
    char *newSequence = malloc(sizeof(char) * (seqlen - numMutations + 1));

    for (int i = numMutations; i < seqlen; i++) {
        newSequence[i - numMutations] = sequence[i];
    }
    newSequence[seqlen - numMutations] = '\0';

    return newSequence;

}

void printMutationsDeletions(
    int *mutations,
    int *deletions,
    size_t length,
    int numMutations
) {

    printf("Mutations (%d) : ", numMutations);
    for (int ix = 0; ix < length; ix++) {
        printf("%d ", mutations[ix]);
    }
    printf("\n");

    printf("Deletions     : ");
    for (int ix = 0; ix < length; ix++) {
        printf("%d ", deletions[ix]);
    }
    printf("\n\n");

}

void getMututationsDeletions(
    const char *reference,
    const char *read
) {

    size_t reflen = strlen(reference);
    int *mutations = calloc(reflen, sizeof(int));
    int *deletions = calloc(reflen, sizeof(int));

    int numMutations = 0;

    // The first pass
    deletions[0]++;
    for (size_t ix = 1; ix < reflen; ix++) {
        if (reference[ix] != read[ix - 1]) {
            mutations[ix]++;
            numMutations++;
        }
    }
    printMutationsDeletions(
        mutations,
        deletions,
        reflen,
        numMutations
    );
    // The second pass
    for (size_t ix = 1; ix < reflen; ix++) {

        // Update the deletions
        deletions[ix-1]--;
        deletions[ix]++;

        // Check if the previously deleted position is now a mutation
        if (reference[ix-1] != read[ix-1]) {
            mutations[ix-1]++;
            numMutations++;
        }
        // Check if the current deleted position was previously a mutation
        if (reference[ix] != read[ix] && mutations[ix] > 0) {
            mutations[ix]--;
            numMutations--;
        }
        printMutationsDeletions(
            mutations,
            deletions,
            reflen,
            numMutations
        );
    }


}


int main(int argc, char **argv) {

    srand(time(NULL));
    char *sequence = randomSequence(177);
    char *mut_sequence = randomMutation(sequence, 1);

    getMututationsDeletions(
        sequence,
        mut_sequence
    );

}
