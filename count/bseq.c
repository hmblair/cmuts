#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdbool.h>

#define A 0x0
#define C 0x1
#define G 0x2
#define T 0x3
#define INVALID 0x4

#define EOH (uint32_t)-1
#define HEADER_UNIT_LENGTH 2
#define HEADER_END_LENGTH 1

#define MAX_LINE_LENGTH 1024
#define BYTES(x) ((x+3)/4)
#define MIN(x,y) ((x < y) ? x : y)
#define MAX(x,y) ((x > y) ? x : y)

char encodeNucleotide(char nucleotide) {
    switch (nucleotide) {
        case 'A': return A;
        case 'C': return C;
        case 'G': return G;
        case 'T': return T;
        default: return INVALID;
    }
}


char* encodeDNASequence(
    const char* sequence,
    size_t* outputSize
) {

    size_t len = strlen(sequence);
    size_t numBytes = BYTES(len);
    char* binarySequence = malloc(numBytes);
    if (binarySequence == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return NULL;
    }

    memset(binarySequence, 0, numBytes);

    for (size_t i = 0; i < len; i++) {

        char encodedNuc = encodeNucleotide(sequence[i]);
        if (encodedNuc == INVALID) {
            free(binarySequence); // Free memory on invalid nucleotide
            fprintf(stderr, "Invalid nucleotide: %c\n", sequence[i]);
            return NULL;
        }
        binarySequence[i / 4] |= (encodedNuc << ((3 - (i % 4)) * 2));

    }

    *outputSize = numBytes;
    return binarySequence;

}


char* decodeDNASequence(
    const char* binarySequence,
    size_t sequenceLength
) {

    size_t binaryLength = BYTES(sequenceLength);
    size_t binaryTerminator = sequenceLength % 4;
    if (binaryTerminator == 0) {
        binaryTerminator = 4;
    }

    char *dnaSequence = malloc(sequenceLength+ 1);
    if (dnaSequence == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return NULL;
    }

    size_t index = 0;
    for (size_t i = 0; i < binaryLength - 1; i++) {

        for (int j = 0; j < 4; j++) {
            char nucleotide = (binarySequence[i] >> ((3 - j) * 2)) & 0x3;
            switch (nucleotide) {
                case A: dnaSequence[index++] = 'A'; break;
                case C: dnaSequence[index++] = 'C'; break;
                case G: dnaSequence[index++] = 'G'; break;
                case T: dnaSequence[index++] = 'T'; break;
            }
        }

    }

    // Deal with the last byte separately, as we may not use
    // the entire byte
    for (int j = 0; j < binaryTerminator; j++) {
        char nucleotide = (binarySequence[binaryLength-1] >> ((3 - j) * 2)) & 0x3;
        switch (nucleotide) {
            case A: dnaSequence[index++] = 'A'; break;
            case C: dnaSequence[index++] = 'C'; break;
            case G: dnaSequence[index++] = 'G'; break;
            case T: dnaSequence[index++] = 'T'; break;
        }
    }

    dnaSequence[index] = '\0';
    return dnaSequence;

}


typedef struct {
    char *header;
    ssize_t headerSize;
    char *sequence;
    ssize_t sequenceSize;
    uint32_t referenceNumber;
    FILE *file;
} Record;


void initRecord(
    Record *record,
    const char *filename
) {

    record->file = fopen(filename, "r");
    if (!record->file) {
        perror("Error opening FASTA file: ");
        return;
    }

    record->header = NULL;
    record->sequence = NULL; 

    record->referenceNumber = -1;

    return;

}


void freeRecord(
    Record record
) {

    free(record.header);
    free(record.sequence);
    fclose(record.file);

}


bool nextRecord(
    Record *record
) {

    int atEOF = fseek(record->file, 1, SEEK_CUR);
    if (atEOF < 0) {
        return false;
    }

    size_t bufferSize;
    record->headerSize = getline(
        &record->header,
        &bufferSize,
        record->file
    ) - 1;
    if (record->headerSize < 0) {
        return false;
    }
    if (record->headerSize > 0 && record->header[record->headerSize] == '\n') {
        record->header[record->headerSize] = '\0';
    }

    record->sequenceSize = getline(
        &record->sequence,
        &bufferSize,
        record->file
    ) - 1;
    if (record->sequenceSize < 0) {
        return false;
    }
    if (record->sequenceSize > 0 && record->sequence[record->sequenceSize] == '\n') {
        record->sequence[record->sequenceSize] = '\0';
    }

    (record->referenceNumber)++;

    return true;

}


static inline bool readHeaderUnit(
    FILE *binaryFile,
    uint32_t *headerUnit
) {

    size_t itemsRead = fread(
        headerUnit,
        sizeof(uint32_t),
        HEADER_UNIT_LENGTH,
        binaryFile
    );
    if (itemsRead < HEADER_UNIT_LENGTH) {
        return false;
    }

    if (headerUnit[0] == EOH) {
        return false;
    }
    return true;

}

static inline void writeHeaderUnit(
    FILE *binaryFile,
    uint32_t *headerUnit
) {

    if (headerUnit[1] > 0) {
        fwrite(
            headerUnit,
            sizeof(uint32_t),
            HEADER_UNIT_LENGTH,
            binaryFile 
        );

        headerUnit[0] = 0;
        headerUnit[1] = 0;
    }

}

static inline void endHeader(
    FILE *binaryFile
) {

    uint32_t eoh = EOH;
    fwrite(
        &eoh,
        sizeof(uint32_t),
        HEADER_END_LENGTH,
        binaryFile 
    );

}

void writeHeader(
    const char *infilename,
    FILE *binaryFile
) {

    Record record;
    initRecord(
        &record,
        infilename
    );

    uint32_t headerUnit[HEADER_UNIT_LENGTH] = {0, 0};

    while (nextRecord(&record)) {

        // If we have hit a sequence of a new length, then
        // write the current header unit to file
        if (record.sequenceSize != headerUnit[0]) {
            writeHeaderUnit(binaryFile, headerUnit);
        }

        headerUnit[0] = record.sequenceSize;
        headerUnit[1]++;

    }

    // Write any final header units to file
    writeHeaderUnit(binaryFile, headerUnit);
    // End the header with the EOH value
    endHeader(binaryFile);
    // Free the record we created
    freeRecord(record);

}


void writeBinarySequences(
    const char *infilename,
    FILE *outputFile
) {

    Record record;
    initRecord(
        &record,
        infilename
    );


    size_t outputSize;
    while (nextRecord(&record)) {

        char* binarySequence = encodeDNASequence(
            record.sequence,
            &outputSize
        );

        if (!binarySequence) {
            printf("The binary sequence failed to be encoded.\n");
            exit(1);
        }

        size_t nwritten = fwrite(
            binarySequence,
            sizeof(char),
            outputSize,
            outputFile
        );

        if (nwritten != outputSize) {
            printf("The binary sequence failed to be written.\n");
            exit(1);
        }

        free(binarySequence);

    }

    freeRecord(record);

}


void getOffset(
    FILE *binaryFile,
    uint32_t sequenceNumber,
    uint32_t *sequenceLength,
    uint32_t *binaryOffset
) {

    uint32_t remainingOffset = sequenceNumber;
    uint32_t headerUnit[2] = {0,0};
    fseek(binaryFile, 0, SEEK_SET);

    do {

        bool status = readHeaderUnit(binaryFile, headerUnit);
        if (!status) {
            fprintf(stderr, "The sequence number %" PRIu32 " cannot be found in the file.\n", sequenceNumber);
            exit(EXIT_FAILURE);
        }

        // Get the number of bytes and the number of sequences
        // indicated by this header unit
        uint32_t nbytes = BYTES(headerUnit[0]);
        uint32_t nseqs = MIN(remainingOffset, headerUnit[1] - 1);
        remainingOffset -= nseqs;
        *binaryOffset += nseqs * nbytes;
        *sequenceLength = headerUnit[0];

        // Take care of this portion of the header
        *binaryOffset += 2 * sizeof(uint32_t);

    } while (remainingOffset > 0);

    // Keep going until we reach the end of the header
    while (readHeaderUnit(binaryFile, headerUnit)) {
        *binaryOffset += 2 * sizeof(uint32_t);
    }
    // Account for the EOH character
    *binaryOffset += sizeof(uint32_t);

}

const char *readBinarySequence(
    FILE *binaryFile,
    uint32_t sequenceNumber
) {

    uint32_t sequenceLength = 0;
    uint32_t binaryOffset = 0;
    getOffset(
        binaryFile,
        sequenceNumber,
        &sequenceLength,
        &binaryOffset
    );

    int numBytes = BYTES(sequenceLength);
    char *binaryBuffer = malloc(numBytes * sizeof(char));
    fseek(binaryFile, binaryOffset, SEEK_SET);

    size_t numBytesRead = fread(
        binaryBuffer,
        sizeof(char),
        numBytes,
        binaryFile 
    );
    if (numBytesRead != numBytes) {
        fprintf(stderr, "Error reading binary data\n");
        exit(1);
    }

    const char *sequence = decodeDNASequence(
        binaryBuffer,
        sequenceLength 
    );

    free(binaryBuffer);
    return sequence;

}


int toBinary(
    const char *infilename,
    const char *outfilename
) {

    // Remove the output file
    remove(outfilename);

    // Open the output file 
    FILE *binaryFile = fopen(outfilename, "rb");
    // Write the header to the file
    writeHeader(infilename, binaryFile);
    // Write the sequences to the file
    writeBinarySequences(infilename, binaryFile);
    
    return EXIT_SUCCESS;
}

int bseq_main(int argc, char **argv) {

    if (argc != 3) {
        printf("usage: tobinary FASTA BINARY\n");
        exit(EXIT_FAILURE);
    }

    const char *infilename = argv[1];
    const char *outfilename = argv[2];

    return toBinary(infilename, outfilename);

}
