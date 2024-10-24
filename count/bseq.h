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
);

void freeRecord(
    Record record
);

bool nextRecord(
    Record *record
);



char* encodeDNASequence(
    const char* sequence,
    size_t* outputSize
);

char* decodeDNASequence(
    const char* binarySequence,
    size_t sequenceLength
);



void writeHeader(
    const char *infilename,
    FILE *binaryFile
);

void writeBinarySequences(
    const char *infilename,
    FILE *outputFile
);

void getOffset(
    FILE *binaryFile,
    uint32_t sequenceNumber,
    uint32_t *sequenceLength,
    uint32_t *binaryOffset
);

const char *readBinarySequence(
    FILE *binaryFile,
    uint32_t sequenceNumber
);
