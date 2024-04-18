#ifndef CONSTANTS_H
#define CONSTANTS_H

#define K 64

#define MAX_LIST_SIZE 10000
#define TRIALS 5e4

// --------------------------- rate 1/2 into two rate 1/1
// ---------------------------

// 3345: 11011100101 x^10 + x^9 + x^7 + x^6 + x^5 + x^2 + 1
// One can divide it into: (x^3 + x^2 + 1)*(x^7 + x^3 + 1)
// x^7 + x^3 + 1 is (binary)10001001, (octal)211, (decimal)137
// x^3 + x^2 + 1 is (binary)1101, (octal)15, (decimal)13

// 3613: 11110001011 x^10 + x^9 + x^8 + x^7 + x^3 + x + 1
// One can divide it into: (x^6 + x^5 + x^2 + x + 1)*(x^4 + x^2 + 1)
// x^6 + x^5 + x^2 + x + 1 is (binary)1100111, (octal)147, (decimal)103
// x^4 + x^2 + 1 is (binary)10101, (octal)25, (decimal)21
// #define RDEN 2
// #define A 3345
// #define B 3613

// #define V_1 7
// #define CC_1 211
// #define CRC_1 13
// #define CRC_LENGTH_1 4

// #define V_2 6
// #define CC_2 147
// #define CRC_2 21
// #define CRC_LENGTH_2 5
// #define PUNCTURE 0
// #define PUNC_RATE 2

// // below are the null values
// #define AB 3325
// #define AC 3237
// #define DC 2711

// #define C 63
// #define D 67

// #define CRC_A 37
// #define CRC_C 51

// #define PUNC_1 0
// #define PUNC_3 5

// --------------------------- rate 1/2 into two rate 1/1
// ---------------------------

// --------------------------- rate 1/3 ---------------------------
// #define RDEN 3
// #define AB 3325
// #define AC 3237
// #define DC 2711
// #define A 45
// #define B 61
// #define C 63
// #define D 67

// #define CRC_A 37
// #define CRC_C 51

// #define PUNCTURE 0

// --------------------------- rate 1/3 ---------------------------
// --------------------------- puncture 1/2 ---------------------------

// #define V 10
// #define HALF_V 5
// #define CONSTRAINT_LENGTH (HALF_V+1)

// #define PUNC_1 0
// #define PUNC_3 5
// #define AB 2267
// #define AC 3631
// #define DC 2353
// #define A 51
// #define B 57
// #define C 61
// #define D 73

// #define CRC_A 41
// #define CRC_C 49

// #define PUNCTURE 1
// #define PUNC_RATE 2
// #define BLOCK_SIZE (PUNC_RATE * (K + V))

// // below are null values
// #define RDEN 3
// --------------------------- puncture 1/2 ---------------------------

// --------------------------- puncture 1/3, V=16 ---------------------------
// #define V 16
// #define HALF_V 8
// #define CONSTRAINT_LENGTH (HALF_V+1)

// #define PUNC_1 0
// #define PUNC_3 5
// #define AB 262221
// #define AC 333523
// #define DC 270507
// #define A 513
// #define B 447
// #define C 711
// #define D 777

// #define CRC_A 331
// #define CRC_C 457

// #define PUNCTURE 1
// #define PUNC_RATE 2
// #define BLOCK_SIZE (PUNC_RATE * (K + V))

// // below are null values
// #define RDEN 3

// --------------------------- puncture 1/3, V=16 ---------------------------

// --------------------------- puncture 1/3, V=14 ---------------------------
#define V 14
#define HALF_V 7
#define CONSTRAINT_LENGTH (HALF_V+1)

#define PUNC_1 0
#define PUNC_3 5
#define AB 51233
#define AC 44157
#define DC 76315
#define A 267
#define B 225
#define C 251
#define D 305

#define CRC_A 183
#define CRC_C 169

#define PUNCTURE 1
#define PUNC_RATE 2
#define BLOCK_SIZE (PUNC_RATE * (K + V))

// below are null values
#define RDEN 3

// --------------------------- puncture 1/3, V=14 ---------------------------

#endif