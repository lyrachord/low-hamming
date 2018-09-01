#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <gmp.h>
#include <chrono>

// Number of bits of each factor
#define BITS 160

// Hamming weight
#define H 7

// Block size
#define K 10

// Max weight per block
#define MAX_WEIGHT 4



void inverse(mpz_t, mpz_t, mpz_t);
int guess_next_block(int, mpz_t);
