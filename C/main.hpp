#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <gmp.h>
#include <chrono>
#include <fstream>
#include <algorithm>

using namespace std;

// Number of bits of each factor
#define BITS 512

// One bit "1" in each PROB bits
#define PROB 16

// Block size
#define K 16

// Max weight per block
#define MAX_WEIGHT 5

// How many candidates to test per list
#define SCOPE 256


void generate_N();
void inverse(mpz_t, mpz_t, mpz_t);
void guess_block(int);
int guess_next_block(int, mpz_t, mpz_t);
int hamming(string);
string format(string);