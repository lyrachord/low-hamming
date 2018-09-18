#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <gmp.h>
#include <chrono>
#include <fstream>
#include <algorithm>
#include <string>
#include "parameters.hpp"

using namespace std;

void generate_N(mpz_t, mpz_t, mpz_t);
bool test_keys(mpz_t, mpz_t);
void inverse(mpz_t, mpz_t, mpz_t);

int hamming(string);
void state_sanity();
void print_keys(char *, char *, char *, mpz_t, mpz_t, mpz_t, bool);
string format(string, bool);
void next_candidate(mpz_t, mpz_t);

void carryless_multiply(mpz_t, mpz_t);

void array_bits(mpz_t, int[]);

