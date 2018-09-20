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

void generate_N_weight(mpz_t, mpz_t, mpz_t);
void generate_N_distance(mpz_t, mpz_t, mpz_t);

bool test_keys_weight(mpz_t, mpz_t);
bool test_keys_distance(mpz_t, mpz_t);

bool test_candidate(mpz_t);

void inverse(mpz_t, mpz_t, mpz_t);

int hamming(string);
void state_sanity_weight(mpz_t, mpz_t, mpz_t);
void state_sanity_distance(mpz_t, mpz_t, mpz_t);
void print_keys(char *, char *, char *, mpz_t, mpz_t, mpz_t, bool);
string format(string, bool);
void next_candidate_weight(mpz_t, mpz_t);
void next_candidate_distance(mpz_t, mpz_t);

void carryless_multiply(mpz_t, mpz_t);

void array_bits(mpz_t, int[]);

