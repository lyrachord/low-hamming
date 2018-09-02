#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <gmp.h>
#include <chrono>
#include <fstream>
#include <algorithm>

using namespace std;

void generate_N();
bool test_keys(mpz_t, mpz_t);

void inverse(mpz_t, mpz_t, mpz_t);
void guess_block(int);
int guess_next_block(int, mpz_t, mpz_t);
int hamming(string);

void state_sanity();
void print_keys();



string format(string);

