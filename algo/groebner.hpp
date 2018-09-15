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

void state_sanity();
void print_keys();

void carryless_multiply(mpz_t, mpz_t);

void array_bits(mpz_t, int[]);

string format(string);

