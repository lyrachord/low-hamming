#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <gmp.h>
#include <fstream>


using namespace std;

int generate_naive();
int generate_naive_2();

// Prints the list of X with n more bits set.
void set_bits(mpz_t X, int n);

string format(string);
