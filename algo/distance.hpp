#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <gmp.h>
#include <chrono>
#include <fstream>
#include <algorithm>
#include "tools.hpp"
#include "parameters.hpp"

using namespace std;

void guess_block(int);
int guess_next_block(int, mpz_t, mpz_t);