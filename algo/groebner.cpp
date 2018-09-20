#include "groebner.hpp"
#include "parameters.hpp"
using namespace std;
using namespace std::chrono;


char * p_str, * q_str, * n_str ;

mpz_t P, Q, N;

bool FORMAT = false;

int main(int argc, char** argv){

  srand(time(NULL));
  // Fetch input. 
  
  if(argc == 2)
  {
    if((string)argv[1] == "-f"){
      FORMAT = true;
    }  
  }

  int big = 512;
  int small = 512;

  gmp_randstate_t state;
  gmp_randinit_default (state);
  
  mpz_t x,X, aux;

  mpz_inits(x,X, aux, NULL);
  

  mpz_t modulus;
  mpz_init(modulus);
  mpz_ui_pow_ui(modulus, 2, big);



  mpz_t small_modulus;
  mpz_init(small_modulus);
  mpz_ui_pow_ui(small_modulus, 2, small);


  mpz_urandomb(X, state, big);
  mpz_urandomb(aux, state, small);

  mpz_mod(x, X, small_modulus);

  int iterations = 1;

  high_resolution_clock::time_point t1a = high_resolution_clock::now();
  
  // A: Invert in R_modulus, multiply and reduce
  for(int i=0; i < iterations; i++)
  {
    gmp_printf("X = %Zd\n", X);
    mpz_invert(X, X, modulus);
    gmp_printf("X inverted = %Zd\n", X);
    mpz_mul(X, X, X);
    gmp_printf("X squared = %Zd\n", X);
    mpz_mod(X, X, modulus);    
    gmp_printf("X reduced = %Zd\n", X);
  }
  
  high_resolution_clock::time_point t2a = high_resolution_clock::now();
  auto duration_a = duration_cast<microseconds>( t2a - t1a ).count();

  cout << "Completed in "<< duration_a <<" microseconds." << endl;


  high_resolution_clock::time_point t1b = high_resolution_clock::now();
  
  // B: Multiply, substract, invert
  for(int i=0; i < iterations; i++)
  {
    gmp_printf("x = %Zd\n", x);
    mpz_mul(x,aux,x);
    gmp_printf("x squared = %Zd\n", x);
    mpz_sub(x,x,aux);
    gmp_printf("x minus = %Zd\n", x);
    mpz_invert(x,x, small_modulus);
    gmp_printf("x inverted = %Zd\n", x);

  }
  high_resolution_clock::time_point t2d = high_resolution_clock::now();
  auto duration_b = duration_cast<microseconds>( t2d - t1b ).count();

  cout << "Completed in "<< duration_b <<" microseconds." << endl;
  


  return 0;

}
