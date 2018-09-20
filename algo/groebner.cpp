#include "groebner.hpp"
#include "parameters.hpp"
using namespace std;
using namespace std::chrono;


char * p_str, * q_str, * n_str ;

mpz_t P, Q, N;

bool FORMAT = true;

int main(int argc, char** argv){

  srand(time(NULL));
  // Fetch input. 
  
  if(argc == 2)
  {
    if((string)argv[1] == "-f"){
      FORMAT = false;
    }  
  }

  gmp_randstate_t state;
  gmp_randinit_default (state);
  
  mpz_t X;
  mpz_init(X);

  mpz_t X0;
  mpz_init(X0);

  mpz_urandomb(X, state, BITS);
  while(mpz_probab_prime_p (X, 50) > 0){
    mpz_urandomb(X, state, BITS);
    gmp_printf("%Zd\n", X);
  }
  mpz_set(X0,X);

  int h0 = mpz_popcount(X);
  int h = h0;
  int cont = 1;
  while(h>3*h0/4)
  {
    h = mpz_popcount(X);
    cout << h << endl;
    //string x_str = mpz_get_str(NULL, 2, X);
    //gmp_printf("%Zd,", X);
    //cout << format(x_str, FORMAT);
    mpz_add(X, X, X0);
    cont ++;
  }
  cout << cont << ", "<< h0 << ", " << h <<endl;
  
  return 0;

}
