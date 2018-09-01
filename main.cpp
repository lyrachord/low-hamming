#include "main.hpp"

using namespace std;
using namespace std::chrono;


char * n_str ;

mpz_t P, Q, N;


int main(int argc, char** argv){

  mpz_inits(P, Q,N, NULL);
  
  mpz_set_str(P, "\
    10000000  \
    10000000  \
    10000000  \
    10000000  \
    10000000  \
    10000000  \
    10000000  \
    10000001", 2);

  mpz_set_str(Q, "\
    10000001   \
    00000001   \
    00000001   \
    00000001   \
    00000001   \
    00000001   \
    00000001   \
    00000001", 2);
  
  mpz_mul(N, P, Q);


  n_str = mpz_get_str(n_str, 2, N);

  gmp_printf("p = %Zd\nq = %Zd\nn = %Zd\n", P, Q, N);
  cout << "n[2] = " << n_str << endl;



  mpz_t guess_p, guess_q, invp, modulus;
  mpz_inits(guess_p, guess_q, invp, modulus, NULL);
  mpz_ui_pow_ui(modulus, 2, 8);

  mpz_set_str(guess_p, "00000011", 2);

  gmp_printf("Guess_p: %Zd\n", guess_p);


  inverse(invp, guess_p, modulus);

  mpz_mul(guess_q, invp, N);

  mpz_mod(guess_q, guess_q, modulus);
  gmp_printf("Guess_q: %Zd\n", guess_q);

  char * guess_q_str = mpz_get_str(NULL, 2, guess_q);

  cout << "Guess_q[2]: " << guess_q_str << endl;
  

  return 0;

}


void inverse (mpz_t inv, mpz_t x, mpz_t modulus){


  mpz_t a, b, u, x_aux;

  mpz_inits(a, b, u, x_aux, NULL);
  mpz_set(x_aux, x);

  mpz_set_ui(a, 0);
  mpz_init_set(b, modulus);
  mpz_init_set_ui(u, 1);

  mpz_t q, r;
  mpz_inits(q, r, NULL);

  mpz_t old_x, old_a, old_b, old_u;
  mpz_inits(old_x, old_a, old_b, old_u);

  while(mpz_cmp_ui(x_aux,0)>0){
    //gmp_printf("x = %Zd, a = %Zd, b = %Zd, u = %Zd\n",x_aux,a,b,u);

    // Save values for multiple affectation
    mpz_set(old_x, x_aux);
    mpz_set(old_a, a);
    mpz_set(old_b, b);
    mpz_set(old_u, u);

    mpz_fdiv_q(q, old_b, old_x);
    
    mpz_mod(x_aux, old_b, old_x);
    mpz_set(a, old_u);
    mpz_set(b, old_x);
    mpz_mul(u, q, old_u);
    mpz_neg(u, u);
    mpz_add(u, old_a, u);

  }

  if(mpz_cmp_ui(b,1)==0){
    mpz_mod(inv, a, modulus);
    return;
  }

  cout << "Error in division (not coprime)" << endl;
  return;
}

// Extended Euclides
/*
//https://stackoverflow.com/questions/14736514/optimising-code-for-modular-arithmetic

function inverse(x, m)
    a, b, u := 0, m, 1
    while x > 0
        q, r := divide(b, x)
        x, a, b, u := b % x, u, x, a - q * u
    if b == 1 return a % m
    error "must be coprime"
*/
