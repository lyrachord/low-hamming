#include "main.hpp"

using namespace std;
// Number of bits of each factor


char * p_str, * q_str, * n_str ;

mpz_t P, Q, N;


int main(int argc, char** argv){

  //srand(time(NULL));
  mpz_inits(P, Q, N, NULL);
  
  mpz_set_ui(P,1);
  mpz_set_ui(Q,1);

  mpz_t aux;
  mpz_init(aux);
  for(int i = 1; i < BITS-1; i++ )
  {
    if(rand()%10==0){
      mpz_ui_pow_ui(aux,2,i);
      mpz_add(P,P,aux);
    }
    if(rand()%10==0){
      mpz_ui_pow_ui(aux,2,i);
      mpz_add(Q,Q,aux);
    }
  }
  mpz_ui_pow_ui(aux,2,BITS-1);
  mpz_add(P,P,aux);
  mpz_add(Q,Q,aux);


  mpz_mul(N, P, Q);


  p_str = mpz_get_str(p_str, 2, P);
  q_str = mpz_get_str(q_str, 2, Q);
  n_str = mpz_get_str(n_str, 2, N);

  gmp_printf("p = %Zd\nq = %Zd\nn = %Zd\n", P, Q, N);
  cout << "p[2] = " << p_str << endl;
  cout << "q[2] = " << q_str << endl;
  cout << "n[2] = " << n_str << endl;

  mpz_t guess_p, guess_q;
  mpz_inits(guess_p, guess_q, NULL);

  int block = 1;
  mpz_mod_ui(guess_p, P, pow(2,block*K));
  mpz_mod_ui(guess_q, Q, pow(2,block*K));
  //mpz_set_ui(guess_p,0);
  //mpz_set_ui(guess_q,0);

  char * guess_p_str = mpz_get_str(NULL, 2, guess_p);
  cout << guess_p_str << endl;
  
  guess_next_block(block, guess_p);
  return 0;

}



int guess_next_block(int block, mpz_t previous_p){


  mpz_t guess_p, guess_q, invp, modulus, small_modulus;
  mpz_inits(guess_p, guess_q, invp, modulus, small_modulus, NULL);
  mpz_ui_pow_ui(modulus, 2, K*(block+1));
  mpz_ui_pow_ui(small_modulus, 2, K*(block));

  // Read guess_p from candidates list
  mpz_set_str(guess_p, "00000000000", 2);
  mpz_mul(guess_p,guess_p,small_modulus);
  mpz_add(guess_p,guess_p,previous_p);
  gmp_printf("Guess_p: %Zd\n", guess_p);
  char * guess_p_str = mpz_get_str(NULL, 2, guess_p);
  cout << guess_p_str << endl;
  // Compute inverse
  inverse(invp, guess_p, modulus);

  // Compute equivalent for q
  mpz_mul(guess_q, invp, N);

  // Look bits of q
  // CONTINUE HERE
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
