#include "weight.hpp"

#define COUT std::cout.width(K);std::cout<<

using namespace std;
using namespace std::chrono;


string CURRENT_DIR = "/home/botvinnik/coding/low_hamming/";


char * p_str, * q_str, * n_str ;

mpz_t P, Q, N;

bool FORMAT = true;

int main(int argc, char** argv){

  std::cout.fill( '0' );

  // Fetch input. 
  //char *p;
  //int block = (int) strtol(argv[1],&p,10);
  if(argc == 2)
  {
    if((string)argv[1] == "-f"){
      FORMAT = false;
    }  
  }
  int block = 0;

  // Generate keys
  srand(time(NULL));
  mpz_inits(P, Q, N, NULL);
  generate_N_weight(P, Q, N);
  int tries = 1;
  while( test_keys_weight(P,Q)==false ){
    generate_N_weight(P, Q, N);
    tries++;
  }

  print_keys(p_str, q_str, n_str, P, Q, N, FORMAT);
  state_sanity_weight(P, Q, N);
  cout << "Found keys after "<< tries << " tries" << endl << endl;

  // Start guessing 

  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  
  for(int i = block; i< (int)BITS/K; i++){
    cout << "Block " << i+1 << ": ";//<<endl;
    guess_block(i);
  }
  //guess_block(3);
  
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>( t2 - t1 ).count();

  cout << "Completed in "<< duration <<" microseconds." << endl;
  return 0;

}


void guess_block(int block){

  mpz_t guess_p, guess_q, real_p, small_modulus, modulus;
  mpz_inits(guess_p, guess_q, real_p, small_modulus, modulus, NULL);

  mpz_ui_pow_ui(modulus, 2, K*(block+1));
  mpz_ui_pow_ui(small_modulus, 2, K*(block));
  mpz_mod(guess_p, P, small_modulus);
  mpz_mod(real_p, P, modulus);
  mpz_mod(guess_q, Q, small_modulus);

  
  //old_guess_next_block(block, guess_p, real_p, guess_q);
  guess_next_block(block, guess_p, real_p, guess_q);
}

int guess_next_block(int block, mpz_t previous_p, mpz_t real_p, mpz_t previous_q){

  cout << endl;

  if(block == 0){
    old_guess_next_block(block, previous_p, real_p, previous_q);
    return 0;
  }
  string success="";

  int cont = 0;
  mpz_t invp, Ndiv, next_p, small_modulus, big_modulus;
  mpz_inits(invp, Ndiv, next_p, small_modulus, big_modulus, NULL);
  mpz_ui_pow_ui(small_modulus, 2, K);
  mpz_ui_pow_ui(big_modulus, 2, K*block);
  
  // Store the real string of bits of p in this range, in order to compare later
  mpz_t block_real_p;
  mpz_init(block_real_p);
  mpz_tdiv_r_2exp(block_real_p, real_p, K*(block+1));
  mpz_tdiv_q_2exp(block_real_p, real_p, K*block);


  mpz_t block_guess_p, block_guess_q;
  mpz_inits(block_guess_p, block_guess_q, NULL);

  // Make a guess on p
  mpz_set_ui(block_guess_p, 0);
  

  mpz_t aux;
  mpz_init(aux);

  // We have 
  // guess_block_q = 
  //     previous_p^{-1}( (N-previous_p*previous_q)/big_mod - previous_q*guess_block_p) mod 2^{small modulus}

  // Compute inverse
  mpz_invert(invp, previous_p, small_modulus);

  // Compute (N- previous_p * previous_q) / big modulus mod small_modulus
  mpz_mul(aux, previous_p, previous_q);
  mpz_sub(Ndiv, N, aux);
  mpz_tdiv_q_2exp(Ndiv, Ndiv, K*block);
  mpz_mod(Ndiv, Ndiv, small_modulus);


  while ( cont < SCOPE )
  {
    mpz_mul(aux, previous_q, block_guess_p);
    mpz_mod(aux, aux, small_modulus);
    mpz_sub(aux, Ndiv, aux);
    mpz_mul(block_guess_q, aux, invp);
    mpz_mod(block_guess_q, block_guess_q, small_modulus);
    
    if(mpz_popcount(block_guess_q) <= MAX_WEIGHT){

      if(mpz_cmp(block_guess_p,block_real_p)==0){
          //cout << "\033[1;33m" << format(guess_p_str) << "\033[0m";
        success += "\033[1;32mFound in position "+to_string(cont+1)+"\033[0m\n"; 
      }
    /*cout << endl << "\033[1;34m";
    COUT  guess_p_str << " - "; COUT  guess_q_str;
    cout << "\033[0m  :  " <<  mpz_popcount(guess_q) << endl;*/
      cont ++ ;
    }

    
    next_candidate_weight(block_guess_p, block_guess_p);
    if (mpz_popcount(block_guess_p)>MAX_WEIGHT){
      break;
    } 
  }

  cout << cont << " candidates"<<endl;


  if(success == ""){
    success = "\033[1;31mNot found\033[0m\n";
  }
  cout << success;
  return cont;
}


int old_guess_next_block(int block, mpz_t previous_p, mpz_t real_p, mpz_t previous_q){

  string success="";

  int cont = 0;
  mpz_t guess_p, guess_q, invp, modulus, small_modulus, Nmod, next_p;
  mpz_inits(guess_p, guess_q, invp, modulus, small_modulus, Nmod, next_p, NULL);
  mpz_ui_pow_ui(modulus, 2, K*(block+1));
  mpz_ui_pow_ui(small_modulus, 2, K*(block));
  mpz_mod(Nmod, N, modulus);
  string guess_q_str;
  string guess_p_str;

  // Make a guess
  mpz_t block_guess_p;

  if(block == 0){
    // First block
    mpz_init_set_ui(block_guess_p, 1);
  }
  else{
    mpz_init_set_ui(block_guess_p, 0);
  }

  int scope = SCOPE;
  if (block==0) {scope = SCOPE*SCOPE;}
  while ( cont < scope )
  {

    // guess_p is of the form "guess + 2^{small_mod}*previous_p"
    mpz_mul(guess_p, block_guess_p, small_modulus);
    mpz_add(guess_p, guess_p, previous_p);

    // Compute inverse

    mpz_invert(invp, guess_p, modulus);

    // Compute corresponding guess_q
    mpz_mul(guess_q, invp, Nmod);
    mpz_mod(guess_q, guess_q, modulus);

    // Compute pertinent parts of q: (quotient, not rest of) guess_q / small_mod
    mpz_tdiv_q_2exp(guess_q, guess_q, K*(block));

    // Look bits of q in the correct range
    guess_q_str = mpz_get_str(NULL, 2, guess_q);
    guess_p_str = mpz_get_str(NULL, 2, block_guess_p);

    if(mpz_popcount(guess_q) <= MAX_WEIGHT){
      /*cout << "\033[1;34m";
      COUT  guess_p_str << " - "; COUT  guess_q_str;
      cout << "\033[0m  :  " <<  mpz_popcount(guess_q) << endl;
    */
      if(mpz_cmp(guess_p,real_p)==0){
          //cout << "\033[1;33m" << format(guess_p_str) << "\033[0m";
        success += "\033[1;32mFound in position "+to_string(cont+1)+"\033[0m\n"; 
      }
    /*cout << endl << "\033[1;34m";
    COUT  guess_p_str << " - "; COUT  guess_q_str;
    cout << "\033[0m  :  " <<  mpz_popcount(guess_q) << endl;*/
      cont ++ ;
    }

    
    next_candidate_weight(block_guess_p, block_guess_p);
    if (mpz_popcount(block_guess_p)>MAX_WEIGHT){
      break;
    } 
  }

  cout << cont << " candidates"<<endl;


  if(success == ""){
    success = "\033[1;31mNot found\033[0m\n";
  }
  cout << success;
  return cont;
}

