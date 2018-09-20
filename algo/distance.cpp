#include "distance.hpp"

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
  generate_N_distance(P, Q, N);
  int tries = 1;
  while( test_keys_distance(P,Q)==false ){
    generate_N_distance(P, Q, N);
    tries++;
  }

  print_keys(p_str, q_str, n_str, P, Q, N, FORMAT);
  state_sanity_distance(P, Q, N);
  cout << "Found keys after "<< tries << " tries" << endl << endl;

  // Start guessing 

  mpz_t p2;
  mpz_init(p2);
  mpz_mul(p2,P,P);
  string p2_str = mpz_get_str(NULL, 2, p2);
  n_str = mpz_get_str(NULL, 2, N);

  cout << format(p2_str, FORMAT) << endl;
  cout << format(n_str, FORMAT) << endl;

  cout << mpz_hamdist(p2,N) << endl;


  return 0;
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
  
  guess_next_block(block, guess_p, real_p);

}

int guess_next_block(int block, mpz_t previous_p, mpz_t real_p){

  string success="";

  mpz_t limit;
  mpz_init(limit);
  mpz_ui_pow_ui(limit, 2, K);

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

  while ( cont < SCOPE )
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

    if(mpz_hamdist(block_guess_p, guess_q) <= MAX_WEIGHT){
      /*cout << "\033[1;34m";
      COUT  guess_p_str << " - "; COUT  guess_q_str;
      cout << "\033[0m  :  " <<  mpz_popcount(guess_q) << endl;*/
    
      if(mpz_cmp(guess_p,real_p)==0){
          //cout << "\033[1;33m" << format(guess_p_str) << "\033[0m";
        success += "\033[1;32mFound in position "+to_string(cont+1)+"\033[0m\n"; 
      }
    /*cout << endl << "\033[1;34m";
    COUT  guess_p_str << " - "; COUT  guess_q_str;
    cout << "\033[0m  :  " <<  mpz_popcount(guess_q) << endl;*/
      cont ++ ;
    }

    
    next_candidate_distance(block_guess_p, block_guess_p);
    if(mpz_cmp(block_guess_p, limit) >= 0){
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
