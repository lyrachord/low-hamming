#include "distance.hpp"

#define COUT std::cout.width(K);std::cout<<

using namespace std;
using namespace std::chrono;


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

  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  
  for(int i = block; i< (int)BITS/K; i++){
    cout << "Block " << i+1 << ": ";//<<endl;
    guess_block(i);
  }
  
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
    guess_first_block(real_p);
    return 0;
  }
  string success="";

  int cont = 0;
  mpz_t invp, Ndiv, small_modulus, big_modulus;
  mpz_inits(invp, Ndiv, small_modulus, big_modulus, NULL);
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
    
    if(mpz_hamdist(block_guess_q, block_guess_p) <= MAX_WEIGHT){


      if(mpz_cmp(block_guess_p,block_real_p)==0){
          //cout << "\033[1;33m" << format(guess_p_str) << "\033[0m";
        success += "\033[1;32mFound in position "+to_string(cont+1)+"\033[0m\n"; 
      }

      /*string guess_q_str = mpz_get_str(NULL, 2, block_guess_q);
      string guess_p_str = mpz_get_str(NULL, 2, block_guess_p);
      cout << endl << "\033[1;34m";
      COUT  guess_p_str << " - " << COUT  guess_q_str;
      cout << "\033[0m  :  ";*/
      cont ++ ;
    }

    mpz_add_ui(block_guess_p, block_guess_p, 1);
    if (mpz_cmp(block_guess_p, small_modulus)>0){
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


int guess_first_block(mpz_t real_p){

  string success="";

  int cont = 0;
  mpz_t invp, modulus, Nmod;
  mpz_inits(invp, modulus, Nmod, NULL);
  mpz_ui_pow_ui(modulus, 2, K);
  mpz_mod(Nmod, N, modulus);
  string guess_q_str;
  string guess_p_str;

  mpz_t block_guess_p, block_guess_q;

  // Make a guess  
  mpz_init_set_ui(block_guess_p, 1);
  mpz_init(block_guess_q);

  while ( cont < SCOPE )
  {

    // Compute inverse
    mpz_invert(invp, block_guess_p, modulus);

    // Compute corresponding guess_q
    mpz_mul(block_guess_q, invp, Nmod);
    mpz_mod(block_guess_q, block_guess_q, modulus);



    if(mpz_hamdist(block_guess_q, block_guess_p) <= MAX_WEIGHT){
      /*cout << "\033[1;34m";
      COUT  guess_p_str << " - "; COUT  guess_q_str;
      cout << "\033[0m  :  " <<  mpz_popcount(guess_q) << endl;
    */
      if(mpz_cmp(block_guess_p,real_p)==0){
          //cout << "\033[1;33m" << format(guess_p_str) << "\033[0m";
        success += "\033[1;32mFound in position "+to_string(cont+1)+"\033[0m\n"; 
      }
      // Look bits of q in the correct range, for display only
      /*string guess_q_str = mpz_get_str(NULL, 2, block_guess_q);
      string guess_p_str = mpz_get_str(NULL, 2, block_guess_p);
      cout << endl << "\033[1;34m";
      COUT  guess_p_str << " - " << COUT  guess_q_str;
      cout << "\033[0m  :  " <<endl;*/
      cont ++ ;
    }

    mpz_add_ui(block_guess_p, block_guess_p, 2);
    if (mpz_cmp(block_guess_p, modulus)>0){
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
