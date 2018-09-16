#include "main.hpp"


using namespace std;
using namespace std::chrono;


string CURRENT_DIR = "/home/botvinnik/coding/low_hamming/";


char * p_str, * q_str, * n_str ;

mpz_t P, Q, N;

bool FORMAT = false;

int main(int argc, char** argv){

  // Fetch input. 
  //char *p;
  //int block = (int) strtol(argv[1],&p,10);
  if(argc == 2)
  {
    if((string)argv[1] == "-f"){
      FORMAT = true;
    }  
  }
  int block = 0;
  
  // Generate keys
  srand(time(NULL));
  mpz_inits(P, Q, N, NULL);
  generate_N(P, Q, N);
  int tries = 1;
  while( test_keys(P,Q)==false ){
    generate_N(P, Q, N);
    tries++;
  }
  print_keys(p_str, q_str, n_str, P, Q, N, FORMAT);
  state_sanity();

  cout << "Found keys after "<< tries << " tries" << endl << endl;

  // Start guessing 

  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  
  for(int i = block; i< (int)BITS/K; i++){
    cout << "Block " << i+1 << ":";
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
  
  guess_next_block(block, guess_p, real_p);

}

int guess_next_block(int block, mpz_t previous_p, mpz_t real_p){

  string success="";

  int cont = 0;
  mpz_t guess_p, guess_q, invp, modulus, small_modulus, Nmod;
  mpz_inits(guess_p, guess_q, invp, modulus, small_modulus, Nmod, NULL);
  mpz_ui_pow_ui(modulus, 2, K*(block+1));
  mpz_ui_pow_ui(small_modulus, 2, K*(block));
  mpz_mod(Nmod, N, modulus);

  // Read guess_p from candidates list
  string filename = CURRENT_DIR + "precomp/"+to_string(K)+"-"+to_string(MAX_WEIGHT);
  ifstream candidates_file (filename);
  string line, guess_p_str, guess_q_str, guess_q_substr;

  if (candidates_file.is_open())
  {
    while ( getline (candidates_file,line) && cont < SCOPE )
    {
      mpz_set_str(guess_p, line.c_str(), 2);

      if(block == 0){
      // First block: if even add 1
        if(mpz_even_p (guess_p)){
          mpz_add_ui(guess_p, guess_p, 1);
        }
      }
      else{
      // guess_p is of the form "guess + 2^{small_mod}*previous_p"
        mpz_mul(guess_p, guess_p, small_modulus);
        mpz_add(guess_p, guess_p, previous_p);
      }
      
      guess_p_str = mpz_get_str(NULL, 2, guess_p);

      // Compute inverse
      inverse(invp, guess_p, modulus);
      
      // Compute corresponding guess_q
      mpz_mul(guess_q, invp, Nmod);
      
      mpz_mod(guess_q, guess_q, modulus);
      
      // Look bits of q in the correct range
      guess_q_str = mpz_get_str(NULL, 2, guess_q);
      guess_q_substr = guess_q_str.substr(1,K-1);
      if(hamming(guess_q_substr) <= MAX_WEIGHT){
        if(mpz_cmp(guess_p,real_p)==0){
          //cout << "\033[1;33m" << format(guess_p_str) << "\033[0m";
          success += "\033[1;32mFound in position "+to_string(cont)+"\033[0m\n"; 
        }
        
        cont ++ ;
      }
    }
    candidates_file.close();
  }
  else{
    cout << "! Error opening " << filename << endl;
  }
  if(success == ""){
    success = "\033[1;31mNot found\033[0m\n";
  }
  cout << success;
  return cont;
}
