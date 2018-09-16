#include "groebner.hpp"
#include "parameters.hpp"
using namespace std;
using namespace std::chrono;


char * p_str, * q_str, * n_str ;

mpz_t P, Q, N;

bool FORMAT = false;

int main(int argc, char** argv){

  // Fetch input. 
  
  if(argc == 2)
  {
    if((string)argv[1] == "-f"){
      FORMAT = true;
    }  
  }
  
  // Generate keys
  srand(time(NULL));
  int tries = 1;
  mpz_inits(P, Q, N, NULL);
  generate_N(P, Q, N);

  while( test_keys(P,Q)==false ){
    generate_N(P, Q, N);
    tries++;
  }
  cout << "Found keys after "<< tries << " tries" << endl << endl;

  print_keys(p_str, q_str, n_str, P, Q, N, FORMAT);
  state_sanity();


  carryless_multiply(P,Q);

  return 0;

}
