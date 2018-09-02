#include "main.hpp"
#include "parameters.hpp"
using namespace std;

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
  generate_N();
  print_keys();
  state_sanity();


  for(int i = block; i< (int)BITS/K; i++){
    cout << "Block " << i+1 << ":";
    guess_block(i);
  }
  

  return 0;

}

void state_sanity(){

  cout << "log(p) = " << BITS << endl;
  cout << "block_size = " << K << endl;
  cout << "number_of_blocks = " << (int) BITS/K << endl;
  cout << "list_scope = " << SCOPE << endl;

  // Algorithms runs in scope^(number of blocks).
  // Display log2() of that number. Should be less than 64.

  mpz_t sanity;
  mpz_init(sanity);
  mpz_ui_pow_ui(sanity, SCOPE, (int)BITS/K);
    

  cout << endl << "\033[1;31mRuntime: 2^"<< mpz_sizeinbase (sanity, 2) << "\033[0m" << endl << endl;


}


void guess_block(int block){

    mpz_t guess_p, guess_q, real_p, small_modulus, modulus;
  mpz_inits(guess_p, guess_q, real_p, small_modulus, modulus, NULL);
  

  mpz_ui_pow_ui(modulus, 2, K*(block+1));
  mpz_ui_pow_ui(small_modulus, 2, K*(block));
  mpz_mod(guess_p, P, small_modulus);
  mpz_mod(real_p, P, modulus);


  
  //mpz_set_ui(guess_p,0);
  

  
  guess_next_block(block, guess_p, real_p);
  //cout << "candidates: " << candidates << endl;

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
      guess_q_substr = guess_q_str.substr(0,K);
      if(hamming(guess_q_substr) <= MAX_WEIGHT){
        if(mpz_cmp(guess_p,real_p)==0){
          //cout << "\033[1;33m" << format(guess_p_str) << "\033[0m";
          
          success += "\033[1;32mFound in position "+to_string(cont)+"\033[0m\n"; 
        }
        else{
          //cout << "Guess_p[2]:" << format(guess_p_str);
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

void generate_N(){

  mpz_inits(P, Q, N, NULL);

  //Sample P,Q  
  mpz_set_ui(P,1);
  mpz_set_ui(Q,1);

  mpz_t aux;
  mpz_init(aux);
  for(int i = 1; i < BITS-1; i++ )
  {
    if(rand()%PROB==0){
      mpz_ui_pow_ui(aux,2,i);
      mpz_add(P,P,aux);
    }
    if(rand()%PROB==0){
      mpz_ui_pow_ui(aux,2,i);
      mpz_add(Q,Q,aux);
    }
  }
  mpz_ui_pow_ui(aux,2,BITS-1);
  mpz_add(P,P,aux);
  mpz_add(Q,Q,aux);
  
  mpz_mul(N, P, Q);



}

void print_keys(){

  p_str = mpz_get_str(p_str, 2, P);
  q_str = mpz_get_str(q_str, 2, Q);
  n_str = mpz_get_str(n_str, 2, N);

  gmp_printf("p = %Zd\nq = %Zd\nn = %Zd\n", P, Q, N);
  cout << "p[2] = " << format(p_str) << endl;
  cout << "q[2] = " << format(q_str) << endl;
  cout << "n[2] = " << format(n_str) << endl;
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
  mpz_inits(old_x, old_a, old_b, old_u, NULL);
  

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

int hamming(string s){
  return count(s.begin(), s.end(), '1');
}

string format(string ss){

  if (FORMAT == false) return ss+"\n";
  string s = ss;
  int prepend = s.length()%K;
  if(prepend != 0){
    for(int i = 0; i < K-prepend; i++){
      s = "0"+s;
    }
  }

  string ret = "\n";
  for(int i = 0; i<=(int) (s.length()-1)/K; i++){

    ret = ret + s.substr(i*K,K)+"\n";
  }
  return ret;
}