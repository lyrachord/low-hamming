#include "tools.hpp"
using namespace std;
using namespace std::chrono;

void state_sanity_weight(mpz_t P,mpz_t Q,mpz_t N){

  int blocks = (int) BITS/K;

  cout << "log(p) = " << BITS << endl;
  cout << "hamming(p) = " << mpz_popcount(P) << endl;
  cout << "hamming(q) = " << mpz_popcount(Q) << endl;
  cout << "block_size = " << K << endl;
  cout << "number_of_blocks = " << blocks << endl;
  cout << "list_scope = " << SCOPE << endl;

  // Algorithms runs in scope^(number of blocks).
  // Display log2() of that number. Should be less than 64.

  mpz_t sanity;
  mpz_init(sanity);
  mpz_ui_pow_ui(sanity, SCOPE, (int)BITS/K);


  cout << endl << "\033[1;31mTotal Runtime: 2^"<< mpz_sizeinbase (sanity, 2) << "\033[0m" << endl << endl;
}

void state_sanity_distance(mpz_t P,mpz_t Q,mpz_t N){

  int blocks = (int) BITS/K;

  cout << "log(p) = " << BITS << endl;
  cout << "hamming_dist(p, q) = " << mpz_hamdist(P,Q) << endl;
  cout << "block_size = " << K << endl;
  cout << "number_of_blocks = " << blocks << endl;
  cout << "hamming = " << blocks*MAX_WEIGHT << endl;
  cout << "list_scope = " << SCOPE << endl;

  // Algorithms runs in scope^(number of blocks).
  // Display log2() of that number. Should be less than 64.

  mpz_t sanity;
  mpz_init(sanity);
  mpz_ui_pow_ui(sanity, SCOPE, (int)BITS/K);


  cout << endl << "\033[1;31mTotal Runtime: 2^"<< mpz_sizeinbase (sanity, 2) << "\033[0m" << endl << endl;
}


void generate_N_weight(mpz_t P, mpz_t Q, mpz_t N){

  //Sample P,Q  
  mpz_set_ui(P,1);
  mpz_set_ui(Q,1);

  mpz_t aux;
  mpz_init(aux);

  int index_p, index_q;

  for(int i = 1; i < SET_BITS-1; i++ )
  {
    index_p = 0;
    index_q = 0;

    while(mpz_tstbit(P, index_p)==1){
      index_p = rand()%(BITS-1) + 1;
    }
    while(mpz_tstbit(Q, index_q)==1){
      index_q = rand()%(BITS-1) + 1;
    }
    mpz_setbit(P, index_p);
    mpz_setbit(Q, index_q);
  }

  mpz_setbit(P, BITS-1);
  mpz_setbit(Q, BITS-1);
  
  mpz_mul(N, P, Q);
}



void generate_N_weight_poisson(mpz_t P, mpz_t Q, mpz_t N){

  //Sample P,Q  
  mpz_set_ui(P,1);
  mpz_set_ui(Q,1);

  mpz_t aux;
  mpz_init(aux);

  for(int i = 1; i < BITS-2; i++ )
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

void generate_N_distance(mpz_t P, mpz_t Q, mpz_t N){

  gmp_randstate_t state;
  gmp_randinit_default (state);
  //Sample P 
  mpz_urandomb(P, state, BITS);

  mpz_setbit (P, 0);
  mpz_setbit (P, BITS-1);

  mpz_set(Q,P);

  for(int i = 1; i < BITS-2; i++ )
  {
    if(rand()%PROB==0){
      mpz_combit(Q,i);   
    }
  }
  mpz_mul(N, P, Q);
}

void print_keys(
  char * p_str, char * q_str, char * n_str,
  mpz_t P, mpz_t Q, mpz_t N, bool FORMAT){

  p_str = mpz_get_str(p_str, 2, P);
  q_str = mpz_get_str(q_str, 2, Q);
  n_str = mpz_get_str(n_str, 2, N);

  gmp_printf("p = %Zd\nq = %Zd\nn = %Zd\n", P, Q, N);
  cout << "p[2] = " << format(p_str, FORMAT) << endl;
  cout << "q[2] = " << format(q_str, FORMAT) << endl;
  cout << "n[2] = " << format(n_str, FORMAT) << endl;
}

void inverse(mpz_t inv, mpz_t x, mpz_t modulus){
  mpz_invert(inv, x, modulus);
}

void inverse_SE (mpz_t inv, mpz_t x, mpz_t modulus){


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

int hamming(string s){
  return count(s.begin(), s.end(), '1');
}


int hamming_distance(string s1, string s2){

  int diff_bits = 0;
  for(unsigned int i=0; i < s1.length(); i++){
    if(s1[i] != s2[i]){
      diff_bits ++;
    }
  }
  return diff_bits;
}

string format(string ss, bool FORMAT){

  if (FORMAT == false) return ss+"\n";
  string s = ss;
  int prepend = s.length()%K;
  if(prepend != 0){
    for(int i = 0; i < K-prepend; i++){
      s = "0"+s;
    }
  }

  string ret = "\n";
  int block_number;
  for(int i = 0; i<=(int) (s.length()-1)/K; i++){

    block_number = (int) (s.length()-1)/K - i + 1;
    ret = ret +  s.substr(i*K,K)+ "  /"+ to_string(block_number)+"\n";
  }
  return ret;
}

bool test_keys_weight(mpz_t p, mpz_t q){

  string prepend = "";
  for(int i = 0; i< K- BITS % K; i++){
    prepend += "0";  
  }

  string p_str = prepend + mpz_get_str(NULL, 2, p);
  string q_str = prepend + mpz_get_str(NULL, 2, q);

  string limb_p;
  string limb_q;

  for(int i = 0; i <= (int) BITS/K; i++){
    limb_p = p_str.substr(i*K,K);
    limb_q = q_str.substr(i*K,K);
    if(hamming(limb_p) > MAX_WEIGHT || hamming(limb_q) > MAX_WEIGHT){
      return false;
    }

  }
  return true;
}

bool test_keys_distance(mpz_t p, mpz_t q){

  string prepend = "";
  for(int i = 0; i< K- BITS % K; i++){
    prepend += "0";  
  }

  string p_str = prepend + mpz_get_str(NULL, 2, p);
  string q_str = prepend + mpz_get_str(NULL, 2, q);

  string limb_p;
  string limb_q;

  for(int i = 0; i <= (int) BITS/K; i++){
    limb_p = p_str.substr(i*K,K);
    limb_q = q_str.substr(i*K,K);
    if(hamming_distance(limb_p, limb_q) > MAX_WEIGHT){
      return false;
    }
  }
  return true;
}

void carryless_multiply(mpz_t x, mpz_t y){

  cout << "Carryless multiplication of p,q" << endl;
  
  int bitsx[BITS];
  int bitsy[BITS];
  array_bits(x, bitsx);
  array_bits(y, bitsy);

  int res[2*BITS];  

  for(int i = 0; i < 2*BITS; i++)
  { 
    res[i] = 0;
    for(int j = 0; j <= i; j++)
    {
      if(j < BITS && i-j < BITS) {res[i] = res[i]+bitsx[j]*bitsy[i-j];}
    }
  }


  
  for(int i=0; i<2*BITS-2; i++)
  {
    cout << res[i];
  }
  cout << res[2*BITS-2] <<  endl;
  
} 

void array_bits(mpz_t x, int res[BITS]){

  mpz_t bits_aux[BITS];
  for(int i=0; i < BITS; i++){
    mpz_init(bits_aux[i]);
  }

  mpz_t power_of_two;
  mpz_init_set_ui(power_of_two, 2);
  mpz_t aux, S;
  mpz_inits(aux, S, NULL);

  for(int i = 0; i < BITS; i++ ){
    mpz_mod(aux, x, power_of_two);
    mpz_sub(aux, aux, S);
    mpz_add(S, S, aux);    
    mpz_set(bits_aux[i],aux);
    mpz_mul_ui(power_of_two,power_of_two, 2);
  }

  mpz_set_ui(power_of_two,1);
  for(int i = 0; i < BITS; i++ ){
    mpz_div(bits_aux[i], bits_aux[i], power_of_two);
    mpz_mul_ui(power_of_two,power_of_two, 2);

    res[i] = mpz_get_ui(bits_aux[i]);
  }

}

void next_candidate_weight(mpz_t x, mpz_t next){


  mpz_t power_of_two, limit;
  mpz_inits(power_of_two, limit, NULL);

  // If x is 0
  if(mpz_popcount(x) == 0){
    mpz_set_ui(next, 1);
    return;
  }
  
  mpz_t result;
  mpz_init_set_ui(result, 0);


  // Add least set bit
  mpz_set(result, x);
  int least_set_bit = mpz_scan1(result, 0);

  mpz_ui_pow_ui(power_of_two, 2, least_set_bit);
  mpz_add(result, result, power_of_two);

  int hx = mpz_popcount(x);
  int hresult = mpz_popcount(result);
  
  mpz_ui_pow_ui(limit, 2, K);  

  // If limit is reached
  if(mpz_tstbit(result, K)==1)
  { 
    mpz_set_ui(result, 0);
    for (int i = 0; i < hx+1; i++)
    {
      mpz_setbit(result, i);
    }
    mpz_set(next, result);
    return;    
  }

  if(hresult == hx){

    mpz_set(next, result);
    return;
  }
  else{
    int lost_bits = hx - hresult;
    for (int i = 0; i < lost_bits; i++)
    {
      mpz_setbit(result, i);
    }
    mpz_set(next, result);
  }
}

void next_candidate_distance(mpz_t x, mpz_t next){

  mpz_add_ui(next, x, 1);
}