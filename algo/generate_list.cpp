#include "parameters.hpp"
#include "generate_list.hpp"
#include <bitset>
using namespace std;


int main(){


  generate_naive_2();

  return 0;
}


// Still very stupid
int generate_naive_2(){

  int cont = 0;
  mpz_t X;
  mpz_init_set_ui(X,0);

  // Hamming weight 0

  for(int n = 0; n <= MAX_WEIGHT; n ++)
  {
    set_bits(X,n);
  }
  return cont;
}

void set_bits(mpz_t X, int n){

  string Xstr = format(mpz_get_str(NULL, 2, X));
  mpz_t pow, Xaux;
  mpz_inits(pow, Xaux, NULL);

  if(n == 0){
    cout << format(Xstr) << endl;
    return;
  }

  for (int i = 0; i < K; i++){
    mpz_ui_pow_ui(pow, 2, i);
    
    if(Xstr[K-1-i]== '0'){
      mpz_add(Xaux, X, pow);
      set_bits(Xaux,n-1);
    }
    
    mpz_set(Xaux,X);
  }

}


int generate_naive(){

  int cont = 0;

  // first add 0 to the list
  cout << bitset<K>(0) << endl;
  for( unsigned int weight = 1; weight <= MAX_WEIGHT; weight++){
    for(int i = 1; i < pow(2,K); i++){
      bitset<K> binary = bitset<K>(i);
      if(binary.count() == weight)
      {
        cout << binary <<endl;
        cont ++ ;
      }
    }
  }
  return cont;
}

string format(string ss){

 
  string s = ss;
  int prepend = s.length()%K;
  if(prepend != 0){
    for(int i = 0; i < K-prepend; i++){
      s = "0"+s;
    }
  }

  string ret = "";
  for(int i = 0; i<=(int) (s.length()-1)/K; i++){

    ret = ret + s.substr(i*K,K);
  }
  return ret;
}


