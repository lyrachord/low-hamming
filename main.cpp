#include "main.hpp"

using namespace std;
using namespace std::chrono;

int main(int argc, char** argv){

  srand(time(NULL));
  gmp_randstate_t state;
  gmp_randinit_mt (state);
  gmp_randseed_ui (state, rand());

  cout << "Hello" << endl;


  int z;
  for(z=24; z<25; z++)
  {
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    multi(state, z);
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<seconds>( t2 - t1 ).count();

    cout << "Hadamard - Multiplied 2^" << z << " = " << pow(2,z) << " polynomials in Z_{2^"<<BITS<<"}/(X^"<<N<<"+1) " << endl;
    cout << "in " << duration << " seconds " << endl;

  }
  
  return 0;

}
void multi(gmp_randstate_t state, int pols)
{


  mpz_t a,b,res;
  mpz_inits (a,b,res,NULL);
  mpz_urandomb(a, state, BITS);
  mpz_urandomb(b, state, BITS);
  //gmp_printf("a = %Zd\n\nb = %Zd\n\n",a,b);
  
  // Multiply 2**pols polynomials
  int j=0;
  long long int quantity = pow(2,pols);
  for(j=0; j < quantity; j++){
    // Multiply N integers of BITS bits
    int i;
    for(i=0; i < N; i++ ){
      mpz_mul(res,a,b);
    }
  }
  
}

