#include "main.hpp"
#include <bitset>
using namespace std;


int main(){

  int cont = 0;


  for( unsigned int weight = 1; weight <= MAX_WEIGHT; weight++)

  for(int i = 1; i < pow(2,K); i++){
    bitset<K> binary = bitset<K>(i);
    if(binary.count() == weight)
    {
      cout << binary <<endl;
      cont ++ ;
    }
  }

  return 0;
}