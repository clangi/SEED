#include "../src/rnd_namespace.h"

using namespace std;

int main(int argc,char *argv[]){
  int my_seed = 11;
  double my_max = 500.0;
  double aa;

  // rnd_gen::set_rng_seed(my_seed);
  // for(int i = 0; i < 101; i++){
  //    aa = rnd_gen::get_uniform(0.0,my_max);
  // }

  rnd_gen::set_rng_seed(10);
  for (int i = 0; i < 10; i++){
    cout << rnd_gen::get_uniform_int0(10) << "\n";
  }

  for(int i = 0; i < 1; i++){
    cout << rnd_gen::get_uniform(0.0,my_max) << "\n";
  }

  for (int i = 0; i < 10; i++){
    cout << rnd_gen::get_uniform_int0(10) << "\n";
  }

}
