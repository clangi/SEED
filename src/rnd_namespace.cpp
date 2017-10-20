#include "rnd_namespace.h"

namespace rnd_gen {

  void set_rng_seed(int s){
    rnd_gen::rng.seed(s);
  }

  double get_uniform(double min, double max){
    boost::random::uniform_real_distribution<double> dist(min, max);
    return dist(rnd_gen::rng);
  }

  int get_uniform_int0(int max){
    boost::random::uniform_int_distribution<> dist(0, max);
    return dist(rnd_gen::rng);
  }

  int get_uniform_int(int min, int max){
    boost::random::uniform_int_distribution<> dist(min, max);
    return dist(rnd_gen::rng);
  }


}
