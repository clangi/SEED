#include "rnd_namespace.h"

namespace rnd_gen {

  void set_rng_seed(int s){
    rnd_gen::rng.seed(s);
  }

  double get_uniform(double min, double max){
    boost::random::uniform_real_distribution<double> dist(min, max);
    return dist(rnd_gen::rng);
  }

}
