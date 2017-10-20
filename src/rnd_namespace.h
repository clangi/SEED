#ifndef _RANDOM_NAMESPACE_H
#define _RANDOM_NAMESPACE_H
#include "boost/random.hpp"

/* globals and functions for random number generation */
namespace rnd_gen {
  boost::random::mt19937 rng;
  void set_rng_seed();
  double get_uniform(double min, double max);
  int get_uniform_int0(int max);
  int get_uniform_int(int min, int max);
  bool toss_coin();
}

 #include "rnd_namespace.cpp"

#endif
