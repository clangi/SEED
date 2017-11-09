/* parameter class:
   for the moment we put inside here only the MC parameters.
   It would be nice for the future to completely migrate all the
   params in this class */
#ifndef PARAMETER_H
#define PARAMETER_H

class Parameter // This is much more like a POD struct than a class
{
  public:
  Parameter();  // default constructor
  ~Parameter(); // default destructor

  char do_mc;
  double mc_temp;
  double mc_max_tran_step;
  double mc_max_rot_step;
  double mc_rot_freq;
  int mc_niter;
  int mc_el_freq;
  double sa_alpha;
  int mc_rand_seed;
};

#endif
