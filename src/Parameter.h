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
  double mc_max_tran_step_fine; // for finer trans move
  double mc_tran_freq; // freq for tran move
  double mc_tran_fine_freq; // rel. freq for finer trans move

  double mc_max_rot_step;
  double mc_max_rot_step_fine; // for finer rot move
  double mc_rot_freq; // freq for rot move (tran_freq = 1 - mc_rot_freq)
  double mc_rot_fine_freq; // rel. freq for finer rot move

  int mc_niter;
  int mc_el_freq;
  double sa_alpha;
  int mc_rand_seed;
};

#endif
