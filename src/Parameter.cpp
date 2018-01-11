#include "Parameter.h"

Parameter::Parameter() : do_mc('n'), mc_temp(0.0), mc_max_tran_step(0.0),
  mc_max_tran_step_fine(0.0),
  mc_tran_freq(0.5),
  mc_tran_fine_freq(0.5),
  mc_max_rot_step(0.0),
  mc_max_rot_step_fine(0.0),
  mc_rot_freq(0.5),
  mc_rot_fine_freq(0.5),
  mc_niter_out(0),
  mc_niter_in(0),
  sa_alpha(1.0),
  mc_rand_seed(-1) { } // default constructor

Parameter::~Parameter() { } // default destructor
