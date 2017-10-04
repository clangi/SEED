#include "Parameter.h"

Parameter::Parameter() : do_mc('n'), mc_temp(0.0), mc_max_tran_step(0.0),
  mc_max_rot_step(0.0), mc_niter(0) { } // default constructor

Parameter::~Parameter() { } // default destructor
