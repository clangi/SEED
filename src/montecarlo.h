#ifndef _MONTECARLO_H
#define _MONTECARLO_H

#include "quaternion.h"
#include "Parameter.h"

void rot_move(double **RoSFCo, int FrAtNu, Parameter& seed_par);
void trans_move(double **coord, int num_at, Parameter& par);

#endif
