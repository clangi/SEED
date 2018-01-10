#ifndef _MONTECARLO_H
#define _MONTECARLO_H

#include "quaternion.h"
#include "Parameter.h"

void rot_move(double **RoSFCo, int FrAtNu, double max_step);
void trans_move(double **coord, int num_at, double max_step);

#endif
