#include "montecarlo.h"
#include "rnd_namespace.h"
#include "funct.h"

void rot_move(double** RoSFCo, int FrAtNu, Parameter & seed_par) {
  double mc_rot_step;
  double ct[4]; // centroid -> centre of rotation
  mc_rot_step = rnd_gen::get_uniform(0, seed_par.mc_max_rot_step);
  Quaternion<double> q = Quaternion<double>::unitRandom(mc_rot_step);

  CentroidCalc(RoSFCo, FrAtNu, ct);

  for (int i = 1; i < FrAtNu; i++){
    q.quatConjugateVecRef(RoSFCo[i], ct);
  }
}
