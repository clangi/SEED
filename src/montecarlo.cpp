#include "montecarlo.h"
#include "rnd_namespace.h"
#include "funct.h"

void rot_move(double** RoSFCo, int FrAtNu, double max_step) {
  double mc_rot_step;
  double ct[4]; // centroid -> centre of rotation
  mc_rot_step = rnd_gen::get_uniform(0.0, max_step);
  Quaternion<double> q = Quaternion<double>::unitRandom(mc_rot_step);

  CentroidCalc(RoSFCo, FrAtNu, ct);

  for (int i = 1; i <= FrAtNu; i++){
    q.quatConjugateVecRef(RoSFCo[i], ct);
  }
}

void trans_move(double **coord, int num_at, double max_step) {
  double t[4];
  for (int i = 1; i <= 3; i++) {
    t[i] = rnd_gen::get_uniform(-max_step, max_step);
  }
  for (int i = 1; i <= num_at; i++){
    coord[i][1] = coord[i][1] + t[1];
    coord[i][2] = coord[i][2] + t[2];
    coord[i][3] = coord[i][3] + t[3];
  }
  //std::cout << "trans_move: " << t[1] << " " << t[2] << " " << t[3] << std::endl;
}
