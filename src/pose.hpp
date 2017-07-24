#ifndef _POSE_H
#define _POSE_H

namespace Seed {
  class Pose {
  private:
    int ReAt_idx, FrAt_idx;
    double RotAxis[4];
    double RotAngle;
    double **Coord;
  public:
    Pose();
    Pose(int ReAt_idx_, int FrAt_idx_, double RotAxis_[4], double RotAngle,
          double **Coord_);
    ~Pose();
    double *getAxis();
    double getAngle();
    double getCoord(int at, int xyz);
    int getFrAtIdx();

    void print();
  };
}

#endif
