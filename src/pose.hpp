#ifndef _POSE_H
#define _POSE_H

#include <string>

namespace Seed {
  class Pose {
  private:
    int Fr_nu, ReAt_idx, FrAt_idx;
    double RotAxis[4];
    double RotAngle;
    double **Coord;
  public:
    Pose();
    Pose(int Fr_nu_, int ReAt_idx_, int FrAt_idx_, double RotAxis_[4], double RotAngle,
          double **Coord_);
    ~Pose();
    double *getAxis();
    double getAngle() const;
    double getCoord(int at, int xyz) const;
    double *getAtom(int idx) const;
    double *getAtom(std::string) const;
    int getFrAtIdx() const;
    int getReAtIdx() const;
    int getFr_nu() const;

    void print();
  };

  struct PoseEnergy
  {
    double VW, In, Dr, Df, Tot;
    void calcTot(){
      Tot = VW + In + Dr + Df;
    };
  };

}

#endif
