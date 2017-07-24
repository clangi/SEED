//#include <stdlib.h>
//#include <stdio.h>
#include <iostream>
#include "pose.hpp"
#include "funct.h"

namespace Seed {
  // Default constructor
  Pose::Pose(): ReAt_idx(0), FrAt_idx(0), RotAngle(0.0) {
    for (int i = 1; i <= 3; i++)
      RotAxis[i] = 0.0;
    Coord = nullptr;
  }
  // Constructor
  Pose::Pose(int ReAt_idx_, int FrAt_idx_, double RotAxis_[4], double RotAngle_,
        double **Coord_): ReAt_idx(ReAt_idx_), FrAt_idx(FrAt_idx_),
                             RotAngle(RotAngle_), Coord(Coord_) {
    for (int i = 1; i <= 3; i++)
      RotAxis[i] = RotAxis_[i];
    NormVe(&RotAxis[1], &RotAxis[2], &RotAxis[3]);
  }
  // destructor
  Pose::~Pose(void){ }

  // Properties
  double *Pose::getAxis(){
    return RotAxis;
  }
  double Pose::getAngle(){
    return RotAngle;
  }
  double Pose::getCoord(int at, int xyz){
    return Coord[at][xyz];
  }
  int Pose::getFrAtIdx(){
    return FrAt_idx;
  }

  void Pose::print(){
    std::cout << "ReAt_idx = " << ReAt_idx << std::endl;
    std::cout << "FrAt_idx = " << FrAt_idx << std::endl;
    std::cout << "RotAngle = " << RotAngle << std::endl;
    std::cout << "First atom: " << Coord[1][1] << " " << Coord[1][2] << " " <<
         Coord[1][3] << std::endl; // !add check for nullptr
  }
}
