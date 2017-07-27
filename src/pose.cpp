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
  Pose::Pose(int Fr_nu_, int ReAt_idx_, int FrAt_idx_, double RotAxis_[4], double RotAngle_,
        double **Coord_): Fr_nu(Fr_nu_), ReAt_idx(ReAt_idx_), FrAt_idx(FrAt_idx_),
                             RotAngle(RotAngle_), Coord(Coord_) {
    for (int i = 1; i <= 3; i++)
      RotAxis[i] = RotAxis_[i]; // This may be probably computed directly hacing ReAt_idx and FrAt_idx !!!
    NormVe(&RotAxis[1], &RotAxis[2], &RotAxis[3]);
  }
  // destructor
  Pose::~Pose(void){ }

  // Properties
  double *Pose::getAxis(){
    return RotAxis;
  }
  double Pose::getAngle() const {
    return RotAngle;
  }
  double Pose::getCoord(int at, int xyz) const {
    return Coord[at][xyz];
  }
  double *Pose::getAtom(int idx) const {
    return Coord[idx]; // need to be checked!
  }
  double *Pose::getAtom(std::string idx) const {
    if (idx == "FragDA"){
      return Coord[FrAt_idx];
    }
    else {
      std::cerr << "Error in pose atom number request" << std::endl;
      exit(12);
    }
  }
  int Pose::getFrAtIdx() const {
    return FrAt_idx;
  }
  int Pose::getReAtIdx() const {
    return ReAt_idx;
  }
  int Pose::getFr_nu() const {
    return Fr_nu;
  }

  void Pose::print(){
    std::cout << "ReAt_idx = " << ReAt_idx << std::endl;
    std::cout << "FrAt_idx = " << FrAt_idx << std::endl;
    std::cout << "RotAngle = " << RotAngle << std::endl;
    std::cout << "First atom: " << Coord[1][1] << " " << Coord[1][2] << " " <<
         Coord[1][3] << std::endl; // !add check for nullptr
  }
}
