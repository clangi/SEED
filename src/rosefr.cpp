/*
*    This file is part of SEED.
*
*    Copyright (C) 2017, Caflisch Lab, University of Zurich
*
*    SEED is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    SEED is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <stdio.h>
#include "funct.h"
#include <quaternion.h> //quaternion class definition

//This may be improved as it seems to be a bottleneck in the computation. clangini
void RoSeFr(int ReCuVe,int *RexxAt,double **ReCoor,int FrCuVe,int *FrxxAt,
            int FrAtNu,double **SeFrCo,double AnglRo,double **RoSFCo)
/* This function rotates the seeded fragment around the axis built with
   the concerned atom of the seeded fragment and the concerned atom of the
   receptor with an angle AnglRo given in radian (in the polar case, the axis
   is acceptor (or donor) of the seeded fragment - donor (or acceptor) of the
   receptor) :
   ReCuVe  receptor current vector
   FrCuVe  fragment current vector
   RotaVe  rotation vector */
{
  int i;
  double RotaVe[4];

  RotaVe[1]=ReCoor[RexxAt[ReCuVe]][1]-SeFrCo[FrxxAt[FrCuVe]][1];
  RotaVe[2]=ReCoor[RexxAt[ReCuVe]][2]-SeFrCo[FrxxAt[FrCuVe]][2];
  RotaVe[3]=ReCoor[RexxAt[ReCuVe]][3]-SeFrCo[FrxxAt[FrCuVe]][3];

  for (i=1;i<=FrAtNu;i++) {
    RoArVe(SeFrCo[i][1]-SeFrCo[FrxxAt[FrCuVe]][1],
           SeFrCo[i][2]-SeFrCo[FrxxAt[FrCuVe]][2],
           SeFrCo[i][3]-SeFrCo[FrxxAt[FrCuVe]][3],
           RotaVe[1],RotaVe[2],RotaVe[3],AnglRo,&RoSFCo[i][1],&RoSFCo[i][2],
           &RoSFCo[i][3]);
    RoSFCo[i][1]=RoSFCo[i][1]+SeFrCo[FrxxAt[FrCuVe]][1];
    RoSFCo[i][2]=RoSFCo[i][2]+SeFrCo[FrxxAt[FrCuVe]][2];
    RoSFCo[i][3]=RoSFCo[i][3]+SeFrCo[FrxxAt[FrCuVe]][3];
  }

}

//clangini START
// void RoSeFrQ(int ReCuVe,int *RexxAt,double **ReCoor,int FrCuVe,int *FrxxAt,
//             int FrAtNu,double **SeFrCo,double AnglRo,double **RoSFCo)
// /* This function is equivalent to RoSeFr but uses a quaternion to perform
//    rotations */
// {
//   int i;
//   double RotaVe[4];
//   Quaternion<double> q; //This calls default constructor
//
//   RotaVe[1]=ReCoor[RexxAt[ReCuVe]][1]-SeFrCo[FrxxAt[FrCuVe]][1];
//   RotaVe[2]=ReCoor[RexxAt[ReCuVe]][2]-SeFrCo[FrxxAt[FrCuVe]][2];
//   RotaVe[3]=ReCoor[RexxAt[ReCuVe]][3]-SeFrCo[FrxxAt[FrCuVe]][3];
//   NormVe(&RotaVe[1],&RotaVe[2],&RotaVe[3]); //NOTE the axis should be normalized!
//   q.fromAngleAxis(AnglRo,RotaVe);
//   for (i=1;i<=FrAtNu;i++) {
//     //Need to add back the
//     RoSFCo[i][1]=RoSFCo[i][1]+SeFrCo[FrxxAt[FrCuVe]][1];
//     RoSFCo[i][2]=RoSFCo[i][2]+SeFrCo[FrxxAt[FrCuVe]][2];
//     RoSFCo[i][3]=RoSFCo[i][3]+SeFrCo[FrxxAt[FrCuVe]][3];
//   }
//
// }
