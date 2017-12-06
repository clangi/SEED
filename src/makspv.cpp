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
#include <math.h>
#include "funct.h"

void MakSpV(int DAType,int DoAcAt,int HyOrZe,double SphAng,int SphPoN,
            int *DANumb,int *ReDATy_L,int *ReDAAt_L,int *ReHydN_L,
            double **ReVeCo_L)
/* This function makes the acceptor or donor vectors on the sphere part,
   the sphere top vector excepted (already made in the calling function) :
  DAType  donor or acceptor type (0 for donor,1 for acceptor)
  DoAcAt  current donor or acceptor atom on which one puts the vectors
  HyOrZe  hydrogen number in the vector or 0 if no hydrogen involved
  NumAnT  total number of theta angles
  NumAnP  number of phi angles for each theta angle
  IntPar  integer part
  ReaPar  real part
  TheSte  theta step
  PhiSte  phi step
  TopVec  top vector number (top of the sphere part)
  TranVe  translation vector
  RotAng  rotation angle
  RotaAx  rotation axis */
{
  int i,j,NumAnT,NumAnP,TopVec;
  double IntPar,ReaPar,PiT180,TheSte,PhiSte,TranVe[4],RotAng,RotaAx[4];
  double DIntPar;

  //PiT180=3.1415927/180; clangini
  PiT180=M_PI/180; //clangini
  TopVec=*DANumb;

  if (DAType==0) {
/* Translation given by the translation of the z axis origin towards the top
   vector origin */
    TranVe[1]=ReVeCo_L[TopVec][1];
    TranVe[2]=ReVeCo_L[TopVec][2];
    TranVe[3]=ReVeCo_L[TopVec][3];
/* Rotation given by the rotation of the translated z axis extremity towards
   the top vector extremity. The rotation vector is perpendicular to the plane
   formed by the translated z axis and the top vector. */
    RotAng=PlaAng(ReVeCo_L[TopVec][1],ReVeCo_L[TopVec][2],ReVeCo_L[TopVec][3],
                  ReVeCo_L[TopVec][4],ReVeCo_L[TopVec][5],ReVeCo_L[TopVec][6],
                  0+TranVe[1],0+TranVe[2],1+TranVe[3]);
    VectPr(0,0,1,ReVeCo_L[TopVec][4]-ReVeCo_L[TopVec][1],ReVeCo_L[TopVec][5]-
           ReVeCo_L[TopVec][2],ReVeCo_L[TopVec][6]-ReVeCo_L[TopVec][3],
           &RotaAx[1],&RotaAx[2],&RotaAx[3]);
  }
  else {
/* Translation given by the translation of the z axis origin towards the top
   vector extremity */
    TranVe[1]=ReVeCo_L[TopVec][4];
    TranVe[2]=ReVeCo_L[TopVec][5];
    TranVe[3]=ReVeCo_L[TopVec][6];
/* Rotation given by the rotation of the translated z axis extremity towards
   the top vector origin. The rotation vector is perpendicular to the plane
   formed by the translated z axis and the top vector. */
    RotAng=PlaAng(ReVeCo_L[TopVec][4],ReVeCo_L[TopVec][5],ReVeCo_L[TopVec][6],
                  ReVeCo_L[TopVec][1],ReVeCo_L[TopVec][2],ReVeCo_L[TopVec][3],
                  0+TranVe[1],0+TranVe[2],1+TranVe[3]);
    VectPr(0,0,1,ReVeCo_L[TopVec][1]-ReVeCo_L[TopVec][4],ReVeCo_L[TopVec][2]-
           ReVeCo_L[TopVec][5],ReVeCo_L[TopVec][3]-ReVeCo_L[TopVec][6],
           &RotaAx[1],&RotaAx[2],&RotaAx[3]);
  }

/* Compute the number of theta angles */
  NumAnT=ffloor(sqrtf(SphPoN*SphAng*PiT180/(2*(1-cosf(SphAng*PiT180)))));
  TheSte=SphAng*PiT180/NumAnT;

  for (i=1;i<=NumAnT;i++) {

/* Compute the number of phi angles for each theta angle */
    ReaPar=modff(2*NumAnT*sinf(TheSte*i),&DIntPar);
    IntPar=DIntPar;
    if (ReaPar<=0.5)
      NumAnP=IntPar;
    else
      NumAnP=IntPar+1;

    if (NumAnP) {
      PhiSte=2*3.1415927/NumAnP;

      for (j=1;j<=NumAnP;j++) {
        *DANumb=*DANumb+1;
        ReDATy_L[*DANumb]=DAType;
        ReDAAt_L[*DANumb]=DoAcAt;
        ReHydN_L[*DANumb]=HyOrZe;

        if (DAType==0) {
          ReVeCo_L[*DANumb][1]=ReVeCo_L[TopVec][1];
          ReVeCo_L[*DANumb][2]=ReVeCo_L[TopVec][2];
          ReVeCo_L[*DANumb][3]=ReVeCo_L[TopVec][3];

          ReVeCo_L[*DANumb][4]=1.0*sinf(i*TheSte)*cosf(j*PhiSte);
          ReVeCo_L[*DANumb][5]=1.0*sinf(i*TheSte)*sinf(j*PhiSte);
          ReVeCo_L[*DANumb][6]=1.0*cosf(i*TheSte);

          RoArVe(ReVeCo_L[*DANumb][4],ReVeCo_L[*DANumb][5],ReVeCo_L[*DANumb][6],
                 RotaAx[1],RotaAx[2],RotaAx[3],RotAng,&ReVeCo_L[*DANumb][4],
                 &ReVeCo_L[*DANumb][5],&ReVeCo_L[*DANumb][6]);

          ReVeCo_L[*DANumb][4]=ReVeCo_L[*DANumb][4]+TranVe[1];
          ReVeCo_L[*DANumb][5]=ReVeCo_L[*DANumb][5]+TranVe[2];
          ReVeCo_L[*DANumb][6]=ReVeCo_L[*DANumb][6]+TranVe[3];
        }
        else {
          ReVeCo_L[*DANumb][4]=ReVeCo_L[TopVec][4];
          ReVeCo_L[*DANumb][5]=ReVeCo_L[TopVec][5];
          ReVeCo_L[*DANumb][6]=ReVeCo_L[TopVec][6];

          ReVeCo_L[*DANumb][1]=1.0*sinf(i*TheSte)*cosf(j*PhiSte);
          ReVeCo_L[*DANumb][2]=1.0*sinf(i*TheSte)*sinf(j*PhiSte);
          ReVeCo_L[*DANumb][3]=1.0*cosf(i*TheSte);

          RoArVe(ReVeCo_L[*DANumb][1],ReVeCo_L[*DANumb][2],ReVeCo_L[*DANumb][3],
                 RotaAx[1],RotaAx[2],RotaAx[3],RotAng,&ReVeCo_L[*DANumb][1],
                 &ReVeCo_L[*DANumb][2],&ReVeCo_L[*DANumb][3]);

          ReVeCo_L[*DANumb][1]=ReVeCo_L[*DANumb][1]+TranVe[1];
          ReVeCo_L[*DANumb][2]=ReVeCo_L[*DANumb][2]+TranVe[2];
          ReVeCo_L[*DANumb][3]=ReVeCo_L[*DANumb][3]+TranVe[3];
        }

      }
    }

  }

}
