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
#include <string.h>
#include <math.h>
#include "nrutil.h"
#include "funct.h"

void FindFrSym(int FrAtNu,double **FrCoor,char **FrAtTy,int *UndisAt_fr)
/* This function finds the undistinguishable fragment atoms (or symmetries)
   by means of rotation :
   UndisAt_fr  undistinguishable fragment atoms (0 if the atom is not taken
               into account because of symmetry, 1 if it is taken into account)
   GeomCen_fr  geometrical center of the fragment
   MidPoint  coordinates of the middle point between the two considered atoms
   VectRot  rotation vector
   FrCoor_rot  coordinates of the rotated fragment
   FrAtTy_rot  atom types of the rotated fragment
   FoundSim  1 if an atom of the rotated fragment was found similar to the
               current original fragment atom, 0 if not
   TotNumSimAt  total number of similar atoms between the original and rotated
                fragments
   PairIdAt  pairs of identical atoms (1->yes,0->no)
   NormVect  norm of the vector
   SquaredDis  squared distance
   AnglArr  array containing the various rotation angle values */
{
  int i,j,k,k1,k2,FoundSim,TotNumSimAt,**PairIdAt;
  double **FrCoor_rot,GeomCen_fr[4],MidPoint[4],VectRot[4],PiT180,SquaredDis,
        NormVect,AnglArr[12];
  char **FrAtTy_rot;

  //PiT180=3.1415927/180; clangini
  PiT180 = M_PI/180; //clangini

/* Allocate memory */
  FrCoor_rot=dmatrix(1,FrAtNu,1,3);
  FrAtTy_rot=cmatrix(1,FrAtNu,1,7);
  PairIdAt=imatrix(1,FrAtNu,1,FrAtNu);

/* Initialization of FrAtTy_rot and UndisAt_fr */
  for (i=1;i<=FrAtNu;i++) {
    strcpy(&FrAtTy_rot[i][1],&FrAtTy[i][1]);
    UndisAt_fr[i]=1;
  }

/* Initialization of AnglArr */
  AnglArr[1]=180.0;
  AnglArr[2]=120.0;
  AnglArr[3]=60.0;
  AnglArr[4]=90.0;
  AnglArr[5]=72.0;
  AnglArr[6]=144.0;
  AnglArr[7]=216.0;
  AnglArr[8]=240.0;
  AnglArr[9]=270.0;
  AnglArr[10]=288.0;
  AnglArr[11]=300.0;

/* Compute the geometrical center of the fragment */
  GeomCen_fr[1]=0.0;
  GeomCen_fr[2]=0.0;
  GeomCen_fr[3]=0.0;
  for (i=1;i<=FrAtNu;i++) {
    GeomCen_fr[1]=GeomCen_fr[1]+FrCoor[i][1];
    GeomCen_fr[2]=GeomCen_fr[2]+FrCoor[i][2];
    GeomCen_fr[3]=GeomCen_fr[3]+FrCoor[i][3];
  }
  GeomCen_fr[1]=GeomCen_fr[1]/FrAtNu;
  GeomCen_fr[2]=GeomCen_fr[2]/FrAtNu;
  GeomCen_fr[3]=GeomCen_fr[3]/FrAtNu;

/* --------------------------------------------------------------------- */
/* First algorithm. Rotation through the middle point between atom pairs */
/* --------------------------------------------------------------------- */

/* Loop over the fragment atom pairs. Each pair (for which the atom types are
   the same) is taken into account only once */
  for (i=1;i<=FrAtNu;i++) {
    for (j=1;j<=FrAtNu;j++) {

      if (!(strcmp(&FrAtTy[i][1],&FrAtTy[j][1]))) {

/* Compute the middle point between the two considered atoms */
        MidPoint[1]=0.5*(FrCoor[i][1]+FrCoor[j][1]);
        MidPoint[2]=0.5*(FrCoor[i][2]+FrCoor[j][2]);
        MidPoint[3]=0.5*(FrCoor[i][3]+FrCoor[j][3]);

/* Compute the rotation vector which joins the middle point with the
   geometrical center */
        VectRot[1]=GeomCen_fr[1]-MidPoint[1];
        VectRot[2]=GeomCen_fr[2]-MidPoint[2];
        VectRot[3]=GeomCen_fr[3]-MidPoint[3];

/* Continue only if the norm of VectRot is larger than 0.01 */
        NormVect=VeNorm(VectRot[1],VectRot[2],VectRot[3]);

        if (NormVect>0.01) {

/* Compute the rotated image of the fragment around VectRot with an angle of
   180 */
          for (k=1;k<=FrAtNu;k++) {
            RoArVe(FrCoor[k][1]-MidPoint[1],FrCoor[k][2]-MidPoint[2],
                   FrCoor[k][3]-MidPoint[3],VectRot[1],VectRot[2],VectRot[3],
                   /*PiT180*180*/M_PI,&FrCoor_rot[k][1],&FrCoor_rot[k][2],
                   &FrCoor_rot[k][3]);
            FrCoor_rot[k][1]=FrCoor_rot[k][1]+MidPoint[1];
            FrCoor_rot[k][2]=FrCoor_rot[k][2]+MidPoint[2];
            FrCoor_rot[k][3]=FrCoor_rot[k][3]+MidPoint[3];
          }

/* Check if the two fragments (the original and the rotated) are the same */
          TotNumSimAt=0;
          for (k1=1;k1<=FrAtNu;k1++) {
            for (k2=1;k2<=FrAtNu;k2++) {
              PairIdAt[k1][k2]=0;
            }
          }

          for (k1=1;k1<=FrAtNu;k1++) {           /* original fragment */

            FoundSim=0;

            for (k2=1;k2<=FrAtNu;k2++) {         /* rotated fragment */

              if (!(FoundSim)) {

                if (!(strcmp(&FrAtTy[k1][1],&FrAtTy_rot[k2][1]))) {

                  SquaredDis=DistSq(FrCoor[k1][1],FrCoor[k1][2],FrCoor[k1][3],
                                    FrCoor_rot[k2][1],FrCoor_rot[k2][2],
                                    FrCoor_rot[k2][3]);

                  if (SquaredDis<=0.001) {
                    FoundSim=1;
                    TotNumSimAt=TotNumSimAt+1;
                    PairIdAt[k1][k2]=1;
                  }

                }

              }

            }

          }

          if (TotNumSimAt==FrAtNu) {
            for (k1=1;k1<=FrAtNu;k1++) {
              for (k2=1;k2<=FrAtNu;k2++) {
                if (PairIdAt[k1][k2]) {
                  if ((UndisAt_fr[k1])&&(UndisAt_fr[k2])&&(k1!=k2)) {
                    UndisAt_fr[k1]=0;
                  }
                }
              }
            }
          }

        }

      }

    }
  }

/* ----------------------------------------------------- */
/* Second algorithm. Rotation through each fragment atom */
/* ----------------------------------------------------- */

/* Loop over the fragment atoms */
  for (i=1;i<=FrAtNu;i++) {

/* Compute the rotation vector which joins the fragment atom with the
   geometrical center */
    VectRot[1]=GeomCen_fr[1]-FrCoor[i][1];
    VectRot[2]=GeomCen_fr[2]-FrCoor[i][2];
    VectRot[3]=GeomCen_fr[3]-FrCoor[i][3];

/* Continue only if the norm of VectRot is larger than 0.01 */
    NormVect=VeNorm(VectRot[1],VectRot[2],VectRot[3]);

    if (NormVect>0.01) {

/* Loop over the various rotation angle values */
      for (j=1;j<=11;j++) {

/* Compute the rotated image of the fragment around VectRot with an angle of
   AnglArr[j] */
        for (k=1;k<=FrAtNu;k++) {
          RoArVe(FrCoor[k][1]-FrCoor[i][1],FrCoor[k][2]-FrCoor[i][2],
                 FrCoor[k][3]-FrCoor[i][3],VectRot[1],VectRot[2],VectRot[3],
                 PiT180*AnglArr[j],&FrCoor_rot[k][1],&FrCoor_rot[k][2],
                 &FrCoor_rot[k][3]);
          FrCoor_rot[k][1]=FrCoor_rot[k][1]+FrCoor[i][1];
          FrCoor_rot[k][2]=FrCoor_rot[k][2]+FrCoor[i][2];
          FrCoor_rot[k][3]=FrCoor_rot[k][3]+FrCoor[i][3];
        }

/* Check if the two fragments (the original and the rotated) are the same */
        TotNumSimAt=0;
        for (k1=1;k1<=FrAtNu;k1++) {
          for (k2=1;k2<=FrAtNu;k2++) {
            PairIdAt[k1][k2]=0;
          }
        }

        for (k1=1;k1<=FrAtNu;k1++) {           /* original fragment */

          FoundSim=0;

          for (k2=1;k2<=FrAtNu;k2++) {         /* rotated fragment */

            if (!(FoundSim)) {

              if (!(strcmp(&FrAtTy[k1][1],&FrAtTy_rot[k2][1]))) {

                SquaredDis=DistSq(FrCoor[k1][1],FrCoor[k1][2],FrCoor[k1][3],
                                  FrCoor_rot[k2][1],FrCoor_rot[k2][2],
                                  FrCoor_rot[k2][3]);

                if (SquaredDis<=0.001) {
                  FoundSim=1;
                  TotNumSimAt=TotNumSimAt+1;
                  PairIdAt[k1][k2]=1;
                }

              }

            }

          }

        }

        if (TotNumSimAt==FrAtNu) {
          for (k1=1;k1<=FrAtNu;k1++) {
            for (k2=1;k2<=FrAtNu;k2++) {
              if (PairIdAt[k1][k2]) {
                if ((UndisAt_fr[k1])&&(UndisAt_fr[k2])&&(k1!=k2)) {
                  UndisAt_fr[k1]=0;
                }
              }
            }
          }
        }

      }

    }

  }

/* Desallocate memory */
  free_dmatrix(FrCoor_rot,1,FrAtNu,1,3);
  free_cmatrix(FrAtTy_rot,1,FrAtNu,1,7);
  free_imatrix(PairIdAt,1,FrAtNu,1,FrAtNu);

}
