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
#include "nrutil.h"
#include "funct.h"

void ReLAIC(int ReAtNu,double **ReCoor,double LaVdWR,double *ReMaxC,double *ReMinC,
            int *CubNum,int ****CubFAI,int ****CubLAI,int **CubLiA)
/* This function constructs the list of receptor atoms which are in the cubes
   of a grid for the checking of bumps in CoNuBu :
   CubNum  the number of cubes along each direction (1->x,2->y,3->z)
   CubFAI  index of the first atoms in each cube
   CubLAI  index of the last atoms in each cube
   CubLiA  list of atoms in each cube
   LiAInd  index of the atom list
   CubAtX  cube number along x in which each receptor atom is
   CubAtY  cube number along y in which each receptor atom is
   CubAtZ  cube number along z in which each receptor atom is
   NuAtoC  number of atoms in the current cube */
{
  int i,j,k,l,LiAInd,*CubAtX,*CubAtY,*CubAtZ,*CubLiA_L,NuAtoC,
      ***CubFAI_L,***CubLAI_L;

/* Find the number of cubes along x, y and z axis */
  for (i=1;i<=3;i++)
    CubNum[i]=ffloor((ReMaxC[i]-ReMinC[i])/(2*LaVdWR))+1;

  CubFAI_L=i3tensor(1,CubNum[1],1,CubNum[2],1,CubNum[3]);
  CubLAI_L=i3tensor(1,CubNum[1],1,CubNum[2],1,CubNum[3]);
  CubLiA_L=ivector(1,ReAtNu+10);
  *CubFAI=CubFAI_L;
  *CubLAI=CubLAI_L;
  *CubLiA=CubLiA_L;

/* Find to which cube belongs each atom of the receptor */
  CubAtX=ivector(1,ReAtNu);
  CubAtY=ivector(1,ReAtNu);
  CubAtZ=ivector(1,ReAtNu);
  for (i=1;i<=ReAtNu;i++) {
    CubAtX[i]=ffloor((ReCoor[i][1]-ReMinC[1])/(2*LaVdWR))+1;
    CubAtY[i]=ffloor((ReCoor[i][2]-ReMinC[2])/(2*LaVdWR))+1;
    CubAtZ[i]=ffloor((ReCoor[i][3]-ReMinC[3])/(2*LaVdWR))+1;
  }

/* Find the list of atoms which are in each cube and the indexes for the first
   and last atoms */
  LiAInd=1;

  for (i=1;i<=CubNum[1];i++) {
    for (j=1;j<=CubNum[2];j++) {
      for (k=1;k<=CubNum[3];k++) {

        CubFAI_L[i][j][k]=LiAInd;
        NuAtoC=0;

        for (l=1;l<=ReAtNu;l++) {

          if ((CubAtX[l]==i)&&(CubAtY[l]==j)&&(CubAtZ[l]==k)) {
            NuAtoC=NuAtoC+1;
            CubLiA_L[LiAInd]=l;
            LiAInd=LiAInd+1;
          }

        }

        if (NuAtoC==0) {
          CubFAI_L[i][j][k]=0;
          CubLAI_L[i][j][k]=0;
        }
        else
          CubLAI_L[i][j][k]=LiAInd-1;

      }
    }
  }

  free_ivector(CubAtX,1,ReAtNu);
  free_ivector(CubAtY,1,ReAtNu);
  free_ivector(CubAtZ,1,ReAtNu);

}
