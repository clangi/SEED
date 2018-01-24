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

void VdWRaM(int NumbAT,double *VdWRad,double *LaVdWR)
/* This function finds the largest van der Waals radius :
   LaVdWR  the largest van der Waals radius */
{
  int i;

  *LaVdWR=VdWRad[1];
  for (i=2;i<=NumbAT;i++) {
    if (VdWRad[i]>*LaVdWR)
      *LaVdWR=VdWRad[i];
  }

}



void ReMMCo(int ReAtNu,double **ReCoor,double *ReMaxC,double *ReMinC)
/* This function finds the maximal and minimal coordinates of the receptor :
   ReMaxC  the maximal coordinates of the receptor (1->x,2->y,3->z)
   ReMinC  the minimal coordinates of the receptor (1->x,2->y,3->z) */
{
  int i,j;

/* Initialization */
  for (i=1;i<=3;i++) {
    ReMaxC[i]=ReCoor[1][i];
    ReMinC[i]=ReCoor[1][i];
  }

/* Find the maximal and minimal coordinates of the receptor */
  for (i=2;i<=ReAtNu;i++) {
    for (j=1;j<=3;j++) {
      if (ReCoor[i][j]>ReMaxC[j])
        ReMaxC[j]=ReCoor[i][j];
      else if (ReCoor[i][j]<ReMinC[j])
        ReMinC[j]=ReCoor[i][j];
    }
  }

}



void BSMMCo(double **ReCoor,int BSAtNu,int *BSAtLi,double *BSMaxC,double *BSMinC)
/* This function finds the maximal and minimal coordinates of the binding
   site of the receptor :
   BSMaxC  the maximal coordinates of the binding site (1->x,2->y,3->z)
   BSMinC  the minimal coordinates of the binding site (1->x,2->y,3->z) */
{
  int i,j;

/* Initialization */
  for (i=1;i<=3;i++) {
    BSMaxC[i]=ReCoor[BSAtLi[1]][i];
    BSMinC[i]=ReCoor[BSAtLi[1]][i];
  }

/* Find the maximal and minimal coordinates of the binding site of the
   receptor */
  for (i=1;i<=BSAtNu;i++) {
    for (j=1;j<=3;j++) {
      if (ReCoor[BSAtLi[i]][j]>BSMaxC[j])
        BSMaxC[j]=ReCoor[BSAtLi[i]][j];
      else if (ReCoor[BSAtLi[i]][j]<BSMinC[j])
        BSMinC[j]=ReCoor[BSAtLi[i]][j];
    }
  }

}
