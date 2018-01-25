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

void Simila(int FrNu_1,int FrNu_2,int FrAtNu,int *FrAtEl_nu,double ***FrCoPo,
            double **SimWei,double SimExp,double *FraSim)
/* This function computes the not normalized similarity number between two
   fragments with the GSEAL method :
   SimWei  similarity weight factors for GSEAL
   SimExp  similarity exponential factor for GSEAL
   FraSim  not normalized similarity number between the two fragments
   FrNu_1  number of the first fragment
   FrNu_2  number of the second fragment */
{
  double SimTot,Dista_sq;
  int i1,i2;

/* One fragment with the other fragment */
  SimTot=0.0;
  for (i1=1;i1<=FrAtNu;i1++) {
    for (i2=1;i2<=FrAtNu;i2++) {

      if(SimWei[FrAtEl_nu[i1]][FrAtEl_nu[i2]] != 0)
	{
	  Dista_sq=DistSq(FrCoPo[FrNu_1][i1][1],FrCoPo[FrNu_1][i1][2],
			  FrCoPo[FrNu_1][i1][3],FrCoPo[FrNu_2][i2][1],
			  FrCoPo[FrNu_2][i2][2],FrCoPo[FrNu_2][i2][3]);
	  SimTot=SimTot+SimWei[FrAtEl_nu[i1]][FrAtEl_nu[i2]]*
	    expf(-SimExp*Dista_sq);
	}
    }
  }

  *FraSim=SimTot;

}
