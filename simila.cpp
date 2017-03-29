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
