#include <stdio.h>
#include <math.h>
#include "funct.h"

float CoFrEn(int FrAtNu,float *FrPaCh,float **RoSFCo,float *BSMinC,
             float CoGrIn,float CoGrSi,int *CoGPoN,float ***CoGrRP,
             FILE *FPaOut,
	     int *print)
/* This function evaluates the coulombic interaction energy between the current
   fragment and the receptor whose effect is precalculated in a grid :
   GrCo_i  "integer coordinate" along the x axis of the grid point just "below"
           the x coordinate of the current fragment atom
   GrCo_j  "integer coordinate" along the y axis of the grid point just "below"
           the y coordinate of the current fragment atom
   GrCo_k  "integer coordinate" along the z axis of the grid point just "below"
           the z coordinate of the current fragment atom
   Weig_t  weight factor for the computation of the coulombic energy
   Weig_u  weight factor for the computation of the coulombic energy
   Weig_v  weight factor for the computation of the coulombic energy
   CoEnEv_gr  evaluation of the coulombic interaction energy on a grid */
{
  int i,GrCo_i,GrCo_j,GrCo_k;
  float Weig_t,Weig_u,Weig_v,CoEnEv_gr;

  CoEnEv_gr=0.0;

  for (i=1;i<=FrAtNu;i++) {

    GrCo_i=ffloor((RoSFCo[i][1]-(BSMinC[1]-CoGrIn))/CoGrSi)+1;
    GrCo_j=ffloor((RoSFCo[i][2]-(BSMinC[2]-CoGrIn))/CoGrSi)+1;
    GrCo_k=ffloor((RoSFCo[i][3]-(BSMinC[3]-CoGrIn))/CoGrSi)+1;

    if ((GrCo_i>=1)&&(GrCo_i<CoGPoN[1])&&(GrCo_j>=1)&&(GrCo_j<CoGPoN[2])&&
        (GrCo_k>=1)&&(GrCo_k<CoGPoN[3])) {

      Weig_t=(RoSFCo[i][1]-(BSMinC[1]-CoGrIn+(GrCo_i-1)*CoGrSi))/CoGrSi;
      Weig_u=(RoSFCo[i][2]-(BSMinC[2]-CoGrIn+(GrCo_j-1)*CoGrSi))/CoGrSi;
      Weig_v=(RoSFCo[i][3]-(BSMinC[3]-CoGrIn+(GrCo_k-1)*CoGrSi))/CoGrSi;

      CoEnEv_gr=CoEnEv_gr + FrPaCh[i] * (
             (1-Weig_t)*(1-Weig_u)*(1-Weig_v)*CoGrRP[GrCo_i][GrCo_j][GrCo_k] +
             Weig_t*(1-Weig_u)*(1-Weig_v)*CoGrRP[GrCo_i+1][GrCo_j][GrCo_k] +
             Weig_t*Weig_u*(1-Weig_v)*CoGrRP[GrCo_i+1][GrCo_j+1][GrCo_k] +
             (1-Weig_t)*Weig_u*(1-Weig_v)*CoGrRP[GrCo_i][GrCo_j+1][GrCo_k] +
             (1-Weig_t)*(1-Weig_u)*Weig_v*CoGrRP[GrCo_i][GrCo_j][GrCo_k+1] +
             Weig_t*(1-Weig_u)*Weig_v*CoGrRP[GrCo_i+1][GrCo_j][GrCo_k+1] +
             Weig_t*Weig_u*Weig_v*CoGrRP[GrCo_i+1][GrCo_j+1][GrCo_k+1] +
             (1-Weig_t)*Weig_u*Weig_v*CoGrRP[GrCo_i][GrCo_j+1][GrCo_k+1] ) ;

    }

#ifndef NOWARN
    else if((*print)==0) {
      ++(*print); /* avoid printing error msg more than once */
      fprintf(FPaOut,"WARNING at least one fragment is not completely in the ");
      /* fprintf(FPaOut,"WARNING One fragment is not completely in the "); */
      fprintf(FPaOut,"coulombic energy grid\n\n");
      return 1e6; /* New -> give penalty */
    }
#else
    {
      return 1e6; /* New -> give penalty */
    }
#endif

  }

  return CoEnEv_gr;

}
