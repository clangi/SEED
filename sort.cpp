#include <stdio.h>
#include <math.h>
#include "funct.h"
#include <float.h> /*DEBUGCHECK*/ //clangini

void Sort(int N,int *IndArr,double *SorArr)
/* This function sorts an array of doubles SorArr and gives also the new order
   in the array IndArr. The result goes from the smallest to the largest. */
{
  int m,k,i,j,l,nn,tempin,loopva;
  double aln2i,lognb2,temp;

  aln2i=1.0/0.69314718;
  lognb2= (double) ((int) (logf(N)*aln2i+5.e-7));

  m=N;
  for (nn=1;nn<=lognb2;nn++) {
    m=m/2;
    k=N-m;
    for (j=1;j<=k;j++) {
      i=j;
      loopva=1;
      while (loopva) {
        loopva=0;
        l=i+m;

/* 	//DEBUG CHECK */
	if(isnan(SorArr[i]) || SorArr[i] != SorArr[i] )
	{
	    SorArr[i]=FLT_MAX;
	    printf("NAN sort() 1 : %f\n",SorArr[i]);
	}
	if(SorArr[l]!=SorArr[l] || isnan(SorArr[l]))
	{
	    SorArr[l]=FLT_MAX;
	    printf("NAN sort() 2 : %f\n",SorArr[l]);
	}
/* 	if(isnan(IndArr[i]) || IndArr[i] != IndArr[i] ) */
/* 	{ */
/* 	    printf("NAN sort() 3 : %lf\n",IndArr[i]); */
/* 	    IndArr[i]=FLT_MAX; */
/* 	} */
/* 	if(IndArr[l]!=IndArr[l] || isnan(IndArr[l])) */
/* 	{ */
/* 	    printf("NAN sort() 4 : %lf\n",IndArr[l]); */
/* 	    IndArr[l]=FLT_MAX;	     */
/* 	} */

        if (SorArr[l]<SorArr[i]) {
          temp=SorArr[i];
          SorArr[i]=SorArr[l];
          SorArr[l]=temp;
          tempin=IndArr[i];
          IndArr[i]=IndArr[l];
          IndArr[l]=tempin;
          i=i-m;
          if (i>=1)
            loopva=1;
        }
      }
    }
  }

}
