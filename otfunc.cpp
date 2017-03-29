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
