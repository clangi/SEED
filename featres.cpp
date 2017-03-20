#include <stdio.h>
#include "nrutil.h"
#include "funct.h"

void FeatRes(int ReAtNu,float **ReCoor,int ReReNu,int *ReResN,float *RePaCh,
             int *FiAtRes,int *LaAtRes,int *AtReprRes,float *TotChaRes,
             FILE *FPaOut)
/* This function determines the first and last atoms and computes the total 
   charge for each residue. It also finds the residue-representative atom 
   which is the closest to the geometrical (for uncharged residue) or the 
   charge (for charged residue) center :
   FiAtRes  first atom of the residue
   LaAtRes  last atom of the residue
   TotChaRes  total charge of the residue
   AtReprRes  residue-representative atom (the closest to the geometrical 
              (for uncharged residue) or the charge (for charged residue) 
              center)
   GeChCent  geometrical or charge center
   DiCentAto  squared distances between the geometrical or charge center and 
              the atoms of the current residue
   IndexArray  array of indexes
   AtoResNum  number of atoms of the current residue */
{
  int i,j,k,UsVal1,*IndexArray,AtoResNum;
  float GeChCent[4],*DiCentAto;

/* Determine the first and last atoms for each residue */
  UsVal1=ReResN[1];
  j=1;
  FiAtRes[j]=1;
  for (i=1;i<=ReAtNu;i++) {
    if (ReResN[i]!=UsVal1) {
      LaAtRes[j]=i-1;
      j=j+1;
      FiAtRes[j]=i;
      UsVal1=ReResN[i];
    }
  }
  LaAtRes[j]=ReAtNu;

/* Compute the total charge of each residue */
  for (i=1;i<=ReReNu;i++) {
    TotChaRes[i]=0.0;
    for (j=FiAtRes[i];j<=LaAtRes[i];j++)
      TotChaRes[i]=TotChaRes[i]+RePaCh[j];
  }

/* Find the residue atom which is the closest to the geometrical (for uncharged 
   residue) or the charge (for charged residue) center */
  for (i=1;i<=ReReNu;i++) {
    
    AtoResNum=LaAtRes[i]-FiAtRes[i]+1;
    GeChCent[1]=0.0;
    GeChCent[2]=0.0;
    GeChCent[3]=0.0;

    if ((TotChaRes[i]>(-0.1))&&(TotChaRes[i]<0.1)) {
      for (j=FiAtRes[i];j<=LaAtRes[i];j++) {
        GeChCent[1]=GeChCent[1]+ReCoor[j][1];
        GeChCent[2]=GeChCent[2]+ReCoor[j][2];
        GeChCent[3]=GeChCent[3]+ReCoor[j][3];
      }
      GeChCent[1]=GeChCent[1]/((float) AtoResNum);
      GeChCent[2]=GeChCent[2]/((float) AtoResNum);
      GeChCent[3]=GeChCent[3]/((float) AtoResNum);
    }
    else {
      for (j=FiAtRes[i];j<=LaAtRes[i];j++) {
        GeChCent[1]=GeChCent[1]+RePaCh[j]*ReCoor[j][1];
        GeChCent[2]=GeChCent[2]+RePaCh[j]*ReCoor[j][2];
        GeChCent[3]=GeChCent[3]+RePaCh[j]*ReCoor[j][3];
      }
      GeChCent[1]=GeChCent[1]/TotChaRes[i];
      GeChCent[2]=GeChCent[2]/TotChaRes[i];
      GeChCent[3]=GeChCent[3]/TotChaRes[i];
    }

    DiCentAto=vector(1,AtoResNum);
    IndexArray=ivector(1,AtoResNum);
    k=0;
    for (j=FiAtRes[i];j<=LaAtRes[i];j++) {
      k=k+1;
      IndexArray[k]=k;
      DiCentAto[k]=DistSq(GeChCent[1],GeChCent[2],GeChCent[3],
                          ReCoor[j][1],ReCoor[j][2],ReCoor[j][3]);
    }

    Sort(AtoResNum,IndexArray,DiCentAto);
    AtReprRes[i]=FiAtRes[i]+IndexArray[1]-1;

    free_vector(DiCentAto,1,AtoResNum);
    free_ivector(IndexArray,1,AtoResNum);

  }    

/* Write in the output file */
  fprintf(FPaOut,"Residue  First  Last  Tot_charge  Residue-representative\n");
  for (i=1;i<=ReReNu;i++)
    fprintf(FPaOut,"%5d  %6d  %6d  %10.4f  %6d\n",i,FiAtRes[i],LaAtRes[i],
            TotChaRes[i],AtReprRes[i]);
  fprintf(FPaOut,"\n");  

}
