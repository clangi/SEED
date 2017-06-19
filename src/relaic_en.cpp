#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "funct.h"

void ReLAIC_en(int ReReNu,double **ReCoor,double GrSiCu_en,double *ReMinC,
               int *CubNum_en,int ***CubFAI_en,int ***CubLAI_en,int *CubLiA_en,
               int *AtReprRes,FILE *FPaOut)
/* This function constructs the list of receptor residues which are in the
   cubes of a grid for the energy evaluation. The position of each residue is
   considered to be the position of its residue-representative atom :
   CubNum_en  the number of cubes along each direction (1->x,2->y,3->z)
   CubFAI_en  index of the first residues in each cube
   CubLAI_en  index of the last residues in each cube
   CubLiA_en  list of residues in each cube
   LiAInd  index of the residue list
   CubAtX  cube number along x in which each receptor residue is
   CubAtY  cube number along y in which each receptor residue is
   CubAtZ  cube number along z in which each receptor residue is
   NuAtoC  number of residues in the current cube */
{
  int i,j,k,l,LiAInd,*CubAtX,*CubAtY,*CubAtZ,NuAtoC;

/* Find to which cube belongs each residue of the receptor */
  CubAtX=ivector(1,ReReNu);
  CubAtY=ivector(1,ReReNu);
  CubAtZ=ivector(1,ReReNu);
  for (i=1;i<=ReReNu;i++) {
    CubAtX[i]=ffloor((ReCoor[AtReprRes[i]][1]-ReMinC[1])/GrSiCu_en)+1;
    CubAtY[i]=ffloor((ReCoor[AtReprRes[i]][2]-ReMinC[2])/GrSiCu_en)+1;
    CubAtZ[i]=ffloor((ReCoor[AtReprRes[i]][3]-ReMinC[3])/GrSiCu_en)+1;
  }

/* Find the list of residues which are in each cube and the indexes for the
   first and last residues */
  LiAInd=1;

  for (i=1;i<=CubNum_en[1];i++) {
    for (j=1;j<=CubNum_en[2];j++) {
      for (k=1;k<=CubNum_en[3];k++) {

        CubFAI_en[i][j][k]=LiAInd;
        NuAtoC=0;

        for (l=1;l<=ReReNu;l++) {

          if ((CubAtX[l]==i)&&(CubAtY[l]==j)&&(CubAtZ[l]==k)) {
            NuAtoC=NuAtoC+1;
            CubLiA_en[LiAInd]=l;
            LiAInd=LiAInd+1;
          }

        }

        if (NuAtoC==0) {
          CubFAI_en[i][j][k]=0;
          CubLAI_en[i][j][k]=0;
        }
        else
          CubLAI_en[i][j][k]=LiAInd-1;

      }
    }
  }

/* Write a message in case of incoherence */
  if ((LiAInd-1)!=ReReNu) {
    fprintf(FPaOut,"\nWARNING The distribution of the residue-representative ");
    fprintf(FPaOut,"atoms on the grid is \n");
    fprintf(FPaOut,"        problematic\n\n");
  }

  free_ivector(CubAtX,1,ReReNu);
  free_ivector(CubAtY,1,ReReNu);
  free_ivector(CubAtZ,1,ReReNu);

}
