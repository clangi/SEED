#include <stdio.h>
#include "funct.h"

void ListChRes(int BSResN,int *BSReNu,int ReReNu,double **ReCoor,int *AtReprRes,
               double *TotChaRes,double PsSpRa,FILE *FPaOut,int *NuChResEn,
               int *LiChResEn)
/* This function finds which charged residues will always be taken into account
   for the evaluation of the energy :
   LiChResEn  list of charged residues kept for the energy evaluation
   NuChResEn  number of charged residues kept for the energy evaluation */
{
  int i,j,UsVal1;
  double Dista1,Cutoff1;

  Cutoff1=(PsSpRa+0.3*PsSpRa)*(PsSpRa+0.3*PsSpRa);

/* Loop over the charged residues */
  *NuChResEn=0;
  for (i=1;i<=ReReNu;i++) {

    if ((TotChaRes[i]<=(-0.1))||(TotChaRes[i]>=0.1)) {

/* Loop over the residues of the binding site. The current charged residue is
   kept if at least one of the distances between that residue and the ones of
   the binding site is smaller than a cutoff. The distance is computed between
   the residue-representative atoms */
      UsVal1=1;
      for (j=1;j<=BSResN;j++) {
        if (UsVal1) {
          Dista1=DistSq(ReCoor[AtReprRes[i]][1],ReCoor[AtReprRes[i]][2],
                        ReCoor[AtReprRes[i]][3],ReCoor[AtReprRes[BSReNu[j]]][1],
                        ReCoor[AtReprRes[BSReNu[j]]][2],
                        ReCoor[AtReprRes[BSReNu[j]]][3]);
          if (Dista1<=Cutoff1) {
            *NuChResEn=*NuChResEn+1;
            LiChResEn[*NuChResEn]=i;
            UsVal1=0;
          }
        }
      }

    }

  }

/* Write in the output file */
  fprintf(FPaOut,"The following %d charged residues are kept ",*NuChResEn);
  fprintf(FPaOut,"for the energy evaluation\n");
  for (i=1;i<=*NuChResEn;i++)
    fprintf(FPaOut,"%5d\n",LiChResEn[i]);
  fprintf(FPaOut,"\n");

}
