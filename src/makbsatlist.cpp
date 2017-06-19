#include <stdio.h>
#include "nrutil.h"

void MakBSAtList(int ReAtNu,int *ReResN,int BSResN,int *BSReNu,int *BSAtNu,
                 int **BSAtLi)
/* This function makes the list of atoms which are in the binding site :
   BSAtNu  number of atoms in the binding site
   BSAtLi  list of the binding site atoms */
{
  int i,j,BSAtLi_co;

/* Count the number of atoms which are in the binding site */
  *BSAtNu=0;
  for (i=1;i<=BSResN;i++) {
    for (j=1;j<=ReAtNu;j++) {
      if (ReResN[j]==BSReNu[i]) 
        *BSAtNu=*BSAtNu+1;
    }
  }

/* Make the list of the binding site atoms */
  *BSAtLi=ivector(1,*BSAtNu);
  BSAtLi_co=0;
  for (i=1;i<=BSResN;i++) {
    for (j=1;j<=ReAtNu;j++) {
      if (ReResN[j]==BSReNu[i]) {
        BSAtLi_co=BSAtLi_co+1;
        (*BSAtLi)[BSAtLi_co]=j;
      }
    }
  }

}
