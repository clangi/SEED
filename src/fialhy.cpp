#include <stdio.h>

void FiAlHy(int FrAtNu,int FrBdNu,int *FrAtEl_nu,int **FrBdAr,int *AliHyd,
            int *NonAliHy)
/* This function finds the aliphatic hydrogens (i.e. the hydrogens attached 
   to carbon) :
   AliHyd  0 for aliphatic hydrogens and 1 for other hydrogens and 
           other atoms 
   NonAliHy  number of non aliphatic-hydrogen atoms */
{
  int i,j;

  for (i=1;i<=FrAtNu;i++) {
    AliHyd[i]=1;
    if (FrAtEl_nu[i]==1) {
      for (j=1;j<=FrBdNu;j++) {
        if (((FrBdAr[j][1]==i) && (FrAtEl_nu[FrBdAr[j][2]]==6)) ||
            ((FrBdAr[j][2]==i) && (FrAtEl_nu[FrBdAr[j][1]]==6)))
          AliHyd[i]=0;
      }
    }  
  }

/* Count the number of non aliphatic-hydrogen atoms */
  *NonAliHy=0;
  for (i=1;i<=FrAtNu;i++) 
    *NonAliHy=*NonAliHy+AliHyd[i];

}
