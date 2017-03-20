#include <stdio.h>
#include <string.h>
#include "funct.h"
#include "nrutil.h"

#ifndef _STRLENGTH
#define _STRLENGTH 500
#endif

void ReReFi_mol2(char *RecFil,int *ReAtNu,int *ReBdNu,int *ReReNu,
                 char ***ReAtEl,float ***ReCoor,char ***ReAtTy,int **ReResN,
                 float **RePaCh,int ***ReBdAr)
/* This function reads the data of the receptor file (RecFil) in the mol2
   format :
   ReAtNu  number of receptor atoms
   ReBdNu  number of receptor bonds
   ReReNu  number of receptor residues
   ReAtEl  receptor atoms elements
   ReCoor  receptor atoms coordinates
   ReAtTy  receptor atoms types
   ReResN  receptor residues numbers
   RePaCh  receptor partial charges
   ReBdAr  receptor bonds array */
{
  FILE *FilePa;
  char StrLin[_STRLENGTH],**ReAtEl_L,**ReAtTy_L,dummy2[10];
  float **ReCoor_L,*RePaCh_L;
  int i,dummy1,*ReResN_L;

  FilePa=fopen(RecFil,"r");

/* Read ReAtNu ReBdNu ReReNu */
  fgets_wrapper(StrLin,_STRLENGTH,FilePa);
  for (i=0;strncmp(StrLin,"@<TRIPOS>MOLECULE",17);i++)
    fgets_wrapper(StrLin,_STRLENGTH,FilePa);
  fgets_wrapper(StrLin,_STRLENGTH,FilePa);
  fgets_wrapper(StrLin,_STRLENGTH,FilePa);
  sscanf(StrLin,"%d%d%d",ReAtNu,ReBdNu,ReReNu);

/* Read ReAtEl ReCoor ReAtTy ReResN RePaCh */
  ReAtEl_L=cmatrix(1,*ReAtNu,1,5);
  ReCoor_L=matrix(1,*ReAtNu,1,3);
  ReAtTy_L=cmatrix(1,*ReAtNu,1,7);
  ReResN_L=ivector(1,*ReAtNu);
  RePaCh_L=vector(1,*ReAtNu);
  *ReAtEl=ReAtEl_L;
  *ReCoor=ReCoor_L;
  *ReAtTy=ReAtTy_L;
  *ReResN=ReResN_L;
  *RePaCh=RePaCh_L;
  fgets_wrapper(StrLin,_STRLENGTH,FilePa);
  for (i=0;strncmp(StrLin,"@<TRIPOS>ATOM",13);i++)
    fgets_wrapper(StrLin,_STRLENGTH,FilePa);
  for (i=1;i<=*ReAtNu;i++) {
    fgets_wrapper(StrLin,_STRLENGTH,FilePa);
    sscanf(StrLin,"%d%s%f%f%f%s%d%s%f",&dummy1,&ReAtEl_L[i][1],
           &ReCoor_L[i][1],&ReCoor_L[i][2],&ReCoor_L[i][3],&ReAtTy_L[i][1],
           &ReResN_L[i],dummy2,&RePaCh_L[i]);
  }

/* Read ReBdAr */
  *ReBdAr=imatrix(1,*ReBdNu,1,2);
  fgets_wrapper(StrLin,_STRLENGTH,FilePa);
  for (i=0;strncmp(StrLin,"@<TRIPOS>BOND",13);i++)
    fgets_wrapper(StrLin,_STRLENGTH,FilePa);
  for (i=1;i<=*ReBdNu;i++) {
    fgets_wrapper(StrLin,_STRLENGTH,FilePa);
    sscanf(StrLin,"%d%d%d%s",&dummy1,&((*ReBdAr)[i][1]),&((*ReBdAr)[i][2]),
           dummy2);
  }

  fclose(FilePa);
}
