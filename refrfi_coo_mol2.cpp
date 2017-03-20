#include <stdio.h>
#include <string.h>
#include "funct.h"

#ifndef _STRLENGTH
#define _STRLENGTH 500
#endif
/* This function is no longer needed clangini*/
void ReFrFi_coo_mol2(int CurFra,int Ind_num_cn,char **FrFiNa,int FrAtNu,
                     float **FrCoor)
/* This function reads the coordinates of the current conformation of the
   current fragment type in the mol2 file :
   Ind_num_cn  index of the current conformation number
   AtIndex_in  index of the initial atom */
{

  FILE *FilePa;
  char StrLin[_STRLENGTH],dummy2[10];
  int i,AtIndex_in,dummy1;

  FilePa=fopen(FrFiNa[CurFra]+1,"r");

/* Compute initial atom index from which the coordinates have to be read */
  AtIndex_in=FrAtNu*(Ind_num_cn-1)+1;

/* Go up to @<TRIPOS>ATOM */
  fgets_wrapper(StrLin,_STRLENGTH,FilePa);
  for (i=0;strncmp(StrLin,"@<TRIPOS>ATOM",13);i++)
    fgets_wrapper(StrLin,_STRLENGTH,FilePa);

/* Skip some lines if necessary (if it is not the first conformation) */
  for (i=1;i<=(AtIndex_in-1);i++)
    fgets_wrapper(StrLin,_STRLENGTH,FilePa);

/* Read the coordinates of the current conformation Ind_num_cn of the current
   fragment type CurFra */
  for (i=1;i<=FrAtNu;i++) {
    fgets_wrapper(StrLin,_STRLENGTH,FilePa);
    sscanf(StrLin,"%d%s%f%f%f",&dummy1,dummy2,&FrCoor[i][1],&FrCoor[i][2],
           &FrCoor[i][3]);
  }

  fclose(FilePa);

}
