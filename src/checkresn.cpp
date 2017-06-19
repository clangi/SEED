#include <stdio.h>
#include <string.h>
#include "nrutil.h"
#include <math.h> /* dey for gcc 4.0 compatibility */
#include "funct.h"
#include <stdlib.h> /* for exit-fct. */

#ifndef _STRLENGTH
#define _STRLENGTH 500
#endif

void CheckRESN(char **FrFiNa,int FragNu,char **ResN_fr)
/* This function checks the existence of the keyword RESN in the fragment
   files. If it is not found, it is automatically assigned :
   ChkVal        0 if RESN exists, 1 if not
   ResN_fr       name of the fragment residue
   SelCha        characters that are selected to build RESN
   CountUnkResn  counter of unknown RESN */
/* Modification 06.08.2004 NM */
{

  FILE *FilePa;
  char StrLin[_STRLENGTH],dummy1[10],dummy2[10],SelChar[40];
  int i,*ChkVal,CountUnkResn,dig0,dig1,dig2,dig3;

  ChkVal=ivector(1,FragNu);

  for (i=1;i<=FragNu;i++)
  {

    FilePa=fopen(FrFiNa[i]+1,"r");

    ChkVal[i]=1;

    (void) fgets_wrapper(StrLin,_STRLENGTH,FilePa);
    while (strncmp(StrLin,"@<TRIPOS>MOLECULE",17))
    {
      if (!strncmp(StrLin,"# RESN",6))
      {
        sscanf(StrLin,"%s%s%s",dummy1,dummy2,ResN_fr[i]+1);
        ChkVal[i]=0;
      }
      fgets_wrapper(StrLin,_STRLENGTH,FilePa);

      if(FilePa==NULL || feof(FilePa) ) /*this should not happen ever*/
      {
	  printf("WARNING could not find @<TRIPOS>MOLECULE-tag in file %s, exiting !\n",FrFiNa[i]+1);
/*	  fclose(FilePa); */
	  exit (12);
/*	  return 0; */
      }


    }

    fclose(FilePa);

  }


/* automatic assignment of RESN if not found */
  sprintf(SelChar,"ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789\n");
  CountUnkResn=0;
  for (i=1;i<=FragNu;i++)
  {

/* if RESN not found in fragment file */
    if (ChkVal[i])
    {

/* decomposition of CountUnkResn in base 36 (AB...YZ01...89):
   CountUnkResn = dig0*1 + dig1*36 + dig2*1296 + dig3*46656 */
      dig3 = ffloor( ((float)CountUnkResn) / 46656 );
      dig2 = ffloor( ( ((float)CountUnkResn) -
                       ((float)dig3)*46656 ) / 1296 );
      dig1 = ffloor( ( ((float)CountUnkResn) -
                       ((float)dig3)*46656 - ((float)dig2)*1296 )
                     / 36 );
      dig0 = ffloor( ((float)CountUnkResn) -
                     ((float)dig3)*46656 - ((float)dig2)*1296 -
                     ((float)dig1)*36 );

      sprintf(ResN_fr[i]+1,"%c%c%c%c%s",SelChar[dig3],
              SelChar[dig2],SelChar[dig1],SelChar[dig0],"\0");

      CountUnkResn++;

    }

  }


  free_ivector(ChkVal,1,FragNu);

}
