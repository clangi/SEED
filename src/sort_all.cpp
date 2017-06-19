#include <stdio.h>
#include <string.h>
#include "funct.h"
#include "nrutil.h"



void Sort_all(int FragNu,int *SFWrNu_ar,char **ResN_fr,char *Solv_typ)
// REMOVED -> substituted by Sort(...) clangini
/* This function collects the energies of all the fragment types and sort them
   with respect to the total energy. Then the best energies are written in a
   file :
   SFWrNu_to  total number of kept fragment positions during the program run
   FraTyp  fragment type
   FraNumb  fragment numbering
   EneCol_?  energy column number ? (VW_ps,El_inter,...)
   StSeek  the character string one looks for in the output file
   StLeng  length of this above character string
   ToNuLi  total number of lines read that contain information
   EneCol_co  copy of the total energy column
   IndCol  column containing the index
   NuLiWr  number of lines to be written */
{
  FILE *FilePa;
  char StrLin[150],StSeek[150];
  int SFWrNu_to,i,*FraTyp,*FraNumb,j,StLeng,ToNuLi,*IndCol,NuLiWr;
  double *EneCol_1,*EneCol_2,*EneCol_3,*EneCol_4,*EneCol_5,
        *EneCol_6,*EneCol_7,*EneCol_8,*EneCol_9,*EneCol_10,
        *EneCol_co;

/* Compute the total number of kept fragment positions during the program run */
  SFWrNu_to=0;
  for (i=1;i<=FragNu;i++)
    SFWrNu_to=SFWrNu_to+SFWrNu_ar[i];

/* Allocate memory */
  FraTyp=ivector(1,SFWrNu_to);
  FraNumb=ivector(1,SFWrNu_to);
  IndCol=ivector(1,SFWrNu_to);
  EneCol_1=dvector(1,SFWrNu_to);
  EneCol_2=dvector(1,SFWrNu_to);
  EneCol_3=dvector(1,SFWrNu_to);
  EneCol_4=dvector(1,SFWrNu_to);
  EneCol_5=dvector(1,SFWrNu_to);
  EneCol_6=dvector(1,SFWrNu_to);
  EneCol_7=dvector(1,SFWrNu_to);
  EneCol_8=dvector(1,SFWrNu_to);
  EneCol_9=dvector(1,SFWrNu_to);
  EneCol_10=dvector(1,SFWrNu_to);
  EneCol_co=dvector(1,SFWrNu_to);

/* Read the needed informations of all the fragment types */
  FilePa=fopen("./outputs/do_not_touch_ener","r");

  ToNuLi=0;

  for (i=1;i<=FragNu;i++) {

    if (SFWrNu_ar[i]!=0) {

      sprintf(StSeek,"Frag_type %d",i);
      StLeng=strlen(StSeek);

      fgets_wrapper(StrLin,150,FilePa);
      while (strncmp(StrLin,StSeek,StLeng))
        fgets_wrapper(StrLin,150,FilePa);

      for (j=1;j<=SFWrNu_ar[i];j++) {
        fgets_wrapper(StrLin,150,FilePa);
        ToNuLi=ToNuLi+1;
        sscanf(StrLin,"%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&FraNumb[ToNuLi],
               &EneCol_1[ToNuLi],&EneCol_2[ToNuLi],&EneCol_3[ToNuLi],
               &EneCol_4[ToNuLi],&EneCol_5[ToNuLi],&EneCol_6[ToNuLi],
               &EneCol_7[ToNuLi],&EneCol_8[ToNuLi],&EneCol_9[ToNuLi],
               &EneCol_10[ToNuLi]);
        FraTyp[ToNuLi]=i;
      }

    }

  }

  fclose(FilePa);

/* Sorting the lines with respect to the total energy */
  for (i=1;i<=SFWrNu_to;i++) {

    if ((Solv_typ[0]=='s')||(Solv_typ[0]=='b'))
      EneCol_co[i]=EneCol_10[i];
    else if (Solv_typ[0]=='f')
      EneCol_co[i]=EneCol_9[i];

    IndCol[i]=i;


  }



  Sort(SFWrNu_to,IndCol,EneCol_co);

/* Write the first best energies of all the fragment types */
  FilePa=fopen("./outputs/bestener_all_frag","w");

  fprintf(FilePa,"               First best energies of ");
  fprintf(FilePa,"all the fragment types\n\n");
  fprintf(FilePa,"Total number of fragment positions : %d\n\n",SFWrNu_to);
  fprintf(FilePa,"Num Fr_typ Fr_nam Fr_num   VW_f      VW_s      In_f      In_s");
  fprintf(FilePa,"      Dr_f      Dr_s      Df_f      Df_s      To_f      To_s\n");

  if (SFWrNu_to>20000)
    NuLiWr=20000;
  else
    NuLiWr=SFWrNu_to;

  for (i=1;i<=NuLiWr;i++)
    fprintf(FilePa,"%5d%4d%5s%7d%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n",
            i,FraTyp[IndCol[i]],&ResN_fr[FraTyp[IndCol[i]]][1],
            FraNumb[IndCol[i]],EneCol_1[IndCol[i]],
            EneCol_2[IndCol[i]],EneCol_3[IndCol[i]],EneCol_4[IndCol[i]],
            EneCol_5[IndCol[i]],EneCol_6[IndCol[i]],EneCol_7[IndCol[i]],
            EneCol_8[IndCol[i]],EneCol_9[IndCol[i]],EneCol_10[IndCol[i]]);

  fclose(FilePa);

/* Free the dynamically allocated memory */
  free_ivector(FraTyp,1,SFWrNu_to);
  free_ivector(FraNumb,1,SFWrNu_to);
  free_ivector(IndCol,1,SFWrNu_to);
  free_dvector(EneCol_1,1,SFWrNu_to);
  free_dvector(EneCol_2,1,SFWrNu_to);
  free_dvector(EneCol_3,1,SFWrNu_to);
  free_dvector(EneCol_4,1,SFWrNu_to);
  free_dvector(EneCol_5,1,SFWrNu_to);
  free_dvector(EneCol_6,1,SFWrNu_to);
  free_dvector(EneCol_7,1,SFWrNu_to);
  free_dvector(EneCol_8,1,SFWrNu_to);
  free_dvector(EneCol_9,1,SFWrNu_to);
  free_dvector(EneCol_10,1,SFWrNu_to);
  free_dvector(EneCol_co,1,SFWrNu_to);

}
