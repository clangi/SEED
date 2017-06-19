#include <stdio.h>

#ifndef _STRLENGTH
#define _STRLENGTH 500
#endif

void write_charmm(int CurFra,int SFWrNu,int FrAtNu,int *Index_ro,
                  double *Coulo_ro,double *Vande_ro,double *TotEn_ro,
                  char **FrAtTy,double ***FrCoPo,char **ResN_fr,
                  char **FrFiNa_out)
/* This function makes a unique charmm format file :
   WriPat  the path of the file in which one writes the charmm format
   FrAtTy  fragment atoms types
   FraNam  fragment name
   SegIde  identification string for the charmm format
   ToNuAt  incrementation up to the total number of atoms
   FrCoPo  coordinates to be written
   LiFrNu  limit fragment number as output
   FolLin  starting number of the following lines to be written */
{
  FILE *FilePa_1;
  char WriPat[_STRLENGTH],*FraNam,SegIde[10];
  int i,j,ToNuAt,LiFrNu,FolLin;


/* Create WriPat and write the first lines */
  sprintf(WriPat,"%s%s%s","./outputs/",&FrFiNa_out[CurFra][1],".chm\0");
  FilePa_1=fopen(WriPat,"w");
  fprintf(FilePa_1,"* seed crds\n");
  fprintf(FilePa_1,"*\n");

  FraNam=&ResN_fr[CurFra][1];

  ToNuAt=0;

/* Decide the number of fragment which will be written in the output file */
  if (SFWrNu>999)
    LiFrNu=999;
  else
    LiFrNu=SFWrNu;
  fprintf(FilePa_1,"%5d\n",LiFrNu*FrAtNu);

/* Decide the starting number of the following lines to be written */
  if (FrAtNu>=5)
    FolLin=6;
  else
    FolLin=4;


  for (i=1;i<=LiFrNu;i++) {

/* Write the "first" lines */
    sprintf(SegIde,"%d%s",i,"A\0");

    ToNuAt=ToNuAt+1;
    fprintf(FilePa_1,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
            ToNuAt,i,FraNam,&FrAtTy[1][1],FrCoPo[Index_ro[i]][1][1],
            FrCoPo[Index_ro[i]][1][2],FrCoPo[Index_ro[i]][1][3],SegIde,i,
            TotEn_ro[Index_ro[i]]);

    ToNuAt=ToNuAt+1;
    fprintf(FilePa_1,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
            ToNuAt,i,FraNam,&FrAtTy[2][1],FrCoPo[Index_ro[i]][2][1],
            FrCoPo[Index_ro[i]][2][2],FrCoPo[Index_ro[i]][2][3],SegIde,i,
            Vande_ro[Index_ro[i]]);

    ToNuAt=ToNuAt+1;
    fprintf(FilePa_1,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
            ToNuAt,i,FraNam,&FrAtTy[3][1],FrCoPo[Index_ro[i]][3][1],
            FrCoPo[Index_ro[i]][3][2],FrCoPo[Index_ro[i]][3][3],SegIde,i,
            Coulo_ro[Index_ro[i]]);

    if (FrAtNu>=5) {
      ToNuAt=ToNuAt+1;
      fprintf(FilePa_1,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
              ToNuAt,i,FraNam,&FrAtTy[4][1],FrCoPo[Index_ro[i]][4][1],
              FrCoPo[Index_ro[i]][4][2],FrCoPo[Index_ro[i]][4][3],SegIde,i,0.0);

      ToNuAt=ToNuAt+1;
      fprintf(FilePa_1,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
              ToNuAt,i,FraNam,&FrAtTy[5][1],FrCoPo[Index_ro[i]][5][1],
              FrCoPo[Index_ro[i]][5][2],FrCoPo[Index_ro[i]][5][3],SegIde,i,
              TotEn_ro[Index_ro[i]]);
    }

/* Write the following lines */
    for (j=FolLin;j<=FrAtNu;j++) {
      ToNuAt=ToNuAt+1;
      fprintf(FilePa_1,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
              ToNuAt,i,FraNam,&FrAtTy[j][1],FrCoPo[Index_ro[i]][j][1],
              FrCoPo[Index_ro[i]][j][2],FrCoPo[Index_ro[i]][j][3],SegIde,i,0.0);
    }

  }

  fclose(FilePa_1);

}
