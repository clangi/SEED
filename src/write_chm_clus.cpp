#include <stdio.h>

#ifndef _STRLENGTH
#define _STRLENGTH 500
#endif

void write_chm_clus(int CurFra,int SFWrNu,int FrAtNu,int *AliHyd,
                    double *Coulo_ro,double *Vande_ro,double *TotEn_ro,
                    char **FrAtTy,double ***FrCoPo,int *ClusLi_sd,
                    int *ClusLi_sd_01,int NonAliHy,char **ResN_fr,
                    char **FrFiNa_out)
/* This function makes a charmm format file containing the representatives of
   the clusters after the second GSEAL :
   WriPat  the path of the file in which one writes the charmm format
   FrAtTy  fragment atoms types
   FraNam  fragment name
   SegIde  identification string for the charmm format
   ToNuAt  incrementation up to the total number of atoms
   FrCoPo  coordinates to be written
   NumRep  number of cluster representatives
   LiFrNu  limit fragment number as output
   ToNuCR  incrementation up to the total number of cluster representatives to
           be written
   InAtFr  incrementation upon the number of atoms in each fragment */
{
  FILE *FilePa_1;
  char WriPat[_STRLENGTH],*FraNam,SegIde[10];
  int ToNuAt,NumRep,LiFrNu,ToNuCR,InAtFr,i1,i2;


/* Create WriPat and write the first lines */
  sprintf(WriPat,"%s%s%s","./outputs/",&FrFiNa_out[CurFra][1],"_clus.chm\0");
  FilePa_1=fopen(WriPat,"w");
  fprintf(FilePa_1,"* seed crds\n");
  fprintf(FilePa_1,"*\n");

  FraNam=&ResN_fr[CurFra][1];

/* Compute the number of cluster representatives */
  NumRep=0;
  for (i1=1;i1<=SFWrNu;i1++)
    NumRep=NumRep+ClusLi_sd_01[i1];

/* Decide the number of fragment which will be written in the output file
   and write the total number of written atoms */
  if (NumRep>999)
    LiFrNu=999;
  else
    LiFrNu=NumRep;
  fprintf(FilePa_1,"%5d\n",LiFrNu*NonAliHy);

  ToNuAt=0;
  ToNuCR=0;

  for (i1=1;i1<=SFWrNu;i1++) {

    if ((ClusLi_sd_01[i1])&&(ToNuCR<999)) {

      ToNuCR=ToNuCR+1;

/* Write the "first" lines */
      sprintf(SegIde,"%d%s",ToNuCR,"A\0");

      InAtFr=1;

      while (!(AliHyd[InAtFr]))
        InAtFr=InAtFr+1;
      ToNuAt=ToNuAt+1;
      fprintf(FilePa_1,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
              ToNuAt,ToNuCR,FraNam,&FrAtTy[InAtFr][1],
              FrCoPo[ClusLi_sd[i1]][InAtFr][1],FrCoPo[ClusLi_sd[i1]][InAtFr][2],
              FrCoPo[ClusLi_sd[i1]][InAtFr][3],SegIde,ToNuCR,
              TotEn_ro[ClusLi_sd[i1]]);

      InAtFr=InAtFr+1;
      while (!(AliHyd[InAtFr]))
        InAtFr=InAtFr+1;
      ToNuAt=ToNuAt+1;
      fprintf(FilePa_1,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
              ToNuAt,ToNuCR,FraNam,&FrAtTy[InAtFr][1],
              FrCoPo[ClusLi_sd[i1]][InAtFr][1],FrCoPo[ClusLi_sd[i1]][InAtFr][2],
              FrCoPo[ClusLi_sd[i1]][InAtFr][3],SegIde,ToNuCR,
              Vande_ro[ClusLi_sd[i1]]);

      InAtFr=InAtFr+1;
      while (!(AliHyd[InAtFr]))
        InAtFr=InAtFr+1;
      ToNuAt=ToNuAt+1;
      fprintf(FilePa_1,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
              ToNuAt,ToNuCR,FraNam,&FrAtTy[InAtFr][1],
              FrCoPo[ClusLi_sd[i1]][InAtFr][1],FrCoPo[ClusLi_sd[i1]][InAtFr][2],
              FrCoPo[ClusLi_sd[i1]][InAtFr][3],SegIde,ToNuCR,
              Coulo_ro[ClusLi_sd[i1]]);

      if (NonAliHy>=5) {
        InAtFr=InAtFr+1;
        while (!(AliHyd[InAtFr]))
          InAtFr=InAtFr+1;
        ToNuAt=ToNuAt+1;
        fprintf(FilePa_1,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
              ToNuAt,ToNuCR,FraNam,&FrAtTy[InAtFr][1],
              FrCoPo[ClusLi_sd[i1]][InAtFr][1],FrCoPo[ClusLi_sd[i1]][InAtFr][2],
              FrCoPo[ClusLi_sd[i1]][InAtFr][3],SegIde,ToNuCR,0.0);

        InAtFr=InAtFr+1;
        while (!(AliHyd[InAtFr]))
          InAtFr=InAtFr+1;
        ToNuAt=ToNuAt+1;
        fprintf(FilePa_1,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
              ToNuAt,ToNuCR,FraNam,&FrAtTy[InAtFr][1],
              FrCoPo[ClusLi_sd[i1]][InAtFr][1],FrCoPo[ClusLi_sd[i1]][InAtFr][2],
              FrCoPo[ClusLi_sd[i1]][InAtFr][3],SegIde,ToNuCR,
              TotEn_ro[ClusLi_sd[i1]]);
      }

/* Write the following lines */
      for (i2=InAtFr+1;i2<=FrAtNu;i2++) {
        if (AliHyd[i2]) {
          ToNuAt=ToNuAt+1;
          fprintf(FilePa_1,
                  "%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
                  ToNuAt,ToNuCR,FraNam,&FrAtTy[i2][1],
                  FrCoPo[ClusLi_sd[i1]][i2][1],FrCoPo[ClusLi_sd[i1]][i2][2],
                  FrCoPo[ClusLi_sd[i1]][i2][3],SegIde,ToNuCR,0.0);
        }
      }

    }

  }

  fclose(FilePa_1);

}



void write_chm_clus_reduc(int CurFra,int SFWrNu,int FrAtNu,int *AliHyd,
                          double *Coulo_ro,double *Vande_ro,double *TotEn_ro,
                          char **FrAtTy,double ***FrCoPo,int *ClusLi_sd,
                          int *ClusLi_sd_01_reduc,int NonAliHy,char **ResN_fr,
                          char **FrFiNa_out)
/* This function makes a charmm format file containing the representatives of
   the clusters after the second GSEAL but keeping only a reduced number of
   representatives :
   WriPat  the path of the file in which one writes the charmm format
   FrAtTy  fragment atoms types
   FraNam  fragment name
   SegIde  identification string for the charmm format
   ToNuAt  incrementation up to the total number of atoms
   FrCoPo  coordinates to be written
   NumRep  number of cluster representatives
   LiFrNu  limit fragment number as output
   ToNuCR  incrementation up to the total number of cluster representatives to
           be written
   InAtFr  incrementation upon the number of atoms in each fragment */
{
  FILE *FilePa_1;
  char WriPat[_STRLENGTH],*FraNam,SegIde[10];
  int ToNuAt,NumRep,LiFrNu,ToNuCR,InAtFr,i1,i2;


/* Create WriPat and write the first lines */
  sprintf(WriPat,"%s%s%s","./outputs/",&FrFiNa_out[CurFra][1],
          "_clus_reduc.chm\0");
  FilePa_1=fopen(WriPat,"w");
  fprintf(FilePa_1,"* seed crds\n");
  fprintf(FilePa_1,"*\n");

  FraNam=&ResN_fr[CurFra][1];

/* Compute the number of cluster representatives */
  NumRep=0;
  for (i1=1;i1<=SFWrNu;i1++)
    NumRep=NumRep+ClusLi_sd_01_reduc[i1];

/* Decide the number of fragment which will be written in the output file
   and write the total number of written atoms */
  if (NumRep>999)
    LiFrNu=999;
  else
    LiFrNu=NumRep;
  fprintf(FilePa_1,"%5d\n",LiFrNu*NonAliHy);

  ToNuAt=0;
  ToNuCR=0;

  for (i1=1;i1<=SFWrNu;i1++) {

    if ((ClusLi_sd_01_reduc[i1])&&(ToNuCR<999)) {

      ToNuCR=ToNuCR+1;

/* Write the "first" lines */
      sprintf(SegIde,"%d%s",ToNuCR,"A\0");

      InAtFr=1;

      while (!(AliHyd[InAtFr]))
        InAtFr=InAtFr+1;
      ToNuAt=ToNuAt+1;
      fprintf(FilePa_1,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
              ToNuAt,ToNuCR,FraNam,&FrAtTy[InAtFr][1],
              FrCoPo[ClusLi_sd[i1]][InAtFr][1],FrCoPo[ClusLi_sd[i1]][InAtFr][2],
              FrCoPo[ClusLi_sd[i1]][InAtFr][3],SegIde,ToNuCR,
              TotEn_ro[ClusLi_sd[i1]]);

      InAtFr=InAtFr+1;
      while (!(AliHyd[InAtFr]))
        InAtFr=InAtFr+1;
      ToNuAt=ToNuAt+1;
      fprintf(FilePa_1,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
              ToNuAt,ToNuCR,FraNam,&FrAtTy[InAtFr][1],
              FrCoPo[ClusLi_sd[i1]][InAtFr][1],FrCoPo[ClusLi_sd[i1]][InAtFr][2],
              FrCoPo[ClusLi_sd[i1]][InAtFr][3],SegIde,ToNuCR,
              Vande_ro[ClusLi_sd[i1]]);

      InAtFr=InAtFr+1;
      while (!(AliHyd[InAtFr]))
        InAtFr=InAtFr+1;
      ToNuAt=ToNuAt+1;
      fprintf(FilePa_1,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
              ToNuAt,ToNuCR,FraNam,&FrAtTy[InAtFr][1],
              FrCoPo[ClusLi_sd[i1]][InAtFr][1],FrCoPo[ClusLi_sd[i1]][InAtFr][2],
              FrCoPo[ClusLi_sd[i1]][InAtFr][3],SegIde,ToNuCR,
              Coulo_ro[ClusLi_sd[i1]]);

      if (NonAliHy>=5) {
        InAtFr=InAtFr+1;
        while (!(AliHyd[InAtFr]))
          InAtFr=InAtFr+1;
        ToNuAt=ToNuAt+1;
        fprintf(FilePa_1,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
              ToNuAt,ToNuCR,FraNam,&FrAtTy[InAtFr][1],
              FrCoPo[ClusLi_sd[i1]][InAtFr][1],FrCoPo[ClusLi_sd[i1]][InAtFr][2],
              FrCoPo[ClusLi_sd[i1]][InAtFr][3],SegIde,ToNuCR,0.0);

        InAtFr=InAtFr+1;
        while (!(AliHyd[InAtFr]))
          InAtFr=InAtFr+1;
        ToNuAt=ToNuAt+1;
        fprintf(FilePa_1,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
              ToNuAt,ToNuCR,FraNam,&FrAtTy[InAtFr][1],
              FrCoPo[ClusLi_sd[i1]][InAtFr][1],FrCoPo[ClusLi_sd[i1]][InAtFr][2],
              FrCoPo[ClusLi_sd[i1]][InAtFr][3],SegIde,ToNuCR,
              TotEn_ro[ClusLi_sd[i1]]);
      }

/* Write the following lines */
      for (i2=InAtFr+1;i2<=FrAtNu;i2++) {
        if (AliHyd[i2]) {
          ToNuAt=ToNuAt+1;
          fprintf(FilePa_1,
                  "%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
                  ToNuAt,ToNuCR,FraNam,&FrAtTy[i2][1],
                  FrCoPo[ClusLi_sd[i1]][i2][1],FrCoPo[ClusLi_sd[i1]][i2][2],
                  FrCoPo[ClusLi_sd[i1]][i2][3],SegIde,ToNuCR,0.0);
        }
      }

    }

  }

  fclose(FilePa_1);

}



void write_chm_clus_pproc(int CurFra,int NuPosSdCl,int FrAtNu,int *AliHyd,
                          double *In_s_ro,double *VW_s_ro,double *To_s_ro,
                          char **FrAtTy,double ***FrCoPo,int *FrPosAr_sort,
                          int NonAliHy,char **ResN_fr,char **FrFiNa_out)
/* This function makes a charmm format file containing the postprocessed
   positions (of the conserved second clusters) :
   WriPat  the path of the file in which one writes the charmm format
   FrAtTy  fragment atoms types
   FraNam  fragment name
   SegIde  identification string for the charmm format
   ToNuAt  incrementation up to the total number of atoms
   FrCoPo  coordinates to be written
   NumRep  number of cluster representatives
   LiFrNu  limit fragment number as output
   ToNuCR  incrementation up to the total number of cluster representatives to
           be written
   InAtFr  incrementation upon the number of atoms in each fragment */
{
  FILE *FilePa_1;
  char WriPat[_STRLENGTH],*FraNam,SegIde[10];
  int ToNuAt,NumRep,LiFrNu,ToNuCR,InAtFr,i1,i2;


/* Create WriPat and write the first lines */
  sprintf(WriPat,"%s%s%s","./outputs/",&FrFiNa_out[CurFra][1],
          "_clus_pproc.chm\0");
  FilePa_1=fopen(WriPat,"w");
  fprintf(FilePa_1,"* seed crds\n");
  fprintf(FilePa_1,"*\n");

  FraNam=&ResN_fr[CurFra][1];

/* Number of cluster representatives */
  NumRep=NuPosSdCl;

/* Decide the number of fragment which will be written in the output file
   and write the total number of written atoms */
  if (NumRep>999)
    LiFrNu=999;
  else
    LiFrNu=NumRep;
  fprintf(FilePa_1,"%5d\n",LiFrNu*NonAliHy);

  ToNuAt=0;
  ToNuCR=0;

  for (i1=1;i1<=NuPosSdCl;i1++) {

    if (ToNuCR<999) {

      ToNuCR=ToNuCR+1;

/* Write the "first" lines */
      sprintf(SegIde,"%d%s",ToNuCR,"A\0");

      InAtFr=1;

      while (!(AliHyd[InAtFr]))
        InAtFr=InAtFr+1;
      ToNuAt=ToNuAt+1;
      fprintf(FilePa_1,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
              ToNuAt,ToNuCR,FraNam,&FrAtTy[InAtFr][1],
              FrCoPo[FrPosAr_sort[i1]][InAtFr][1],
              FrCoPo[FrPosAr_sort[i1]][InAtFr][2],
              FrCoPo[FrPosAr_sort[i1]][InAtFr][3],SegIde,ToNuCR,
              To_s_ro[FrPosAr_sort[i1]]);

      InAtFr=InAtFr+1;
      while (!(AliHyd[InAtFr]))
        InAtFr=InAtFr+1;
      ToNuAt=ToNuAt+1;
      fprintf(FilePa_1,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
              ToNuAt,ToNuCR,FraNam,&FrAtTy[InAtFr][1],
              FrCoPo[FrPosAr_sort[i1]][InAtFr][1],
              FrCoPo[FrPosAr_sort[i1]][InAtFr][2],
              FrCoPo[FrPosAr_sort[i1]][InAtFr][3],SegIde,ToNuCR,
              VW_s_ro[FrPosAr_sort[i1]]);

      InAtFr=InAtFr+1;
      while (!(AliHyd[InAtFr]))
        InAtFr=InAtFr+1;
      ToNuAt=ToNuAt+1;
      fprintf(FilePa_1,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
              ToNuAt,ToNuCR,FraNam,&FrAtTy[InAtFr][1],
              FrCoPo[FrPosAr_sort[i1]][InAtFr][1],
              FrCoPo[FrPosAr_sort[i1]][InAtFr][2],
              FrCoPo[FrPosAr_sort[i1]][InAtFr][3],SegIde,ToNuCR,
              In_s_ro[FrPosAr_sort[i1]]);

      if (NonAliHy>=5) {
        InAtFr=InAtFr+1;
        while (!(AliHyd[InAtFr]))
          InAtFr=InAtFr+1;
        ToNuAt=ToNuAt+1;
        fprintf(FilePa_1,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
              ToNuAt,ToNuCR,FraNam,&FrAtTy[InAtFr][1],
              FrCoPo[FrPosAr_sort[i1]][InAtFr][1],
              FrCoPo[FrPosAr_sort[i1]][InAtFr][2],
              FrCoPo[FrPosAr_sort[i1]][InAtFr][3],SegIde,ToNuCR,0.0);

        InAtFr=InAtFr+1;
        while (!(AliHyd[InAtFr]))
          InAtFr=InAtFr+1;
        ToNuAt=ToNuAt+1;
        fprintf(FilePa_1,"%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
              ToNuAt,ToNuCR,FraNam,&FrAtTy[InAtFr][1],
              FrCoPo[FrPosAr_sort[i1]][InAtFr][1],
              FrCoPo[FrPosAr_sort[i1]][InAtFr][2],
              FrCoPo[FrPosAr_sort[i1]][InAtFr][3],SegIde,ToNuCR,
              To_s_ro[FrPosAr_sort[i1]]);
      }

/* Write the following lines */
      for (i2=InAtFr+1;i2<=FrAtNu;i2++) {
        if (AliHyd[i2]) {
          ToNuAt=ToNuAt+1;
          fprintf(FilePa_1,
                  "%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
                  ToNuAt,ToNuCR,FraNam,&FrAtTy[i2][1],
                  FrCoPo[FrPosAr_sort[i1]][i2][1],
                  FrCoPo[FrPosAr_sort[i1]][i2][2],
                  FrCoPo[FrPosAr_sort[i1]][i2][3],SegIde,ToNuCR,0.0);
        }
      }

    }

  }

  fclose(FilePa_1);

}



void write_chm_clus_pprocr(int CurFra,int NuSdClKe,int FrAtNu,char **FrAtTy,
                           double ***FrCoPo,char **ResN_fr,char **FrFiNa_out,
                           int *ReprSdClAr)
/* This function makes a charmm format file containing the postprocessed
   representative positions (of the conserved second clusters) :
   WriPat  the path of the file in which one writes the charmm format
   FrAtTy  fragment atoms types
   FraNam  fragment name
   SegIde  identification string for the charmm format
   ToNuAt  incrementation up to the total number of atoms
   FrCoPo  coordinates to be written
   NumRep  number of cluster representatives
   LiFrNu  limit fragment number as output
   ToNuCR  incrementation up to the total number of cluster representatives to
           be written */
{
  FILE *FilePa_1;
  char WriPat[_STRLENGTH],*FraNam,SegIde[10];
  int ToNuAt,NumRep,LiFrNu,ToNuCR,i1,i2;


/* Create WriPat and write the first lines */
  sprintf(WriPat,"%s%s%s","./outputs/",&FrFiNa_out[CurFra][1],
          "_clus_pprocr.chm\0");
  FilePa_1=fopen(WriPat,"w");
  fprintf(FilePa_1,"* seed crds\n");
  fprintf(FilePa_1,"*\n");

  FraNam=&ResN_fr[CurFra][1];

/* Number of cluster representatives */
  NumRep=NuSdClKe;

/* Decide the number of fragment which will be written in the output file
   and write the total number of written atoms */
  if (NumRep>999)
    LiFrNu=999;
  else
    LiFrNu=NumRep;
  fprintf(FilePa_1,"%5d\n",LiFrNu*FrAtNu);

  ToNuAt=0;
  ToNuCR=0;

  for (i1=1;i1<=NuSdClKe;i1++) {

    if (ToNuCR<999) {

      ToNuCR=ToNuCR+1;

      sprintf(SegIde,"%d%s",ToNuCR,"A\0");

/* Write the following lines */
      for (i2=1;i2<=FrAtNu;i2++) {
        ToNuAt=ToNuAt+1;
        fprintf(FilePa_1,
                "%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n",
                ToNuAt,ToNuCR,FraNam,&FrAtTy[i2][1],
                FrCoPo[ReprSdClAr[i1]][i2][1],
                FrCoPo[ReprSdClAr[i1]][i2][2],
                FrCoPo[ReprSdClAr[i1]][i2][3],SegIde,ToNuCR,0.0);
      }

    }

  }

  fclose(FilePa_1);

}
