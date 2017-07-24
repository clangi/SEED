#include <stdio.h>
#include <vector>
#include "nrutil.h"
#include "funct.h"

void FrAcDo(int FrAtNu,int *FrAtEl_nu,double **FrCoor,int *HybFrAt,int FrBdNu,
            int **FrBdAr,int *FrAcNu,int *FrDoNu,int **FrDATy,int **FrDAAt,
            int **FrHydN,double ***FrVeCo,FILE *FPaOut,char *OutFil)
/* This function constructs the list of the donor and acceptor vectors for the
   current fragment :
   FrAcNu  number of acceptor vectors in the fragment
   FrDoNu  number of donor vectors in the fragment
   FrDATy  fragment donor acceptor type (0 for donor,1 for acceptor)
   FrDAAt  fragment donor or acceptor atom numbers
   FrHydN  fragment hydrogen number in the current vector (0 if no hydrogen
           involved)
   FrVeCo  coordinates of the vector whose orientation goes from the donor
           to the acceptor (for the fragment)
   LiToNu  total number of linked atoms
   LiHyNu  number of linked hydrogen atoms
   LiAtom  array of linked atoms
   AcAtNu  acceptor atom number
   DoAtNu  donor atom number
   HelpA1  help atom number 1
   HelpA2  help atom number 2
   HyAtNu  hydrogen atom number
   HyAtN2  hydrogen atom number 2
   HyAtN3  hydrogen atom number 3
   DANumb  donor and acceptor number
   Angl    angle */
{
/*  FILE *FilePa;*/

  int i,j,k,*LiToNu,*LiHyNu,**LiAtom,AcAtNu,DoAtNu,HelpA1,HelpA2,HyAtNu,
      HyAtN2,HyAtN3,DANumb;
  double PiT180,Angl;

  HyAtNu = -1;
  HyAtN2 = -1;
  HyAtN3 = -1;

  //PiT180=3.1415927/180; // clangini
  PiT180=M_PI/180; //clangini

/* Find the total number of linked atoms, the linked atoms and the number of
   linked hydrogen atoms for each atom of the fragment */
  LiToNu=ivector(1,FrAtNu);
  LiHyNu=ivector(1,FrAtNu);
  LiAtom=imatrix(1,FrAtNu,1,10);

  for (i=1;i<=FrAtNu;i++) {
    LiToNu[i]=0;
    LiHyNu[i]=0;
    for (j=1;j<=FrBdNu;j++) {
      if (FrBdAr[j][1]==i) {
        LiToNu[i]=LiToNu[i]+1;
        LiAtom[i][LiToNu[i]]=FrBdAr[j][2];
        if (FrAtEl_nu[FrBdAr[j][2]]==1) LiHyNu[i]=LiHyNu[i]+1;
      }
      else if (FrBdAr[j][2]==i) {
        LiToNu[i]=LiToNu[i]+1;
        LiAtom[i][LiToNu[i]]=FrBdAr[j][1];
        if (FrAtEl_nu[FrBdAr[j][1]]==1) LiHyNu[i]=LiHyNu[i]+1;
      }
    }

  }

/* Determine the number of acceptor and donor vectors in the fragment */
  *FrAcNu=0;
  *FrDoNu=0;

  for (i=1;i<=FrAtNu;i++) {

    if (FrAtEl_nu[i]==7) {
      if ((LiToNu[i]==2)&&(LiHyNu[i]==0)) {
        *FrAcNu=*FrAcNu+1;
        fprintf(FPaOut,"Atom N %5d     1 acceptor vector\n",i);
      }
      else if ((LiToNu[i]==3)&&(LiHyNu[i]==1)) {
        *FrDoNu=*FrDoNu+1;
        fprintf(FPaOut,"Atom N %5d     1 donor vector\n",i);
      }
      else if ((LiToNu[i]==3)&&(LiHyNu[i]==2)) {
        *FrDoNu=*FrDoNu+2;
        fprintf(FPaOut,"Atom N %5d     2 donor vectors\n",i);
      }
      else if ((LiToNu[i]==4)&&(LiHyNu[i]==1)) {
        *FrDoNu=*FrDoNu+1;
        fprintf(FPaOut,"Atom N %5d     1 donor vector\n",i);
      }
      else if ((LiToNu[i]==4)&&(LiHyNu[i]==2)) {
        *FrDoNu=*FrDoNu+2;
        fprintf(FPaOut,"Atom N %5d     2 donor vectors\n",i);
      }
      else if ((LiToNu[i]==4)&&(LiHyNu[i]==3)) {
        *FrDoNu=*FrDoNu+3;
        fprintf(FPaOut,"Atom N %5d     3 donor vectors\n",i);
      }
    }
    else if (FrAtEl_nu[i]==8) {
      if ((LiToNu[i]==1)&&(LiHyNu[i]==0)) {
        *FrAcNu=*FrAcNu+2;
        fprintf(FPaOut,"Atom O %5d     2 acceptor vectors\n",i);
      }
      else if ((LiToNu[i]==2)&&(LiHyNu[i]==0)) {
        if ((HybFrAt[LiAtom[i][1]]==3)||(HybFrAt[LiAtom[i][2]]==3)) {
          *FrAcNu=*FrAcNu+2;
          fprintf(FPaOut,"Atom O %5d     2 acceptor vectors\n",i);
        }
        else {
          *FrAcNu=*FrAcNu+1;
          fprintf(FPaOut,"Atom O %5d     1 acceptor vector\n",i);
        }
      }
      else if ((LiToNu[i]==2)&&(LiHyNu[i]==1)) {
        if ((HybFrAt[LiAtom[i][1]]==3)||(HybFrAt[LiAtom[i][2]]==3)) {
          *FrDoNu=*FrDoNu+1;
          *FrAcNu=*FrAcNu+2;
          fprintf(FPaOut,"Atom O %5d     1 donor and 2 acceptor vectors\n",i);
        }
        else {
          *FrDoNu=*FrDoNu+1;
          *FrAcNu=*FrAcNu+1;
          fprintf(FPaOut,"Atom O %5d     1 donor and 1 acceptor vectors\n",i);
        }
      }
      /* dey added section in order to be able to dock water */
      else if ((LiToNu[i]==2)&&(LiHyNu[i]==2)) {
        if ((HybFrAt[LiAtom[i][1]]==3)||(HybFrAt[LiAtom[i][2]]==3)) {
          *FrDoNu=*FrDoNu+2;
          *FrAcNu=*FrAcNu+2;
          fprintf(FPaOut,"Atom O %5d     2 donor and 2 acceptor vectors\n",i);
        }
        else {
          *FrDoNu=*FrDoNu+1;
          *FrAcNu=*FrAcNu+1;
          fprintf(FPaOut,"Atom O %5d     1 donor and 1 acceptor vectors\n",i);
        }
      }
    }
    else if (FrAtEl_nu[i]==16) {
      if ((LiToNu[i]==2)&&(LiHyNu[i]==1)) {
        if ((HybFrAt[LiAtom[i][1]]==3)||(HybFrAt[LiAtom[i][2]]==3)) {
          *FrDoNu=*FrDoNu+1;
          *FrAcNu=*FrAcNu+2;
          fprintf(FPaOut,"Atom S %5d     1 donor and 2 acceptor vectors\n",i);
        }
      }
    }

  }

  if (*FrAcNu+*FrDoNu) fprintf(FPaOut,"\n");

/* Allocate the memory */
  *FrDATy=ivector(1,*FrAcNu+*FrDoNu);
  *FrDAAt=ivector(1,*FrAcNu+*FrDoNu);
  *FrHydN=ivector(1,*FrAcNu+*FrDoNu);
  *FrVeCo=dmatrix(1,*FrAcNu+*FrDoNu,1,6);

/* Determine the vectors for the donors and acceptors */
  DANumb=0;

  for (i=1;i<=FrAtNu;i++) {

/* Case where the fragment atom is a N */
    if (FrAtEl_nu[i]==7) {

      if ((LiToNu[i]==2)&&(LiHyNu[i]==0)) {
        AcAtNu=i;
        HelpA1=LiAtom[i][1];
        HelpA2=LiAtom[i][2];

        DANumb=DANumb+1;
        (*FrDATy)[DANumb]=1;
        (*FrDAAt)[DANumb]=AcAtNu;
        (*FrHydN)[DANumb]=0;
        for (k=4;k<=6;k++)
          (*FrVeCo)[DANumb][k]=FrCoor[AcAtNu][k-3];
        Angl=PlaAng(FrCoor[AcAtNu][1],FrCoor[AcAtNu][2],FrCoor[AcAtNu][3],
                    FrCoor[HelpA1][1],FrCoor[HelpA1][2],FrCoor[HelpA1][3],
                    FrCoor[HelpA2][1],FrCoor[HelpA2][2],FrCoor[HelpA2][3]);
        RotPla(FrCoor[HelpA2][1],FrCoor[HelpA2][2],FrCoor[HelpA2][3],
               FrCoor[AcAtNu][1],FrCoor[AcAtNu][2],FrCoor[AcAtNu][3],
               FrCoor[HelpA1][1],FrCoor[HelpA1][2],FrCoor[HelpA1][3],
               /*(PiT180*180*2-Angl)/2 clangini*/ (M_PI*2-Angl)/2,
               &(*FrVeCo)[DANumb][1],
               &(*FrVeCo)[DANumb][2],&(*FrVeCo)[DANumb][3]);
        PoCoVe(FrCoor[AcAtNu][1],FrCoor[AcAtNu][2],FrCoor[AcAtNu][3],
               (*FrVeCo)[DANumb][1],(*FrVeCo)[DANumb][2],(*FrVeCo)[DANumb][3],
               1.0,&(*FrVeCo)[DANumb][1],&(*FrVeCo)[DANumb][2],
               &(*FrVeCo)[DANumb][3]);
      }

      else if ((LiToNu[i]==3)&&(LiHyNu[i]==1)) {
        DoAtNu=i;
        for (j=1;j<=3;j++) {
          if (FrAtEl_nu[LiAtom[i][j]]==1) HyAtNu=LiAtom[i][j];
        }

        DANumb=DANumb+1;
        (*FrDATy)[DANumb]=0;
        (*FrDAAt)[DANumb]=DoAtNu;
        (*FrHydN)[DANumb]=HyAtNu;
        for (k=1;k<=3;k++)
          (*FrVeCo)[DANumb][k]=FrCoor[DoAtNu][k];
        PoCoVe(FrCoor[DoAtNu][1],FrCoor[DoAtNu][2],FrCoor[DoAtNu][3],
               FrCoor[HyAtNu][1],FrCoor[HyAtNu][2],FrCoor[HyAtNu][3],
               1.0,&(*FrVeCo)[DANumb][4],&(*FrVeCo)[DANumb][5],
               &(*FrVeCo)[DANumb][6]);
      }

      else if ((LiToNu[i]==3)&&(LiHyNu[i]==2)) {
        DoAtNu=i;
        HyAtNu=0;
        for (j=1;j<=3;j++) {
          if (FrAtEl_nu[LiAtom[i][j]]==1) {
            if (!HyAtNu)
              HyAtNu=LiAtom[i][j];
            else
              HyAtN2=LiAtom[i][j];
          }
        }

        DANumb=DANumb+1;
        (*FrDATy)[DANumb]=0;
        (*FrDAAt)[DANumb]=DoAtNu;
        (*FrHydN)[DANumb]=HyAtNu;
        for (k=1;k<=3;k++)
          (*FrVeCo)[DANumb][k]=FrCoor[DoAtNu][k];
        PoCoVe(FrCoor[DoAtNu][1],FrCoor[DoAtNu][2],FrCoor[DoAtNu][3],
               FrCoor[HyAtNu][1],FrCoor[HyAtNu][2],FrCoor[HyAtNu][3],
               1.0,&(*FrVeCo)[DANumb][4],&(*FrVeCo)[DANumb][5],
               &(*FrVeCo)[DANumb][6]);

        DANumb=DANumb+1;
        (*FrDATy)[DANumb]=0;
        (*FrDAAt)[DANumb]=DoAtNu;
        (*FrHydN)[DANumb]=HyAtN2;
        for (k=1;k<=3;k++)
          (*FrVeCo)[DANumb][k]=FrCoor[DoAtNu][k];
        PoCoVe(FrCoor[DoAtNu][1],FrCoor[DoAtNu][2],FrCoor[DoAtNu][3],
               FrCoor[HyAtN2][1],FrCoor[HyAtN2][2],FrCoor[HyAtN2][3],
               1.0,&(*FrVeCo)[DANumb][4],&(*FrVeCo)[DANumb][5],
               &(*FrVeCo)[DANumb][6]);
      }

      else if ((LiToNu[i]==4)&&(LiHyNu[i]==1)) {
        DoAtNu=i;
        for (j=1;j<=4;j++) {
          if (FrAtEl_nu[LiAtom[i][j]]==1) HyAtNu=LiAtom[i][j];
        }

        DANumb=DANumb+1;
        (*FrDATy)[DANumb]=0;
        (*FrDAAt)[DANumb]=DoAtNu;
        (*FrHydN)[DANumb]=HyAtNu;
        for (k=1;k<=3;k++)
          (*FrVeCo)[DANumb][k]=FrCoor[DoAtNu][k];
        PoCoVe(FrCoor[DoAtNu][1],FrCoor[DoAtNu][2],FrCoor[DoAtNu][3],
               FrCoor[HyAtNu][1],FrCoor[HyAtNu][2],FrCoor[HyAtNu][3],
               1.0,&(*FrVeCo)[DANumb][4],&(*FrVeCo)[DANumb][5],
               &(*FrVeCo)[DANumb][6]);
      }

      else if ((LiToNu[i]==4)&&(LiHyNu[i]==2)) {
        DoAtNu=i;
        HyAtNu=0;
        for (j=1;j<=4;j++) {
          if (FrAtEl_nu[LiAtom[i][j]]==1) {
            if (!HyAtNu)
              HyAtNu=LiAtom[i][j];
            else
              HyAtN2=LiAtom[i][j];
          }
        }

        DANumb=DANumb+1;
        (*FrDATy)[DANumb]=0;
        (*FrDAAt)[DANumb]=DoAtNu;
        (*FrHydN)[DANumb]=HyAtNu;
        for (k=1;k<=3;k++)
          (*FrVeCo)[DANumb][k]=FrCoor[DoAtNu][k];
        PoCoVe(FrCoor[DoAtNu][1],FrCoor[DoAtNu][2],FrCoor[DoAtNu][3],
               FrCoor[HyAtNu][1],FrCoor[HyAtNu][2],FrCoor[HyAtNu][3],
               1.0,&(*FrVeCo)[DANumb][4],&(*FrVeCo)[DANumb][5],
               &(*FrVeCo)[DANumb][6]);

        DANumb=DANumb+1;
        (*FrDATy)[DANumb]=0;
        (*FrDAAt)[DANumb]=DoAtNu;
        (*FrHydN)[DANumb]=HyAtN2;
        for (k=1;k<=3;k++)
          (*FrVeCo)[DANumb][k]=FrCoor[DoAtNu][k];
        PoCoVe(FrCoor[DoAtNu][1],FrCoor[DoAtNu][2],FrCoor[DoAtNu][3],
               FrCoor[HyAtN2][1],FrCoor[HyAtN2][2],FrCoor[HyAtN2][3],
               1.0,&(*FrVeCo)[DANumb][4],&(*FrVeCo)[DANumb][5],
               &(*FrVeCo)[DANumb][6]);
      }

      else if ((LiToNu[i]==4)&&(LiHyNu[i]==3)) {
        DoAtNu=i;
        HyAtNu=0;
        HyAtN2=0;
        for (j=1;j<=4;j++) {
          if (FrAtEl_nu[LiAtom[i][j]]==1) {
            if (!HyAtNu)
              HyAtNu=LiAtom[i][j];
            else if (!HyAtN2)
              HyAtN2=LiAtom[i][j];
            else
              HyAtN3=LiAtom[i][j];
          }
        }

        DANumb=DANumb+1;
        (*FrDATy)[DANumb]=0;
        (*FrDAAt)[DANumb]=DoAtNu;
        (*FrHydN)[DANumb]=HyAtNu;
        for (k=1;k<=3;k++)
          (*FrVeCo)[DANumb][k]=FrCoor[DoAtNu][k];
        PoCoVe(FrCoor[DoAtNu][1],FrCoor[DoAtNu][2],FrCoor[DoAtNu][3],
               FrCoor[HyAtNu][1],FrCoor[HyAtNu][2],FrCoor[HyAtNu][3],
               1.0,&(*FrVeCo)[DANumb][4],&(*FrVeCo)[DANumb][5],
               &(*FrVeCo)[DANumb][6]);

        DANumb=DANumb+1;
        (*FrDATy)[DANumb]=0;
        (*FrDAAt)[DANumb]=DoAtNu;
        (*FrHydN)[DANumb]=HyAtN2;
        for (k=1;k<=3;k++)
          (*FrVeCo)[DANumb][k]=FrCoor[DoAtNu][k];
        PoCoVe(FrCoor[DoAtNu][1],FrCoor[DoAtNu][2],FrCoor[DoAtNu][3],
               FrCoor[HyAtN2][1],FrCoor[HyAtN2][2],FrCoor[HyAtN2][3],
               1.0,&(*FrVeCo)[DANumb][4],&(*FrVeCo)[DANumb][5],
               &(*FrVeCo)[DANumb][6]);

        DANumb=DANumb+1;
        (*FrDATy)[DANumb]=0;
        (*FrDAAt)[DANumb]=DoAtNu;
        (*FrHydN)[DANumb]=HyAtN3;
        for (k=1;k<=3;k++)
          (*FrVeCo)[DANumb][k]=FrCoor[DoAtNu][k];
        PoCoVe(FrCoor[DoAtNu][1],FrCoor[DoAtNu][2],FrCoor[DoAtNu][3],
               FrCoor[HyAtN3][1],FrCoor[HyAtN3][2],FrCoor[HyAtN3][3],
               1.0,&(*FrVeCo)[DANumb][4],&(*FrVeCo)[DANumb][5],
               &(*FrVeCo)[DANumb][6]);
      }

    }

/* Case where the fragment atom is a O */
    else if (FrAtEl_nu[i]==8) {

      if ((LiToNu[i]==1)&&(LiHyNu[i]==0)) {
        AcAtNu=i;
        HelpA1=LiAtom[i][1];
        HelpA2=0;
        for (j=1;((j<=FrBdNu)&&(!HelpA2));j++) {
          if ((FrBdAr[j][1]==HelpA1)&&(FrBdAr[j][2]!=i))
            HelpA2=FrBdAr[j][2];
          else if ((FrBdAr[j][2]==HelpA1)&&(FrBdAr[j][1]!=i))
            HelpA2=FrBdAr[j][1];
        }

        if (!HelpA2) {
          fprintf(FPaOut,"WARNING Was not able to find a non-O linked ");
          fprintf(FPaOut,"atom for atom %d\n\n",HelpA1);
          fclose(FPaOut);
          FPaOut=fopen(OutFil,"a");
        }

        DANumb=DANumb+1;
        (*FrDATy)[DANumb]=1;
        (*FrDAAt)[DANumb]=AcAtNu;
        (*FrHydN)[DANumb]=0;
        for (k=4;k<=6;k++)
          (*FrVeCo)[DANumb][k]=FrCoor[AcAtNu][k-3];
        RotPla(FrCoor[HelpA2][1],FrCoor[HelpA2][2],FrCoor[HelpA2][3],
               FrCoor[AcAtNu][1],FrCoor[AcAtNu][2],FrCoor[AcAtNu][3],
               FrCoor[HelpA1][1],FrCoor[HelpA1][2],FrCoor[HelpA1][3],
               135*PiT180,&(*FrVeCo)[DANumb][1],&(*FrVeCo)[DANumb][2],
               &(*FrVeCo)[DANumb][3]);
        PoCoVe(FrCoor[AcAtNu][1],FrCoor[AcAtNu][2],FrCoor[AcAtNu][3],
               (*FrVeCo)[DANumb][1],(*FrVeCo)[DANumb][2],(*FrVeCo)[DANumb][3],
               1.0,&(*FrVeCo)[DANumb][1],&(*FrVeCo)[DANumb][2],
               &(*FrVeCo)[DANumb][3]);

        DANumb=DANumb+1;
        (*FrDATy)[DANumb]=1;
        (*FrDAAt)[DANumb]=AcAtNu;
        (*FrHydN)[DANumb]=0;
        for (k=4;k<=6;k++)
          (*FrVeCo)[DANumb][k]=FrCoor[AcAtNu][k-3];
        RotPla(FrCoor[HelpA2][1],FrCoor[HelpA2][2],FrCoor[HelpA2][3],
               FrCoor[AcAtNu][1],FrCoor[AcAtNu][2],FrCoor[AcAtNu][3],
               FrCoor[HelpA1][1],FrCoor[HelpA1][2],FrCoor[HelpA1][3],
               225*PiT180,&(*FrVeCo)[DANumb][1],&(*FrVeCo)[DANumb][2],
               &(*FrVeCo)[DANumb][3]);
        PoCoVe(FrCoor[AcAtNu][1],FrCoor[AcAtNu][2],FrCoor[AcAtNu][3],
               (*FrVeCo)[DANumb][1],(*FrVeCo)[DANumb][2],(*FrVeCo)[DANumb][3],
               1.0,&(*FrVeCo)[DANumb][1],&(*FrVeCo)[DANumb][2],
               &(*FrVeCo)[DANumb][3]);
      }

      else if ((LiToNu[i]==2)&&(LiHyNu[i]==0)) {
        AcAtNu=i;
        HelpA1=LiAtom[i][1];
        HelpA2=LiAtom[i][2];
        if ((HybFrAt[HelpA1]==3)||(HybFrAt[HelpA2]==3)) {

          DANumb=DANumb+1;
          (*FrDATy)[DANumb]=1;
          (*FrDAAt)[DANumb]=AcAtNu;
          (*FrHydN)[DANumb]=0;
          for (k=4;k<=6;k++)
            (*FrVeCo)[DANumb][k]=FrCoor[AcAtNu][k-3];
          RoArVe(FrCoor[HelpA2][1]-FrCoor[AcAtNu][1],FrCoor[HelpA2][2]-
                 FrCoor[AcAtNu][2],FrCoor[HelpA2][3]-FrCoor[AcAtNu][3],
                 FrCoor[HelpA1][1]-FrCoor[AcAtNu][1],FrCoor[HelpA1][2]-
                 FrCoor[AcAtNu][2],FrCoor[HelpA1][3]-FrCoor[AcAtNu][3],
                 120*PiT180,&(*FrVeCo)[DANumb][1],&(*FrVeCo)[DANumb][2],
                 &(*FrVeCo)[DANumb][3]);
          (*FrVeCo)[DANumb][1]=(*FrVeCo)[DANumb][1]+FrCoor[AcAtNu][1];
          (*FrVeCo)[DANumb][2]=(*FrVeCo)[DANumb][2]+FrCoor[AcAtNu][2];
          (*FrVeCo)[DANumb][3]=(*FrVeCo)[DANumb][3]+FrCoor[AcAtNu][3];
          PoCoVe(FrCoor[AcAtNu][1],FrCoor[AcAtNu][2],FrCoor[AcAtNu][3],
                 (*FrVeCo)[DANumb][1],(*FrVeCo)[DANumb][2],(*FrVeCo)[DANumb][3],
                 1.0,&(*FrVeCo)[DANumb][1],&(*FrVeCo)[DANumb][2],
                 &(*FrVeCo)[DANumb][3]);

          DANumb=DANumb+1;
          (*FrDATy)[DANumb]=1;
          (*FrDAAt)[DANumb]=AcAtNu;
          (*FrHydN)[DANumb]=0;
          for (k=4;k<=6;k++)
            (*FrVeCo)[DANumb][k]=FrCoor[AcAtNu][k-3];
          RoArVe(FrCoor[HelpA2][1]-FrCoor[AcAtNu][1],FrCoor[HelpA2][2]-
                 FrCoor[AcAtNu][2],FrCoor[HelpA2][3]-FrCoor[AcAtNu][3],
                 FrCoor[HelpA1][1]-FrCoor[AcAtNu][1],FrCoor[HelpA1][2]-
                 FrCoor[AcAtNu][2],FrCoor[HelpA1][3]-FrCoor[AcAtNu][3],
                 240*PiT180,&(*FrVeCo)[DANumb][1],&(*FrVeCo)[DANumb][2],
                 &(*FrVeCo)[DANumb][3]);
          (*FrVeCo)[DANumb][1]=(*FrVeCo)[DANumb][1]+FrCoor[AcAtNu][1];
          (*FrVeCo)[DANumb][2]=(*FrVeCo)[DANumb][2]+FrCoor[AcAtNu][2];
          (*FrVeCo)[DANumb][3]=(*FrVeCo)[DANumb][3]+FrCoor[AcAtNu][3];
          PoCoVe(FrCoor[AcAtNu][1],FrCoor[AcAtNu][2],FrCoor[AcAtNu][3],
                 (*FrVeCo)[DANumb][1],(*FrVeCo)[DANumb][2],(*FrVeCo)[DANumb][3],
                 1.0,&(*FrVeCo)[DANumb][1],&(*FrVeCo)[DANumb][2],
                 &(*FrVeCo)[DANumb][3]);
        }
        else {
          DANumb=DANumb+1;
          (*FrDATy)[DANumb]=1;
          (*FrDAAt)[DANumb]=AcAtNu;
          (*FrHydN)[DANumb]=0;
          for (k=4;k<=6;k++)
            (*FrVeCo)[DANumb][k]=FrCoor[AcAtNu][k-3];
          Angl=PlaAng(FrCoor[AcAtNu][1],FrCoor[AcAtNu][2],FrCoor[AcAtNu][3],
                      FrCoor[HelpA1][1],FrCoor[HelpA1][2],FrCoor[HelpA1][3],
                      FrCoor[HelpA2][1],FrCoor[HelpA2][2],FrCoor[HelpA2][3]);
          RotPla(FrCoor[HelpA2][1],FrCoor[HelpA2][2],FrCoor[HelpA2][3],
                 FrCoor[AcAtNu][1],FrCoor[AcAtNu][2],FrCoor[AcAtNu][3],
                 FrCoor[HelpA1][1],FrCoor[HelpA1][2],FrCoor[HelpA1][3],
                 /*(PiT180*180*2-Angl)/2 clangini*/
                 (M_PI*2-Angl)/2,&(*FrVeCo)[DANumb][1],
                 &(*FrVeCo)[DANumb][2],&(*FrVeCo)[DANumb][3]);
          PoCoVe(FrCoor[AcAtNu][1],FrCoor[AcAtNu][2],FrCoor[AcAtNu][3],
                 (*FrVeCo)[DANumb][1],(*FrVeCo)[DANumb][2],(*FrVeCo)[DANumb][3],
                 1.0,&(*FrVeCo)[DANumb][1],&(*FrVeCo)[DANumb][2],
                 &(*FrVeCo)[DANumb][3]);
        }
      }

      else if ((LiToNu[i]==2)&&(LiHyNu[i]>=1)) { /* Dey check -> was ==, water could not be docked */
        DoAtNu=i;
        AcAtNu=i;
        if (FrAtEl_nu[LiAtom[i][1]]==1) {
          HyAtNu=LiAtom[i][1];
          HelpA1=LiAtom[i][2];
        }
        else {
          HyAtNu=LiAtom[i][2];
          HelpA1=LiAtom[i][1];
        }
        if (HybFrAt[HelpA1]==3) {

          DANumb=DANumb+1;
          (*FrDATy)[DANumb]=0;
          (*FrDAAt)[DANumb]=DoAtNu;
          (*FrHydN)[DANumb]=HyAtNu;
          for (k=1;k<=3;k++)
            (*FrVeCo)[DANumb][k]=FrCoor[DoAtNu][k];
          PoCoVe(FrCoor[DoAtNu][1],FrCoor[DoAtNu][2],FrCoor[DoAtNu][3],
                 FrCoor[HyAtNu][1],FrCoor[HyAtNu][2],FrCoor[HyAtNu][3],
                 1.0,&(*FrVeCo)[DANumb][4],&(*FrVeCo)[DANumb][5],
                 &(*FrVeCo)[DANumb][6]);

          DANumb=DANumb+1;
          (*FrDATy)[DANumb]=1;
          (*FrDAAt)[DANumb]=AcAtNu;
          (*FrHydN)[DANumb]=0;
          for (k=4;k<=6;k++)
            (*FrVeCo)[DANumb][k]=FrCoor[AcAtNu][k-3];
          RoArVe(FrCoor[HyAtNu][1]-FrCoor[AcAtNu][1],FrCoor[HyAtNu][2]-
                 FrCoor[AcAtNu][2],FrCoor[HyAtNu][3]-FrCoor[AcAtNu][3],
                 FrCoor[HelpA1][1]-FrCoor[AcAtNu][1],FrCoor[HelpA1][2]-
                 FrCoor[AcAtNu][2],FrCoor[HelpA1][3]-FrCoor[AcAtNu][3],
                 120*PiT180,&(*FrVeCo)[DANumb][1],&(*FrVeCo)[DANumb][2],
                 &(*FrVeCo)[DANumb][3]);
          (*FrVeCo)[DANumb][1]=(*FrVeCo)[DANumb][1]+FrCoor[AcAtNu][1];
          (*FrVeCo)[DANumb][2]=(*FrVeCo)[DANumb][2]+FrCoor[AcAtNu][2];
          (*FrVeCo)[DANumb][3]=(*FrVeCo)[DANumb][3]+FrCoor[AcAtNu][3];
          PoCoVe(FrCoor[AcAtNu][1],FrCoor[AcAtNu][2],FrCoor[AcAtNu][3],
                 (*FrVeCo)[DANumb][1],(*FrVeCo)[DANumb][2],(*FrVeCo)[DANumb][3],
                 1.0,&(*FrVeCo)[DANumb][1],&(*FrVeCo)[DANumb][2],
                 &(*FrVeCo)[DANumb][3]);

          DANumb=DANumb+1;
          (*FrDATy)[DANumb]=1;
          (*FrDAAt)[DANumb]=AcAtNu;
          (*FrHydN)[DANumb]=0;
          for (k=4;k<=6;k++)
            (*FrVeCo)[DANumb][k]=FrCoor[AcAtNu][k-3];
          RoArVe(FrCoor[HyAtNu][1]-FrCoor[AcAtNu][1],FrCoor[HyAtNu][2]-
                 FrCoor[AcAtNu][2],FrCoor[HyAtNu][3]-FrCoor[AcAtNu][3],
                 FrCoor[HelpA1][1]-FrCoor[AcAtNu][1],FrCoor[HelpA1][2]-
                 FrCoor[AcAtNu][2],FrCoor[HelpA1][3]-FrCoor[AcAtNu][3],
                 240*PiT180,&(*FrVeCo)[DANumb][1],&(*FrVeCo)[DANumb][2],
                 &(*FrVeCo)[DANumb][3]);
          (*FrVeCo)[DANumb][1]=(*FrVeCo)[DANumb][1]+FrCoor[AcAtNu][1];
          (*FrVeCo)[DANumb][2]=(*FrVeCo)[DANumb][2]+FrCoor[AcAtNu][2];
          (*FrVeCo)[DANumb][3]=(*FrVeCo)[DANumb][3]+FrCoor[AcAtNu][3];
          PoCoVe(FrCoor[AcAtNu][1],FrCoor[AcAtNu][2],FrCoor[AcAtNu][3],
                 (*FrVeCo)[DANumb][1],(*FrVeCo)[DANumb][2],(*FrVeCo)[DANumb][3],
                 1.0,&(*FrVeCo)[DANumb][1],&(*FrVeCo)[DANumb][2],
                 &(*FrVeCo)[DANumb][3]);
        }
        else {

          DANumb=DANumb+1;
          (*FrDATy)[DANumb]=1;
          (*FrDAAt)[DANumb]=AcAtNu;
          (*FrHydN)[DANumb]=0;
          for (k=4;k<=6;k++)
            (*FrVeCo)[DANumb][k]=FrCoor[AcAtNu][k-3];
          Angl=PlaAng(FrCoor[AcAtNu][1],FrCoor[AcAtNu][2],FrCoor[AcAtNu][3],
                      FrCoor[HelpA1][1],FrCoor[HelpA1][2],FrCoor[HelpA1][3],
                      FrCoor[HyAtNu][1],FrCoor[HyAtNu][2],FrCoor[HyAtNu][3]);
          RotPla(FrCoor[HyAtNu][1],FrCoor[HyAtNu][2],FrCoor[HyAtNu][3],
                 FrCoor[AcAtNu][1],FrCoor[AcAtNu][2],FrCoor[AcAtNu][3],
                 FrCoor[HelpA1][1],FrCoor[HelpA1][2],FrCoor[HelpA1][3],
                 /*(PiT180*180*2-Angl)/2 clangini*/
                 (M_PI*2-Angl)/2,&(*FrVeCo)[DANumb][1],
                 &(*FrVeCo)[DANumb][2],&(*FrVeCo)[DANumb][3]);
          PoCoVe(FrCoor[AcAtNu][1],FrCoor[AcAtNu][2],FrCoor[AcAtNu][3],
                 (*FrVeCo)[DANumb][1],(*FrVeCo)[DANumb][2],(*FrVeCo)[DANumb][3],
                 1.0,&(*FrVeCo)[DANumb][1],&(*FrVeCo)[DANumb][2],
                 &(*FrVeCo)[DANumb][3]);

          DANumb=DANumb+1;
          (*FrDATy)[DANumb]=0;
          (*FrDAAt)[DANumb]=DoAtNu;
          (*FrHydN)[DANumb]=HyAtNu;
          for (k=1;k<=3;k++)
            (*FrVeCo)[DANumb][k]=FrCoor[DoAtNu][k];
          PoCoVe(FrCoor[DoAtNu][1],FrCoor[DoAtNu][2],FrCoor[DoAtNu][3],
                 FrCoor[HyAtNu][1],FrCoor[HyAtNu][2],FrCoor[HyAtNu][3],
                 1.0,&(*FrVeCo)[DANumb][4],&(*FrVeCo)[DANumb][5],
                 &(*FrVeCo)[DANumb][6]);
        }
      }

    }

/* Case where the fragment atom is a S */
    else if (FrAtEl_nu[i]==16) {

      if ((LiToNu[i]==2)&&(LiHyNu[i]==1)) {
        DoAtNu=i;
        AcAtNu=i;
        if (FrAtEl_nu[LiAtom[i][1]]==1) {
          HyAtNu=LiAtom[i][1];
          HelpA1=LiAtom[i][2];
        }
        else {
          HyAtNu=LiAtom[i][2];
          HelpA1=LiAtom[i][1];
        }
        if (HybFrAt[HelpA1]==3) {

          DANumb=DANumb+1;
          (*FrDATy)[DANumb]=0;
          (*FrDAAt)[DANumb]=DoAtNu;
          (*FrHydN)[DANumb]=HyAtNu;
          for (k=1;k<=3;k++)
            (*FrVeCo)[DANumb][k]=FrCoor[DoAtNu][k];
          PoCoVe(FrCoor[DoAtNu][1],FrCoor[DoAtNu][2],FrCoor[DoAtNu][3],
                 FrCoor[HyAtNu][1],FrCoor[HyAtNu][2],FrCoor[HyAtNu][3],
                 1.0,&(*FrVeCo)[DANumb][4],&(*FrVeCo)[DANumb][5],
                 &(*FrVeCo)[DANumb][6]);

          DANumb=DANumb+1;
          (*FrDATy)[DANumb]=1;
          (*FrDAAt)[DANumb]=AcAtNu;
          (*FrHydN)[DANumb]=0;
          for (k=4;k<=6;k++)
            (*FrVeCo)[DANumb][k]=FrCoor[AcAtNu][k-3];
          RoArVe(FrCoor[HyAtNu][1]-FrCoor[AcAtNu][1],FrCoor[HyAtNu][2]-
                 FrCoor[AcAtNu][2],FrCoor[HyAtNu][3]-FrCoor[AcAtNu][3],
                 FrCoor[HelpA1][1]-FrCoor[AcAtNu][1],FrCoor[HelpA1][2]-
                 FrCoor[AcAtNu][2],FrCoor[HelpA1][3]-FrCoor[AcAtNu][3],
                 120*PiT180,&(*FrVeCo)[DANumb][1],&(*FrVeCo)[DANumb][2],
                 &(*FrVeCo)[DANumb][3]);
          (*FrVeCo)[DANumb][1]=(*FrVeCo)[DANumb][1]+FrCoor[AcAtNu][1];
          (*FrVeCo)[DANumb][2]=(*FrVeCo)[DANumb][2]+FrCoor[AcAtNu][2];
          (*FrVeCo)[DANumb][3]=(*FrVeCo)[DANumb][3]+FrCoor[AcAtNu][3];
          PoCoVe(FrCoor[AcAtNu][1],FrCoor[AcAtNu][2],FrCoor[AcAtNu][3],
                 (*FrVeCo)[DANumb][1],(*FrVeCo)[DANumb][2],(*FrVeCo)[DANumb][3],
                 1.0,&(*FrVeCo)[DANumb][1],&(*FrVeCo)[DANumb][2],
                 &(*FrVeCo)[DANumb][3]);

          DANumb=DANumb+1;
          (*FrDATy)[DANumb]=1;
          (*FrDAAt)[DANumb]=AcAtNu;
          (*FrHydN)[DANumb]=0;
          for (k=4;k<=6;k++)
            (*FrVeCo)[DANumb][k]=FrCoor[AcAtNu][k-3];
          RoArVe(FrCoor[HyAtNu][1]-FrCoor[AcAtNu][1],FrCoor[HyAtNu][2]-
                 FrCoor[AcAtNu][2],FrCoor[HyAtNu][3]-FrCoor[AcAtNu][3],
                 FrCoor[HelpA1][1]-FrCoor[AcAtNu][1],FrCoor[HelpA1][2]-
                 FrCoor[AcAtNu][2],FrCoor[HelpA1][3]-FrCoor[AcAtNu][3],
                 240*PiT180,&(*FrVeCo)[DANumb][1],&(*FrVeCo)[DANumb][2],
                 &(*FrVeCo)[DANumb][3]);
          (*FrVeCo)[DANumb][1]=(*FrVeCo)[DANumb][1]+FrCoor[AcAtNu][1];
          (*FrVeCo)[DANumb][2]=(*FrVeCo)[DANumb][2]+FrCoor[AcAtNu][2];
          (*FrVeCo)[DANumb][3]=(*FrVeCo)[DANumb][3]+FrCoor[AcAtNu][3];
          PoCoVe(FrCoor[AcAtNu][1],FrCoor[AcAtNu][2],FrCoor[AcAtNu][3],
                 (*FrVeCo)[DANumb][1],(*FrVeCo)[DANumb][2],(*FrVeCo)[DANumb][3],
                 1.0,&(*FrVeCo)[DANumb][1],&(*FrVeCo)[DANumb][2],
                 &(*FrVeCo)[DANumb][3]);
        }
      }

    }

  }

  if ((*FrAcNu+*FrDoNu)!=DANumb) {
    fprintf(FPaOut,"WARNING The expected total number of vectors is not ");
    fprintf(FPaOut,"satisfied\n\n");
/* remove::     fprintf(FPaOut,"%d %d %d\n",DANumb,*FrAcNu,*FrDoNu); */
  }


/*
  FilePa=fopen("./outputs/vectfra.mol2","w");
  fprintf(FilePa,"# TRIPOS MOL2 file generated by WITNOTP\n");
  fprintf(FilePa,"\n");
  fprintf(FilePa,"@<TRIPOS>MOLECULE\n");
  fprintf(FilePa,"vectfra\n");
  fprintf(FilePa,"%d 0 0 0 0\n",2*DANumb);
  fprintf(FilePa,"****\n");
  fprintf(FilePa,"USER_CHARGES\n");
  fprintf(FilePa,"INVALID_CHARGES\n");
  fprintf(FilePa,"\n");
  fprintf(FilePa,"@<TRIPOS>ATOM\n");

  for (i=1;i<=DANumb;i++) {
      fprintf(FilePa,"%-5d %-3d%11.4f%11.4f%11.4f ***** 1    ****  0.000\n",
              i*2-1,i*2-1,(*FrVeCo)[i][1],(*FrVeCo)[i][2],(*FrVeCo)[i][3]);
      fprintf(FilePa,"%-5d %-3d%11.4f%11.4f%11.4f ***** 1    ****  0.000\n",
              i*2,i*2,(*FrVeCo)[i][4],(*FrVeCo)[i][5],(*FrVeCo)[i][6]);
  }
  fclose(FilePa);
*/


  free_ivector(LiToNu,1,FrAtNu);
  free_ivector(LiHyNu,1,FrAtNu);
  free_imatrix(LiAtom,1,FrAtNu,1,10);



}
