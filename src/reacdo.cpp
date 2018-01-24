/*
*    This file is part of SEED.
*
*    Copyright (C) 2017, Caflisch Lab, University of Zurich
*
*    SEED is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    SEED is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "funct.h"
#include <stdlib.h> /* for exit-fct.*/

#ifdef ENABLE_MPI
  #include <mpi.h>
  #ifndef MASTERRANK
    #define MASTERRANK 0
  #endif
#endif

void ReAcDo(int BSAtNu,int *BSAtLi,int *ReAtEl_nu,double **ReCoor,
            int *HybReAt,int ReBdNu,int **ReBdAr,double SphAng,int SphPoN,
            int *ReAcNu,int *ReDoNu,int **ReDATy,int **ReDAAt,int **ReHydN,
            double ***ReVeCo,FILE *FPaOut,char *OutFil,int BSMeNu,int *BSMeAN,
            double **BSMeVE)
/* This function constructs the list of the donor and acceptor vectors for
   the binding site of the receptor. Then each of them is considered as the
   vector pointing towards the top of the part of the sphere on which other
   vectors are constructed (by the function MakSpV) :
   ReAcNu  number of acceptor vectors in the binding site of the receptor
   ReDoNu  number of donor vectors in the binding site of the receptor
   ReDATy  receptor donor acceptor type (0 for donor,1 for acceptor)
   ReDAAt  receptor donor or acceptor atom numbers
   ReHydN  receptor hydrogen number in the current vector (0 if no hydrogen
           involved)
   ReVeCo  coordinates of the vector whose orientation goes from the donor
           to the acceptor (for the binding site of the receptor)
   LiToNu  total number of linked atoms
   LiHyNu  number of linked hydrogen atoms
   LiAtom  array of linked atoms
   DoAtNu  donor atom number
   HyAtNu  hydrogen atom number
   HyAtN2  hydrogen atom number 2
   HyAtN3  hydrogen atom number 3
   AcAtNu  acceptor atom number
   HelpA1  help atom number 1
   HelpA2  help atom number 2
   DANumb  donor and acceptor number
   PoSpAN  actual number of points on the part of the sphere
   NumAnT  total number of theta angles
   NumAnP  number of phi angles for each theta angle
   IntPar  integer part
   ReaPar  real part
   Angl    angle
   diff    used to check if metal coordination sites are defined correctly
   NuOccu  number of occurrences */
{
  FILE *FilePa;
  int i,j,k,*LiToNu,*LiHyNu,**LiAtom,AcAtNu,DoAtNu,HelpA1,HelpA2,HyAtNu,
      HyAtN2,HyAtN3,DANumb,PoSpAN,NumAnT,NumAnP,NuOccu;
  double PiT180,Angl,ReaPar,IntPar,diff;
  double DIntPar;


  HyAtNu = -1;
  HyAtN2 = -1;
  HyAtN3 = -1;

  //PiT180=3.1415927/180;//clangini
  PiT180=M_PI/180;//clangini

/* Find the total number of linked atoms, the linked atoms and the number of
   linked hydrogen atoms for each atom of the receptor binding site */
  LiToNu=ivector(1,BSAtNu);
  LiHyNu=ivector(1,BSAtNu);
  LiAtom=imatrix(1,BSAtNu,1,10);

  for (i=1;i<=BSAtNu;i++) {
    LiToNu[i]=0;
    LiHyNu[i]=0;
    for (j=1;j<=ReBdNu;j++) {
      if (ReBdAr[j][1]==BSAtLi[i]) {
        LiToNu[i]=LiToNu[i]+1;
        LiAtom[i][LiToNu[i]]=ReBdAr[j][2];
        if (ReAtEl_nu[ReBdAr[j][2]]==1) LiHyNu[i]=LiHyNu[i]+1;
      }
      else if (ReBdAr[j][2]==BSAtLi[i]) {
        LiToNu[i]=LiToNu[i]+1;
        LiAtom[i][LiToNu[i]]=ReBdAr[j][1];
        if (ReAtEl_nu[ReBdAr[j][1]]==1) LiHyNu[i]=LiHyNu[i]+1;
      }
    }
  }

/* Count the actual number of points on the part of the sphere */
  PoSpAN=0;
  NumAnT=floor(sqrt(SphPoN*SphAng*PiT180/(2*(1-cos(SphAng*PiT180)))));
  for (i=1;i<=NumAnT;i++) {
    //ReaPar=modff(2*NumAnT*sinf((SphAng*PiT180/NumAnT)*i),&DIntPar);
    ReaPar=modf(2*NumAnT*sin((SphAng*PiT180/NumAnT)*i),&DIntPar);
    IntPar=DIntPar;
    if (ReaPar<=0.5)
      NumAnP=IntPar;
    else
      NumAnP=IntPar+1;
    PoSpAN=PoSpAN+NumAnP;
  }
/* Add the point on the top of the sphere part */
  PoSpAN=PoSpAN+1;

  fprintf(FPaOut,"Actual number of points on the sphere part :\n%d\n\n",PoSpAN);

/* Determine the number of acceptor and donor vectors in the receptor binding
   site */
  *ReAcNu=0;
  *ReDoNu=0;

  for (i=1;i<=BSAtNu;i++) {

    if (ReAtEl_nu[BSAtLi[i]]==7) {
      if ((LiToNu[i]==2)&&(LiHyNu[i]==0)) {
        *ReAcNu=*ReAcNu+1;
        fprintf(FPaOut,"Atom N %5d     1 acceptor vector\n",BSAtLi[i]);
      }
      else if ((LiToNu[i]==3)&&(LiHyNu[i]==1)) {
        *ReDoNu=*ReDoNu+1;
        fprintf(FPaOut,"Atom N %5d     1 donor vector\n",BSAtLi[i]);
      }
      else if ((LiToNu[i]==3)&&(LiHyNu[i]==2)) {
        *ReDoNu=*ReDoNu+2;
        fprintf(FPaOut,"Atom N %5d     2 donor vectors\n",BSAtLi[i]);
      }
      else if ((LiToNu[i]==4)&&(LiHyNu[i]==1)) {
        *ReDoNu=*ReDoNu+1;
        fprintf(FPaOut,"Atom N %5d     1 donor vector\n",BSAtLi[i]);
      }
      else if ((LiToNu[i]==4)&&(LiHyNu[i]==2)) {
        *ReDoNu=*ReDoNu+2;
        fprintf(FPaOut,"Atom N %5d     2 donor vectors\n",BSAtLi[i]);
      }
      else if ((LiToNu[i]==4)&&(LiHyNu[i]==3)) {
        *ReDoNu=*ReDoNu+3;
        fprintf(FPaOut,"Atom N %5d     3 donor vectors\n",BSAtLi[i]);
      }
    }
    else if (ReAtEl_nu[BSAtLi[i]]==8) {
      if ((LiToNu[i]==1)&&(LiHyNu[i]==0)) {
        *ReAcNu=*ReAcNu+2;
        fprintf(FPaOut,"Atom O %5d     2 acceptor vectors\n",BSAtLi[i]);
      }
      else if ((LiToNu[i]==2)&&(LiHyNu[i]==0)) {
        if ((HybReAt[LiAtom[i][1]]==3)||(HybReAt[LiAtom[i][2]]==3)) {
          *ReAcNu=*ReAcNu+2;
          fprintf(FPaOut,"Atom O %5d     2 acceptor vectors\n",BSAtLi[i]);
        }
        else {
          *ReAcNu=*ReAcNu+1;
          fprintf(FPaOut,"Atom O %5d     1 acceptor vector\n",BSAtLi[i]);
        }
      }
      else if ((LiToNu[i]==2)&&(LiHyNu[i]==1)) {
        if ((HybReAt[LiAtom[i][1]]==3)||(HybReAt[LiAtom[i][2]]==3)) {
          *ReDoNu=*ReDoNu+1;
          *ReAcNu=*ReAcNu+2;
          fprintf(FPaOut,"Atom O %5d     1 donor and 2 acceptor vectors\n",
                  BSAtLi[i]);
        }
        else {
          *ReDoNu=*ReDoNu+1;
          *ReAcNu=*ReAcNu+1;
          fprintf(FPaOut,"Atom O %5d     1 donor and 1 acceptor vectors\n",
                  BSAtLi[i]);
        }
      }
      else if ((LiToNu[i]==2)&&(LiHyNu[i]==2)) {
        *ReAcNu=*ReAcNu+2;
        *ReDoNu=*ReDoNu+2;
        fprintf(FPaOut,"Atom O %5d     2 acceptor and 2 donor vectors\n",
                BSAtLi[i]);
      }
    }
    else if (ReAtEl_nu[BSAtLi[i]]==16) {
      if ((LiToNu[i]==2)&&(LiHyNu[i]==1)) {
        if ((HybReAt[LiAtom[i][1]]==3)||(HybReAt[LiAtom[i][2]]==3)) {
          *ReDoNu=*ReDoNu+1;
          *ReAcNu=*ReAcNu+2;
          fprintf(FPaOut,"Atom S %5d     1 donor and 2 acceptor vectors\n",
                  BSAtLi[i]);
        }
      }
    }
    else if ((ReAtEl_nu[BSAtLi[i]]==25)||(ReAtEl_nu[BSAtLi[i]]==30)||
             (ReAtEl_nu[BSAtLi[i]]==13)||(ReAtEl_nu[BSAtLi[i]]==20)||
             (ReAtEl_nu[BSAtLi[i]]==29)||(ReAtEl_nu[BSAtLi[i]]==26)||
             (ReAtEl_nu[BSAtLi[i]]==19)||(ReAtEl_nu[BSAtLi[i]]==12)||
             (ReAtEl_nu[BSAtLi[i]]==14)) {
      NuOccu=0;
      for (j=1;j<=BSMeNu;j++) {
        if (BSMeAN[j]==BSAtLi[i]) NuOccu=NuOccu+1;
      }
      *ReDoNu=*ReDoNu+NuOccu;
      fprintf(FPaOut,"Atom Me %5d    %d donor vector(s)\n",BSAtLi[i],NuOccu);
    }

    fflush(stdout);/*new */
  }

  *ReAcNu=(*ReAcNu)*PoSpAN;
  *ReDoNu=(*ReDoNu)*PoSpAN;

  if (*ReAcNu+*ReDoNu) fprintf(FPaOut,"\n");

/* Allocate the memory */
  *ReDATy=ivector(1,*ReAcNu+*ReDoNu);
  *ReDAAt=ivector(1,*ReAcNu+*ReDoNu);
  *ReHydN=ivector(1,*ReAcNu+*ReDoNu);
  *ReVeCo=dmatrix(1,*ReAcNu+*ReDoNu,1,6);

/* **************************************************
   Determine the vectors for the donors and acceptors
   ************************************************** */
  DANumb=0;

  for (i=1;i<=BSAtNu;i++) {

/* Case where the receptor binding site atom is a N */
    if (ReAtEl_nu[BSAtLi[i]]==7) {

      if ((LiToNu[i]==2)&&(LiHyNu[i]==0)) {
        AcAtNu=BSAtLi[i];
        HelpA1=LiAtom[i][1];
        HelpA2=LiAtom[i][2];

        DANumb=DANumb+1;
        (*ReDATy)[DANumb]=1;
        (*ReDAAt)[DANumb]=AcAtNu;
        (*ReHydN)[DANumb]=0;
        for (k=4;k<=6;k++)
          (*ReVeCo)[DANumb][k]=ReCoor[AcAtNu][k-3];
        Angl=PlaAng(ReCoor[AcAtNu][1],ReCoor[AcAtNu][2],ReCoor[AcAtNu][3],
                    ReCoor[HelpA1][1],ReCoor[HelpA1][2],ReCoor[HelpA1][3],
                    ReCoor[HelpA2][1],ReCoor[HelpA2][2],ReCoor[HelpA2][3]);
        RotPla(ReCoor[HelpA2][1],ReCoor[HelpA2][2],ReCoor[HelpA2][3],
               ReCoor[AcAtNu][1],ReCoor[AcAtNu][2],ReCoor[AcAtNu][3],
               ReCoor[HelpA1][1],ReCoor[HelpA1][2],ReCoor[HelpA1][3],
               /*(PiT180*180*2-Angl)/2 clangini*/
               (M_PI*2-Angl)/2,&(*ReVeCo)[DANumb][1],
               &(*ReVeCo)[DANumb][2],&(*ReVeCo)[DANumb][3]);
        PoCoVe(ReCoor[AcAtNu][1],ReCoor[AcAtNu][2],ReCoor[AcAtNu][3],
               (*ReVeCo)[DANumb][1],(*ReVeCo)[DANumb][2],(*ReVeCo)[DANumb][3],
               1.0,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
               &(*ReVeCo)[DANumb][3]);
        MakSpV(1,AcAtNu,0,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
               *ReHydN,*ReVeCo);
      }

      else if ((LiToNu[i]==3)&&(LiHyNu[i]==1)) {
        DoAtNu=BSAtLi[i];
        for (j=1;j<=3;j++) {
          if (ReAtEl_nu[LiAtom[i][j]]==1) HyAtNu=LiAtom[i][j];
        }

        DANumb=DANumb+1;
        (*ReDATy)[DANumb]=0;
        (*ReDAAt)[DANumb]=DoAtNu;
        (*ReHydN)[DANumb]=HyAtNu;
        for (k=1;k<=3;k++)
          (*ReVeCo)[DANumb][k]=ReCoor[DoAtNu][k];
        PoCoVe(ReCoor[DoAtNu][1],ReCoor[DoAtNu][2],ReCoor[DoAtNu][3],
               ReCoor[HyAtNu][1],ReCoor[HyAtNu][2],ReCoor[HyAtNu][3],
               1.0,&(*ReVeCo)[DANumb][4],&(*ReVeCo)[DANumb][5],
               &(*ReVeCo)[DANumb][6]);
        MakSpV(0,DoAtNu,HyAtNu,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
               *ReHydN,*ReVeCo);
      }

      else if ((LiToNu[i]==3)&&(LiHyNu[i]==2)) {
        DoAtNu=BSAtLi[i];
        HyAtNu=0;
        for (j=1;j<=3;j++) {
          if (ReAtEl_nu[LiAtom[i][j]]==1) {
            if (!HyAtNu)
              HyAtNu=LiAtom[i][j];
            else
              HyAtN2=LiAtom[i][j];
          }
        }

        DANumb=DANumb+1;
        (*ReDATy)[DANumb]=0;
        (*ReDAAt)[DANumb]=DoAtNu;
        (*ReHydN)[DANumb]=HyAtNu;
        for (k=1;k<=3;k++)
          (*ReVeCo)[DANumb][k]=ReCoor[DoAtNu][k];
        PoCoVe(ReCoor[DoAtNu][1],ReCoor[DoAtNu][2],ReCoor[DoAtNu][3],
               ReCoor[HyAtNu][1],ReCoor[HyAtNu][2],ReCoor[HyAtNu][3],
               1.0,&(*ReVeCo)[DANumb][4],&(*ReVeCo)[DANumb][5],
               &(*ReVeCo)[DANumb][6]);
        MakSpV(0,DoAtNu,HyAtNu,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
               *ReHydN,*ReVeCo);

        DANumb=DANumb+1;
        (*ReDATy)[DANumb]=0;
        (*ReDAAt)[DANumb]=DoAtNu;
        (*ReHydN)[DANumb]=HyAtN2;
        for (k=1;k<=3;k++)
          (*ReVeCo)[DANumb][k]=ReCoor[DoAtNu][k];
        PoCoVe(ReCoor[DoAtNu][1],ReCoor[DoAtNu][2],ReCoor[DoAtNu][3],
               ReCoor[HyAtN2][1],ReCoor[HyAtN2][2],ReCoor[HyAtN2][3],
               1.0,&(*ReVeCo)[DANumb][4],&(*ReVeCo)[DANumb][5],
               &(*ReVeCo)[DANumb][6]);
        MakSpV(0,DoAtNu,HyAtN2,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
               *ReHydN,*ReVeCo);
      }

      else if ((LiToNu[i]==4)&&(LiHyNu[i]==1)) {
        DoAtNu=BSAtLi[i];
        for (j=1;j<=4;j++) {
          if (ReAtEl_nu[LiAtom[i][j]]==1) HyAtNu=LiAtom[i][j];
        }

        DANumb=DANumb+1;
        (*ReDATy)[DANumb]=0;
        (*ReDAAt)[DANumb]=DoAtNu;
        (*ReHydN)[DANumb]=HyAtNu;
        for (k=1;k<=3;k++)
          (*ReVeCo)[DANumb][k]=ReCoor[DoAtNu][k];
        PoCoVe(ReCoor[DoAtNu][1],ReCoor[DoAtNu][2],ReCoor[DoAtNu][3],
               ReCoor[HyAtNu][1],ReCoor[HyAtNu][2],ReCoor[HyAtNu][3],
               1.0,&(*ReVeCo)[DANumb][4],&(*ReVeCo)[DANumb][5],
               &(*ReVeCo)[DANumb][6]);
        MakSpV(0,DoAtNu,HyAtNu,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
               *ReHydN,*ReVeCo);
      }

      else if ((LiToNu[i]==4)&&(LiHyNu[i]==2)) {
        DoAtNu=BSAtLi[i];
        HyAtNu=0;
        for (j=1;j<=4;j++) {
          if (ReAtEl_nu[LiAtom[i][j]]==1) {
            if (!HyAtNu)
              HyAtNu=LiAtom[i][j];
            else
              HyAtN2=LiAtom[i][j];
          }
        }

        DANumb=DANumb+1;
        (*ReDATy)[DANumb]=0;
        (*ReDAAt)[DANumb]=DoAtNu;
        (*ReHydN)[DANumb]=HyAtNu;
        for (k=1;k<=3;k++)
          (*ReVeCo)[DANumb][k]=ReCoor[DoAtNu][k];
        PoCoVe(ReCoor[DoAtNu][1],ReCoor[DoAtNu][2],ReCoor[DoAtNu][3],
               ReCoor[HyAtNu][1],ReCoor[HyAtNu][2],ReCoor[HyAtNu][3],
               1.0,&(*ReVeCo)[DANumb][4],&(*ReVeCo)[DANumb][5],
               &(*ReVeCo)[DANumb][6]);
        MakSpV(0,DoAtNu,HyAtNu,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
               *ReHydN,*ReVeCo);

        DANumb=DANumb+1;
        (*ReDATy)[DANumb]=0;
        (*ReDAAt)[DANumb]=DoAtNu;
        (*ReHydN)[DANumb]=HyAtN2;
        for (k=1;k<=3;k++)
          (*ReVeCo)[DANumb][k]=ReCoor[DoAtNu][k];
        PoCoVe(ReCoor[DoAtNu][1],ReCoor[DoAtNu][2],ReCoor[DoAtNu][3],
               ReCoor[HyAtN2][1],ReCoor[HyAtN2][2],ReCoor[HyAtN2][3],
               1.0,&(*ReVeCo)[DANumb][4],&(*ReVeCo)[DANumb][5],
               &(*ReVeCo)[DANumb][6]);
        MakSpV(0,DoAtNu,HyAtN2,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
               *ReHydN,*ReVeCo);
      }

      else if ((LiToNu[i]==4)&&(LiHyNu[i]==3)) {
        DoAtNu=BSAtLi[i];
        HyAtNu=0;
        HyAtN2=0;
        for (j=1;j<=4;j++) {
          if (ReAtEl_nu[LiAtom[i][j]]==1) {
            if (!HyAtNu)
              HyAtNu=LiAtom[i][j];
            else if (!HyAtN2)
              HyAtN2=LiAtom[i][j];
            else
              HyAtN3=LiAtom[i][j];
          }
        }

        DANumb=DANumb+1;
        (*ReDATy)[DANumb]=0;
        (*ReDAAt)[DANumb]=DoAtNu;
        (*ReHydN)[DANumb]=HyAtNu;
        for (k=1;k<=3;k++)
          (*ReVeCo)[DANumb][k]=ReCoor[DoAtNu][k];
        PoCoVe(ReCoor[DoAtNu][1],ReCoor[DoAtNu][2],ReCoor[DoAtNu][3],
               ReCoor[HyAtNu][1],ReCoor[HyAtNu][2],ReCoor[HyAtNu][3],
               1.0,&(*ReVeCo)[DANumb][4],&(*ReVeCo)[DANumb][5],
               &(*ReVeCo)[DANumb][6]);
        MakSpV(0,DoAtNu,HyAtNu,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
               *ReHydN,*ReVeCo);

        DANumb=DANumb+1;
        (*ReDATy)[DANumb]=0;
        (*ReDAAt)[DANumb]=DoAtNu;
        (*ReHydN)[DANumb]=HyAtN2;
        for (k=1;k<=3;k++)
          (*ReVeCo)[DANumb][k]=ReCoor[DoAtNu][k];
        PoCoVe(ReCoor[DoAtNu][1],ReCoor[DoAtNu][2],ReCoor[DoAtNu][3],
               ReCoor[HyAtN2][1],ReCoor[HyAtN2][2],ReCoor[HyAtN2][3],
               1.0,&(*ReVeCo)[DANumb][4],&(*ReVeCo)[DANumb][5],
               &(*ReVeCo)[DANumb][6]);
        MakSpV(0,DoAtNu,HyAtN2,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
               *ReHydN,*ReVeCo);

        DANumb=DANumb+1;
        (*ReDATy)[DANumb]=0;
        (*ReDAAt)[DANumb]=DoAtNu;
        (*ReHydN)[DANumb]=HyAtN3;
        for (k=1;k<=3;k++)
          (*ReVeCo)[DANumb][k]=ReCoor[DoAtNu][k];
        PoCoVe(ReCoor[DoAtNu][1],ReCoor[DoAtNu][2],ReCoor[DoAtNu][3],
               ReCoor[HyAtN3][1],ReCoor[HyAtN3][2],ReCoor[HyAtN3][3],
               1.0,&(*ReVeCo)[DANumb][4],&(*ReVeCo)[DANumb][5],
               &(*ReVeCo)[DANumb][6]);
        MakSpV(0,DoAtNu,HyAtN3,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
               *ReHydN,*ReVeCo);
      }

    }

/* Case where the receptor binding site atom is a O */
    else if (ReAtEl_nu[BSAtLi[i]]==8) {

      if ((LiToNu[i]==1)&&(LiHyNu[i]==0)) {
        AcAtNu=BSAtLi[i];
        HelpA1=LiAtom[i][1];
        HelpA2=0;
        for (j=1;((j<=ReBdNu)&&(!HelpA2));j++) {
          if ((ReBdAr[j][1]==HelpA1)&&(ReBdAr[j][2]!=BSAtLi[i]))
            HelpA2=ReBdAr[j][2];
          else if ((ReBdAr[j][2]==HelpA1)&&(ReBdAr[j][1]!=BSAtLi[i]))
            HelpA2=ReBdAr[j][1];
        }

        if (!HelpA2) {
          fprintf(FPaOut,"WARNING Was not able to find a non-O linked ");
          fprintf(FPaOut,"atom for atom %d\n\n",HelpA1);
          fclose(FPaOut);
          FPaOut=fopen(OutFil,"a");
        }

        DANumb=DANumb+1;
        (*ReDATy)[DANumb]=1;
        (*ReDAAt)[DANumb]=AcAtNu;
        (*ReHydN)[DANumb]=0;
        for (k=4;k<=6;k++)
          (*ReVeCo)[DANumb][k]=ReCoor[AcAtNu][k-3];
        RotPla(ReCoor[HelpA2][1],ReCoor[HelpA2][2],ReCoor[HelpA2][3],
               ReCoor[AcAtNu][1],ReCoor[AcAtNu][2],ReCoor[AcAtNu][3],
               ReCoor[HelpA1][1],ReCoor[HelpA1][2],ReCoor[HelpA1][3],
               135*PiT180,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
               &(*ReVeCo)[DANumb][3]);
        PoCoVe(ReCoor[AcAtNu][1],ReCoor[AcAtNu][2],ReCoor[AcAtNu][3],
               (*ReVeCo)[DANumb][1],(*ReVeCo)[DANumb][2],(*ReVeCo)[DANumb][3],
               1.0,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
               &(*ReVeCo)[DANumb][3]);
        MakSpV(1,AcAtNu,0,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
               *ReHydN,*ReVeCo);

        DANumb=DANumb+1;
        (*ReDATy)[DANumb]=1;
        (*ReDAAt)[DANumb]=AcAtNu;
        (*ReHydN)[DANumb]=0;
        for (k=4;k<=6;k++)
          (*ReVeCo)[DANumb][k]=ReCoor[AcAtNu][k-3];
        RotPla(ReCoor[HelpA2][1],ReCoor[HelpA2][2],ReCoor[HelpA2][3],
               ReCoor[AcAtNu][1],ReCoor[AcAtNu][2],ReCoor[AcAtNu][3],
               ReCoor[HelpA1][1],ReCoor[HelpA1][2],ReCoor[HelpA1][3],
               225*PiT180,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
               &(*ReVeCo)[DANumb][3]);
        PoCoVe(ReCoor[AcAtNu][1],ReCoor[AcAtNu][2],ReCoor[AcAtNu][3],
               (*ReVeCo)[DANumb][1],(*ReVeCo)[DANumb][2],(*ReVeCo)[DANumb][3],
               1.0,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
               &(*ReVeCo)[DANumb][3]);
        MakSpV(1,AcAtNu,0,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
               *ReHydN,*ReVeCo);
      }

      else if ((LiToNu[i]==2)&&(LiHyNu[i]==0)) {
        AcAtNu=BSAtLi[i];
        HelpA1=LiAtom[i][1];
        HelpA2=LiAtom[i][2];
        if ((HybReAt[HelpA1]==3)||(HybReAt[HelpA2]==3)) {

          DANumb=DANumb+1;
          (*ReDATy)[DANumb]=1;
          (*ReDAAt)[DANumb]=AcAtNu;
          (*ReHydN)[DANumb]=0;
          for (k=4;k<=6;k++)
            (*ReVeCo)[DANumb][k]=ReCoor[AcAtNu][k-3];
          RoArVe(ReCoor[HelpA2][1]-ReCoor[AcAtNu][1],ReCoor[HelpA2][2]-
                 ReCoor[AcAtNu][2],ReCoor[HelpA2][3]-ReCoor[AcAtNu][3],
                 ReCoor[HelpA1][1]-ReCoor[AcAtNu][1],ReCoor[HelpA1][2]-
                 ReCoor[AcAtNu][2],ReCoor[HelpA1][3]-ReCoor[AcAtNu][3],
                 120*PiT180,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
                 &(*ReVeCo)[DANumb][3]);
          (*ReVeCo)[DANumb][1]=(*ReVeCo)[DANumb][1]+ReCoor[AcAtNu][1];
          (*ReVeCo)[DANumb][2]=(*ReVeCo)[DANumb][2]+ReCoor[AcAtNu][2];
          (*ReVeCo)[DANumb][3]=(*ReVeCo)[DANumb][3]+ReCoor[AcAtNu][3];
          PoCoVe(ReCoor[AcAtNu][1],ReCoor[AcAtNu][2],ReCoor[AcAtNu][3],
                 (*ReVeCo)[DANumb][1],(*ReVeCo)[DANumb][2],(*ReVeCo)[DANumb][3],
                 1.0,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
                 &(*ReVeCo)[DANumb][3]);
          MakSpV(1,AcAtNu,0,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
                 *ReHydN,*ReVeCo);

          DANumb=DANumb+1;
          (*ReDATy)[DANumb]=1;
          (*ReDAAt)[DANumb]=AcAtNu;
          (*ReHydN)[DANumb]=0;
          for (k=4;k<=6;k++)
            (*ReVeCo)[DANumb][k]=ReCoor[AcAtNu][k-3];
          RoArVe(ReCoor[HelpA2][1]-ReCoor[AcAtNu][1],ReCoor[HelpA2][2]-
                 ReCoor[AcAtNu][2],ReCoor[HelpA2][3]-ReCoor[AcAtNu][3],
                 ReCoor[HelpA1][1]-ReCoor[AcAtNu][1],ReCoor[HelpA1][2]-
                 ReCoor[AcAtNu][2],ReCoor[HelpA1][3]-ReCoor[AcAtNu][3],
                 240*PiT180,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
                 &(*ReVeCo)[DANumb][3]);
          (*ReVeCo)[DANumb][1]=(*ReVeCo)[DANumb][1]+ReCoor[AcAtNu][1];
          (*ReVeCo)[DANumb][2]=(*ReVeCo)[DANumb][2]+ReCoor[AcAtNu][2];
          (*ReVeCo)[DANumb][3]=(*ReVeCo)[DANumb][3]+ReCoor[AcAtNu][3];
          PoCoVe(ReCoor[AcAtNu][1],ReCoor[AcAtNu][2],ReCoor[AcAtNu][3],
                 (*ReVeCo)[DANumb][1],(*ReVeCo)[DANumb][2],(*ReVeCo)[DANumb][3],
                 1.0,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
                 &(*ReVeCo)[DANumb][3]);
          MakSpV(1,AcAtNu,0,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
                 *ReHydN,*ReVeCo);
        }
        else {
          DANumb=DANumb+1;
          (*ReDATy)[DANumb]=1;
          (*ReDAAt)[DANumb]=AcAtNu;
          (*ReHydN)[DANumb]=0;
          for (k=4;k<=6;k++)
            (*ReVeCo)[DANumb][k]=ReCoor[AcAtNu][k-3];
          Angl=PlaAng(ReCoor[AcAtNu][1],ReCoor[AcAtNu][2],ReCoor[AcAtNu][3],
                      ReCoor[HelpA1][1],ReCoor[HelpA1][2],ReCoor[HelpA1][3],
                      ReCoor[HelpA2][1],ReCoor[HelpA2][2],ReCoor[HelpA2][3]);
          RotPla(ReCoor[HelpA2][1],ReCoor[HelpA2][2],ReCoor[HelpA2][3],
                 ReCoor[AcAtNu][1],ReCoor[AcAtNu][2],ReCoor[AcAtNu][3],
                 ReCoor[HelpA1][1],ReCoor[HelpA1][2],ReCoor[HelpA1][3],
                 /*(PiT180*180*2-Angl)/2 clangini*/
                 (M_PI*2-Angl)/2,&(*ReVeCo)[DANumb][1],
                 &(*ReVeCo)[DANumb][2],&(*ReVeCo)[DANumb][3]);
          PoCoVe(ReCoor[AcAtNu][1],ReCoor[AcAtNu][2],ReCoor[AcAtNu][3],
                 (*ReVeCo)[DANumb][1],(*ReVeCo)[DANumb][2],(*ReVeCo)[DANumb][3],
                 1.0,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
                 &(*ReVeCo)[DANumb][3]);
          MakSpV(1,AcAtNu,0,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
                 *ReHydN,*ReVeCo);
        }
      }

      else if ((LiToNu[i]==2)&&(LiHyNu[i]==1)) {
        DoAtNu=BSAtLi[i];
        AcAtNu=BSAtLi[i];
        if (ReAtEl_nu[LiAtom[i][1]]==1) {
          HyAtNu=LiAtom[i][1];
          HelpA1=LiAtom[i][2];
        }
        else {
          HyAtNu=LiAtom[i][2];
          HelpA1=LiAtom[i][1];
        }
        if (HybReAt[HelpA1]==3) {

          DANumb=DANumb+1;
          (*ReDATy)[DANumb]=0;
          (*ReDAAt)[DANumb]=DoAtNu;
          (*ReHydN)[DANumb]=HyAtNu;
          for (k=1;k<=3;k++)
            (*ReVeCo)[DANumb][k]=ReCoor[DoAtNu][k];
          PoCoVe(ReCoor[DoAtNu][1],ReCoor[DoAtNu][2],ReCoor[DoAtNu][3],
                 ReCoor[HyAtNu][1],ReCoor[HyAtNu][2],ReCoor[HyAtNu][3],
                 1.0,&(*ReVeCo)[DANumb][4],&(*ReVeCo)[DANumb][5],
                 &(*ReVeCo)[DANumb][6]);
          MakSpV(0,DoAtNu,HyAtNu,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
                 *ReHydN,*ReVeCo);

          DANumb=DANumb+1;
          (*ReDATy)[DANumb]=1;
          (*ReDAAt)[DANumb]=AcAtNu;
          (*ReHydN)[DANumb]=0;
          for (k=4;k<=6;k++)
            (*ReVeCo)[DANumb][k]=ReCoor[AcAtNu][k-3];
          RoArVe(ReCoor[HyAtNu][1]-ReCoor[AcAtNu][1],ReCoor[HyAtNu][2]-
                 ReCoor[AcAtNu][2],ReCoor[HyAtNu][3]-ReCoor[AcAtNu][3],
                 ReCoor[HelpA1][1]-ReCoor[AcAtNu][1],ReCoor[HelpA1][2]-
                 ReCoor[AcAtNu][2],ReCoor[HelpA1][3]-ReCoor[AcAtNu][3],
                 120*PiT180,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
                 &(*ReVeCo)[DANumb][3]);
          (*ReVeCo)[DANumb][1]=(*ReVeCo)[DANumb][1]+ReCoor[AcAtNu][1];
          (*ReVeCo)[DANumb][2]=(*ReVeCo)[DANumb][2]+ReCoor[AcAtNu][2];
          (*ReVeCo)[DANumb][3]=(*ReVeCo)[DANumb][3]+ReCoor[AcAtNu][3];
          PoCoVe(ReCoor[AcAtNu][1],ReCoor[AcAtNu][2],ReCoor[AcAtNu][3],
                 (*ReVeCo)[DANumb][1],(*ReVeCo)[DANumb][2],(*ReVeCo)[DANumb][3],
                 1.0,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
                 &(*ReVeCo)[DANumb][3]);
          MakSpV(1,AcAtNu,0,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
                 *ReHydN,*ReVeCo);

          DANumb=DANumb+1;
          (*ReDATy)[DANumb]=1;
          (*ReDAAt)[DANumb]=AcAtNu;
          (*ReHydN)[DANumb]=0;
          for (k=4;k<=6;k++)
            (*ReVeCo)[DANumb][k]=ReCoor[AcAtNu][k-3];
          RoArVe(ReCoor[HyAtNu][1]-ReCoor[AcAtNu][1],ReCoor[HyAtNu][2]-
                 ReCoor[AcAtNu][2],ReCoor[HyAtNu][3]-ReCoor[AcAtNu][3],
                 ReCoor[HelpA1][1]-ReCoor[AcAtNu][1],ReCoor[HelpA1][2]-
                 ReCoor[AcAtNu][2],ReCoor[HelpA1][3]-ReCoor[AcAtNu][3],
                 240*PiT180,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
                 &(*ReVeCo)[DANumb][3]);
          (*ReVeCo)[DANumb][1]=(*ReVeCo)[DANumb][1]+ReCoor[AcAtNu][1];
          (*ReVeCo)[DANumb][2]=(*ReVeCo)[DANumb][2]+ReCoor[AcAtNu][2];
          (*ReVeCo)[DANumb][3]=(*ReVeCo)[DANumb][3]+ReCoor[AcAtNu][3];
          PoCoVe(ReCoor[AcAtNu][1],ReCoor[AcAtNu][2],ReCoor[AcAtNu][3],
                 (*ReVeCo)[DANumb][1],(*ReVeCo)[DANumb][2],(*ReVeCo)[DANumb][3],
                 1.0,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
                 &(*ReVeCo)[DANumb][3]);
          MakSpV(1,AcAtNu,0,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
                 *ReHydN,*ReVeCo);
        }
        else {

          DANumb=DANumb+1;
          (*ReDATy)[DANumb]=1;
          (*ReDAAt)[DANumb]=AcAtNu;
          (*ReHydN)[DANumb]=0;
          for (k=4;k<=6;k++)
            (*ReVeCo)[DANumb][k]=ReCoor[AcAtNu][k-3];
          Angl=PlaAng(ReCoor[AcAtNu][1],ReCoor[AcAtNu][2],ReCoor[AcAtNu][3],
                      ReCoor[HelpA1][1],ReCoor[HelpA1][2],ReCoor[HelpA1][3],
                      ReCoor[HyAtNu][1],ReCoor[HyAtNu][2],ReCoor[HyAtNu][3]);
          RotPla(ReCoor[HyAtNu][1],ReCoor[HyAtNu][2],ReCoor[HyAtNu][3],
                 ReCoor[AcAtNu][1],ReCoor[AcAtNu][2],ReCoor[AcAtNu][3],
                 ReCoor[HelpA1][1],ReCoor[HelpA1][2],ReCoor[HelpA1][3],
                 /*(PiT180*180*2-Angl)/2 clangini*/
                 (M_PI*2-Angl)/2,&(*ReVeCo)[DANumb][1],
                 &(*ReVeCo)[DANumb][2],&(*ReVeCo)[DANumb][3]);
          PoCoVe(ReCoor[AcAtNu][1],ReCoor[AcAtNu][2],ReCoor[AcAtNu][3],
                 (*ReVeCo)[DANumb][1],(*ReVeCo)[DANumb][2],(*ReVeCo)[DANumb][3],
                 1.0,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
                 &(*ReVeCo)[DANumb][3]);
          MakSpV(1,AcAtNu,0,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
                 *ReHydN,*ReVeCo);

          DANumb=DANumb+1;
          (*ReDATy)[DANumb]=0;
          (*ReDAAt)[DANumb]=DoAtNu;
          (*ReHydN)[DANumb]=HyAtNu;
          for (k=1;k<=3;k++)
            (*ReVeCo)[DANumb][k]=ReCoor[DoAtNu][k];
          PoCoVe(ReCoor[DoAtNu][1],ReCoor[DoAtNu][2],ReCoor[DoAtNu][3],
                 ReCoor[HyAtNu][1],ReCoor[HyAtNu][2],ReCoor[HyAtNu][3],
                 1.0,&(*ReVeCo)[DANumb][4],&(*ReVeCo)[DANumb][5],
                 &(*ReVeCo)[DANumb][6]);
          MakSpV(0,DoAtNu,HyAtNu,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
                 *ReHydN,*ReVeCo);
        }
      }

      else if ((LiToNu[i]==2)&&(LiHyNu[i]==2)) {
        DoAtNu=BSAtLi[i];
        AcAtNu=BSAtLi[i];
        HyAtNu=LiAtom[i][1];
        HyAtN2=LiAtom[i][2];
        HelpA1=LiAtom[i][2];

        DANumb=DANumb+1;
        (*ReDATy)[DANumb]=0;
        (*ReDAAt)[DANumb]=DoAtNu;
        (*ReHydN)[DANumb]=HyAtNu;
        for (k=1;k<=3;k++)
          (*ReVeCo)[DANumb][k]=ReCoor[DoAtNu][k];
        PoCoVe(ReCoor[DoAtNu][1],ReCoor[DoAtNu][2],ReCoor[DoAtNu][3],
               ReCoor[HyAtNu][1],ReCoor[HyAtNu][2],ReCoor[HyAtNu][3],
               1.0,&(*ReVeCo)[DANumb][4],&(*ReVeCo)[DANumb][5],
               &(*ReVeCo)[DANumb][6]);
        MakSpV(0,DoAtNu,HyAtNu,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
               *ReHydN,*ReVeCo);

        DANumb=DANumb+1;
        (*ReDATy)[DANumb]=0;
        (*ReDAAt)[DANumb]=DoAtNu;
        (*ReHydN)[DANumb]=HyAtN2;
        for (k=1;k<=3;k++)
          (*ReVeCo)[DANumb][k]=ReCoor[DoAtNu][k];
        PoCoVe(ReCoor[DoAtNu][1],ReCoor[DoAtNu][2],ReCoor[DoAtNu][3],
               ReCoor[HyAtN2][1],ReCoor[HyAtN2][2],ReCoor[HyAtN2][3],
               1.0,&(*ReVeCo)[DANumb][4],&(*ReVeCo)[DANumb][5],
               &(*ReVeCo)[DANumb][6]);
        MakSpV(0,DoAtNu,HyAtN2,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
               *ReHydN,*ReVeCo);

        DANumb=DANumb+1;
        (*ReDATy)[DANumb]=1;
        (*ReDAAt)[DANumb]=AcAtNu;
        (*ReHydN)[DANumb]=0;
        for (k=4;k<=6;k++)
          (*ReVeCo)[DANumb][k]=ReCoor[AcAtNu][k-3];
        RoArVe(ReCoor[HyAtNu][1]-ReCoor[AcAtNu][1],ReCoor[HyAtNu][2]-
               ReCoor[AcAtNu][2],ReCoor[HyAtNu][3]-ReCoor[AcAtNu][3],
               ReCoor[HelpA1][1]-ReCoor[AcAtNu][1],ReCoor[HelpA1][2]-
               ReCoor[AcAtNu][2],ReCoor[HelpA1][3]-ReCoor[AcAtNu][3],
               120*PiT180,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
               &(*ReVeCo)[DANumb][3]);
        (*ReVeCo)[DANumb][1]=(*ReVeCo)[DANumb][1]+ReCoor[AcAtNu][1];
        (*ReVeCo)[DANumb][2]=(*ReVeCo)[DANumb][2]+ReCoor[AcAtNu][2];
        (*ReVeCo)[DANumb][3]=(*ReVeCo)[DANumb][3]+ReCoor[AcAtNu][3];
        PoCoVe(ReCoor[AcAtNu][1],ReCoor[AcAtNu][2],ReCoor[AcAtNu][3],
               (*ReVeCo)[DANumb][1],(*ReVeCo)[DANumb][2],(*ReVeCo)[DANumb][3],
               1.0,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
               &(*ReVeCo)[DANumb][3]);
        MakSpV(1,AcAtNu,0,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
               *ReHydN,*ReVeCo);

        DANumb=DANumb+1;
        (*ReDATy)[DANumb]=1;
        (*ReDAAt)[DANumb]=AcAtNu;
        (*ReHydN)[DANumb]=0;
        for (k=4;k<=6;k++)
          (*ReVeCo)[DANumb][k]=ReCoor[AcAtNu][k-3];
        RoArVe(ReCoor[HyAtNu][1]-ReCoor[AcAtNu][1],ReCoor[HyAtNu][2]-
               ReCoor[AcAtNu][2],ReCoor[HyAtNu][3]-ReCoor[AcAtNu][3],
               ReCoor[HelpA1][1]-ReCoor[AcAtNu][1],ReCoor[HelpA1][2]-
               ReCoor[AcAtNu][2],ReCoor[HelpA1][3]-ReCoor[AcAtNu][3],
               240*PiT180,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
               &(*ReVeCo)[DANumb][3]);
        (*ReVeCo)[DANumb][1]=(*ReVeCo)[DANumb][1]+ReCoor[AcAtNu][1];
        (*ReVeCo)[DANumb][2]=(*ReVeCo)[DANumb][2]+ReCoor[AcAtNu][2];
        (*ReVeCo)[DANumb][3]=(*ReVeCo)[DANumb][3]+ReCoor[AcAtNu][3];
        PoCoVe(ReCoor[AcAtNu][1],ReCoor[AcAtNu][2],ReCoor[AcAtNu][3],
               (*ReVeCo)[DANumb][1],(*ReVeCo)[DANumb][2],(*ReVeCo)[DANumb][3],
               1.0,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
               &(*ReVeCo)[DANumb][3]);
        MakSpV(1,AcAtNu,0,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
               *ReHydN,*ReVeCo);
      }

    }

/* Case where the receptor binding site atom is a S */
    else if (ReAtEl_nu[BSAtLi[i]]==16) {

      if ((LiToNu[i]==2)&&(LiHyNu[i]==1)) {
        DoAtNu=BSAtLi[i];
        AcAtNu=BSAtLi[i];
        if (ReAtEl_nu[LiAtom[i][1]]==1) {
          HyAtNu=LiAtom[i][1];
          HelpA1=LiAtom[i][2];
        }
        else {
          HyAtNu=LiAtom[i][2];
          HelpA1=LiAtom[i][1];
        }
        if (HybReAt[HelpA1]==3) {

          DANumb=DANumb+1;
          (*ReDATy)[DANumb]=0;
          (*ReDAAt)[DANumb]=DoAtNu;
          (*ReHydN)[DANumb]=HyAtNu;
          for (k=1;k<=3;k++)
            (*ReVeCo)[DANumb][k]=ReCoor[DoAtNu][k];
          PoCoVe(ReCoor[DoAtNu][1],ReCoor[DoAtNu][2],ReCoor[DoAtNu][3],
                 ReCoor[HyAtNu][1],ReCoor[HyAtNu][2],ReCoor[HyAtNu][3],
                 1.0,&(*ReVeCo)[DANumb][4],&(*ReVeCo)[DANumb][5],
                 &(*ReVeCo)[DANumb][6]);
          MakSpV(0,DoAtNu,HyAtNu,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
                 *ReHydN,*ReVeCo);

          DANumb=DANumb+1;
          (*ReDATy)[DANumb]=1;
          (*ReDAAt)[DANumb]=AcAtNu;
          (*ReHydN)[DANumb]=0;
          for (k=4;k<=6;k++)
            (*ReVeCo)[DANumb][k]=ReCoor[AcAtNu][k-3];
          RoArVe(ReCoor[HyAtNu][1]-ReCoor[AcAtNu][1],ReCoor[HyAtNu][2]-
                 ReCoor[AcAtNu][2],ReCoor[HyAtNu][3]-ReCoor[AcAtNu][3],
                 ReCoor[HelpA1][1]-ReCoor[AcAtNu][1],ReCoor[HelpA1][2]-
                 ReCoor[AcAtNu][2],ReCoor[HelpA1][3]-ReCoor[AcAtNu][3],
                 120*PiT180,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
                 &(*ReVeCo)[DANumb][3]);
          (*ReVeCo)[DANumb][1]=(*ReVeCo)[DANumb][1]+ReCoor[AcAtNu][1];
          (*ReVeCo)[DANumb][2]=(*ReVeCo)[DANumb][2]+ReCoor[AcAtNu][2];
          (*ReVeCo)[DANumb][3]=(*ReVeCo)[DANumb][3]+ReCoor[AcAtNu][3];
          PoCoVe(ReCoor[AcAtNu][1],ReCoor[AcAtNu][2],ReCoor[AcAtNu][3],
                 (*ReVeCo)[DANumb][1],(*ReVeCo)[DANumb][2],(*ReVeCo)[DANumb][3],
                 1.0,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
                 &(*ReVeCo)[DANumb][3]);
          MakSpV(1,AcAtNu,0,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
                 *ReHydN,*ReVeCo);

          DANumb=DANumb+1;
          (*ReDATy)[DANumb]=1;
          (*ReDAAt)[DANumb]=AcAtNu;
          (*ReHydN)[DANumb]=0;
          for (k=4;k<=6;k++)
            (*ReVeCo)[DANumb][k]=ReCoor[AcAtNu][k-3];
          RoArVe(ReCoor[HyAtNu][1]-ReCoor[AcAtNu][1],ReCoor[HyAtNu][2]-
                 ReCoor[AcAtNu][2],ReCoor[HyAtNu][3]-ReCoor[AcAtNu][3],
                 ReCoor[HelpA1][1]-ReCoor[AcAtNu][1],ReCoor[HelpA1][2]-
                 ReCoor[AcAtNu][2],ReCoor[HelpA1][3]-ReCoor[AcAtNu][3],
                 240*PiT180,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
                 &(*ReVeCo)[DANumb][3]);
          (*ReVeCo)[DANumb][1]=(*ReVeCo)[DANumb][1]+ReCoor[AcAtNu][1];
          (*ReVeCo)[DANumb][2]=(*ReVeCo)[DANumb][2]+ReCoor[AcAtNu][2];
          (*ReVeCo)[DANumb][3]=(*ReVeCo)[DANumb][3]+ReCoor[AcAtNu][3];
          PoCoVe(ReCoor[AcAtNu][1],ReCoor[AcAtNu][2],ReCoor[AcAtNu][3],
                 (*ReVeCo)[DANumb][1],(*ReVeCo)[DANumb][2],(*ReVeCo)[DANumb][3],
                 1.0,&(*ReVeCo)[DANumb][1],&(*ReVeCo)[DANumb][2],
                 &(*ReVeCo)[DANumb][3]);
          MakSpV(1,AcAtNu,0,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,
                 *ReHydN,*ReVeCo);
        }
      }

    }

/* Case where the receptor binding site atom is a metal */
    else if ((ReAtEl_nu[BSAtLi[i]]==25)||(ReAtEl_nu[BSAtLi[i]]==30)||
             (ReAtEl_nu[BSAtLi[i]]==13)||(ReAtEl_nu[BSAtLi[i]]==20)||
             (ReAtEl_nu[BSAtLi[i]]==29)||(ReAtEl_nu[BSAtLi[i]]==26)||
             (ReAtEl_nu[BSAtLi[i]]==19)||(ReAtEl_nu[BSAtLi[i]]==12)||
             (ReAtEl_nu[BSAtLi[i]]==14)) {

      for (j=1;j<=BSMeNu;j++) {
        if (BSMeAN[j]==BSAtLi[i]) {
          DoAtNu=BSAtLi[i];
          DANumb=DANumb+1;
          (*ReDATy)[DANumb]=0;
          (*ReDAAt)[DANumb]=DoAtNu;
          (*ReHydN)[DANumb]=0;
	  diff = 0;
          for (k=1;k<=3;k++) {
            (*ReVeCo)[DANumb][k]=ReCoor[DoAtNu][k];
	    diff += fabs(ReCoor[DoAtNu][k] - BSMeVE[j][k]);
	  }
	  if(diff<1e-6) /* Coordinate site defined same as metal coordinates */
	    {
	      fprintf(stderr,"WARNING coordinates of metal coordination site for atom %d are the same as metal coordinates!",BSMeAN[j]);
	      fprintf(stderr,"Cannot determine polar vectors, exiting...\n");

	      exit(13);
	    }
          PoCoVe(ReCoor[DoAtNu][1],ReCoor[DoAtNu][2],ReCoor[DoAtNu][3],
                 BSMeVE[j][1],BSMeVE[j][2],BSMeVE[j][3],1.0,
                 &(*ReVeCo)[DANumb][4],&(*ReVeCo)[DANumb][5],
                 &(*ReVeCo)[DANumb][6]);
          MakSpV(0,DoAtNu,0,SphAng,SphPoN,&DANumb,*ReDATy,*ReDAAt,*ReHydN,
                 *ReVeCo);
        }
      }

    }

  }

  if ((*ReAcNu+*ReDoNu)!=DANumb) {
    fprintf(FPaOut,"WARNING The expected total number of vectors is not ");
    fprintf(FPaOut,"satisfied\n\n");
  }


  #ifdef ENABLE_MPI
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  if(myrank == MASTERRANK){
  #endif
  FilePa=fopen("./outputs/polar_rec.mol2","w");
  fprintf(FilePa,"# TRIPOS MOL2 file generated by WITNOTP\n");
  fprintf(FilePa,"\n");
  fprintf(FilePa,"@<TRIPOS>MOLECULE\n");
  fprintf(FilePa,"polar_rec\n");
  fprintf(FilePa,"%d 0 0 0 0\n",DANumb);
  fprintf(FilePa,"****\n");
  fprintf(FilePa,"USER_CHARGES\n");
  fprintf(FilePa,"INVALID_CHARGES\n");
  fprintf(FilePa,"\n");
  fprintf(FilePa,"@<TRIPOS>ATOM\n");

  for (i=1;i<=DANumb;i++) {
/*
      fprintf(FilePa,"%-5d %-3d%11.4f%11.4f%11.4f ***** 1    ****  0.000\n",
              i*2-1,i*2-1,(*ReVeCo)[i][1],(*ReVeCo)[i][2],(*ReVeCo)[i][3]);
      fprintf(FilePa,"%-5d %-3d%11.4f%11.4f%11.4f ***** 1    ****  0.000\n",
              i*2,i*2,(*ReVeCo)[i][4],(*ReVeCo)[i][5],(*ReVeCo)[i][6]);
*/
    if ((*ReDATy)[i]==0)
      fprintf(FilePa,"%-5d %-3d%11.4f%11.4f%11.4f ***** 1    ****  0.000\n",
              i,i,(*ReVeCo)[i][4],(*ReVeCo)[i][5],(*ReVeCo)[i][6]);
    else
      fprintf(FilePa,"%-5d %-3d%11.4f%11.4f%11.4f ***** 1    ****  0.000\n",
              i,i,(*ReVeCo)[i][1],(*ReVeCo)[i][2],(*ReVeCo)[i][3]);
  }
  fclose(FilePa);
  #ifdef ENABLE_MPI
  }
  #endif



  free_ivector(LiToNu,1,BSAtNu);
  free_ivector(LiHyNu,1,BSAtNu);
  free_imatrix(LiAtom,1,BSAtNu,1,10);

}
