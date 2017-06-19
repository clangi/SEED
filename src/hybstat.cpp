#include <stdio.h>
#include <math.h>
#include "funct.h"

void HybStAt(int xxAtNu,int *xxAtEl_nu,double **xxCoor,int xxBdNu,int **xxBdAr,
             int *HybxxAt,FILE *FPaOut)
/* This function finds the hybridization state of the carbon atoms :
   HybxxAt  hybridization state (2 -> sp2, 3 -> sp3, 0 -> no hybridization
                                 state found)
   NumbLA  number of linked atoms
   LiAtAr  array of linked atoms
   Distan  distance between two atoms
   AngPla  planar angle
   TriPro  triple scalar product */
{
  int i,j,NumbLA,LiAtAr[10];
  double Distan,AngPla,TriPro,PiT180;

  //PiT180=3.1415927/180; clangini
  PiT180=M_PI/180; // clangini

/* Initialization */
  for (i=1;i<=xxAtNu;i++)
    HybxxAt[i]=0;

/* Loop over the atoms */
  for (i=1;i<=xxAtNu;i++) {

/* One tries to find the hybridization state only for the carbon atoms */
    if (xxAtEl_nu[i]==6) {

/* Find and count the number of atoms which are linked to the current atom */
      NumbLA=0;
      for (j=1;j<=xxBdNu;j++) {
        if (xxBdAr[j][1]==i) {
          NumbLA=NumbLA+1;
          LiAtAr[NumbLA]=xxBdAr[j][2];
        }
        else if (xxBdAr[j][2]==i) {
          NumbLA=NumbLA+1;
          LiAtAr[NumbLA]=xxBdAr[j][1];
        }
      }

/* Find the hybridization state */
      if (NumbLA==4) {
        HybxxAt[i]=3;
      }
      else if (NumbLA==3) {
        TriPro=TrProd(xxCoor[i][1],xxCoor[i][2],xxCoor[i][3],
                xxCoor[LiAtAr[1]][1],xxCoor[LiAtAr[1]][2],xxCoor[LiAtAr[1]][3],
                xxCoor[LiAtAr[2]][1],xxCoor[LiAtAr[2]][2],xxCoor[LiAtAr[2]][3],
                xxCoor[LiAtAr[3]][1],xxCoor[LiAtAr[3]][2],xxCoor[LiAtAr[3]][3]);
        if (TriPro<0.3)
          HybxxAt[i]=2;
        else
          HybxxAt[i]=3;
      }
      else if (NumbLA==2) {
        AngPla=PlaAng(xxCoor[i][1],xxCoor[i][2],xxCoor[i][3],
                xxCoor[LiAtAr[1]][1],xxCoor[LiAtAr[1]][2],xxCoor[LiAtAr[1]][3],
                xxCoor[LiAtAr[2]][1],xxCoor[LiAtAr[2]][2],xxCoor[LiAtAr[2]][3]);
        if ((AngPla>(115*PiT180))&&(AngPla<(125*PiT180)))
          HybxxAt[i]=2;
        else
          HybxxAt[i]=3;
      }
      else if (NumbLA==1) {
        Distan=DistSq(xxCoor[i][1],xxCoor[i][2],xxCoor[i][3],
                xxCoor[LiAtAr[1]][1],xxCoor[LiAtAr[1]][2],xxCoor[LiAtAr[1]][3]);
        Distan=sqrtf(Distan);
        if (xxAtEl_nu[LiAtAr[1]]==6) {
          if (Distan>1.4)
            HybxxAt[i]=3;
          else
            HybxxAt[i]=2;
        }
        else if (xxAtEl_nu[LiAtAr[1]]==8) {
          if (Distan>1.305)
            HybxxAt[i]=3;
          else
            HybxxAt[i]=2;
        }
        else if (xxAtEl_nu[LiAtAr[1]]==7) {
          if (Distan>1.375)
            HybxxAt[i]=3;
          else
            HybxxAt[i]=2;
        }
        else if (xxAtEl_nu[LiAtAr[1]]==16) {
          if (Distan>1.685)
            HybxxAt[i]=3;
          else
            HybxxAt[i]=2;
        }
      }

    }

  }

/* Write the hybridization state for the carbon atoms in the output file */
/*
  fprintf(FPaOut,"Hybridization state for the carbon atoms ");
  fprintf(FPaOut,"(2->sp2,3->sp3,0->not found) :\n");
  for (i=1;i<=xxAtNu;i++) {
    if (xxAtEl_nu[i]==6)
      fprintf(FPaOut,"%5d    %d\n",i,HybxxAt[i]);
  }
*/

/* Write if no hybridization state was found for some carbon atoms */
  for (i=1;i<=xxAtNu;i++) {
    if ((xxAtEl_nu[i]==6)&&(HybxxAt[i]==0)) {
      fprintf(FPaOut,"WARNING No hybridization state found for the ");
      fprintf(FPaOut,"atom %d\n",i);
    }
  }

  fprintf(FPaOut,"\n");

}
