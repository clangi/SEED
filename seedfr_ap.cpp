#include <stdio.h>
#include "funct.h"

void SeedFr_ap(int ReCuVe,double **apol_Vect_re,int *ReApAt,int FrCuVe,
               double **apol_Vect_fr,int *FrApAt,double **FrCoor,double *ReVdWR,
               double *FrVdWR,int FrAtNu,FILE *FPaOut,double **SeFrCo)
/* This function seeds the fragment in the binding site of the receptor in the
   apolar case :
   ReCuVe  receptor current vector
   FrCuVe  fragment current vector
   TranVe  translation vector
   RotAng  rotation angle
   RotaAx  rotation axis
   ApLeng  apolar bond length */
{
  int i;
  double TranVe[4],RotAng,RotaAx[4],ApLeng;

/* Translation given by the translation of the fragment vector origin towards
   the receptor vector origin */
  TranVe[1]=apol_Vect_re[ReCuVe][1]-apol_Vect_fr[FrCuVe][1];
  TranVe[2]=apol_Vect_re[ReCuVe][2]-apol_Vect_fr[FrCuVe][2];
  TranVe[3]=apol_Vect_re[ReCuVe][3]-apol_Vect_fr[FrCuVe][3];

  for (i=1;i<=FrAtNu;i++) {
    SeFrCo[i][1]=FrCoor[i][1]+TranVe[1];
    SeFrCo[i][2]=FrCoor[i][2]+TranVe[2];
    SeFrCo[i][3]=FrCoor[i][3]+TranVe[3];
  }

/* Rotation given by the rotation of the translated fragment vector extremity
   towards the receptor vector extremity. The rotation vector is perpendicular
   to the plane formed by the translated fragment and receptor vectors.
   If the angle is near PI or 0, the vector is built in another way. */

  RotAng=PlaAng(apol_Vect_re[ReCuVe][1],apol_Vect_re[ReCuVe][2],
                apol_Vect_re[ReCuVe][3],apol_Vect_re[ReCuVe][4],
                apol_Vect_re[ReCuVe][5],apol_Vect_re[ReCuVe][6],
                apol_Vect_fr[FrCuVe][4]+TranVe[1],
                apol_Vect_fr[FrCuVe][5]+TranVe[2],
                apol_Vect_fr[FrCuVe][6]+TranVe[3]);

  if ((RotAng>3.141)||(RotAng<0.001)) {
    RotaAx[1]=-(apol_Vect_re[ReCuVe][5]-apol_Vect_re[ReCuVe][2]);
    RotaAx[2]=apol_Vect_re[ReCuVe][4]-apol_Vect_re[ReCuVe][1];
    RotaAx[3]=0.0;
  }
  else {
    VectPr(apol_Vect_fr[FrCuVe][4]+TranVe[1]-apol_Vect_re[ReCuVe][1],
           apol_Vect_fr[FrCuVe][5]+TranVe[2]-apol_Vect_re[ReCuVe][2],
           apol_Vect_fr[FrCuVe][6]+TranVe[3]-apol_Vect_re[ReCuVe][3],
           apol_Vect_re[ReCuVe][4]-apol_Vect_re[ReCuVe][1],
           apol_Vect_re[ReCuVe][5]-apol_Vect_re[ReCuVe][2],
           apol_Vect_re[ReCuVe][6]-apol_Vect_re[ReCuVe][3],
           &RotaAx[1],&RotaAx[2],&RotaAx[3]);
  }

  for (i=1;i<=FrAtNu;i++) {
    RoArVe(SeFrCo[i][1]-apol_Vect_re[ReCuVe][1],
           SeFrCo[i][2]-apol_Vect_re[ReCuVe][2],
           SeFrCo[i][3]-apol_Vect_re[ReCuVe][3],RotaAx[1],RotaAx[2],RotaAx[3],
           RotAng,&SeFrCo[i][1],&SeFrCo[i][2],&SeFrCo[i][3]);
    SeFrCo[i][1]=SeFrCo[i][1]+apol_Vect_re[ReCuVe][1];
    SeFrCo[i][2]=SeFrCo[i][2]+apol_Vect_re[ReCuVe][2];
    SeFrCo[i][3]=SeFrCo[i][3]+apol_Vect_re[ReCuVe][3];
  }

/* Another translation in order to get the desired length between the concerned
   receptor and fragment atoms. One has to take into account that the vectors
   have already a length of 1.0 */

  TranVe[1]=apol_Vect_re[ReCuVe][4]-apol_Vect_re[ReCuVe][1];
  TranVe[2]=apol_Vect_re[ReCuVe][5]-apol_Vect_re[ReCuVe][2];
  TranVe[3]=apol_Vect_re[ReCuVe][6]-apol_Vect_re[ReCuVe][3];

  ApLeng=0.0;
  ApLeng=ReVdWR[ReApAt[ReCuVe]]+FrVdWR[FrApAt[FrCuVe]];

  if (ApLeng<1.0)
    fprintf(FPaOut,"WARNING Apolar bond length was less than 1.0\n\n");

  for (i=1;i<=FrAtNu;i++) {
    SeFrCo[i][1]=SeFrCo[i][1]+(ApLeng-1.0)*TranVe[1];
    SeFrCo[i][2]=SeFrCo[i][2]+(ApLeng-1.0)*TranVe[2];
    SeFrCo[i][3]=SeFrCo[i][3]+(ApLeng-1.0)*TranVe[3];
  }

}
