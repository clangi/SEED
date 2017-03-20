#include <stdio.h>
#include "funct.h"

void SeedFr(int ReCuVe,float **ReVeCo,int FrCuVe,float **FrVeCo,int FrAtNu,
            float **FrCoor,int *ReDATy,float **SeFrCo,int *FrAtoTyp_nu,
            int *FrDAAt,int *ReAtoTyp_nu,int *ReDAAt,float **BLAtTy,
            FILE *FPaOut)
/* This function seeds the fragment in the binding site of the receptor :
   ReCuVe  receptor current vector
   FrCuVe  fragment current vector
   TranVe  translation vector
   RotAng  rotation angle
   RotaAx  rotation axis 
   HBLeng  hydrogen bond length */
{
  int i;
  float TranVe[4],RotAng,RotaAx[4],HBLeng;

/* Translation given by the translation of the fragment donor atom vector 
   towards the receptor donor atom vector */
  TranVe[1]=ReVeCo[ReCuVe][1]-FrVeCo[FrCuVe][1];
  TranVe[2]=ReVeCo[ReCuVe][2]-FrVeCo[FrCuVe][2];
  TranVe[3]=ReVeCo[ReCuVe][3]-FrVeCo[FrCuVe][3];

  for (i=1;i<=FrAtNu;i++) {
    SeFrCo[i][1]=FrCoor[i][1]+TranVe[1];
    SeFrCo[i][2]=FrCoor[i][2]+TranVe[2];
    SeFrCo[i][3]=FrCoor[i][3]+TranVe[3];
  }

/* Rotation given by the rotation of the translated fragment acceptor atom 
   vector towards the receptor acceptor atom vector. The rotation vector is 
   perpendicular to the plane formed by the translated fragment and receptor 
   vectors. If the angle is near PI or 0, the vector is built in another way. */

  RotAng=PlaAng(ReVeCo[ReCuVe][1],ReVeCo[ReCuVe][2],ReVeCo[ReCuVe][3],
                ReVeCo[ReCuVe][4],ReVeCo[ReCuVe][5],ReVeCo[ReCuVe][6],
                FrVeCo[FrCuVe][4]+TranVe[1],FrVeCo[FrCuVe][5]+TranVe[2],
                FrVeCo[FrCuVe][6]+TranVe[3]);

  if ((RotAng>3.141)||(RotAng<0.001)) {
    RotaAx[1]=-(ReVeCo[ReCuVe][5]-ReVeCo[ReCuVe][2]);
    RotaAx[2]=ReVeCo[ReCuVe][4]-ReVeCo[ReCuVe][1];
    RotaAx[3]=0.0;
  }
  else {
    VectPr(FrVeCo[FrCuVe][4]+TranVe[1]-ReVeCo[ReCuVe][1],
           FrVeCo[FrCuVe][5]+TranVe[2]-ReVeCo[ReCuVe][2],
           FrVeCo[FrCuVe][6]+TranVe[3]-ReVeCo[ReCuVe][3],
           ReVeCo[ReCuVe][4]-ReVeCo[ReCuVe][1],ReVeCo[ReCuVe][5]-
           ReVeCo[ReCuVe][2],ReVeCo[ReCuVe][6]-ReVeCo[ReCuVe][3],
           &RotaAx[1],&RotaAx[2],&RotaAx[3]);
  }

  for (i=1;i<=FrAtNu;i++) {
    RoArVe(SeFrCo[i][1]-ReVeCo[ReCuVe][1],SeFrCo[i][2]-ReVeCo[ReCuVe][2],
           SeFrCo[i][3]-ReVeCo[ReCuVe][3],RotaAx[1],RotaAx[2],RotaAx[3],
           RotAng,&SeFrCo[i][1],&SeFrCo[i][2],&SeFrCo[i][3]);
    SeFrCo[i][1]=SeFrCo[i][1]+ReVeCo[ReCuVe][1];
    SeFrCo[i][2]=SeFrCo[i][2]+ReVeCo[ReCuVe][2];
    SeFrCo[i][3]=SeFrCo[i][3]+ReVeCo[ReCuVe][3];
  }

/* Another translation in order to get the desired length between the receptor 
   donor (acceptor) atom and the fragment acceptor (donor) atom. One has to 
   take into account that the vectors have already a length of 1.0 */
  
  if (ReDATy[ReCuVe]==0) {
    TranVe[1]=ReVeCo[ReCuVe][4]-ReVeCo[ReCuVe][1];
    TranVe[2]=ReVeCo[ReCuVe][5]-ReVeCo[ReCuVe][2];
    TranVe[3]=ReVeCo[ReCuVe][6]-ReVeCo[ReCuVe][3];
  }
  else {
    TranVe[1]=ReVeCo[ReCuVe][1]-ReVeCo[ReCuVe][4];
    TranVe[2]=ReVeCo[ReCuVe][2]-ReVeCo[ReCuVe][5];
    TranVe[3]=ReVeCo[ReCuVe][3]-ReVeCo[ReCuVe][6];
  }

  HBLeng=BLAtTy[FrAtoTyp_nu[FrDAAt[FrCuVe]]][ReAtoTyp_nu[ReDAAt[ReCuVe]]];

  for (i=1;i<=FrAtNu;i++) {
    SeFrCo[i][1]=SeFrCo[i][1]+(HBLeng-1.0)*TranVe[1];
    SeFrCo[i][2]=SeFrCo[i][2]+(HBLeng-1.0)*TranVe[2];
    SeFrCo[i][3]=SeFrCo[i][3]+(HBLeng-1.0)*TranVe[3];
  }

}
