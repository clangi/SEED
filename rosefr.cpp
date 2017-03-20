#include <stdio.h>
#include "funct.h"

void RoSeFr(int ReCuVe,int *RexxAt,float **ReCoor,int FrCuVe,int *FrxxAt,
            int FrAtNu,float **SeFrCo,float AnglRo,float **RoSFCo)
/* This function rotates the seeded fragment around the axis built with 
   the concerned atom of the seeded fragment and the concerned atom of the 
   receptor with an angle AnglRo given in radian (in the polar case, the axis 
   is acceptor (or donor) of the seeded fragment - donor (or acceptor) of the 
   receptor) :
   ReCuVe  receptor current vector
   FrCuVe  fragment current vector 
   RotaVe  rotation vector */
{
  int i;
  float RotaVe[4];
 
  RotaVe[1]=ReCoor[RexxAt[ReCuVe]][1]-SeFrCo[FrxxAt[FrCuVe]][1];
  RotaVe[2]=ReCoor[RexxAt[ReCuVe]][2]-SeFrCo[FrxxAt[FrCuVe]][2];
  RotaVe[3]=ReCoor[RexxAt[ReCuVe]][3]-SeFrCo[FrxxAt[FrCuVe]][3];

  for (i=1;i<=FrAtNu;i++) {
    RoArVe(SeFrCo[i][1]-SeFrCo[FrxxAt[FrCuVe]][1],
           SeFrCo[i][2]-SeFrCo[FrxxAt[FrCuVe]][2],
           SeFrCo[i][3]-SeFrCo[FrxxAt[FrCuVe]][3],
           RotaVe[1],RotaVe[2],RotaVe[3],AnglRo,&RoSFCo[i][1],&RoSFCo[i][2],
           &RoSFCo[i][3]);
    RoSFCo[i][1]=RoSFCo[i][1]+SeFrCo[FrxxAt[FrCuVe]][1];
    RoSFCo[i][2]=RoSFCo[i][2]+SeFrCo[FrxxAt[FrCuVe]][2];
    RoSFCo[i][3]=RoSFCo[i][3]+SeFrCo[FrxxAt[FrCuVe]][3];
  }

}
