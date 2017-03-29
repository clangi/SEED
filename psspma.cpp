#include <stdio.h>

void PsSpMa(double PsSpRa,double GrSiCu_en,int PsSpNC,int ***PsSphe)
/* This function constructs the pseudo-sphere for the evaluation of the
   energy. Here, the center is supposed to be (0,0,0) :
  PsSpNC  number of cubes along the radius of the pseudo-sphere for the
          evaluation of the energy
  PsSphe  pseudo-sphere in a cube (1 <-> inner part of the sphere,
                                   0 <-> outer part of the sphere)
  DistVa  squared distance from the origin
  SqRoot3  square root of 3 */
{
  int i,j,k;
  double DistVa,SqRoot3;

  SqRoot3=1.7320508;

  for (i=-PsSpNC;i<=PsSpNC;i++) {
    for (j=-PsSpNC;j<=PsSpNC;j++) {
      for (k=-PsSpNC;k<=PsSpNC;k++) {

        DistVa=GrSiCu_en*GrSiCu_en*(i*i+j*j+k*k);

        if (DistVa<=((PsSpRa+SqRoot3*GrSiCu_en)*(PsSpRa+SqRoot3*GrSiCu_en)))
          PsSphe[i+PsSpNC+1][j+PsSpNC+1][k+PsSpNC+1]=1;
        else
          PsSphe[i+PsSpNC+1][j+PsSpNC+1][k+PsSpNC+1]=0;

      }
    }
  }

}
