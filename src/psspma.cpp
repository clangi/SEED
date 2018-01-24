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
