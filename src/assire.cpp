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
#include <string.h>
#include <math.h>
#include "funct.h"

int AssiRE(int NumbAT,char **AtTyAr,double *VdWRad,double *VdWEne,int xxAtNu,
            char **xxAtTy,double *xxVdWR,double *xxVdWE,FILE *FPaOut,
            int *AtENAr,int *xxAtEl_nu,int *xxAtoTyp_nu)
/* This function assigns the van der Waals radii, energies and atom element
   numbers :
   xxAtNu  number of atoms
   xxAtTy  atoms types
   xxVdWR  van der Waals radii array
   xxVdWE  van der Waals energies array
   xxAtEl_nu  atom element numbers array
   xxAtoTyp_nu  atom type numbering array
   AssiTe  assignation test */
{
  int i,j,AssiTe,MissingTypes;
  MissingTypes = 0;

  for (i=1;i<=xxAtNu;i++) {
    //std::cout << "Atom type "<<i<<": "<<xxAtTy[i]<<std::endl;
    AssiTe=1;

    for (j=1;j<=NumbAT;j++) {
      if (!(strcmp(&xxAtTy[i][1],&AtTyAr[j][1]))) {
        AssiTe=0;
        xxVdWR[i]=VdWRad[j];
        xxVdWE[i]=VdWEne[j];
        xxAtEl_nu[i]=AtENAr[j];
        xxAtoTyp_nu[i]=j;
      }
    }

    if (AssiTe)
      {
	fprintf(FPaOut,"WARNING The following atom type is unknown : %s\n\n",&xxAtTy[i][1]);
	++MissingTypes;
      }
  }
  return (MissingTypes==0);
}



void SqRoEn(int xxAtNu,double *xxVdWE,double *xxVdWE_sr)
/* This function computes the square root of the van der Waals energies :
   xxAtNu     number of atoms
   xxVdWE     van der Waals energies array
   xxVdWE_sr  square root of the van der Waals energies array */
{
  int i;

  for (i=1;i<=xxAtNu;i++)
    xxVdWE_sr[i]=sqrtf(xxVdWE[i]);

}
