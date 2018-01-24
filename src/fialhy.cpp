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

void FiAlHy(int FrAtNu,int FrBdNu,int *FrAtEl_nu,int **FrBdAr,int *AliHyd,
            int *NonAliHy)
/* This function finds the aliphatic hydrogens (i.e. the hydrogens attached 
   to carbon) :
   AliHyd  0 for aliphatic hydrogens and 1 for other hydrogens and 
           other atoms 
   NonAliHy  number of non aliphatic-hydrogen atoms */
{
  int i,j;

  for (i=1;i<=FrAtNu;i++) {
    AliHyd[i]=1;
    if (FrAtEl_nu[i]==1) {
      for (j=1;j<=FrBdNu;j++) {
        if (((FrBdAr[j][1]==i) && (FrAtEl_nu[FrBdAr[j][2]]==6)) ||
            ((FrBdAr[j][2]==i) && (FrAtEl_nu[FrBdAr[j][1]]==6)))
          AliHyd[i]=0;
      }
    }  
  }

/* Count the number of non aliphatic-hydrogen atoms */
  *NonAliHy=0;
  for (i=1;i<=FrAtNu;i++) 
    *NonAliHy=*NonAliHy+AliHyd[i];

}
