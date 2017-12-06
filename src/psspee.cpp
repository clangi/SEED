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

void PsSpEE(int FrAtNu,int ReAtNu,double *ReVdWE_sr,double *FrVdWE_sr,
            double *ReVdWR,double *FrVdWR,double *VWEnEv_ps,double **SDFrRe_ps)
/* This function evaluates the energy using a pseudo-sphere approach as a
   cutoff :
   RFSqDi  squared distance between one atom of the receptor and one atom
           of the fragment
   SumRad  sum of the van der Waals radii */
{
  int i,j;
  double SumRad,SumRad_p6,RFSqDi,RFSqDi_p3;

  *VWEnEv_ps=0;

  for (i=1;i<=FrAtNu;i++) {

    for (j=1;j<=ReAtNu;j++) {

      RFSqDi=SDFrRe_ps[i][j];

      if (RFSqDi>=0.0) {

        if (RFSqDi>(1.e-4)) {

          RFSqDi_p3=RFSqDi*RFSqDi*RFSqDi;
          SumRad=ReVdWR[j]+FrVdWR[i];
          SumRad_p6=SumRad*SumRad*SumRad*SumRad*SumRad*SumRad;

          *VWEnEv_ps=*VWEnEv_ps+ReVdWE_sr[j]*FrVdWE_sr[i]*
                     SumRad_p6*((SumRad_p6/(RFSqDi_p3*RFSqDi_p3))-
                     (2/RFSqDi_p3));
        }

        else
          *VWEnEv_ps=*VWEnEv_ps+(1.e+12);

      }

    }

  }

}
