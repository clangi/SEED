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
#include <iostream>
#include "funct.h"

using namespace std;

double **convert2pp_return(int irl, int irh,int icl, int ich,double *first,
                            double *mat_pp_rows[]){
  int i,ncol,nrow;
  ncol = ich - icl + 1;
  nrow = irh - irl + 1;

  for(i=0;i<nrow;i++){
    mat_pp_rows[i] = first + ncol*i - 1;
  }
  return (mat_pp_rows - 1);
}

void print_pp(int irl, int irh,int icl, int ich,double **mat){
  int i,j;
  for (i = irl;i<=irh; i++){
    for (j = icl;j<=ich; j++){
      cout << mat[i][j] << " ";
    }
    cout << "\n";
  }
}
