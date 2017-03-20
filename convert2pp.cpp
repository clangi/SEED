#include <iostream>
#include "funct.h"

using namespace std;

float **convert2pp_return(int irl, int irh,int icl, int ich,float *first,
                            float *mat_pp_rows[]){
  int i,ncol,nrow;
  ncol = ich - icl + 1;
  nrow = irh - irl + 1;

  for(i=0;i<nrow;i++){
    mat_pp_rows[i] = first + ncol*i - 1;
  }
  return (mat_pp_rows - 1);
}

void print_pp(int irl, int irh,int icl, int ich,float **mat){
  int i,j;
  for (i = irl;i<=irh; i++){
    for (j = icl;j<=ich; j++){
      cout << mat[i][j] << " ";
    }
    cout << "\n";
  }
}
