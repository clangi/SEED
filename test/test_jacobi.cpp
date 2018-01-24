//#include "eigen.h"
#include <iostream>
#include <iomanip>
#include "../nrutil.h"
#include "../eigen.h"

using namespace std;

int main(int argc, char *argv[]){
  double **mymatrix;
  //const double * m_diag;
  int i,j;
  int n = 4;
  mymatrix = dmatrix(1,n,1,n);
  mymatrix[1][1] = 1.0;
  mymatrix[1][2] = 3.5;
  mymatrix[1][3] = 2.3;
  mymatrix[1][4] = 4.1;
  mymatrix[2][2] = 2.0;
  mymatrix[2][3] = 6.25;
  mymatrix[2][4] = 1.5;
  mymatrix[3][3] = 3.0;
  mymatrix[3][4] = 2.06;
  mymatrix[4][4] = 4.0;
  for(i=2;i<=n;i++){
    for(j=1;j < i; j++)
      mymatrix[i][j] = mymatrix[j][i];
  }
  //in R: aaa <- rbind(c(1,3.5,2.3,4.1),c(3.5,2,6.25,1.5),
  //                   c(2.3,6.25,3,2.06),c(4.1,1.5,2.06,4))
  cout << "Starting Matrix:\n";
  for(i=1;i<=n;i++){
    for(j=1;j<=n;j++){
      cout << setw(5) << mymatrix[i][j] << " ";
    }
    cout << "\n";
  }

  Jaco mymatrix_d(mymatrix,n);
  cout << "Eigenvector Matrix:\n";
  //m_diag = mymatrix_d.diagonalized();
  for(i=1;i<=n;i++){
    for(j=1;j<=n;j++){
      //cout << setw(5) << m_diag[i][j] << " ";
      cout << setw(10) << mymatrix_d.eigenvectors()[i][j] << " ";
    }
    cout << "\n";
  }
  cout << "Eigenvalues:\n";
  //m_diag = mymatrix_d.diagonalized();
  for(i=1;i<=n;i++){
    cout << setw(8) << mymatrix_d.eigenvalues()[i] << " ";
    cout << "\n";
  }

  free_dmatrix(mymatrix,1,n,1,n);
}
