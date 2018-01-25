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

/* Eigenvalues routines */
#ifndef _EIGEN_SEED_H
#define _EIGEN_SEED_H

#define MAXITERJAC 50

#include "nrutil.h"
#include <limits>
#include <cmath>

/* Class for Jacobi diagonalization adapted from Numerical Recipes in C++ */
void rotate(double ** const a,const double s,const double tau,const int i,
  const int j, const int k, const int l);

class Jaco
{
  private:
    int n;
    double **a,**v;
    double *eigvals;
    int rot;
    const static double EPSILON;

    void deallocate(void);
    void eigen_sort(void);

  public:
    Jaco(double **aa, int a_size);
    ~Jaco();
    const double * eigenvalues(void);
    const double ** eigenvectors(void);
    //double ** diagonalized(void);
};

/* Adapted from ALGOL code in
   [The Jacobi Method for Real Symmetric Matrices, H. Rutishauser, 1971] */
Jaco::Jaco(double **aa, int a_size){
  int i, j, p, q; //indexes
  double sm, c, s, t, h, g, tau, theta, tresh; // sm = sum of off diagonal moduli
  double *b, *z;
  // Member initialization:
  n = a_size;
  a = dmatrix(1,n,1,n);
  v = dmatrix(1,n,1,n);
  eigvals = dvector(1,n);

  // inizialization
  for(p = 1; p <= n; p++){
    for(q = 1; q <= n; q++){
      a[p][q] = aa[p][q]; // copy aa not to modify the original data
      v[p][q] = 0.0;
    }
    v[p][p] = 1.0;
  }
  rot = 0;
  // Initialization of supporting variables:
  b = dvector(1,n);
  z = dvector(1,n);
  for(p = 1; p <= n; p++){
    b[p] = eigvals[p] = a[p][p];
    z[p] = 0.0;
  }
  // swp:
  for (i=1; i <= MAXITERJAC; i++) {
    sm = 0.0;
    for(p = 1;p <= (n-1); p++) {
      for (q = p+1; q <= n; q++){
        sm += std::abs(a[p][q]);
      }
    }
    // exit condition:
    if (sm == 0.0) {
      free_dvector(b, 1, n);
      free_dvector(z, 1, n);
      eigen_sort(); // sort eigenvalues
      return;
    }
    if (i < 4)
      tresh = 0.2 * sm/(n*n);
    else
      tresh = 0.0;
    for (p = 1; p <= (n-1); p++) {
      for (q = (p + 1); q <= n; q++) {
        g = 100.0 * std::abs(a[p][q]); //Only rotates if diagonal element is small
        if (i > 4 && g <= EPSILON * std::abs(eigvals[p]) &&
            g <= EPSILON * std::abs(eigvals[q])) // As in NR in C++
          a[p][q] = 0.0; // set elem to zero and skip the iteration
        else if (std::abs(a[p][q]) > tresh) { // rotate
          h = eigvals[q] - eigvals[p];
          if (g <= EPSILON*std::abs(h))
            t = (a[p][q])/h;
          else { // computing tan of rotation angle
            theta = 0.5 * h/a[p][q];
            t = 1.0/(std::abs(theta) + std::sqrt(1.0 + theta*theta));
            if (theta < 0.0)
              t = -t;
          }
          c = 1.0 / std::sqrt(1 + t*t);
          s = t * c; // ~angle
          tau = s / (1.0 + c);

          h = t * a[p][q];
          z[p] -= h;
          z[q] += h;
          eigvals[p] -= h;
          eigvals[q] += h;
          // note the a is diagonal so we evaluate
          // elements only in the upper triangular half
          a[p][q] = 0.0;
          for (j = 1; j <= p-1; j++)
            rotate(a, s, tau, j, p, j, q);
          for(j = p+1; j <= q-1; j++)
            rotate(a, s, tau, p, j, j, q);
          for (j = q + 1; j <= n; j++)
            rotate(a, s, tau, p, j, q, j);
          for (j = 1; j <= n; j++)
            rotate(v, s, tau, j, p, j, q);
          ++rot;
        } // end rotate
      }
    }
    for (p = 1; p <= n; p++) {
      b[p] += z[p];
      eigvals[p] = b[p];
      z[p] = 0.0;
    }
  }
  throw("Maximum number of iterations exceeded in Jacobi diagonalization");
}

Jaco::~Jaco(){
  deallocate();
}

void Jaco::deallocate(void){
  free_dmatrix(a,1,n,1,n);
  free_dmatrix(v,1,n,1,n);
  free_dvector(eigvals,1,n);
}

inline void rotate(double ** const mat, const double angle, const double tau,
  const int i, const int j, const int k, const int l){
    double m_ij = mat[i][j]; // temporaries
    double m_kl = mat[k][l];
    mat[i][j] = m_ij - angle*(m_kl + m_ij*tau);
    mat[k][l] = m_kl + angle*(m_ij - m_kl*tau);
  }

// void Jaco::eigsrt(void){//classic implementation of insertion sort
//   int i,j,k;
//   double temp;
//   for (i = 2; i <= n; i++){
//     temp = d[i]
//     j = i - 1;
//     while((j > 0) && (d[j] < temp)){
//       d[j+1] = d[j];
//       j--;
//     }
//     d[j+1] = temp;
//     j++;
//     if (j!=i){
//       for (k=1; k <= n; k++){
//         temp = v[k][j];
//         v[k][j] = v[k][i];
//         v[k][i] = temp;
//       }
//     }
//   }
// }

void Jaco::eigen_sort(void){//straight insertion sort as in Numerical Recipes
  //Sort eigenvalues (and eigenvectors), from larger to smaller.
  int i,j,k;
  double temp;
  for(i = 1; i < n; i++){
    temp = eigvals[i];
    k = i;
    for(j = i+1; j <= n; j++){
      if (eigvals[j] >= temp){
        temp = eigvals[j];
        k = j;
      }
    }
    if (k != i){
      eigvals[k] = eigvals[i];
      eigvals[i] = temp;
      for (j=1;j<=n;j++){
        temp = v[j][i];
        v[j][i] = v[j][k];
        v[j][k] = temp;
      }
    }
  }
}

const double * Jaco::eigenvalues(void){
  const double * d_ptr = &eigvals[0];
  return d_ptr;
}
const double ** Jaco::eigenvectors(void){
  const double ** v_ptr = const_cast<const double **> (&v[0]);
  return v_ptr;
}
// double ** Jaco::diagonalized(void){
//     double ** a_ptr = &a[0];
//     return a_ptr;
// }

const double Jaco::EPSILON = std::numeric_limits<double>::epsilon();



#endif
