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

/* Class for Jacobi diagonalization adapted from Numerical Recipes */
void jaco_rot(double ** const a,const double s,const double tau,const int i,
  const int j, const int k, const int l);

class Jaco
{
  private:
    int n;
    double **a,**v;
    double *eigvals;
    int nrot;
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

Jaco::Jaco(double **aa, int a_size){
  int i, j, ip, iq; //indexes
  double tresh, theta, tau, t, s, h, g, c;
  double sum_off; // sum of off diagonal moduli
  double *b, *z;
  // Member initialization:
  n = a_size;
  a = dmatrix(1,n,1,n);
  v = dmatrix(1,n,1,n);
  eigvals = dvector(1,n);

  // inizialization
  for(i = 1;i <= n;i++){
    for(j=1;j<=n;j++){
      a[i][j] = aa[i][j];
      v[i][j] = 0.0;
    }
    v[i][i] = 1.0;
  }
  nrot = 0;
  // Initialization of supporting variables:
  b = dvector(1,n);
  z = dvector(1,n);
  for(i = 1; i <= n; i++){
    b[i] = eigvals[i] = a[i][i];
    z[i] = 0.0;
  }
  // Iterations:
  for (i=1; i <= MAXITERJAC; i++) {
    sum_off = 0.0;
    for(ip = 1;ip <= (n-1); ip++) {
      for (iq = ip+1; iq <= n; iq++){
        sum_off += std::abs(a[ip][iq]);
      }
    }
    if (sum_off == 0.0) {
      free_dvector(b, 1, n);
      free_dvector(z, 1, n);
      eigen_sort(); // sort eigenvalues
      return;
    }
    if (i < 4)
      tresh = 0.2 * sum_off/(n*n);
    else
      tresh = 0.0;
    for (ip = 1; ip <= (n-1); ip++) {
      for (iq = (ip + 1); iq <= n; iq++) {
        g = 100.0 * std::abs(a[ip][iq]); //Only rotates if diagonal element is small
        if (i > 4 && g <= EPSILON * std::abs(eigvals[ip]) && g <= EPSILON * std::abs(eigvals[iq]))
          a[ip][iq] = 0.0; // skip the iteration
        else if (std::abs(a[ip][iq]) > tresh) {
          h = eigvals[iq] - eigvals[ip];
          if (g <= EPSILON*std::abs(h))
            t = (a[ip][iq])/h;
          else {
            theta = 0.5 * h/(a[ip][iq]);
            t = 1.0/(std::abs(theta) + std::sqrt(1.0 + theta*theta));
            if (theta < 0.0)
              t = -t;
          }
          c = 1.0 / std::sqrt(1 + t*t);
          s = t * c; // ~angle
          tau = s / (1.0 + c);

          h = t * a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          eigvals[ip] -= h;
          eigvals[iq] += h;
          // note the a is diagonal so we evaluate
          // elements only in the upper triangular half
          a[ip][iq] = 0.0;
          for (j = 1; j < ip; j++)
            jaco_rot(a, s, tau, j, ip, j, iq);
          for(j = ip+1; j < iq; j++)
            jaco_rot(a, s, tau, ip, j, j, iq);
          for (j = iq + 1; j <= n; j++)
            jaco_rot(a, s, tau, ip, j, iq, j);
          for (j = 1; j <= n; j++)
            jaco_rot(v, s, tau, j, ip, j, iq);
          ++nrot;
        }
      }
    }
    for (ip = 1; ip <= n; ip++) {
      b[ip] += z[ip];
      eigvals[ip] = b[ip];
      z[ip] = 0.0;
    }
  }
  throw("Too many iterations in routine jacobi");
}

Jaco::~Jaco(){
  deallocate();
}

void Jaco::deallocate(void){
  free_dmatrix(a,1,n,1,n);
  free_dmatrix(v,1,n,1,n);
  free_dvector(eigvals,1,n);
}

inline void jaco_rot(double ** const mat, const double angle, const double tau,
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
