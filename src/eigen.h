/* Implementation of eigenvalues routines */
#ifndef _EIGEN_SEED_H
#define _EIGEN_SEED_H

#define MAXITERJAC 50

#include "nrutil.h"
#include <limits>
#include <cmath>

/* Jacobi diagonalization, as explained in "Numerical Recipes in C++" */
class Jacobi
{
  private:
    int n;
    double **a,**v;
    double *d;
    int nrot;
    const static double EPS;

    void deallocate(void);
    void rot(double ** const a,const double s,const double tau,const int i,
      const int j, const int k, const int l);
    void eigsrt(void);

  public:
    Jacobi(double **aa, int sz);
    ~Jacobi();
    const double * eigenvalues(void);
    const double ** eigenvectors(void);
    //double ** diagonalized(void);
};

Jacobi::Jacobi(double **aa, int sz){
  int i,j,ip,iq;
  double tresh, theta, tau, t, sm, s, h, g, c;
  double *b,*z;
  // Member initialization:
  n = sz;
  a = dmatrix(1,n,1,n);
  v = dmatrix(1,n,1,n);
  d = dvector(1,n);

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
  for(i=1;i<=n;i++){
    b[i] = d[i] = a[i][i];
    z[i] = 0.0;
  }
  // Iterations:
  for (i=1; i <= MAXITERJAC;i++) {
    sm=0.0;
    for(ip=1;ip <= (n-1); ip++) {
      for (iq = ip+1; iq <= n; iq++){
        sm += std::abs(a[ip][iq]);
      }
    }
    if (sm == 0.0) {
      free_dvector(b,1,n);
      free_dvector(z,1,n);
      eigsrt(); // implement this
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh = 0.0;
    for (ip=1;ip<=n-1;ip++) {
      for (iq=ip+1;iq<=n;iq++) {
        g=100.0*std::abs(a[ip][iq]); //Only rotates if diagonal element is small
        if (i > 4 && g <= EPS*std::abs(d[ip]) && g <= EPS*std::abs(d[iq]))
          a[ip][iq]=0.0;
        else if (std::abs(a[ip][iq]) > tresh) {
          h=d[iq]-d[ip];
          if (g <= EPS*std::abs(h))
            t=(a[ip][iq])/h;
          else {
            theta=0.5*h/(a[ip][iq]);
            t=1.0/(std::abs(theta)+std::sqrt(1.0+theta*theta));
            if (theta < 0.0)
              t = -t;
          }
          c=1.0/std::sqrt(1+t*t);
          s=t*c;
          tau=s/(1.0+c);
          h=t*a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[ip][iq]=0.0;
          for (j=1;j<ip;j++)
            rot(a,s,tau,j,ip,j,iq);
          for(j=ip+1;j<iq;j++)
            rot(a,s,tau,ip,j,j,iq);
          for (j=iq+1;j<=n;j++)
            rot(a,s,tau,ip,j,iq,j);
          for (j=1;j<=n;j++)
            rot(v,s,tau,j,ip,j,iq);
          ++nrot;
        }
      }
    }
    for (ip=1;ip<=n;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  throw("Too many iterations in routine jacobi");
}

Jacobi::~Jacobi(){
  deallocate();
}

void Jacobi::deallocate(void){
  free_dmatrix(a,1,n,1,n);
  free_dmatrix(v,1,n,1,n);
  free_dvector(d,1,n);
}

inline void Jacobi::rot(double ** const a, const double s, const double tau,
  const int i, const int j, const int k, const int l){
    double g = a[i][j];
    double h = a[k][l];
    a[i][j] = g - s*(h+g*tau);
    a[k][l] = h+s*(g-h*tau);
  }

// void Jacobi::eigsrt(void){//classic implementation of insertion sort
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

void Jacobi::eigsrt(void){//straight insertion sort as in Numerical Recipes
  //Sort eigenvalues (and eigenvectors), from larger to smaller.
  int i,j,k;
  double temp;
  for(i = 1; i < n; i++){
    temp = d[i];
    k = i;
    for(j = i+1; j <= n; j++){
      if (d[j] >= temp){
        temp = d[j];
        k = j;
      }
    }
    if (k != i){
      d[k] = d[i];
      d[i] = temp;
      for (j=1;j<=n;j++){
        temp = v[j][i];
        v[j][i] = v[j][k];
        v[j][k] = temp;
      }
    }
  }
}

const double * Jacobi::eigenvalues(void){
  const double * d_ptr = &d[0];
  return d_ptr;
}
const double ** Jacobi::eigenvectors(void){
  const double ** v_ptr = const_cast<const double **> (&v[0]);
  return v_ptr;
}
// double ** Jacobi::diagonalized(void){
//     double ** a_ptr = &a[0];
//     return a_ptr;
// }

const double Jacobi::EPS = std::numeric_limits<double>::epsilon();



#endif
