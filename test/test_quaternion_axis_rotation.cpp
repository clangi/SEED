#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <math.h>
//#include "../nrutil.h"
void print_summary();
#include "quaternion_class_test.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif

using namespace std;
typedef vector< Quaternion<float> >::iterator itq;

void print_summary(){
  cout<<"Total number: " << Quaternion<float>::howMany() << endl;
  cout<<"Created: " << Quaternion<float>::howManyCreated() << endl;
  cout<<"Destroyed: " << Quaternion<float>::howManyDestroyed() << endl;
}

void normVec(double *v){
  double v1 = v[1];
  double v2 = v[2];
  double v3 = v[3];
  double norm = sqrt(v1*v1+v2*v2+v3*v3);
  v[1] = v1/norm;
  v[2] = v2/norm;
  v[3] = v3/norm;
}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

// main for testing
int main(int argc,char *argv[]){

  Quaternion<double> q;

  double axis[4], ref[4], point[4];
  for (int i=0; i<=3; i++){
    axis[i] = rand();
    ref[i] = fRand(0,10);
    point[i] = fRand(0,10);
  }
  cout<< "axis:\n" << axis[1] << "  " << axis[2] << "  "
      << axis[3] << endl;
  normVec(axis);
  cout << setprecision(20)
      << "Normalized axis:\n" << axis[1] << "  " << axis[2] << "  "
      << axis[3] << "  Norm:  "
      << sqrt(axis[1]*axis[1]+axis[2]*axis[2]+axis[3]*axis[3]) << endl;
  cout << setprecision(20)
      << "Reference:\n" << ref[1] << "  " << ref[2] << "  "
      << ref[3] << "  Norm:  "
      << sqrt(ref[1]*ref[1]+ref[2]*ref[2]+ref[3]*ref[3]) << endl;
  cout << setprecision(20)
      << "Point to rotate:\n" << point[1] << "  " << point[2] << "  "
      << point[3] << "  Norm:  "
      << sqrt(point[1]*point[1]+point[2]*point[2]+point[3]*point[3]) << endl;

  double angle = M_PI/180.0;
  int prec = 13;
  q.fromAngleAxis(angle, axis);
  cout << q << endl;

  int max_rotation = 360000;
  cout << left << setprecision(prec)
      << "Point: " << setw(24) << point[1] << setw(24) << point[2] << setw(24)
      << point[3] << endl;
  for (int i=1; i<= max_rotation; i++){
    q.quatConjugateVecRef<double>(point,ref[1],ref[2],ref[3]);
    if ((i % 360) == 0){
      cout << setprecision(prec)
          << "Point: " << setw(24) << point[1] << setw(24) << point[2] << setw(24)
          << point[3] << "  iter: " << i << endl;
    }
  }

  // double angle = M_PI*2;
  // int prec = 13;
  // q.fromAngleAxis(angle, axis);
  // cout << q << endl;
  //
  // int max_try = 100000;
  // cout << left << setprecision(prec)
  //     << "Point: " << setw(24) << point[1] << setw(24) << point[2] << setw(24)
  //     << point[3] << endl;
  // for (int i=1; i<= max_try; i++){
  //   q.quatConjugateVecRef<double>(point,ref[1],ref[2],ref[3]);
  //   if ((i % 360) == 0){
  //     cout << setprecision(prec)
  //         << "Point: " << setw(24) << point[1] << setw(24) << point[2] << setw(24)
  //         << point[3] << "  iter: " << i << endl;
  //   }
  // }

}
