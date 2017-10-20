#include <iostream>
#include <vector>
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

// main for testing
int main(int argc,char *argv[]){

  print_summary();
  Quaternion<float> a;
  a.print_quat();
  print_summary();
  Quaternion<float> b(2.3,3.0,4.0,5.0);
  b.print_quat();
  print_summary();
  //Quaternion<float> c(1.0,0.0,0.0,0.0);
  a = Quaternion<float>::randomQuaternion(10);
  a.print_quat();
  print_summary();
  Quaternion<float> c = Quaternion<float>::randomQuaternion(7);
  // Quaternion<float> c;
  // c = Quaternion<float>::randomQuaternion(10);
  c.print_quat();
  print_summary();

  cout << "\n========================\n" << endl;
  Quaternion<float> q1(1,2,3,4);
  Quaternion<float> q2(1,2,3,4);
  q1.print_quat();
  q1.normalize().print_quat();
  q2 = q1.normalize();
  q1.print_quat();
  q2.print_quat();

  q1.norm_inplace();
  q1.print_quat();
}
