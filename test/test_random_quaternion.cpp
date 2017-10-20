#include <iostream>
#include <vector>
#include "math.h"
#include "../src/funct.h"
#include "../src/quaternion.h"


// #ifndef PI
// #define PI 3.14159265358979323846
// #endif


// void print_summary(){
//   cout<<"Total number: " << Quaternion<float>::howMany() << endl;
//   cout<<"Created: " << Quaternion<float>::howManyCreated() << endl;
//   cout<<"Destroyed: " << Quaternion<float>::howManyDestroyed() << endl;
// }

// template <class T>
// void get_angleAxis(Quaternion<T>& q, point& axis, double& angle){
//   angle = acos()
// }

using namespace std;
// main for testing
int main(int argc,char *argv[]){
  Quaternion<double> q1, q2, q3;
  point my_axis;
  double my_angle;
  double mc_step = 5.0;

  q1 = Quaternion<double>::unitRandom();
  // q1.print_quat();
  q1.get_AngleAxis(my_axis, my_angle);
  // cout << q1.norm() << "\n";
  // cout << "Axis: [" << my_axis.x << " " << my_axis.y << " " << my_axis.z << "]\n";
  // cout << "Angle: " << my_angle << " \n";

  // q2 = Quaternion<double>::unitRandom(mc_step);
  // q2.print_quat();
  // cout << q2.norm() << "\n";
  // q2.get_AngleAxis(my_axis, my_angle);
  // cout << "Axis: [" << my_axis.x << " " << my_axis.y << " " << my_axis.z << "]\n";
  // cout << "Angle: " << my_angle << " \n";
  for (int i = 0; i < 10000; i++){
    q1 = Quaternion<double>::unitRandom(mc_step);
    q1.get_AngleAxis(my_axis, my_angle);
    cout << my_axis.x << " " << my_axis.y << " " << my_axis.z << "\n";
  }
}
