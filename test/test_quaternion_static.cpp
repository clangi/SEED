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

  Quaternion<float> a;
  Quaternion<float> b(2.3,3.0,4.0,5.0);
  Quaternion<float> c(1.0,0.0,0.0,0.0);

  print_summary();
  a.print_quat();
  a = b;
  a.print_quat();

  for(int i=0; i < 10; i++){
    Quaternion<float> a(i,i,i,i); // This is constructed and destroyed at each
      //iteration. NOTE that the scope is hiding the quaternion a declared in
      //the main function
    a.print_quat();
  }

  {
    Quaternion<float> a(1,1,1,1);
  }

  a.print_quat();

  print_summary();
  a = Quaternion<float>(1.0,1.1,1.2,1.3); //This creates a temporary
      //which is created and destructed;
  a.print_quat();
  print_summary();

  a = Quaternion<float>(); //This creates a temporary
      //which is created and destructed;
  a.print_quat();
  print_summary();

  vector< Quaternion<float> > my_vector (100,Quaternion<float>());
  print_summary();

  vector< Quaternion<float> > my_vector2 (100);
  print_summary();

  vector< Quaternion<float> > my_vector3;
  print_summary();
  for(int i = 0; i < 100; i++){
    my_vector3.push_back(Quaternion<float>(1,2,3,4));//This is creating/destroying one
     //Quaternion at each iteration!
  }
  print_summary();
  my_vector[2].print_quat();
  my_vector2[5].print_quat();
  my_vector3[6].print_quat();

  for(itq i = my_vector.begin(); i != my_vector.end(); i++){
    *i = Quaternion<float>(5,6,7,8);
  }
  my_vector[10].print_quat();
  print_summary();

  for(itq i = my_vector2.begin(); i != my_vector2.end(); i++){
    (*i).Set(10,11,12,13);
  }
  my_vector2[10].print_quat();
  print_summary();

  float vec[3];
  vec[0] = 1.0; vec[1] = 3.5; vec[2] = 0.5;
  cout << "Original vec: " << vec[0] << " " << vec[1] << " "
       << vec[2] << " " << vec[3] << endl;
  Quaternion<float> q(0,1,2,3);
  q.print_quat();
  print_summary();
  q.norm_inplace();
  q.print_quat();
  cout << "Before rotation =================" << endl;
  print_summary();
  q.quatConjugateVec(vec-1);
  cout << "After rotation  ==================" << endl;
  cout << "Rot vec: " << vec[0] << " " << vec[1] << " "
       << vec[2] << " " << vec[3] << endl;
  q.print_quat();
  print_summary();

  Quaternion<float> qq(0,1,2,3);
  Quaternion<float> pp(qq);
  qq.print_quat();
  pp.print_quat();
  print_summary();
  qq += pp;
  qq.print_quat();
  print_summary();

  cout << "Set test" << endl;
  qq.Set(1,2,3,4).print_quat();
  //qq.print_quat();
  cout << "\nTest from angle-axis\n";
  float angle = PI;
  float axis[3] = {1.0,1.0,1.0};
  print_summary();
  qq.fromAngleAxis(angle,(axis-1));
  qq.print_quat();
  pp.fromAngleAxis(angle,axis[0],axis[1],axis[2]).print_quat();
  print_summary();
}
