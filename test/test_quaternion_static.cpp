#include <iostream>
//#include "../nrutil.h"
#include "quaternion_class_test.h"

using namespace std;

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

  print_summary();
  a = Quaternion<float>(1.0,1.1,1.2,1.3); //This creates a temporary
      //which is created and destructed;
  a.print_quat();
  print_summary();

  a = Quaternion<float>(); //This creates a temporary
      //which is created and destructed;
  a.print_quat();

  print_summary();
}
