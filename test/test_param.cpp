#include <iostream>
#include "../src/Parameter.h"

using namespace std;

int main(int argc,char *argv[]){
  Parameter my_param;
  int a;
  float a_float;

  cout << "Temperature is: " << my_param.mc_temp << endl;
  cout << "a: " << a << " a_float: " << a_float << endl;
}
