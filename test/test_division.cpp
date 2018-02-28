#include <iostream>
#include "../src/Parameter.h"

using namespace std;

int main(int argc,char *argv[]){
  float a, b;
  int result0, result1, result2, result3;
  int div1, div2;

  result0 = 2.5;
  result1 = 5.0/2.0;
  result2 = 5.0/1.0;
  result3 = 5.0/4.0;

  a = 2.3;
  b = 0.4;
  div1 = a/b;
  cout << result0 << " " << result1 << " " << result2 << " " << result3 << endl;
  cout << div1 << endl;
}
