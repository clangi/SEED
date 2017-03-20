#include <iostream>
//#include "../nrutil.h"
#include "../quaternion.h"
// main for testing

int main(int argc,char *argv[]){
  float aaa[5], ref[4];
  for (int i=0; i < 5; i++){
    aaa[i] = i;
  }
  ref[1] = 1.3; 
  ref[2] = 2.4;
  ref[3] = 1.9;
  // Test of the different constructor
  Quaternion<float> a(0,1,0,0);
  Quaternion<float> b(0,0,1,0);
  Quaternion<float> c(2.3,3.0,4.0,5.0);
  Quaternion<float> d(b);

  std::cout<< "prova" << std::endl;
  std::cout<< a <<std::endl;
  std::cout<< b <<std::endl;
  std::cout<< c <<std::endl;
  std::cout<< d <<std::endl;

  //a.print_quat();
  //b = a;

  std::cout<<  a * b << std::endl;
  std::cout<<  -(b * a) << std::endl;
  std::cout<< "prova" << std::endl;
  std::cout<< c.conj() <<std::endl;
  std::cout<< c.norm() <<std::endl;
  std::cout<< c.norm2() << std::endl;
  std::cout<< 2*c <<std::endl;
  std::cout << c.inverse() << std::endl;
  std::cout<< d <<std::endl;
  std::cout<< c/2 <<std::endl;
  std::cout<< c/c <<std::endl;
  std::cout<< (d == b) << std::endl;
  std::cout<< (d != b) << std::endl;

  std::cout<< "prova" << std::endl;
  std::cout<< c <<std::endl;
  std::cout<< c.normalize() <<std::endl;
  std::cout<< c <<std::endl;
  //c = c.normalize();
  std::cout<< c <<std::endl;

  std::cout<< aaa[1]<<", "<<aaa[2]<<", "<<aaa[3]<<std::endl;
  c.normalize().quatConjugateVec(aaa);
  std::cout<<aaa[1]<<", "<<aaa[2]<<", "<<aaa[3]<<std::endl;
  for (int i=0; i < 5; i++){
    aaa[i] = i;
  }
  std::cout<<aaa[1]<<", "<<aaa[2]<<", "<<aaa[3]<<std::endl;
  c.normalize().quatConjugateVecRef(aaa,ref);
  std::cout<<aaa[1]<<", "<<aaa[2]<<", "<<aaa[3]<<std::endl;
  return 1;
}
