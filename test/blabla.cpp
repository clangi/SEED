#include "quaternion.h"
#include "math.h"

template <class T>
Quaternion<T>::Quaternion(){
  w = 0.0;
  x = 0.0;
  y = 0.0;
  z = 0.0;
}

template <class T>
Quaternion<T>::Quaternion(T w_i, T x_i, T y_i, T z_i){
  w = w_i;
  x = x_i;
  y = y_i;
  z = z_i;
}

template <class T>
Quaternion<T>::Quaternion(T *v){
  w = v[1];
  x = v[2];
  y = v[3];
  z = v[4];
}

template <class T>
Quaternion<T>::Quaternion(const Quaternion<T>& q){//Copy constructor
  w = q.w;
  x = q.x;
  y = q.y;
  z = q.z;
}

template <class T>
Quaternion<T>::~Quaternion(){
}
// Overloaded operators
template <class T>
const Quaternion<T>& Quaternion<T>::operator=(const Quaternion<T>& q){
  //assignment operator
  //if (this != &q){
    w = q.w;
    x = q.x;
    y = q.y;
    z = q.z;
  //}
  return *this;
}

template <class T>
const Quaternion<T>& Quaternion<T>::operator+=(const Quaternion<T>& q){
  w += q.w;
  x += q.x;
  y += q.y;
  z += q.z;
  return *this;
}
template <class T>
const Quaternion<T>& Quaternion<T>::operator-=(const Quaternion<T>& q){
  w -= q.w;
  x -= q.x;
  y -= q.y;
  z -= q.z;
  return (*this);
}
template <class T>
Quaternion<T> Quaternion<T>::operator*(const Quaternion<T>& q) const {
  return Quaternion(
    w*q.w - x*q.x - y*q.y - z*q.z,
    w*q.x + x*q.w + y*q.z - z*q.y,
    w*q.y + y*q.w + z*q.x - x*q.z,
    w*q.z + z*q.w + x*q.y - y*q.x);
}
template <class T>
Quaternion<T> Quaternion<T>::operator/(const T s){
  if (s == 0){
    std::cerr<<"Trying to divide quaternion by 0!"<<std::endl;
  }
  return Quaternion(w/s,x/s,y/s,z/s);
}
template <class T>
Quaternion<T> Quaternion<T>::operator/(Quaternion<T>& q){
  return ((*this)*q.inverse());
}

template <class T>
bool Quaternion<T>::operator==(const Quaternion<T>& q) const {
  return (w==q.w && x==q.x && y==q.y && z==q.z);
}
template <class T>
bool Quaternion<T>::operator!=(const Quaternion<T>& q) const {
  return !((*this) == q);
}
// template <class T>
// Quaternion<T> operator+(const Quaternion<T>& lhs,const Quaternion<T>& rhs){
//   return Quaternion<T>(lhs) += rhs;
// }
//
//template <class T>
//std::ostream& operator<<(std::ostream& output, const Quaternion<T>& q){
//  output << "[" << q.w << ", " << "(" << q.x << ", " << q.y << ", " << q.z << ")]";
//  return output;
//}
template<class T>
void Quaternion<T>::print_quat(void) const {
  std::cout<< w <<", "<< x <<", "<< y <<", "<< z <<std::endl;
}

template <class T>
Quaternion<T> Quaternion<T>::conjugate(void){
  return Quaternion(w, -x, -y, -z);
}
template <class T>
T Quaternion<T>::norm2(void){ //Squared norm
  return (w*w + x*x + y*y + z*z);
}
template <class T>
T Quaternion<T>::norm(void){
  return sqrt(norm2());
}
template <class T>
Quaternion<T> Quaternion<T>::inverse(void){
  return conjugate()/norm2();
}

template<class T>
Quaternion<T> Quaternion<T>::normalize(void){
  return (*this)/(this->norm());
}

template<class T>
template<typename V>
void Quaternion<T>::quatRotateVec(V *v){
  //WARNING: it modifies the original vector v!
  Quaternion<T> qv(0, v[1],v[2],v[3]);
  Quaternion<T> qm = (*this)*qv*conjugate();
  v[1] = qm.x;
  v[2] = qm.y;
  v[3] = qm.z;
}
#include <iostream>
//#include "../nrutil.h"
#include "quaternion.h"
// main for testing

int main(int argc,char *argv[]){
  float aaa[5];
  for (int i=0; i < 5; i++){
    aaa[i] = i;
  }
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
  std::cout<< c.conjugate() <<std::endl;
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

  std::cout<< aaa[1]<<", "<<aaa[2]<<", "<<aaa[3]<<", "<<std::endl;
  c.normalize().quatRotateVec(aaa);
  std::cout<<aaa[1]<<", "<<aaa[2]<<", "<<aaa[3]<<", "<<std::endl;
  return 1;
}
