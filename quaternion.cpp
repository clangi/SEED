//#include "quaternion.h"
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
Quaternion<T> Quaternion<T>::conj(void){
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
  return conj()/norm2();
}

template<class T>
Quaternion<T> Quaternion<T>::normalize(void){
  return (*this)/(this->norm());
}
template<class T>
void Quaternion<T>::norm_inplace(void){
  w = w/(this->norm());
  x = x/(this->norm());
  y = y/(this->norm());
  z = z/(this->norm());
}

template<class T>
template<typename V>
void Quaternion<T>::quatConjugateVec(V *v){
  //WARNING: it modifies the original vector v!
  Quaternion<T> qv(0.0, v[1],v[2],v[3]);
  Quaternion<T> qm = (*this)*qv*conj();
  v[1] = qm.x;
  v[2] = qm.y;
  v[3] = qm.z;
}
template<class T>
template<typename V>
void Quaternion<T>::quatConjugateVec(V *v1,V *v2,V *v3){//Overloaded
  Quaternion<T> qv(0.0, (*v1),(*v2),(*v3));
  Quaternion<T> qm = (*this)*qv*conj();
  *v1 = qm.x;
  *v2 = qm.y;
  *v3 = qm.z;
}
template<class T>
template<typename V>
void Quaternion<T>::quatConjugateVecRef(V *v,V *ref){
  //WARNING: it modifies the original vector v!
  Quaternion<T> qv(0.0, v[1]-ref[1],v[2]-ref[2],v[3]-ref[3]);
  Quaternion<T> qm = (*this)*qv*conj();
  v[1] = qm.x + ref[1];
  v[2] = qm.y + ref[2];
  v[3] = qm.z + ref[3];
}
template<class T>
template<typename V>
void Quaternion<T>::quatConjugateVec(V *v1,V *v2,V *v3,V *ref){//Overloaded
  Quaternion<T> qv(0.0, (*v1)-ref[1],(*v2)-ref[2],(*v3)-ref[3]);
  Quaternion<T> qm = (*this)*qv*conj();
  *v1 = qm.x + ref[1];
  *v2 = qm.y + ref[2];
  *v3 = qm.z + ref[3];
}

//Methods for setting quaternion components
template<class T>
const Quaternion<T>& Quaternion<T>::Set(T w_i, T x_i, T y_i, T z_i)
{
  w = w_i;
  x = x_i;
  y = y_i;
  z = z_i;
  return *this;
}

template<class T>
const Quaternion<T>& Quaternion<T>::fromAngleAxis(T angle, const T *axis){
  const float halfAngle = 0.5*angle;
  const float sin_half = sin(halfAngle);
  w = cos(halfAngle);
  x = axis[1]*sin_half;
  y = axis[2]*sin_half;
  z = axis[3]*sin_half;
  return (*this);
}
template<class T>
const Quaternion<T>& Quaternion<T>::fromAngleAxis(T angle, T ax1, T ax2, T ax3){
  const float halfAngle = 0.5*angle;
  const float sin_half = sin(halfAngle);
  w = cos(halfAngle);
  x = ax1*sin_half;
  y = ax2*sin_half;
  z = ax3*sin_half;
  return (*this);
}
