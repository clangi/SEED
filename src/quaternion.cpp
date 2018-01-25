/*
*    This file is part of SEED.
*
*    Copyright (C) 2017, Caflisch Lab, University of Zurich
*
*    SEED is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    SEED is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
//#include "quaternion.h"
#include "math.h"
#include "rnd_namespace.h"

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
  T my_norm = this->norm();
  w = w/my_norm;
  x = x/my_norm;
  y = y/my_norm;
  z = z/my_norm;
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
void Quaternion<T>::quatConjugateVecRef(V *v,V ref1,V ref2,V ref3){
  //WARNING: it modifies the original vector v!
  Quaternion<T> qv(0.0, v[1]-ref1,v[2]-ref2,v[3]-ref3);
  Quaternion<T> qm = (*this)*qv*conj();
  v[1] = qm.x + ref1;
  v[2] = qm.y + ref2;
  v[3] = qm.z + ref3;
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
  const T halfAngle = 0.5*angle;
  const T sin_half = sin(halfAngle);
  w = cos(halfAngle);
  x = axis[1]*sin_half;
  y = axis[2]*sin_half;
  z = axis[3]*sin_half;
  return (*this);
}
template<class T>
const Quaternion<T>& Quaternion<T>::fromAngleAxis(T angle, T ax1, T ax2, T ax3){
  const T halfAngle = 0.5*angle;
  const T sin_half = sin(halfAngle);
  w = cos(halfAngle);
  x = ax1*sin_half;
  y = ax2*sin_half;
  z = ax3*sin_half;
  return (*this);
}

template<class T>
void Quaternion<T>::get_AngleAxis(point& axis, T& theta) {
  T sin_theta_half = sqrt(1 - w*w);
  theta = acos(w)*2;
  axis.x = x/sin_theta_half;
  axis.y = y/sin_theta_half;
  axis.z = z/sin_theta_half;
}


template<class T>
Quaternion<T> Quaternion<T>::unitRandom() {
  T x0, x1, x2, theta1, theta2, r1, r2;
  x0 = rnd_gen::get_uniform(0, 1);
  x1 = rnd_gen::get_uniform(0, 1);
  x2 = rnd_gen::get_uniform(0, 1);

  theta1 = 2*M_PI*x1;
  theta2 = 2*M_PI*x2;

  r1 = sqrt(1 - x0);
  r2 = sqrt(x0);


  return Quaternion(cos(theta2) * r2, sin(theta1) * r1,
                    cos(theta1) * r1, sin(theta2) * r2);

}

template<class T>
Quaternion<T> Quaternion<T>::unitRandom(T mc_rot_step) {
  T theta, u, u1;
  Quaternion<T> qt;
  theta = rnd_gen::get_uniform(0, 2*M_PI);
  u = rnd_gen::get_uniform(-1, 1);

  u1 = sqrt(1 - u*u);

  qt.fromAngleAxis(mc_rot_step, u1 * cos(theta), u1 * sin(theta), u);
  return qt;
}
