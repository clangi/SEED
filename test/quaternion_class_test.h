#ifndef _QUATERNION_TEST_H
#define _QUATERNION_TEST_H

#include <iostream>
#include "boost/random.hpp"
//#include <math.h>
//#include "nrutil.h"

class RndQuatExpr {
  public:
    int seed;
    RndQuatExpr(int seed_): seed(seed_) { }
};

template <class T>
class Quaternion {
  private:
    T w,x,y,z;
    //static members to count objects:
    static int created;
    static int destroyed;
    static int total_number; 
  public:
    // Constructors
    Quaternion(); // Default constructor: creates a zero quaternion.
    Quaternion(T w_i,T x_i,T y_i, T z_i);
    Quaternion(T *v);
    Quaternion(const Quaternion<T>& q); //copy constructor
    Quaternion(const RndQuatExpr& rndq);
    //default destructor
    ~Quaternion();
    //static methods:
    static int howManyCreated();
    static int howManyDestroyed();
    static int howMany();
    //Overloaded operators
    const Quaternion<T>& operator=(const Quaternion<T>& rhs);
    const Quaternion<T>& operator+=(const Quaternion<T>& q);
    const Quaternion<T>& operator-=(const Quaternion<T>& q);
    Quaternion<T> operator*(const Quaternion<T>& q) const;
    Quaternion<T> operator/(Quaternion<T>& q);
    Quaternion<T> operator/(const T s);
    //const Quaternion<T>& operator*=(const Quaternion<T>& q);
    bool operator==(const Quaternion<T>& q) const;
    bool operator!=(const Quaternion<T>& q) const;

    friend Quaternion operator+(const Quaternion& lhs,
            const Quaternion& rhs){
      return Quaternion(lhs) += rhs;
    }
    friend Quaternion operator-(const Quaternion& lhs,
            const Quaternion& rhs){
      return Quaternion(lhs) -= rhs;
    }
    friend Quaternion operator-(const Quaternion& q){
      return Quaternion()-q;
    }
    friend Quaternion operator*(const Quaternion& q, const T s){
      return Quaternion(q.w*s,q.x*s,q.y*s,q.z*s);
    }
    friend Quaternion operator*(const T s, const Quaternion& q){
      return Quaternion(q.w*s,q.x*s,q.y*s,q.z*s);
    }
    //friend Quaternion<T> operator*(const Quaternion<T>& lhs,
    //                                const Quaternion<T>& rhs);
    friend std::ostream& operator<<(std::ostream& output,
                                    const Quaternion& q){
      output << "[" << q.w << ", " << "(" << q.x <<
               ", " << q.y << ", " << q.z << ")]";
      return output;
    }
    void print_quat(void) const;
    Quaternion<T> conj(void);
    T norm2(void); //Squared norm
    T norm(void);
    Quaternion<T> inverse(void);
    Quaternion<T> normalize(void);
    void norm_inplace(void);

    template<typename V>
    void quatConjugateVec(V *v);
    template<typename V>
    void quatConjugateVec(V *v1,V *v2,V *v3);//Overloaded
    template<typename V>
    void quatConjugateVecRef(V *v,V *ref);
    template<typename V>
    void quatConjugateVecRef(V *v,V ref1,V ref2,V ref3);
    template<typename V>
    void quatConjugateVec(V *v1,V *v2,V *v3,V *ref);//Overloaded
    //template<typename V>
    //void quatConjugateVec(V v1,V v2,V v3,V ref1,V ref2,V ref3);//Overloaded
    // Methods for setting private components
    const Quaternion<T>& Set(T, T, T, T);
    const Quaternion<T>& fromAngleAxis(T angle,const T *axis);
    const Quaternion<T>& fromAngleAxis(T angle,T ax1,T ax2,T ax3);

    /* Random Quaternions */
    const setRandom_mc();
    const setRandom();

    static const RndQuatExpr randomQuaternion(int seed);
    const Quaternion<T>& operator=(const RndQuatExpr rhs);

    static Quaternion<T> getRandom_MC_quaternion();
};


#include "quaternion_class_test.cpp"
#endif
