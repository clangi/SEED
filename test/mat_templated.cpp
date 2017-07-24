#include <iostream>
#include <cstdlib>
using namespace std;

#define Dynamic -1

template <typename T, int NROW, int NCOL>
class Mat {
  T data[NROW * NCOL];
  //int n_elem = NROW * NCOL; // only available in c++11
  int n_elem;
public:
  Mat(T fill){
    n_elem = NROW * NCOL;
    for(int i = 0; i < n_elem; i++)
      data[i] = fill;
  };
  void print(){
    for(int i = 0; i < n_elem; i++)
      cout << "data[" << i << "] = " << data[i] << endl;
  };
};

/* Partial template specialization */
template <typename T>
class Mat<T,Dynamic,Dynamic> {
  T *data;
  int nrow;
  int ncol;
  int n_elem;
public:
  Mat(T fill, int nrow_, int ncol_){
    nrow = nrow_;
    ncol = ncol_;
    n_elem = nrow * ncol;
    data = new T[n_elem]; //static_cast< T* >(malloc(sizeof(T)*n_elem));
    for (int i = 0; i < n_elem; i++){
      data[i] = fill;
    }
  };
  ~Mat(){
    delete [] data;
  };
  void print(){
    for(int i = 0; i < n_elem; i++)
      cout << "data[" << i << "] = " << data[i] << endl;
  };
};
/* we cannot specialize a single member of the function. A specialized
member needs to reside in its specialized function. An better approach may be
the overload of the constructor and the destructor */

int main(int argc, char *argv[]){
  int user_value, ncol, nrow;
  int a, b;

  Mat<char,2,3> myMat('a');
  myMat.print();

  //const int a = 3; // with this it works
  //const int b = 5;
  //a = 3;
  //b = 4;
  //const int A = const_cast<const int&>(a); // does not work
  //const int B = const_cast<const int&>(b);
  //cout << A << " " << B << endl;
  //Mat<char,A,B> myMat_2('b');
  //myMat_2.print();

  cout << "Input nrow, ncol, nvalue: " << endl;
  cin >> nrow >> ncol >> user_value;
  cout << "Inserted: " << nrow << " " << ncol << " " << user_value << endl;
  Mat<int, Dynamic, Dynamic> userMat(user_value, nrow, ncol);
  userMat.print();
}
