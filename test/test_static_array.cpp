#include <iostream>
#include <ctime>
#include "stdlib.h"

using namespace std;

// function to generate and retrun random numbers.
int * getRandom(int my_seed) {

   static int  r[10];

   // set the seed
   //srand( (unsigned)time( NULL ) );
   srand(my_seed);
   
   for (int i = 0; i < 10; ++i) {
      r[i] = rand();
      cout << r[i] << endl;
   }

   return r;
}

// main function to call above defined function.
int main () {

   // a pointer to an int.
   int *p, *q;

   p = getRandom(1);
   
   for ( int i = 0; i < 10; i++ ) {
      cout << "*(p + " << i << ") : ";
      cout << *(p + i) << endl;
   }
   
   q = getRandom(2);

   for ( int i = 0; i < 10; i++ ) {
      cout << "*(q + " << i << ") : ";
      cout << *(q + i) << endl;
   }

   return 0;
}
