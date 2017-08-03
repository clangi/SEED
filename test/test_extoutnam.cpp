#define _STRLENGTH 500
#include "../src/extoutnam.cpp"

int main(int argc, char *argv[]){
  char mystring[] = "/users/clangini/etc/cgenff_library.mol2";
  char mystring_out[500];

  ExtOutNam(mystring, mystring_out);
  printf("mystring = %s\n", mystring);
  printf("mystring_out = %s\n", mystring_out);

  ExtOutNam(mystring, mystring_out, 10);
  printf("mystring = %s\n", mystring);
  printf("mystring_out = %s\n", mystring_out);

}
