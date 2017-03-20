#include <iostream>

using namespace std;

void print(int nrow, int ncol, float mat[][3]);
void print_pp(int irl, int irh,int icl, int ich,float **mat);
//void convert2pp(int irl, int irh,int icl, int ich,float *first, float ***mat_pp);
void convert2pp(int irl, int irh,int icl, int ich,float *first,
  float *mat_pp_rows[], float ***mat_pp);
float **convert2pp_return(int irl, int irh,int icl, int ich,float *first,
  float *mat_pp_rows[]);

int main(int argc, char *argv[]){
  int i,j,icl,ich,irl,irh,nrow,ncol;
  nrow = 2;
  ncol = 3;
  irl = 1; irh = irl+nrow-1; icl = 1; ich = icl+ncol-1;
  float mat[nrow][3];
  float *mat_pp_rows[2];
  float **mat_pp, **mat_pp_r;

  cout << "mat_pp_rows "<< mat_pp_rows << endl;
  for (i = 0;i<nrow; i++){
    for (j = 0;j<ncol; j++){
      mat[i][j] = j + i*ncol;
    }
  }
  for (i = 0;i<nrow; i++){
    for (j = 0;j<ncol; j++){
      cout << mat[i][j] << " ";
    }
    cout << "\n";
  }
  cout << "---------" << endl;
  print(nrow,ncol,mat);
  cout << "---------" << endl;
  cout << "mat: "<< mat << endl;
  cout << "&mat: "<< &mat << endl;
  cout << "mat[0]: "<< mat[0] << endl;
  cout << "&mat[0]: "<< &mat[0] << endl;
  cout << "&mat[0][0]: "<< &mat[0][0] << endl;
  cout << "mat_pp: "<< mat_pp << endl;
  //convert2pp(irl, irh,icl,ich,mat[0],&mat_pp);
  convert2pp(irl, irh,icl,ich,mat[0],mat_pp_rows,&mat_pp);
  cout << "Outside function" << endl;
  cout << "mat_pp_rows "<< mat_pp_rows << endl;
  cout << "mat_pp_rows[0] "<< mat_pp_rows[0] << endl;
  cout << "*mat_pp_rows[0] "<< *mat_pp_rows[0] << endl;
  cout << "mat_pp_rows[1] "<< mat_pp_rows[1] << endl;
  cout << "*mat_pp_rows[1] "<< *mat_pp_rows[1] << endl;
  mat_pp = mat_pp_rows - 1;
  cout << "mat[1][1]: "<< mat_pp[1][1] << endl;
  for (i = 1;i<=2; i++){
    for (j = 1;j<=3; j++){
      cout << mat_pp[i][j] << " ";
    }
    cout << "\n";
  }
  cout << "---------" << endl;
  print_pp(irl,irh,icl,ich,mat_pp);
  cout << "---------" << endl;
  cout << "Sizes: " << endl;
  cout << "sizeof(float) "<< sizeof(float) << endl;
  cout << "sizeof(float*) "<< sizeof(float*) << endl;
  cout << "sizeof(mat) "<< sizeof(mat) << endl;
  cout << "sizeof(mat_pp_rows) "<< sizeof(mat_pp_rows) << endl;
  cout << "sizeof(mat_pp) "<< sizeof(mat_pp) << endl;
  cout << "---------" << endl;
  cout << "Convert with return value" << endl;
  mat_pp_r = convert2pp_return(irl, irh,icl,ich,mat[0],mat_pp_rows);
  print_pp(irl,irh,icl,ich,mat_pp_r);

}

void print(int nrow, int ncol,float mat[][3]){
  int i,j;
  for (i = 0;i<nrow; i++){
    for (j = 0;j<ncol; j++){
      cout << mat[i][j] << " ";
    }
    cout << "\n";
  }
}

void print_pp(int irl, int irh,int icl, int ich,float **mat){
  int i,j;
  cout << "mat_pp[2][3]: "<< mat[2][3] << endl;
  for (i = irl;i<=irh; i++){
    for (j = icl;j<=ich; j++){
      cout << mat[i][j] << " ";
    }
    cout << "\n";
  }
}

float **convert2pp_return(int irl, int irh,int icl, int ich,float *first,
  float *mat_pp_rows[]){
  int i,j,ncol,nrow;
  ncol = ich - icl + 1;
  cout << "In Function" << endl;
  cout << "mat_pp_rows "<< mat_pp_rows << endl;

  cout << "&mat_pp_rows[0]"<< &mat_pp_rows[0] << endl;
  cout << "&mat_pp_rows[1]"<< &mat_pp_rows[1] << endl;

  mat_pp_rows[0] = first - 1;
  mat_pp_rows[1] = first + ncol - 1;

  cout << "mat_pp_rows[0] "<< mat_pp_rows[0] << endl;
  cout << "*mat_pp_rows[0] "<< *mat_pp_rows[0] << endl;
  cout << "mat_pp_rows[1] "<< mat_pp_rows[1] << endl;
  cout << "*mat_pp_rows[1] "<< *mat_pp_rows[1] << endl;

  float **mat = mat_pp_rows - 1;

  cout << "mat[1][1]: "<< mat[1][1] << endl;
  for (i = 1;i<=2; i++){
    for (j = 1;j<=3; j++){
      cout << mat[i][j] << " ";
    }
    cout << "\n";
  }
  return mat;
}

void convert2pp(int irl, int irh,int icl, int ich,float *first,
  float *mat_pp_rows[],float ***mat_pp){
  int i,j,ncol,nrow;
  ncol = ich - icl + 1;
  cout << "In Function" << endl;
  cout << "mat_pp_rows "<< mat_pp_rows << endl;

  cout << "&mat_pp_rows[0]"<< &mat_pp_rows[0] << endl;
  cout << "&mat_pp_rows[1]"<< &mat_pp_rows[1] << endl;

  mat_pp_rows[0] = first - 1;
  mat_pp_rows[1] = first + ncol - 1;

  cout << "mat_pp_rows[0] "<< mat_pp_rows[0] << endl;
  cout << "*mat_pp_rows[0] "<< *mat_pp_rows[0] << endl;
  cout << "mat_pp_rows[1] "<< mat_pp_rows[1] << endl;
  cout << "*mat_pp_rows[1] "<< *mat_pp_rows[1] << endl;

  float **mat = mat_pp_rows - 1;

  cout << "mat[1][1]: "<< mat[1][1] << endl;
  for (i = 1;i<=2; i++){
    for (j = 1;j<=3; j++){
      cout << mat[i][j] << " ";
    }
    cout << "\n";
  }
  mat_pp = &mat;
}

// void convert2pp(int irl, int irh,int icl, int ich,float *first, float ***mat){
//   int i,j,ncol,nrow;
//   float **mat_pp;
//   ncol = ich - icl + 1;
//
//   float *mat_pp_rows[2] = {first, first + ncol};
//   cout << "mat_pp_rows[0] "<< mat_pp_rows[0] << endl;
//   cout << "*mat_pp_rows[0] "<< *mat_pp_rows[0] << endl;
//   cout << "mat_pp_rows[1] "<< mat_pp_rows[1] << endl;
//   cout << "*mat_pp_rows[1] "<< *mat_pp_rows[1] << endl;
//   // cout << "ncol: "<< ncol << endl;
//   cout << "first: "<< first << endl;
//   cout << "&first: "<< &first << endl;
//   // cout << "mat_pp: "<< mat_pp << endl;
//   // mat_pp = &first;
//   // cout << "(*mat_pp): "<< (*mat_pp) << endl;
//   // cout << "mat_pp[0]: "<< mat_pp[0] << endl;
//   // cout << "mat_pp[0][0]: "<< mat_pp[0][0] << endl;
//   // mat_pp = &first - 1;
//   // mat_pp[2] = first + 2 + 1;
//   // cout << "(*mat_pp): "<< (*mat_pp) << endl;
//   // cout << "mat_pp[1]: "<< mat_pp[1] << endl;
//   // mat_pp[1] = first;
//   // cout << "mat_pp[1]: "<< mat_pp[1] << endl;
//   // cout << "mat_pp[2]: "<< mat_pp[2] << endl;
//   // cout << "mat_pp[1][0]: "<< mat_pp[1][0] << endl;
//   // cout << "mat_pp[1][2]: "<< mat_pp[1][2] << endl;
//   // cout << "mat_pp[2][1]: "<< mat_pp[2][1] << endl;
//   //
//   // mat_pp = &first - 1;
//   // mat_pp = &first - 1;
//   // cout << "mat_pp[1][1]: "<< mat_pp[1][1] << endl;
//   // //for(i=irl;i<=irh;i++){
//   // //  mat_pp[i] = first + (i-1)*ncol - 1;
//   // //}
//   mat[0][1] = first  - 1;
//   mat[0][2] = first  + ncol -1;
//   cout << "addresses"<<endl;
//   cout << "mat[1][1]: "<< mat[1][1] << endl;
//   cout << "mat[1][2]: "<< mat[1][2] << endl;
//   cout << "mat[1][3]: "<< mat[1][3] << endl;
//   //cout << "mat_pp[2][0]: "<< mat_pp[2][0] << endl;
//   cout << "mat[2][1]: "<< mat[2][1] << endl;
//   cout << "mat[2][2]: "<< mat[2][2] << endl;
//   cout << "mat[2][3]: "<< mat[2][3] << endl;
//
//   cout << "values"<<endl;
//   cout << "*mat[1][1]: "<< *mat[1][1] << endl;
//   cout << "*mat[1][2]: "<< *(mat[1][2]) << endl;
//   cout << "*mat[1][3]: "<< *mat[1][3] << endl;
//   //cout << "mat_pp[2][0]: "<< mat_pp[2][0] << endl;
//   cout << "*mat[2][1]: "<< *mat[2][1] << endl;
//   cout << "*mat[2][2]: "<< *mat[2][2] << endl;
//   cout << "*mat[2][3]: "<< *mat[2][3] << endl;
//
//   mat_pp[1] = first  - 1;
//   mat_pp[2] = first + ncol -1;
//   cout << "---------" << endl;
//   //cout << "mat_pp[1][0]: "<< mat_pp[1][0] << endl;
//   cout << "mat_pp[1][1]: "<< mat_pp[1][1] << endl;
//   cout << "mat_pp[1][2]: "<< mat_pp[1][2] << endl;
//   cout << "mat_pp[1][3]: "<< mat_pp[1][3] << endl;
//   //cout << "mat_pp[2][0]: "<< mat_pp[2][0] << endl;
//   cout << "mat_pp[2][1]: "<< mat_pp[2][1] << endl;
//   cout << "mat_pp[2][2]: "<< mat_pp[2][2] << endl;
//   cout << "mat_pp[2][3]: "<< mat_pp[2][3] << endl;
//   i = 2; j = 3;
//   cout << "aaamat_pp[2][3]: "<< mat_pp[i][j] << endl;
//   for (i = 1;i<=2; i++){
//     for (j = 1;j<=3; j++){
//       cout << mat_pp[i][j] << " ";
//     }
//     cout << "\n";
//   }
//   mat = &mat_pp;
//
//
//
// }
