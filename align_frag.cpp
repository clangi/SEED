#include <math.h>
#include <iostream>
#include "funct.h"
#include "quaternion.h"
#include "eigen.h"
//#include "nrutil.h"

using namespace std;

bool set_align_ref(float **FrCoor,int FrAtNu, float **ref, float **FrAlSet)
  /*##########################################
  Set a three atom reference for fragment alignment. The first atom of
  the fragment is placed in the origin, the second on the x-axis and the
  third (or if collinear, the fourth, fifth, and so on) on the x-z plane.
  ############################################*/

  /*##########################################
  float **FrCoor - Original fragment coordinates
  int FrAtNu - Number of atoms
  float **ref - New coordinates of the three atoms used as reference
  float **FrAlSet - Old coordinates of the three atoms used as reference
  ############################################*/
{
  float ang123,d12,d23;
  float dAB[4], dAC[4], cPr[4];
  //float align_set[4][4];// ref[4][4];
  int i,j, iat;
  /*##########################################
  int i,j,iat - multipurpose indices
  float ang123 - angle (in range 0-180) between atoms 1,2,3
  float d12 - distance between atom 1 and 2
  float d13 - distance between atom 1 and 3
  dAB - difference between coordinates of second and first atom
  dAC - difference between coordinates of third and first atom
  cPr - Cross product between dAB and dAC to check for collinearity
  ############################################*/
  dAB[1] = FrCoor[2][1] - FrCoor[1][1];
  dAB[2] = FrCoor[2][2] - FrCoor[1][2];
  dAB[3] = FrCoor[2][3] - FrCoor[1][3];
  for(iat=3;iat<=FrAtNu;iat++){
    dAC[1] = FrCoor[iat][1] - FrCoor[1][1];
    dAC[2] = FrCoor[iat][2] - FrCoor[1][2];
    dAC[3] = FrCoor[iat][3] - FrCoor[1][3];

    VectPr(dAB[1],dAB[2],dAB[3],dAC[1],dAC[2],dAC[3],
      &cPr[1],&cPr[2],&cPr[3]);
      //std::cout<< "VectPro " << (cPr[1]*cPr[1]+cPr[2]*cPr[2]+cPr[3]*cPr[3])<< std::endl;
    if ((cPr[1]*cPr[1]+cPr[2]*cPr[2]+cPr[3]*cPr[3]) != 0.0){
      for(i=1; i<=3;i++){
        for(j=1; j<=3;j++){
            //align_set[i][j] = FrCoor[i][j];
          ref[i][j] = 0.0;
        }
      }
      d12 = sqrtf(dAB[1]*dAB[1]+dAB[2]*dAB[2]+dAB[3]*dAB[3]);
      ref[2][1] = d12; //y=z=0
      //std::cout << "d12 " << d12<<std::endl;
      d23 = sqrtf(DistSq(FrCoor[2][1],FrCoor[2][2],FrCoor[2][3],
        FrCoor[iat][1],FrCoor[iat][2],FrCoor[iat][3]));
      //std::cout << "d23 " << d23<<std::endl;
      ang123 = PlaAng(FrCoor[2][1],FrCoor[2][2],FrCoor[2][3], //B
          FrCoor[1][1],FrCoor[1][2],FrCoor[1][3],//A
          FrCoor[iat][1],FrCoor[iat][2],FrCoor[iat][3]);//C
      //std::cout << "ang123 " << ang123<<std::endl;
      //ref[3][2] = 0.0; y = 0
      ref[3][1] = d12 - d23*cosf(ang123);
      ref[3][3] = d23*sinf(ang123);
      //cout<<"ref[2][1] = "<<ref[2][1]<<endl;
      //cout<<"ref[3][1] = "<<ref[3][1]<<endl;
      //cout<<"ref[3][3] = "<<ref[3][3]<<endl;
      FrAlSet[1][1] = FrCoor[1][1]; //The first two atoms of the align set are
      FrAlSet[1][2] = FrCoor[1][2]; //always the first two atoms of the
      FrAlSet[1][3] = FrCoor[1][3]; //fragment. The third may change,
      FrAlSet[2][1] = FrCoor[2][1]; //depending on collinearity.
      FrAlSet[2][2] = FrCoor[2][2];
      FrAlSet[2][3] = FrCoor[2][3];
      FrAlSet[3][1] = FrCoor[iat][1];
      FrAlSet[3][2] = FrCoor[iat][2];
      FrAlSet[3][3] = FrCoor[iat][3];
      return true;
    }
  }
  std::cout << "All atoms are colliner. "
        << "Check input for errors" << std::endl;
  return false;
}

void struct_align(float **FrCoor,int FrAtNu, float **FrAlRef,
                float **FrAlSet,int AlAtNu)
  /*########################
  Perform least square alignment using quaternions. It is a wrapper for align_3D
  float **FrCoor - Fragment coords
  int FrAtNu - Number of atoms
  float **FrAlRef - Reference for alignment
  float **FrAlSet - Alignment set
  int AlAtNu - Number of atoms in the alignment set
  ########################*/
{
  Quaternion<double> *qrot;
  float tvec[4], cto[4], ctn[4];
  int iat;
  /*########################
  Quaternion<double> *qrot - Rotation Quaternion
  float tvec[4] - Tranlation vectors
  float cto[4] - Old centroid of the alignment set
  float ctn[4] - New centroid of the alignment set
  ########################*/
  //cout << "before align 3D" << endl;
  align_3D(AlAtNu,FrAlSet,FrAlRef,&qrot,tvec,cto,ctn);
  //cout << "after align 3D" << endl;
  //qrot->print_quat();
  qrot->norm_inplace();
  //qrot->print_quat();

  //cout << "qrot normalized" << endl;
  for(iat=1;iat<=FrAtNu;iat++){
    qrot->quatConjugateVecRef(&FrCoor[iat][0],cto); // Quaternion conjugation
    //with reference vector
    FrCoor[iat][1] += tvec[1];
    FrCoor[iat][2] += tvec[2];
    FrCoor[iat][3] += tvec[3];
  }
  //cout << "after whole molecule align" << endl;
  delete qrot;
}

template<class T>
void align_3D(int n,float **FrAlSet,float **FrAlRef,Quaternion<T> **qrot,
                float *tvec,float *cto,float *ctn)
  /*########################
  Calculate quaternion and translation vector for least square alignment.
  int n - Number of atoms
  float **FrAlSet - Coordinates of set to align
  float **FrAlRef - Coordinates of reference set
  Quaternion<T> **qrot - Alignment quaternion
  float *tvec - Translation vector
  float *cto - Old centroid
  float *ctn - New cetroid
  ########################*/
{
  Quaternion<double> *qrot_L;
  int iat,i,j;
  tvec[1] = 0.0; tvec[2] = 0.0; tvec[3] = 0.0;
  double **mata,**matb,**matn,**ntmp;
  double **mataT; // auxiliary
  const double **eigv;
  //Calculate centroid:
  CentroidCalc(FrAlSet,n,cto);
  CentroidCalc(FrAlRef,n,ctn);
  //tvec
  tvec[1] = ctn[1] - cto[1];
  tvec[2] = ctn[2] - cto[2];
  tvec[3] = ctn[3] - cto[3];
  //Allocation and initialization to zero matrix:
  mata = zero_dmatrix(1,4,1,4);
  matb = zero_dmatrix(1,4,1,4);
  matn = zero_dmatrix(1,4,1,4);
  ntmp = zero_dmatrix(1,4,1,4);
  mataT = zero_dmatrix(1,4,1,4);
  //Calculation of matn
  for(iat=1;iat<=n;iat++){
    //mata
    //x
    mata[1][2] = -(FrAlSet[iat][1]-cto[1]);
    mata[2][1] = (FrAlSet[iat][1]-cto[1]);
    mata[3][4] = (FrAlSet[iat][1]-cto[1]);
    mata[4][3] = -(FrAlSet[iat][1]-cto[1]);
    //y
    mata[1][3] = -(FrAlSet[iat][2]-cto[2]);
    mata[3][1] = (FrAlSet[iat][2]-cto[2]);
    mata[2][4] = -(FrAlSet[iat][2]-cto[2]);
    mata[4][2] = (FrAlSet[iat][2]-cto[2]);
    //z
    mata[1][4] = -(FrAlSet[iat][3]-cto[3]);
    mata[4][1] = (FrAlSet[iat][3]-cto[3]);
    mata[2][3] = (FrAlSet[iat][3]-cto[3]);
    mata[3][2] = -(FrAlSet[iat][3]-cto[3]);
    //matb
    //x
    matb[1][2] = -(FrAlRef[iat][1]-ctn[1]);
    matb[2][1] = (FrAlRef[iat][1]-ctn[1]);
    matb[3][4] = -(FrAlRef[iat][1]-ctn[1]);
    matb[4][3] = (FrAlRef[iat][1]-ctn[1]);
    //y
    matb[1][3] = -(FrAlRef[iat][2]-ctn[2]);
    matb[3][1] = (FrAlRef[iat][2]-ctn[2]);
    matb[2][4] = (FrAlRef[iat][2]-ctn[2]);
    matb[4][2] = -(FrAlRef[iat][2]-ctn[2]);
    //z
    matb[1][4] = -(FrAlRef[iat][3]-ctn[3]);
    matb[4][1] = (FrAlRef[iat][3]-ctn[3]);
    matb[2][3] = -(FrAlRef[iat][3]-ctn[3]);
    matb[3][2] = (FrAlRef[iat][3]-ctn[3]);

    dm_transpose(4,4,mata,mataT);
    dmm_prod(4,4,4,4,mataT,matb,ntmp);
    for(i=1;i<=4;i++){
      for(j=1;j<=4;j++){
        matn[i][j] = matn[i][j] + ntmp[i][j];
      }
    }
  } //end of loop over atoms
  //cout << "before eigen" << endl;
  Jacobi eigmatn(matn,4); // Run Jacobi diagonalization
  //cout << "before eigen 1" << endl;
  eigv = eigmatn.eigenvectors();
  //cout << "before eigen 2" << endl;
  qrot_L = new Quaternion<double>(eigv[1][1],eigv[2][1],eigv[3][1],eigv[4][1]);
  for(i=1;i<=3;i++){
    tvec[i] = ctn[i] - cto[i];
  }
  //cout << "quaternion defined" << endl;
  (*qrot) = qrot_L;
  //qrot_L->print_quat();
  //Cleaning:
  free_dmatrix(mata,1,4,1,4);
  free_dmatrix(matb,1,4,1,4);
  free_dmatrix(matn,1,4,1,4);
  free_dmatrix(ntmp,1,4,1,4);
  free_dmatrix(mataT,1,4,1,4);
}

void CentroidCalc(float **Coord,int nAt,float *centroid)
  /*########################
  Calculate geometric centroid.
  float **Coord - molecule coordinates
  int nAt - number of atoms
  float *centroid - centroid coordinates
  ########################*/
{
  int i;
  centroid[1] = centroid[2] = centroid[3] = 0.0;
  for(i=1;i<=nAt;i++){
    centroid[1] += Coord[i][1];
    centroid[2] += Coord[i][2];
    centroid[3] += Coord[i][3];
  }
  centroid[1] /= nAt;
  centroid[2] /= nAt;
  centroid[3] /= nAt;
}
