#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "funct.h"

/* #ifdef OMP */
/* #include <omp.h> */
/* #endif */

void VWGReP(const int ReAtNu,float **ReCoor,float *ReVdWR_p3,float *ReVdWE_sr,
            const float VWGrIn, const float VWGrSi,float *BSMinC,int *VWGPoN,
            float *ReMaxC,float *ReMinC,
            float ***VWGrRP_at,float ***VWGrRP_re)
/* This function computes the receptor part of the van der Waals interaction 
   on a grid :
   VWGPoN     number of points for the grid of the van der Waals interaction 
              (1->x,2->y,3->z)
   VWGrRP_at  value of the receptor part of the van der Waals attractive 
              interaction on the grid
   VWGrRP_re  value of the receptor part of the van der Waals repulsive  
              interaction on the grid  
   VWGCor     coordinates of one point of the vdW grid (1->x,2->y,3->z)
   GRSqDi     squared distance between two points (of the grid and receptor) 
   SizGr_LRA  size of the grid containing the list of receptor atoms
   CubNum_LRA  the number of cubes along each direction (1->x,2->y,3->z)
   CubFAI_LRA  index of the first atoms in each cube
   CubLAI_LRA  index of the last atoms in each cube
   CubLiA_LRA  list of atoms in each cube
   CubAto_LRA  cube number in which the current vdW grid point is (1->x,2->y,
               3->z) */
{
  int i,j,k,m,n,j1,k1,l1,CubNum_LRA[4],***CubFAI_LRA,***CubLAI_LRA,
      *CubLiA_LRA,CubAto_LRA[4];
  float VWGCor[4],GRSqDi,VWGrRP_at_in,VWGrRP_re_in,TwoPo7,TwoP12,GRSqDi_p3,
        SizGr_LRA,VWGrCutoff;

/* The size of the grid containing the list of receptor atoms is set as half 
   the real size of the grid in SizGr_LRA because of ReLAIC function */
  SizGr_LRA=5.0;
  VWGrCutoff=(2.0*SizGr_LRA)*(2.0*SizGr_LRA);

/* Construct the list of receptor atoms which are in the cubes of a grid */
  ReLAIC(ReAtNu,ReCoor,SizGr_LRA,ReMaxC,ReMinC,CubNum_LRA,&CubFAI_LRA,&CubLAI_LRA,&CubLiA_LRA);

/* For each vdW grid point, the contribution for the vdW potential of the 
   receptor part involves only the receptor atoms which are in the cube to 
   which the current vdW grid point belongs and the surrounding cubes */
  TwoPo7=2*2*2*2*2*2*2;
  TwoP12=2*2*2*2*2*2*2*2*2*2*2*2;

/* #pragma omp parallel default(shared) private(VWGCor,VWGrRP_at_in,VWGrRP_re_in,CubAto_LRA,m,n,GRSqDi,GRSqDi_p3,j1,k1,l1,i,j,k) */

/* #ifdef OMP */
/* #pragma omp parallel for default(none) shared(VWGrRP_re,CubLAI_LRA,BSMinC,ReMinC,ReMaxC,CubNum_LRA,CubFAI_LRA,CubLiA_LRA,SizGr_LRA,ReCoor,VWGrCutoff,ReVdWE_sr,ReVdWR_p3,TwoPo7,TwoP12,VWGrRP_at,VWGPoN) private(VWGCor,VWGrRP_at_in,VWGrRP_re_in,CubAto_LRA,GRSqDi,GRSqDi_p3,m,n,j1,k1,l1,i,j,k) */
/* #endif */
  for (i=1;i<=VWGPoN[1];i++) {

/* #pragma omp critical */
    /* printf("VWPT %d\n",i);   */

    for (j=1;j<=VWGPoN[2];j++) {
#ifdef OMP
#pragma omp parallel for default(none) shared(i,j,VWGrRP_re,CubLAI_LRA,BSMinC,ReMinC,ReMaxC,CubNum_LRA,CubFAI_LRA,CubLiA_LRA,SizGr_LRA,ReCoor,VWGrCutoff,ReVdWE_sr,ReVdWR_p3,TwoPo7,TwoP12,VWGrRP_at,VWGPoN) private(VWGCor,VWGrRP_at_in,VWGrRP_re_in,CubAto_LRA,m,n,GRSqDi,GRSqDi_p3,j1,k1,l1,k)
#endif
      for (k=1;k<=VWGPoN[3];k++) {
/* Coordinates of one point of the van der Waals grid */
        VWGCor[1]=BSMinC[1]-VWGrIn+(i-1)*VWGrSi;
        VWGCor[2]=BSMinC[2]-VWGrIn+(j-1)*VWGrSi;
        VWGCor[3]=BSMinC[3]-VWGrIn+(k-1)*VWGrSi;

/* Initialization of the "potentials" */
        VWGrRP_at_in=0.0;
        VWGrRP_re_in=0.0;

/* Find the cube in which the current vdW grid point is or the nearest cube 
   if this vdW grid point doesn't belong to any cubes of the (list of receptor 
   atoms) grid */
        for (m=1;m<=3;m++) {
          if (VWGCor[m]<ReMinC[m])
            CubAto_LRA[m]=1;
          else if (VWGCor[m]>ReMaxC[m])
            CubAto_LRA[m]=CubNum_LRA[m];
          else
            CubAto_LRA[m]=ffloor((VWGCor[m]-ReMinC[m])/(2*SizGr_LRA))+1;
        }

/* Go into the 27 cubes when possible */
        for (j1=-1;j1<=1;j1++) {
          for (k1=-1;k1<=1;k1++) {
            for (l1=-1;l1<=1;l1++) {

              if (((CubAto_LRA[1]+j1)>=1)&&((CubAto_LRA[1]+j1)<=CubNum_LRA[1])&&
                  ((CubAto_LRA[2]+k1)>=1)&&((CubAto_LRA[2]+k1)<=CubNum_LRA[2])&&
                  ((CubAto_LRA[3]+l1)>=1)&&((CubAto_LRA[3]+l1)<=CubNum_LRA[3])) {

                if (CubFAI_LRA[CubAto_LRA[1]+j1][CubAto_LRA[2]+k1]
                              [CubAto_LRA[3]+l1]!=0) {

                  for (n=CubFAI_LRA[CubAto_LRA[1]+j1][CubAto_LRA[2]+k1]
                                   [CubAto_LRA[3]+l1];
                       n<=CubLAI_LRA[CubAto_LRA[1]+j1][CubAto_LRA[2]+k1]
                                    [CubAto_LRA[3]+l1];n++) {

                    m=CubLiA_LRA[n];

                    GRSqDi=DistSq(VWGCor[1],VWGCor[2],VWGCor[3],
                                  ReCoor[m][1],ReCoor[m][2],ReCoor[m][3]);

                    if (GRSqDi>(1.e-4)) {
                      if (GRSqDi<=VWGrCutoff) {
                        GRSqDi_p3=GRSqDi*GRSqDi*GRSqDi;
                        VWGrRP_at_in=VWGrRP_at_in+ReVdWE_sr[m]*ReVdWR_p3[m]/
                                                  (GRSqDi_p3);
                        VWGrRP_re_in=VWGrRP_re_in+ReVdWE_sr[m]*ReVdWR_p3[m]*
                                             ReVdWR_p3[m]/(GRSqDi_p3*GRSqDi_p3);
                      }            
                    }
                    else {
                      VWGrRP_re_in=VWGrRP_re_in+(1.e+12)/TwoP12;
                    }

                  }

                }

              }

            }
          }
        }

        VWGrRP_at[i][j][k]=-VWGrRP_at_in*TwoPo7;
        VWGrRP_re[i][j][k]=VWGrRP_re_in*TwoP12;
 
      }
    }
  }

  free_i3tensor(CubFAI_LRA,1,CubNum_LRA[1],1,CubNum_LRA[2],1,CubNum_LRA[3]);
  free_i3tensor(CubLAI_LRA,1,CubNum_LRA[1],1,CubNum_LRA[2],1,CubNum_LRA[3]);
  free_ivector(CubLiA_LRA,1,ReAtNu+10);

}
