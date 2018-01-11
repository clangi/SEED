#include <stdio.h>
#include <math.h>
#include "funct.h"
#include "nrutil.h"

void SqDisFrRe_ps(int FrAtNu,double **RoSFCo,double **ReCoor,double *ReMinC,
            double GrSiCu_en,int *CubNum_en,int ***CubFAI_en,int ***CubLAI_en,
            int *CubLiA_en,int PsSpNC,int ***PsSphe,double **SDFrRe_ps,
            int ReAtNu,double PsSpRa,double *RePaCh,int ReReNu,int *AtReprRes,
            int *FiAtRes,int *LaAtRes,double *TotChaRes,int NuChResEn,
            int *LiChResEn,double **SDFrRe_ps_elec,double **ChFrRe_ps_elec)
/* This function computes the squared distances between the fragment atoms and
   the corresponding receptor atoms of the pseudo-sphere procedure, one list
   for the vdW energy and another for the electrostatic interaction energy.
   In the case of the electrostatic interaction energy, it also provides the
   corresponding partial or total charges :
   SDFrRe_ps  array of squared distances between the fragment atoms and the
              corresponding receptor atoms of the pseudo-sphere procedure
              for the evaluation of the vdW energy
   SDFrRe_ps_elec  same as above but for the evaluation of the electrostatic
                   interaction energy
   ChFrRe_ps_elec  receptor partial or total charges for the electrostatic
                   interaction energy
   ResFlag  residue flag (1 if taken into account, 0 if not)
   CubAto  cube number in which the current atom is (1->along x,2->y,3->z)
   RFSqDi  squared distance between one atom of the receptor and one atom
           of the fragment
   PsSpRa_sq  squared pseudo-sphere radius */
{
  int i,j,k,m,n,CubAto[4],p,*ResFlag;
  double RFSqDi,PsSpRa_sq;

  PsSpRa_sq=PsSpRa*PsSpRa;
  ResFlag=ivector(1,ReReNu);

/* Initialization of SDFrRe_ps and SDFrRe_ps_elec */
  for (i=1;i<=FrAtNu;i++) {
    for (j=1;j<=ReAtNu;j++) {
      SDFrRe_ps[i][j]=-1.0;
      SDFrRe_ps_elec[i][j]=-1.0;
    }
  }

/* Initialization of ChFrRe_ps_elec */
  for (i=1;i<=FrAtNu;i++) {
    for (j=1;j<=ReAtNu;j++) {
      ChFrRe_ps_elec[i][j]=RePaCh[j];
    }
  }

/* ---------------------------- */
/* Loop over the fragment atoms */
/* ---------------------------- */
  for (i=1;i<=FrAtNu;i++) {

/* Initialization of ResFlag */
    for (j=1;j<=ReReNu;j++)
      ResFlag[j]=0;

/* Find the cube or pseudo-cube (negative assignation) in which the current
   fragment atom is */
    for (j=1;j<=3;j++)
      CubAto[j]=ffloor((RoSFCo[i][j]-ReMinC[j])/GrSiCu_en)+1;

/* Loop over the cube of the pseudo-sphere */
    for (j=-PsSpNC;j<=PsSpNC;j++) {
      for (k=-PsSpNC;k<=PsSpNC;k++) {
        for (m=-PsSpNC;m<=PsSpNC;m++) {

/* The cube of the pseudo-sphere must have a value of 1 */
/* This means that I should be within the pseudo-sphere. clangini */
          if (PsSphe[j+PsSpNC+1][k+PsSpNC+1][m+PsSpNC+1]) {

/* Check whether the list of residue-representative atoms cube exists */
            if (((CubAto[1]+j)>=1)&&((CubAto[1]+j)<=CubNum_en[1])&&
                ((CubAto[2]+k)>=1)&&((CubAto[2]+k)<=CubNum_en[2])&&
                ((CubAto[3]+m)>=1)&&((CubAto[3]+m)<=CubNum_en[3])) {

/* Check whether the list of residue-representative atoms cube contains at
   least one residue-representative atom */
              if (CubFAI_en[CubAto[1]+j][CubAto[2]+k][CubAto[3]+m]!=0) {

                for (n=CubFAI_en[CubAto[1]+j][CubAto[2]+k][CubAto[3]+m];
                     n<=CubLAI_en[CubAto[1]+j][CubAto[2]+k][CubAto[3]+m];n++) {

                  RFSqDi=DistSq(RoSFCo[i][1],RoSFCo[i][2],RoSFCo[i][3],
                                ReCoor[AtReprRes[CubLiA_en[n]]][1],
                                ReCoor[AtReprRes[CubLiA_en[n]]][2],
                                ReCoor[AtReprRes[CubLiA_en[n]]][3]);

/* Take the residue into account only if this distance (between the current
   fragment atom and the current residue-representative atom) is smaller or
   equal to the pseudo-sphere radius */
/* Residue is taken into account only if it is within the cut-off distance
   from the current fragment atom. clangini */
                  if (RFSqDi<=PsSpRa_sq) {
                    ResFlag[CubLiA_en[n]]=1;
/* Loop over the current residue atoms */
                    for (p=FiAtRes[CubLiA_en[n]];p<=LaAtRes[CubLiA_en[n]];p++) {
                      RFSqDi=DistSq(RoSFCo[i][1],RoSFCo[i][2],RoSFCo[i][3],
                                    ReCoor[p][1],ReCoor[p][2],ReCoor[p][3]);
                      SDFrRe_ps[i][p] = RFSqDi;
                      SDFrRe_ps_elec[i][p] = RFSqDi;
                    }
                  }

                }

              }

            }
          }

        }
      }
    }

/* Determine which charged residues still have to be taken into account as
   monopoles (not explicitly but only the residue-representative atom with
   the total charge of the residue on it). These monopole residues are added
   only to the list for the electrostatic interaction energy */
    for (j=1;j<=NuChResEn;j++) {
      if (!ResFlag[LiChResEn[j]]) {
        RFSqDi=DistSq(RoSFCo[i][1],RoSFCo[i][2],RoSFCo[i][3],
                      ReCoor[AtReprRes[LiChResEn[j]]][1],
                      ReCoor[AtReprRes[LiChResEn[j]]][2],
                      ReCoor[AtReprRes[LiChResEn[j]]][3]);
        SDFrRe_ps_elec[i][AtReprRes[LiChResEn[j]]]=RFSqDi;
        ChFrRe_ps_elec[i][AtReprRes[LiChResEn[j]]]=TotChaRes[LiChResEn[j]];
      }
    }

  }

  free_ivector(ResFlag,1,ReReNu);

}

void SqDisFrRe_ps_vdW(int FrAtNu,double **RoSFCo,double **ReCoor,double *ReMinC,
            double GrSiCu_en,int *CubNum_en,int ***CubFAI_en,int ***CubLAI_en,
            int *CubLiA_en,int PsSpNC,int ***PsSphe,double **SDFrRe_ps,
            int ReAtNu,double PsSpRa,int ReReNu,int *AtReprRes,
            int *FiAtRes,int *LaAtRes)
/* This function is identical to SqDisFrRe_ps() but for the vdW energy only.
   This function computes the squared distances between the fragment atoms and
   the corresponding receptor atoms of the pseudo-sphere procedure, for the vdW
   energy only.
   It will be used in MC run. clangini

   SDFrRe_ps  array of squared distances between the fragment atoms and the
              corresponding receptor atoms of the pseudo-sphere procedure
              for the evaluation of the vdW energy
   ResFlag  residue flag (1 if taken into account, 0 if not)
   CubAto  cube number in which the current atom is (1->along x,2->y,3->z)
   RFSqDi  squared distance between one atom of the receptor and one atom
           of the fragment
   PsSpRa_sq  squared pseudo-sphere radius */
{
  int i,j,k,m,n,CubAto[4],p,*ResFlag;
  double RFSqDi,PsSpRa_sq;

  PsSpRa_sq = PsSpRa*PsSpRa; // squared pseudo-sphere radius
  ResFlag=ivector(1,ReReNu);

/* Initialization of SDFrRe_ps */
  for (i=1;i<=FrAtNu;i++) {
    for (j=1;j<=ReAtNu;j++) {
      SDFrRe_ps[i][j]=-1.0;
    }
  }
/* ---------------------------- */
/* Loop over the fragment atoms */
/* ---------------------------- */
  for (i=1;i<=FrAtNu;i++) {

/* Initialization of ResFlag */
    for (j=1;j<=ReReNu;j++)
      ResFlag[j]=0;

/* Find the cube or pseudo-cube (negative assignation) in which the current
   fragment atom is */
    for (j=1;j<=3;j++)
      CubAto[j]=ffloor((RoSFCo[i][j]-ReMinC[j])/GrSiCu_en)+1;

/* Loop over the cube of the pseudo-sphere */
    for (j=-PsSpNC;j<=PsSpNC;j++) {
      for (k=-PsSpNC;k<=PsSpNC;k++) {
        for (m=-PsSpNC;m<=PsSpNC;m++) {

/* The cube of the pseudo-sphere must have a value of 1 */
/* This means that I should be within the pseudo-sphere. clangini */
          if (PsSphe[j+PsSpNC+1][k+PsSpNC+1][m+PsSpNC+1]) {

/* Check whether the list of residue-representative atoms cube exists */
            if (((CubAto[1]+j)>=1)&&((CubAto[1]+j)<=CubNum_en[1])&&
                ((CubAto[2]+k)>=1)&&((CubAto[2]+k)<=CubNum_en[2])&&
                ((CubAto[3]+m)>=1)&&((CubAto[3]+m)<=CubNum_en[3])) {

/* Check whether the list of residue-representative atoms cube contains at
   least one residue-representative atom */
              if (CubFAI_en[CubAto[1]+j][CubAto[2]+k][CubAto[3]+m]!=0) {

                for (n=CubFAI_en[CubAto[1]+j][CubAto[2]+k][CubAto[3]+m];
                     n<=CubLAI_en[CubAto[1]+j][CubAto[2]+k][CubAto[3]+m];n++) {

                  RFSqDi=DistSq(RoSFCo[i][1],RoSFCo[i][2],RoSFCo[i][3],
                                ReCoor[AtReprRes[CubLiA_en[n]]][1],
                                ReCoor[AtReprRes[CubLiA_en[n]]][2],
                                ReCoor[AtReprRes[CubLiA_en[n]]][3]);

/* Take the residue into account only if this distance (between the current
   fragment atom and the current residue-representative atom) is smaller or
   equal to the pseudo-sphere radius */
/* Residue is taken into account only if it is within the cut-off distance
   from the current fragment atom. clangini */
                  if (RFSqDi <= PsSpRa_sq) {
                    ResFlag[CubLiA_en[n]]=1;
/* Loop over the current residue atoms */
                    for (p = FiAtRes[CubLiA_en[n]];p <= LaAtRes[CubLiA_en[n]]; p++) {
                      RFSqDi = DistSq(RoSFCo[i][1],RoSFCo[i][2],RoSFCo[i][3],
                                     ReCoor[p][1],ReCoor[p][2],ReCoor[p][3]);
                      SDFrRe_ps[i][p] = RFSqDi;
                    }
                  }

                }

              }

            }
          }

        } // loop over cubes of pseudo-sphere: x
      } // loop over cubes of pseudo-sphere: y
    } // loop over cubes of pseudo-sphere: z
  } // loop over frag atom

  free_ivector(ResFlag,1,ReReNu);

}
