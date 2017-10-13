#include <stdio.h>
#include "nrutil.h"
#include <string.h>
#include <stdlib.h> /* for exit-fct.*/
#include "funct.h"

#ifndef _STRLENGTH
#define _STRLENGTH 500
#endif

/*void ReInFi(char *InpFil,char *RecFil,int *BSResN,int **BSReNu,
            int *FragNu,char ***FrFiNa,char *TREFiP,double *SphAng,
            int *SphPoN,int *NuRoAx,double *VdWFaB,double *CoDieV,int *CoDieP,
            double *CoGrIn,double *CoGrSi,char *OutFil,double *BuEvFa,
            double **FrMaEn,double *PsSpRa,double *GrSiCu_en,
            int *FiNuMa,double *GrInSo,double *GrSiSo,double *WaMoRa,
            int *NPtSphere,double *DielWa,double *DielRe,char *ReDesoAlg,
            char *DesoMapAcc,char *DesoMapFile,char *FDexe,char *FDdir,
            double *ReSurfDens_apol,double *PtDensFr,double *Sphere_apol,
            double *NCutapolRatio,double *ScaleDeso,
            double *ScaleVDW,double **SimWei,double *SimExp,
            double *SimCut,double **FrMaEn_sd,double *SimExp_sd,double *SimCut_sd,
            int *BSMeNu,int **BSMeAN,double ***BSMeVE,char *CoGrAcc,
            char *CoGrFile,char *EvalEn,char *Solv_typ,
            char *SpPoCh_opt,double *SpPoCh_cent,double *SpPoCh_rad,
            double *SFDeso_fr,double *SFDeso_re,double *SFVWEn,double *SFIntElec,
            int *NuClusMem,double *RedRPV_rp,double*RedRPV_nkvRatio,double *ScMaBump,
            double *MuFaVdWCoff_ap,int *NuLiEnClus,char ***ApPoChoi,
            double *VWGrIn,double *VWGrSi,double *BumpFaCut,char *VWGrAcc,
            char *VWGrFile,int *MaxPosClus,int *PrintLev,
            int *NumbAT,char ***AtTyAr,int **AtENAr,double **VdWRad,
            double **VdWEne,double ***BLAtTy,int *distrPointBSNumb,
	          double ***distrPointBS,double *angle_rmin,double *angle_rmax,
	          double *mult_fact_rmin,double *mult_fact_rmax,char *EmpCorrB,
            char *gc_opt,int *gc_reprke,double *gc_cutclus,double *gc_endifclus,
            double *gc_weighneg,double *gc_weighpos,int *gc_maxwrite,
            char *write_pproc_opt,char *write_pproc_chm_opt,int *CorrFiNumb)*/
void ReInFi(char *InpFil,char *RecFil,int *BSResN,int **BSReNu,
            char *FrFiNa,char *TREFiP,double *SphAng,
            int *SphPoN,int *NuRoAx,double *VdWFaB,double *CoDieV,int *CoDieP,
            double *CoGrIn,double *CoGrSi,char *OutFil,double *BuEvFa,
            double *FrMaEn,double *PsSpRa,double *GrSiCu_en,
            int *FiNuMa,double *GrInSo,double *GrSiSo,double *WaMoRa,
            int *NPtSphere,double *DielWa,double *DielRe,char *ReDesoAlg,
            char *DesoMapAcc,char *DesoMapFile,char *FDexe,char *FDdir,
            double *ReSurfDens_apol,double *PtDensFr,double *Sphere_apol,
            double *NCutapolRatio/**NCutapol*/,double *ScaleDeso,
            double *ScaleVDW,double **SimWei,double *SimExp,
            double *SimCut,double *FrMaEn_sd,double *SimExp_sd,double *SimCut_sd,
            int *BSMeNu,int **BSMeAN,double ***BSMeVE,char *CoGrAcc,
            char *CoGrFile,char *EvalEn,char *Solv_typ,
            char *SpPoCh_opt,double *SpPoCh_cent,double *SpPoCh_rad,
            double *SFDeso_fr,double *SFDeso_re,double *SFVWEn,double *SFIntElec,
            int *NuClusMem,int *NuPosMem /*clangini*/,double *RedRPV_rp,
            double*RedRPV_nkvRatio,double *ScMaBump, /*int *RedRPV_nkv, dey new */
            double *MuFaVdWCoff_ap,int *NuLiEnClus,char *ApPoChoi,
            double *VWGrIn,double *VWGrSi,double *BumpFaCut,char *VWGrAcc,
            char *VWGrFile,int *MaxPosClus,int *PrintLev,
            int *NumbAT,char ***AtTyAr,int **AtENAr,double **VdWRad,
            double **VdWEne,double ***BLAtTy,int *distrPointBSNumb,
	          double ***distrPointBS,double *angle_rmin,double *angle_rmax,
	          double *mult_fact_rmin,double *mult_fact_rmax,char *EmpCorrB,
            char *gc_opt,int *gc_reprke,double *gc_cutclus,double *gc_endifclus,
            double *gc_weighneg,double *gc_weighpos,int *gc_maxwrite,
            char *write_pproc_opt,char *write_pproc_chm_opt,char *write_best_opt,
            char *write_sumtab_opt,char *write_best_sumtab_opt,double **AtWei,
            Parameter &seed_par)
/* This function reads the data of the input (InpFil) and parameters (TREFiP)
   files :
   OutFil  path of the file containing the output informations
   RecFil  name of the receptor file (mol format)
   BSResN  number of residues in the binding site
   BSReNu  numbers of the residues in the binding site
   FragNu  number of fragments
   FrFiNa  names of the files containing the fragments
   CorrFiNumb  number of corrupt mol2-files
   TREFiP  path of the file which contains the parameters
   SphAng  angle (in degrees) of the part of the sphere
   SphPoN  number of desired points on the sphere
   NuRoAx  number of fragment rotations around each axis
   VdWFaB  factor for the checking of the bumps distances
   CoDieV  dielectric value for the coulombic interaction
   CoDieP  dielectric power for the coulombic interaction (0 or 1,1->r-dielect.)
   CoGrIn  increase of the grid for the coulombic interaction
   CoGrSi  grid size for the coulombic interaction
   BuEvFa  factor of bump evaluation
   FrMaEn  maximal total energy for each type of fragment
   PsSpRa  pseudo-sphere radius as cutoff for the evaluation of the energy
   GrSiCu_en  grid size of the cubes for the evaluation of the energy
   FiNuMa  maximal number of mol2 files to be written for each fragment type
   GrInSo  grid increase for the solvation
   GrSiSo  grid size for the solvation
   WaMoRa  water molecule radius (for the solvation)
   NPtSphere # points per sphere in the generation of the S&R surface
   DielWa  water dielectric (for the solvation)
   DielRe  receptor dielectric (for the solvation)
   ReDesoAlg Algoritm for the evaluation of receptor desolvation (coulomb
             approximation or finite difference calculation)
   DesoMapAcc Desolvation map access (write read or nothing)
   DesoMapFile name of the file to be accessed
   FDexe path of the UHBD executable
   FDdir name of the directory where to temporary write the UHBD outputs
   ReSurfDens_apol Point density (in A^-2) on the desolvation surf of the receptor
   PtDensFr Point density (in A^-2) on the desolvation surf of the fragment
   Nsurfpt_re_Deso # points per sphere to generate SAS for apolar vectors
   Sphere_apol limit water molecule radius (for the desolvation surface)
   ScaleDeso scaling factor multiplying the desolvation of the probe sphere
   ScaleVDW scaling factor multiplying the vdW interactions of the probe sphere
   NCutapol number of desolvation vectors
   NCutapol ratio of desolvation vectors
   SimWei  similarity weight factors for GSEAL
   SimExp  similarity exponential factor for GSEAL
   SimCut  similarity cutoff factor for GSEAL
   FrMaEn_sd  maximal total energy for each type of fragment (second)
   SimExp_sd  similarity exponential factor for GSEAL (second)
   SimCut_sd  similarity cutoff factor for GSEAL (second)
   BSMeNu  total number of vector extremities for the metal atoms of the
           binding site
   BSMeAN  atom numbers for the metal atoms of the binding site
   BSMeVE  coordinates of the vector extremities for the metal atoms of the
           binding site
   CoGrAcc  Coulombic grid access (write read or none)
   CoGrFile  name of the file to be accessed for the Coulombic grid
   EvalEn  compute the energy of the following fragment positions (y,n)
   Solv_typ  type of solvation algorithm: fast or slow (f,s)
   SpPoCh_opt  check if the fragment is in the specified sphere (y,n)
   SpPoCh_cent  center of the sphere for the checking of the fragment position
   SpPoCh_rad  radius of the sphere for the checking of the fragment position
   SFDeso_fr  scaling factor for the fragment desolvation energy
   SFDeso_re  scaling factor for the receptor desolvation energy
   SFVWEn  scaling factor for the vdW energy
   SFIntElec  scaling factor for the electrostatic interaction energy
   NuClusMem  number of cluster members to be written in ?_clus_reduc.chm file
   NuPosMem   number of poses to be written in ?_best_pproc.mol2 file //clangini
   RedRPV_rp  vdW radius of the probe in the reducing of the receptor polar
              vectors
   RedRPV_nkv  maximal number of kept vectors in the reducing of the receptor
               polar vectors -> deprecated
   RedRPV_nkvRatio  Ratio  (kept vectors / (total number of donor & acceptor vectors))
                    to reduce number of polar vectors in binding site
   ScMaBump  scaling factor for computing the maximal number of tolerated
             bumps
   MuFaVdWCoff_ap  multiplicative factor used for computing the van der Waals
                   energy cutoff in the seeding of the apolar fragments
   NuLiEnClus  number of lines to be written in the output file for the
               sorted energies and the two clustering procedures
   ApPoChoi  user choice between apolar or polar seeding for each fragment
             type (a,p,n)
   VWGrIn van der Waals grid increase
   VWGrSi  van der Waals grid size
   BumpFaCut  cutoff (maximal vdW energy) for bump checking with fast energy
              evaluation
   VWGrAcc  van der Waals grid access (write read or none)
   VWGrFile  name of the file to be accessed for the van der Waals grid
   MaxPosClus  maximal number of positions to be clustered
   PrintLev  printing level (0,1->preprocess,2->second clustering)
   NumbAT  number of atom types
   AtTyAr  atom types array
   AtENAr  atom element numbers array
   VdWRad  van der Waals radii
   VdWEne  van der Waals energies
   BLAtTy  bond lengths between atom types in the hydrogen bond making
           distrPointBSNumb  number of points in the binding site used for the
           reduction of polar and apolar vectors using an angle criterion
   distrPointBS  distribution of points in the binding site used for
                 the reduction of polar and apolar vectors using an angle criterion
   angle_rmin  angle cutoff for the reduction of vectors
   angle_rmax  angle cutoff for the reduction of vectors
   mult_fact_rmin  multiplicative factor for the reduction of vectors
   mult_fact_rmax  multiplicative factor for the reduction of vectors
   EmpCorrB  empirical correction term (y,n) to the Coulomb field
             approximation for the accurate screened interaction and
	           fragment desolvation energies */
{
    FILE *FilePa,*FilePa2;/* unused variable : ,*FilChk;*/
  /* char StrLin[_STRLENGTH], **FrFiNa_L,**ApPoChoi_L,**AtTyAr_L; */
  char StrLin[_STRLENGTH], /*FrFiNa_L[_STRLENGTH], ApPoChoi_L[2]*/ **AtTyAr_L;
 /* double *FrMaEn_L,Value_nd,*FrMaEn_sd_L,*VdWRad_L,*VdWEne_L,UsVal;*/
  double Value_nd, *VdWRad_L, *VdWEne_L, UsVal;
  int *BSReNu_L,i,NonDef,AtoEle_1,AtoEle_2,j,Atom1,Atom2,dummy1,UsVal2,k;
  int AtEl; //clangini
  //double AtWei_k; // clangini
  char dummyStr[5]; //clangini
/* ----------------------------------------- */
/* ---------- Read the input file ---------- */
/* ----------------------------------------- */

  FilePa=fopen(InpFil,"r");

  CheckFile(InpFil,'r');
//  std::cout << "file: " << InpFil << std::endl;
/* Path of the file containing the parameters */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%s",TREFiP);

  CheckFile(TREFiP,'r');
//  std::cout << "file: " << TREFiP << std::endl;
/* Path of the file containing the receptor coordinates (mol2 format) */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%s",RecFil);

  CheckFile(RecFil,'r'); /*dey */
//  std::cout << "file: " << RecFil << std::endl;
  if(!CheckMol2File(RecFil))
  {
      fclose(FilePa);
      exit(12);
  }


/* Residues of the receptor which are in the binding site */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%d",BSResN);
  BSReNu_L=ivector(1,*BSResN);
  *BSReNu=BSReNu_L;
  for (i=1;i<=*BSResN;i++) {
    fgets_wrapper(StrLin,_STRLENGTH,FilePa);
    sscanf(StrLin,"%d",&BSReNu_L[i]);
  }

/* Read the number and the coordinates of points in the binding site
   used for the reduction of polar and apolar vectors using an angle
   criterion */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%d",distrPointBSNumb);
  *distrPointBS=dmatrix(1,*distrPointBSNumb,1,3);
  for (i=1;i<=*distrPointBSNumb;i++) {
    fgets_wrapper(StrLin,_STRLENGTH,FilePa);
    sscanf(StrLin,"%lf%lf%lf",&(*distrPointBS)[i][1],&(*distrPointBS)[i][2],
           &(*distrPointBS)[i][3]);
  }

/* Read the atom numbers and vector extremities of the metals which are in the
   binding site */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%d",BSMeNu);
  *BSMeAN=ivector(1,*BSMeNu);
  *BSMeVE=dmatrix(1,*BSMeNu,1,3);
  for (i=1;i<=*BSMeNu;i++) {
    fgets_wrapper(StrLin,_STRLENGTH,FilePa);
    sscanf(StrLin,"%d%lf%lf%lf",&(*BSMeAN)[i],&(*BSMeVE)[i][1],&(*BSMeVE)[i][2],
	   &(*BSMeVE)[i][3]);
  }

/* Select only fragments in sphere (y/n, sphere center, sphere radius) */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%s%lf%lf%lf%lf",SpPoCh_opt,&SpPoCh_cent[1],&SpPoCh_cent[2],
         &SpPoCh_cent[3],SpPoCh_rad);

/* Number of fragments / Compute the energy of these positions (y,n)
   Paths of the files which contain the fragments /
   Apolar,polar seeding or no specification (a,p,n)
   two energy cutoffs (in generating fragments,in the 2nd clustering) */
  SkipComLin(FilePa,StrLin);
  /* clangini 2016 */
  /* sscanf(StrLin,"%d%s",FragNu,EvalEn);*/
  sscanf(StrLin, "%s", EvalEn);
  fgets_wrapper(StrLin,_STRLENGTH,FilePa);
  //sscanf(StrLin, "%s%s%lf%lf", FrFiNa_L, ApPoChoi_L, FrMaEn, FrMaEn_sd);
  sscanf(StrLin, "%s%s%lf%lf", FrFiNa, ApPoChoi, FrMaEn, FrMaEn_sd);
  //*FrFiNa = FrFiNa_L;
  //*ApPoChoi = ApPoChoi_L;

  /*FrFiNa_L=cmatrix(1,*FragNu,1,_STRLENGTH);
  FrMaEn_L=vector(1,*FragNu);
  FrMaEn_sd_L=vector(1,*FragNu);
  ApPoChoi_L=cmatrix(1,*FragNu,1,2);
  *FrFiNa=FrFiNa_L;
  *FrMaEn=FrMaEn_L;
  *FrMaEn_sd=FrMaEn_sd_L;
  *ApPoChoi=ApPoChoi_L; */

  /* *CorrFiNumb=0;
  for (i=1;i<=*FragNu;i++) {
    fgets_wrapper(StrLin,_STRLENGTH,FilePa);
    sscanf(StrLin,"%s%s%lf%lf",FrFiNa_L[i]+1,ApPoChoi_L[i]+1,&FrMaEn_L[i],
                             &FrMaEn_sd_L[i]);

    if(!CheckMol2File( FrFiNa_L[i]+1 ))
    {
	fprintf(stderr,"WARNING, skipping file : %s !\n",FrFiNa_L[i]+1);
	--(*FragNu);
	--i;
	++(*CorrFiNumb);
    }

  }*/

  fclose(FilePa);


/* ---------------------------------------------- */
/* ---------- Read the parameters file ---------- */
/* ---------------------------------------------- */

  FilePa=fopen(TREFiP,"r");

  if(FilePa == NULL ) /*dey new*/
  {
    fprintf (stderr,"WARNING, error occured while trying to open parameter file , exiting !\n");
    fclose(FilePa);
    exit (12);
  }

/* Check file tag in parameter file*/
/*  fgets_wrapper(StrLin,_STRLENGTH,FilePa);
  if(strncmp("#seed.param v3.3.6",StrLin,18)!=0)
  {
      fprintf (stderr,"WARNING, missing tag in parameter file : #seed.param v3.3.6 , exiting !\n");
      exit(12);
  } clangini */
/* Check file tag in parameter file -> updated to v4.0 clangini*/
  fgets_wrapper(StrLin,_STRLENGTH,FilePa);
  if(strncmp("#seed.param v4.0",StrLin,16)!=0)
  {
      fprintf (stderr,"WARNING, missing tag in parameter file : #seed.param v4.0 , exiting !\n");
      exit(12);
  }

/* CLANGINI 2016 Moved some values from seed.inp to seed.par and changed the order of parameters */
/* Dielectric value of the solute (receptor and fragment) */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%lf",DielRe);
  *CoDieV=(*DielRe)*1.43;      /* <-- Coulombic int. : dielectric value */

/* old : Maximal number of conserved vectors : polar case / apolar case */
/* Ratio of conserved vectors to retain : polar case / apolar case */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%lf%lf",RedRPV_nkvRatio,NCutapolRatio);

/* Write *_pproc* files in postprocess mode (y/n) / charmm files */
  //SkipComLin(FilePa,StrLin);
  //sscanf(StrLin,"%s%s",write_pproc_opt,write_pproc_chm_opt);

/* Write *_pproc*.mol2 file and write *_best_pproc*.mol2 file
in postprocess mode (y/n   y/n) */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%s%s",write_pproc_opt,write_best_opt);

/* Write summary tables *_pproc*.dat and *_best_pproc*.dat
in postprocess mode (y/n   y/n) */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%s%s",write_sumtab_opt,write_best_sumtab_opt);

/* Number of cluster members and total number of poses to be
written in output file*/ //clangini
// cluster members for the first clustering!
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%d%d",NuClusMem, NuPosMem);

/* Path of the file containing the output informations */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%s",OutFil);

/* Coulombic grid access (write read or none) / name of the file */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%s%s",CoGrAcc,CoGrFile);

  if(CoGrAcc[0] != 'n') /* new : always calculate from scratch - don't read or write */
    CheckFile(CoGrFile,CoGrAcc[0]); /*dey new*/

/* van der Waals grid access (write read or none) / name of the file */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%s%s",VWGrAcc,VWGrFile);

  if(VWGrAcc[0] != 'n') /* new : always calculate from scratch - don't read or write */
    CheckFile(VWGrFile,VWGrAcc[0]); /*dey new*/

/* D map grid access (write read or none) / name of the file */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%s%s",DesoMapAcc,DesoMapFile);

  if(DesoMapAcc[0] != 'n') /* new : always calculate from scratch - don't read or write */
    CheckFile(DesoMapFile,DesoMapAcc[0]); /*dey new*/

/* Read the scaling factor used for computing the maximal number of tolerated
   bumps, the factor for the checking of the bumps distances and the factor of
   bump evaluation */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%lf%lf%lf",ScMaBump,VdWFaB,BuEvFa);

/* Cutoff (maximal vdW energy) for bump checking with fast energy evaluation */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%lf",BumpFaCut);

/* Read the sphere angle and number of desired points on the part of the
   sphere */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%lf%d",SphAng,SphPoN);

/* Read the number of fragment rotations around each axis */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%d",NuRoAx);

/* Read the parameters used for the reduction of rec polar
   and apolar vectors */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%lf%lf%lf%lf",angle_rmin,angle_rmax,
         mult_fact_rmin,mult_fact_rmax);

/* Reducing of the rec. polar vect. : vdW radius probe */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%lf",RedRPV_rp);

/* Read for the coulombic interaction : dielectric power, grid increase and
   grid size */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%d%lf%lf",CoDieP,CoGrIn,CoGrSi);

/* Read VWGrIn and VWGrSi */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%lf%lf",VWGrIn,VWGrSi);

/* Read the pseudo-sphere radius as cutoff and grid size of the cubes for the
   evaluation of the vdW energy */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%lf%lf",PsSpRa,GrSiCu_en);

/* Multiplic. factor for vdW energy cutoff in seeding apolar fragments */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%lf",MuFaVdWCoff_ap);

/* Maximal number of mol2 files to be written for each fragment type */
/* SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%d",FiNuMa);*/ /* Choice 1 is hard-coded.
                                clangini April 2017
                                simplification of seed.par */
  *FiNuMa = 1; // clangini

/* Solvation algorithm type (f_ast,s_low,b_oth,p_ostprocess) */
/*  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%s",Solv_typ);*/ /* Choice 'p' is hard-coded.
                                    clangini April 2017
                                    simplification of seed.par */

  strcpy(Solv_typ,"p"); //clangini

/* Empirical correction term (y,n) to the Coulomb field approximation
   for the accurate screened interaction and fragment desolvation energies */
/*  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%s",EmpCorrB);*/ /* Choice 'y' is hard-coded.
                                    clangini April 2017
                                    simplification of seed.par */
  strcpy(EmpCorrB, "y"); //clangini

/* Grid increase and grid size for the solvation */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%lf%lf",GrInSo,GrSiSo);

/* Water molecule radius, # points per sphere in the generation of
   the S&R surface, water dielectric for the solvation */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%lf%d%lf",WaMoRa,NPtSphere,DielWa);

/* CLANGINI 2016 We do not use UHBD any more but we keep the setup for retro-compatibility */
/* Solvation : Receptor desolvation algorithm (co or fd) */
/*  SkipComLin(FilePa,StrLin); */
/*  sscanf(StrLin,"%s",ReDesoAlg); */
  sprintf(ReDesoAlg,"%s","co\0");
/* Solvation FD executable and directory */
/*  SkipComLin(FilePa,StrLin); */
/*  sscanf(StrLin,"%s%s",FDexe,FDdir); */
  *FDexe = '\0';
  *FDdir = '\0';
/* CLANGINI 2016 END */

/* Solvation apolar vectors: point densities for rec and frag, probe radius,
                             scaling factors for desolvation and for vdW */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%lf%lf%lf%lf%lf",ReSurfDens_apol,PtDensFr,Sphere_apol,
                                  ScaleDeso,ScaleVDW);

/* Scaling factors for vdW energy, for electrostatic interaction energy and
   for the receptor/fragment desolvation energies */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%lf%lf%lf%lf",SFVWEn,SFIntElec,SFDeso_re,SFDeso_fr);

/*clangini 2016 START*/
/* Similarity weight factors (150 atom elements) for GSEAL */
/*  for (i=1;i<=150;i++) {
    for (j=1;j<=150;j++) {
      if (i==j)
        SimWei[i][j]=1.0;
      else
        SimWei[i][j]=0.0;
    }
  } we want to use atom element 0 for the lone pair
    -> index of SimWei has to start from 0  clangini*/
/* Similarity weight factors (151 atom elements) for GSEAL  clangini*/
  for (i=0;i<=150;i++) {
    for (j=0;j<=150;j++) {
      if (i==j)
        SimWei[i][j]=1.0;
      else
        SimWei[i][j]=0.0;
    }
  }
/*clangini 2016 END*/
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%d",&NonDef);
  if (NonDef) {
    for (i=1;i<=NonDef;i++) {
      fgets_wrapper(StrLin,_STRLENGTH,FilePa);
      sscanf(StrLin,"%d%d%lf",&AtoEle_1,&AtoEle_2,&Value_nd);
      SimWei[AtoEle_1][AtoEle_2]=Value_nd;
      if (AtoEle_1!=AtoEle_2)
        SimWei[AtoEle_2][AtoEle_1]=Value_nd;
    }
  }

/* Similarity exponential and cutoff factors for GSEAL */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%lf%lf",SimExp,SimCut);

/* Similarity exponential and cutoff factors for GSEAL (second) */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%lf%lf",SimExp_sd,SimCut_sd);

/* Maximal number of positions to be clustered */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%d",MaxPosClus);

/* clangini 2016. Hard-code of 'n' option for FFLD, but keep all the variables for future
   modifications (Compare with UHBD options) */
/* build geometrical centers for FFLD */
/*  SkipComLin(FilePa,StrLin); */
/*  sscanf(StrLin,"%s%d%d%lf%lf%lf%lf",gc_opt,gc_maxwrite,gc_reprke,gc_cutclus,
                gc_endifclus,gc_weighneg,gc_weighpos); */
  gc_opt = (char *)"n"; //added explicit cast. clangini
/* clangini 2016 end */

/* Number of lines to be written in the output file for the sorted energies
   and the two clustering procedures;
   Printing level (0,1->preprocess,2->second clustering) */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%d%d",NuLiEnClus,PrintLev);
/* Parameters for MC run */
  SkipComLin(FilePa, StrLin);
  sscanf(StrLin, "%c", &(seed_par.do_mc));

  if (seed_par.do_mc == 'y'){
    SkipComLin(FilePa, StrLin);
    sscanf(StrLin, "%lf", &(seed_par.mc_temp));
    SkipComLin(FilePa, StrLin);
    sscanf(StrLin, "%lf", &(seed_par.mc_max_tran_step));
    SkipComLin(FilePa, StrLin);
    sscanf(StrLin, "%lf", &(seed_par.mc_max_rot_step));
    SkipComLin(FilePa, StrLin);
    sscanf(StrLin, "%d", &(seed_par.mc_niter));
    SkipComLin(FilePa,StrLin);
    sscanf(StrLin, "%d", &(seed_par.mc_rand_seed));
  }
/* CLANGINI 2016 END */

/* Read NumbAT AtTyAr AtENAr VdWRad VdWEne */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%d",NumbAT);
  AtTyAr_L=cmatrix(1,*NumbAT,1,7);
  *AtENAr=ivector(1,*NumbAT);
  VdWRad_L=dvector(1,*NumbAT);
  VdWEne_L=dvector(1,*NumbAT);
  *AtTyAr=AtTyAr_L;
  *VdWRad=VdWRad_L;
  *VdWEne=VdWEne_L;
  for (i=1;i<=*NumbAT;i++) {
    fgets_wrapper(StrLin,_STRLENGTH,FilePa);
    sscanf(StrLin,"%d%s%d%lf%lf",&dummy1,&AtTyAr_L[i][1],&(*AtENAr)[i],
           &VdWRad_L[i],&VdWEne_L[i]);
  }

/* Read and construct BLAtTy */
  SkipComLin(FilePa,StrLin);
  *BLAtTy=dmatrix(1,*NumbAT,1,*NumbAT);
  sscanf(StrLin,"%lf",&UsVal);
  for (i=1;i<=*NumbAT;i++) {
    for (j=1;j<=*NumbAT;j++) {
      (*BLAtTy)[i][j]=UsVal;
    }
  }

  fgets_wrapper(StrLin,_STRLENGTH,FilePa);
  sscanf(StrLin,"%d",&UsVal2);
  for (k=1;k<=UsVal2;k++) {
    fgets_wrapper(StrLin,_STRLENGTH,FilePa);
    sscanf(StrLin,"%d%d%lf",&Atom1,&Atom2,&UsVal);
    for (i=1;i<=*NumbAT;i++) {
      for (j=1;j<=*NumbAT;j++) {
        if ((((*AtENAr)[i]==Atom1)&&((*AtENAr)[j]==Atom2))||
            (((*AtENAr)[i]==Atom2)&&((*AtENAr)[j]==Atom1)))
          (*BLAtTy)[i][j]=UsVal;
      }
    }
  }

  fgets_wrapper(StrLin,_STRLENGTH,FilePa);
  sscanf(StrLin,"%d",&UsVal2);
  for (k=1;k<=UsVal2;k++) {
    fgets_wrapper(StrLin,_STRLENGTH,FilePa);
    sscanf(StrLin,"%d%d%lf",&Atom1,&Atom2,&UsVal);
    (*BLAtTy)[Atom1][Atom2]=UsVal;
    (*BLAtTy)[Atom2][Atom1]=UsVal;
  }

/* Print BLAtTy in a separate file */
  FilePa2=fopen("./outputs/length_hb.gen","w");
  fprintf(FilePa2,"File containing the bond lengths between atom types ");
  fprintf(FilePa2,"in the hydrogen bond \n");
  fprintf(FilePa2,"making procedure\n\n");
  for (i=1;i<=*NumbAT;i++) {
    for (j=i;j<=*NumbAT;j++) {
      fprintf(FilePa2,"%-7s  %-7s    %-11.5f\n",&AtTyAr_L[i][1],
              &AtTyAr_L[j][1],(*BLAtTy)[i][j]);
    }
  }
  fclose(FilePa2);
  /*clangini START Read atomic weights */
  SkipComLin(FilePa,StrLin);
  sscanf(StrLin,"%d",&UsVal2); // Number of elements (without element 0)
  *AtWei=dvector(0,UsVal2);
  for (k=0;k<=UsVal2;k++) {
    fgets_wrapper(StrLin,_STRLENGTH,FilePa);
    sscanf(StrLin,"%s%d%lf",dummyStr,&AtEl,&(*AtWei)[k]);
    if (k != AtEl){
      std::cerr << "Atom is missing in the atomic weight list! Exit program!" << std::endl;
      exit(12);
    }
  }
  /*clangini END*/
  fclose(FilePa);

}



void SkipComLin(FILE *FilePa,char *StrLin)
/* This function skips the comment lines starting with a # */
{
  fgets_wrapper(StrLin,_STRLENGTH,FilePa);
  while (!(strncmp(StrLin,"#",1)))
    fgets_wrapper(StrLin,_STRLENGTH,FilePa);
}

int CheckFile(char*name,char type) /* char * type)*/
/* this function checks the presence of a file */
{
  FILE *FilChk;
  /* char * tmp = (char *) malloc (2); */
  /* tmp[0] = type; */
    /* FilChk= (type == 'r') ? fopen(name,'r') : fopen(name,'w'); */

  if ( type =='r' && (FilChk = fopen ( name, "r" ) ) == NULL )
    {
	fprintf(stderr,"WARNING, could not read file %s, exiting \n",name);
	exit (12);
    }
  else if ( type =='w' && (FilChk = fopen ( name, "w" ) ) == NULL )
    {
      fprintf(stderr,"WARNING, could not write file %s, exiting \n",name);
      exit (12);
    }
  else if (type !='r' && type !='w')
    {
      fprintf(stderr,"WARNING, file access should be either 'r' or 'w' for %s, and not %c\n",name,type);
      exit(12);
    }
  else
    {
      fclose(FilChk);
    }
  /* } */

  /* free(tmp); */

  /* fclose(FilChk); */
  return 1;
}


/* The CheckMol2File() will now be used only to check the receptor file, not the fragment file */
int CheckMol2File(char * name)
{
    /*
      small mol2-file integrity check function
      -> scans for TRIPOS tags
    */

    char StrLin[200];
    FILE *FilChk;
    FilChk=fopen( name,"r");
    if(FilChk == NULL ) /*dey new*/
    {
	fprintf(stderr,"WARNING, could not read file %s ! \n",name);
	return 0;
    }



    while(FilChk!=NULL && !feof(FilChk) )
    {
	fgets_wrapper(StrLin,200,FilChk);

	if(!strncmp("@<TRIPOS>MOLECULE",StrLin,17))
	{
	    break;
	}
    }

    if(FilChk==NULL || feof(FilChk) ||  ferror(FilChk))
    {
	fprintf(stderr,"WARNING, could not find @<TRIPOS>MOLECULE-tag in file %s ! \n",name);
	return 0;
    }


    while(FilChk!=NULL && !feof(FilChk) && !ferror(FilChk)  )
    {
	fgets_wrapper(StrLin,200,FilChk);

	if(!strncmp("@<TRIPOS>ATOM",StrLin,13))
	{
	    break;
	}
    }

    if(FilChk==NULL || feof(FilChk) ||  ferror(FilChk)  ) /* || FilChk==EOF)*/
    {
	fprintf(stderr,"WARNING, could not find @<TRIPOS>ATOM-tag in file %s ! \n",name);
	return 0;
    }

    while(FilChk!=NULL && !feof(FilChk) &&  !ferror(FilChk)  )
    {
	fgets_wrapper(StrLin,200,FilChk);
	if(!strncmp("@<TRIPOS>BOND",StrLin,13))
	{
	    break;
	}
    }

    if(FilChk==NULL || feof(FilChk) ||  ferror(FilChk)  )
    {
	fprintf(stderr,"WARNING, could not find @<TRIPOS>BOND-tag in file %s ! \n",name);
	return 0;
    }

/* clangini 2016: SUBSTRUCTURE-tag is still used for the receptor */
/* clangini 2017: SUBSTRUCTURE-tag is no longer used for the receptor */
  // while(FilChk!=NULL && !feof(FilChk) && !ferror(FilChk)  )
  // {
	//    fgets_wrapper(StrLin,200,FilChk);
  //    if(!strncmp("@<TRIPOS>SUBSTRUCTURE",StrLin,21))
  //    {
  //      break;
  //    }
  //  }
  //
  // if(FilChk==NULL || feof(FilChk) ||  ferror(FilChk)  )
  // {
	//    fprintf(stderr,"WARNING, could not find @<TRIPOS>SUBSTRUCTURE-tag in file %s ! \n",name);
	//     return 0;
  // }
  /* clangini 2017 introduced check for ALT_TYPE */
  // while(FilChk!=NULL && !feof(FilChk) && !ferror(FilChk)  )
  // {
	//    fgets_wrapper(StrLin,200,FilChk);
  //    if(!strncmp("@<TRIPOS>ALT_TYPE",StrLin,17))
  //    {
  //      break;
  //    }
  //  }
  //
  // if(FilChk==NULL || feof(FilChk) ||  ferror(FilChk)  )
  // {
	//    fprintf(stderr,"WARNING, could not find @<TRIPOS>ALT_TYPE-tag in file %s ! \n",name);
	//     return 0;
  // }

    fclose(FilChk);
    return 1;
}
