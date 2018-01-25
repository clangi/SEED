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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"
#include <sys/time.h>
/* CLANGINI 2016 */
//#include <boost/math/constants/constants.hpp>
#include <sys/stat.h> /* to check for directory existence (consider changing using boost/system instead) */
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#ifdef USE_QUATERNION_ROTATION
#include <quaternion.h>
#endif

#if  __cplusplus > 199711L
#include <unordered_map> //C++11 feature. Implemented as hash table
#else
#include <map>
#endif
/* CLANGINI 2016 END */
#include "funct.h"

#include <float.h>


#define _STRLENGTH 500


/*
   allocation block size , avoid "do_not_touch_ener"
 */
#define _ALLOCSIZE 100000

int main(int argc,char *argv[])
  /* Main part of the program :
     FPaOut  pointer for the output file
     SeFrCo  seeded fragment coordinates
     RoSFCo  coordinates of the rotated seeded fragment
     WriPat  path for the output file
     WriNam  name to be written in the output file
     SFWrNu  number of written seeded fragments
     AnglRo  fragment rotation angle
     TREFiP  path of the file which contains the parameters
     ReVdWR  receptor van der Waals radii
     ReVdWE  receptor van der Waals energies
     FrVdWR  fragment van der Waals radii
     FrVdWE  fragment van der Waals energies
     ReVdWE_sr  square root of receptor van der Waals energies
     FrVdWE_sr  square root of fragment van der Waals energies
     FrAcce  fragment accepted or not (1 -> yes, 0 -> no)
     BumpMa  maximum number of tolerated bumps
     NuRoAx  number of fragment rotations around each axis
     RecHyN  receptor hydrogen number in the current vector (0 if no hydrogen)
     FraHyN  fragment hydrogen number in the current vector (0 if no hydrogen)
     FrMaEn  maximal total energy for each type of fragment
     FrNuTo  total number of generated fragments of the current type
     FrNuPB  number of generated fragments of the current type that passed the
     bump checking
     VWEnEv_ps  evaluation of the van der Waals energy (pseudo-sphere approach)
     ToNFrP  total number of fragment positions
     FrCoPo  vector of pointers pointing to matrices of fragment coordinates
     FraEqu  equivalence between the fragments
     ClusNu  number of fragment clusters
     FraSim_no  normalization for the similarity number between two fragments
     Dist_sq  squared distance
     ClusIn  indexes for the cluster list done with GSEAL
     ClusLi  cluster list done with GSEAL
     FraSim_no_sd  normalization for the similarity number between two fragments
     for the second GSEAL
     FraEqu_sd  equivalence between the fragments for the second GSEAL
     PrintVa  variable used in the printing process
     ClusVa_1  variable used in the clustering process
     ClusVa_2  variable used in the clustering process
     ClusVa_3  variable used in the clustering process
     ClusLi_sd  cluster list done with the second GSEAL
     ClusLi_sd_01  cluster indexation for the second GSEAL (1->all the represen-
     tatives,0->in the clusters)
     FilChk  file existence checking
     ReAtEl_nu  atom element numbers array of the receptor
     FrAtEl_nu  atom element numbers array of the fragment
     HybReAt  hybridization state of the receptor atoms (2 -> sp2, 3 -> sp3,
     0 -> no hybridization state found)
     HybFrAt  hybridization state of the fragment atoms (2 -> sp2, 3 -> sp3,
     0 -> no hybridization state found)
     ReApAt  receptor atom which is the origin of the vector used for the
     seeding of apolar fragments
     FrApAt  fragment atom which is the extremity of the vector used for the
     seeding of apolar fragments
     apol_Vect_re  receptor vectors used for the seeding of apolar fragments
     apol_Vect_fr  fragment vectors used for the seeding of apolar fragments
     FrMinC  minimal coordinates of the fragment in each direction (x,y,z)
     FrMaxC  maximal coordinates of the fragment in each direction (x,y,z)
     SpPoCh_bool  fragment position in the specified sphere (1->yes,0->no)
     FrNuSp  number of fragments kept in the specified sphere
     SpPoCh_gc  geometric center of the fragment
     SpPoCh_dist  squared distance between the geometric center of the fragment
     and the center of the specified sphere
     SpPoCh_rad_sq  squared radius of the sphere for the checking of the fragment
     position
     SFWrNu_ar  array of the finally kept number of fragment positions for each
     fragment type
     ClusLi_sd_01_reduc  same comment as for ClusLi_sd_01 but keeping only a
     reduced number of representatives
     ClusVa_4  variable used in the clustering process
     ReAtoTyp_nu  receptor atom type numbering
  FrAtoTyp_nu  fragment atom type numbering
  VdWCoff_ap  van der Waals energy cutoff in the seeding of apolar fragments
  FrNuVdW_ap  number of fragments that passed vdW energy checking (only for
      apolar fragments)
  FileTemp  temporary file used to write the energies during the generation
  of the fragment positions (before the sorting)  -> removed Dey
ClusVa_wr  variable used in the clustering process (for writing)
  ChkExit  1 if any problem (program exits), 0 if not
  FrCoNu  number of conformations of the current fragment type
  Ind_num_cn  index number of the current conformation
  FrCoor_clus  coordinates of a fragment for computing the normalization
  in the clustering procedure
  ConfArr  array of conformations numbers
  FraSim_no_max  normalization (maximum) for the clustering procedure
  ReVdWR_p3  power of 3 of the receptor van der Waals radii
  FrVdWR_p3  power of 3 of the fragment van der Waals radii
  VW_f  van der Waals energy (fast algorithm)
  VW_s  van der Waals energy (slow algorithm)
  In_f  electrostatic interaction energy (fast algorithm)
  In_s  electrostatic interaction energy (slow algorithm)
  Dr_f  desolvation of the receptor (fast algorithm)
  Dr_s  desolvation of the receptor (slow algorithm)
  Df_f  desolvation of the fragment (fast algorithm)
  Df_s  desolvation of the fragment (slow algorithm)
  To_f  total energy (fast algorithm)
To_s  total energy (slow algorithm)
  Index_both  index array for sorting lists if both (SLOW and FAST) methods are used
  TotEn_both  total energy array for sorting lists if both (SLOW and FAST) methods are used
  ClusLi_sd_pproc  cluster indexation for the second GSEAL for the postprocessing
  (2->representative of the first GSEAL,1->representative of the
   second GSEAL(only a reduced number),0->in the clusters)
  NuSdClKe  number of second clusters which are kept (energy of the cluster
      representative below a user-given cutoff)
  NuPosSdCl  total number of kept positions which are members of the kept
second clusters (representative included)
  FrPosAr_pproc  array containing the number of the fragment position for the
  postprocessing
  SdClusAr_pproc  array containing the number of the second cluster to which
  belongs the fragment position for the postprocessing
  TotEnSdClus_pproc  array containing the total energy (slow method) of the
  second cluster (kept positions) for the postprocessing
  IntVar1  integer variable used as an index
  IntVar2  integer variable used as an index
  Index_pproc  array of indexes for the postprocessing (sorting step)
  FrPosAr_sort  array containing the number of the fragment position for the
    postprocessing after sorting
  SdClusAr_sort  array containing the number of the second cluster to which
    belongs the fragment position for the postprocessing after
    sorting (new order)
  FlagAr  flag array
  ReprSdClAr  array containing the representatives of the conserved second
    clusters after postprocessing
  CluIndex_sort  array to go from cluster index number (SdClusAr_pproc) to cluster
    index number after inter- and intra-cluster sorting (SdClusAr_sort)
  SFWrNu_init  initial value of total SFWrNu
  gc_opt  option for the calculation of the FFLD geometrical centers
  write_pproc_opt  option for writing *_pproc* files
  write_best_opt  option for writing *_best_pproc* files
  MaxGrIn  maximal grid increase
ChkInGrid  check if ligand lies in grid(s) (1->yes,0->no)
  vdwErr,clbErr print "out of grid"-warnings only once per fragment
TotFra fragment counter (both sane and failed fragments). For the sane only, CurFra will be used. (clangini)
  */

{


  FILE *FPaOut,*FilChk,*FilePa;
  char *InpFil,RecFil[_STRLENGTH],/* *FrFiNa*/ FrFiNa[_STRLENGTH],FragNa[_STRLENGTH],**FrAtEl,**FrAtTy,**FrSyAtTy,
       **FrBdTy,**ReAtEl,**ReAtTy,WriPat[_STRLENGTH],/*WriNam[_STRLENGTH],clangini*/
       TREFiP[_STRLENGTH],**AtTyAr,OutFil[_STRLENGTH],FrSubN[20],FrSubC[20],StrLin[_STRLENGTH],
       CoGrAcc[2],CoGrFile[_STRLENGTH],EvalEn[2],SpPoCh_opt[2],**SubNa,
       FrFiNa_out[_STRLENGTH],/* *ApPoChoi,*/ ApPoChoi[2], VWGrAcc[2],VWGrFile[_STRLENGTH],EmpCorrB[2],
       gc_opt[2],write_pproc_opt[2],write_pproc_chm_opt[2], TabFil[_STRLENGTH]/*clangini*/,
       write_sumtab_opt[2], write_best_opt[2], write_best_sumtab_opt[2];
  int i,*BSReNu,BSResN/*,FragNu cl*/,FrAtNu,FrBdNu,CurFra,**FrBdAr,*AliHyd,ReAtNu/*,TotFra clangini*/,
      ReBdNu,ReReNu,*ReResN,ReAcNu,ReDoNu,*ReDATy,*ReDAAt,*ReHydN,FrAcNu,
      FrDoNu,*FrDATy,*FrDAAt,*FrHydN,j,k,SFWrNu,l,NumbAT,CubNum[4],*CubLiA,
      FrAcce,BumpMa,SphPoN,NuRoAx,RecHyN,FraHyN,CoDieP,CoGPoN[4],VWGPoN[4],
      FrNuTo,FrNuPB,***CubFAI,***CubLAI,CubNum_en[4],***CubFAI_en,***CubLAI_en,
      *CubLiA_en,PsSpNC,***PsSphe,*Index_ro,ToNFrP,FiNuMa,i1,i2,i3,
      *FrAtEl_nu,*FraEqu,ClusNu,*ClusIn,*ClusLi,ClusVa_1,ClusVa_2,*FraEqu_sd,
      PrintVa,ClusVa_3,*ClusLi_sd,*ClusLi_sd_01,*AtENAr,*ReAtEl_nu,BSAtNu,
      *BSAtLi,**ReBdAr,*HybReAt,*HybFrAt,BSMeNu,*BSMeAN,NonAliHy,Nx,Ny,Nz,
      SpPoCh_bool,FrNuSp,/* *SFWrNu_ar, clangini*/SFWrNu_ar,*UndisAt_fr,NuClusMem,
      NuPosMem /*clangini*/,ClusVa_4,
      *ClusLi_sd_01_reduc,/*unused variable RedRPV_nkv,*/ *PolVect_rec_01,*ReAtoTyp_nu,
      *FrAtoTyp_nu,FrNuVdW_ap,NuLiEnClus,ClusVa_wr,*FiAtRes,*LaAtRes,
      *AtReprRes,*LiChResEn,NuChResEn,ChkExit,FrCoNu,Ind_num_cn,*ConfArr,
      /*unused variables *Index_both, :LoopTe,*/*ClusLi_sd_pproc,NuSdClKe,NuPosSdCl,
      *FrPosAr_pproc,
       *SdClusAr_pproc,IntVar1,IntVar2,*Index_pproc,*FrPosAr_sort,
       *SdClusAr_sort,*FlagAr,*ReprSdClAr,MaxPosClus,SFWrNu_init,PrintLev,
       ToNFrP_ap,ij,distrPointBSNumb,gc_reprke,gc_maxwrite,ChkInGrid; /*,CorrFiNumb; clangini*/

      double **FrCoor,*FrPaCh,**ReCoor,*RePaCh,**ReVeCo,**FrVeCo,**SeFrCo,
      /*AnglRo,*/**RoSFCo,*VdWRad,*VdWEne,*ReVdWR,*ReVdWE,*FrVdWR,*FrVdWE,
      LaVdWR,ReMaxC[4],ReMinC[4],SphAng,VdWFaB,BSMaxC[4],BSMinC[4],CoDieV,
      CoGrIn,CoGrSi,***CoGrRP,BuEvFa,*ReVdWE_sr,*FrVdWE_sr,
      FrMaEn,PsSpRa,GrSiCu_en,VWEnEv_ps,
      *Coulo_ro,*Vande_ro,*TotEn_ro,*TotEn_ro_cp,***FrCoPo,**SimWei,SimExp,
      SimCut,FraSim,*FraSim_no,Dist_sq,FrMaEn_sd,SimExp_sd,SimCut_sd,
      *FraSim_no_sd,**BSMeVE,
      SpPoCh_cent[4],SpPoCh_rad,SpPoCh_gc[4],SpPoCh_dist,
      SpPoCh_rad_sq,SFDeso_fr,**SDFrRe_ps,SFDeso_re,SFVWEn,SFIntElec,
      RedRPV_rp,**BLAtTy,ScMaBump,MuFaVdWCoff_ap,VdWCoff_ap,*TotChaRes,
      **SDFrRe_ps_elec,**ChFrRe_ps_elec,**FrCoor_clus,FraSim_no_max,
      VW_f,VW_s,In_f,In_s,Dr_f,Dr_s,Df_f,Df_s,To_f,To_s,
      *VW_f_ro,*VW_s_ro,*In_f_ro,*In_s_ro,*Dr_f_ro,*Dr_s_ro,*Df_f_ro,
      *Df_s_ro,*To_f_ro,*To_s_ro,***VWGrRP_at,***VWGrRP_re,
      *ReVdWR_p3,*FrVdWR_p3,VWGrIn,VWGrSi,BumpFaCut/*unused variable :,*TotEn_both*/,
      *TotEnSdClus_pproc,**distrPointBS,angle_rmin,angle_rmax,
      mult_fact_rmin,mult_fact_rmax,gc_cutclus,gc_endifclus,gc_weighneg,
      gc_weighpos,MaxGrIn,RedRPV_nkvRatio,NCutapolRatio/*newdey*/;
  /* struct tms timevar; */
  struct timeval time_0,time_1,time_2,time_3,time_4,time_5,time_6,time_7,time_8,time_9,time_10,
                 time_11,time_12,time_13,time_14;
  /* long time_1,time_2,time_3,time_4,time_5,time_6,time_7,time_8,time_9,time_10, */
  /*      time_11,time_12,time_13,time_14; */
  time_t runtime;

  char DesoMapAcc[2],DesoMapFile[_STRLENGTH],RecFilPDB[_STRLENGTH],FDexe[_STRLENGTH],FDdir[_STRLENGTH],
       ***GridMat,ReDesoAlg[3],Solv_typ[2],FD[3],CO[3];
  int iat,NPtSphere,NGridx,NGridy,NGridz,*nsurf_re,*nsurf_fr,*pointsrf_fr,
      Nsurfpt_fr,*pointsrf_re,*iatsrf_fr,*iatsrf_re,Nsurfpt_re,
      *nsurf_fr_apol,*pointsrf_fr_apol,Nsurfpt_fr_apol,*iatsrf_fr_apol,
      *nsurf_fr_deso,*pointsrf_fr_deso,Nsurfpt_fr_deso,*iatsrf_fr_deso,
      *iatsrf_re_deso,*nsurf_re_deso,*pointsrf_re_deso,Nsurfpt_re_deso,
      /*unused variable :Nsurfpt_BS,*/isph,nxminBS,nyminBS,nzminBS,nxmaxBS,nymaxBS,nzmaxBS,
      Napol_Vect_re,NPtSphereMax,Nsurfpt_re_apol,*ReApAt,*FrApAt,
      Nsurfpt_re_apol_BS,/*NCutapol,*/NPtSphereMax_Fr,***nlist_re,
      **list_re,MaxNeigh,nstep_re,midg_re,nn,
      vdwErr,clbErr;
  double FrMinC[4],FrMaxC[4],**Frdist2,**apol_Vect_re,**apol_Vect_fr,ReRmax;
  double *ReRad,*ReRadOut,*ReRad2,*ReRadOut2,*FrRad,*FrRadOut,*FrRad2,
         *ReSelfVol,***DeltaPrDeso,GrSiSo,GrInSo,WaMoRa,Sphere_apol,
         DielWa,DielRe,*ReSurf_apol,*ReSurf_deso,*FrSurf_deso,*FrRadOut2,
         *vect,Tr[4],U1[4][4],U2[4][4],ReDesoElec,FrDesoElec,
         Kelec,Ksolv,UnitVol,ReFrIntElec,
         *XGrid,*YGrid,*ZGrid,ReSurfDens_apol,PtDensFr,Rmin,
         pi4,ScaleDeso,ScaleVDW,FrSolvEn,SurfDens_deso,
         rslop,rscale_re,FrRmax,corr_re_desoco,corr_re_desofd,corr_fr_deso,
         corr_scrint,corr_fast_deso,*ReSelfVol_corrB;
  struct point Min,Max,*surfpt_re,*surfpt_fr,*surfpt_fr_apol,
               *surfpt_fr_deso_orig,*surfpt_fr_deso,*surfpt_ex,*surfpt_re_deso,
               cbmid_re;
  /* CLANGINI 2016 */
  std::vector<int> PosToRem;
  struct stat DirExist;
  char TabLin[_STRLENGTH]; //clangini
  std::ofstream TabOutStream; //clangini
  int HeAtCo = 0; // HeAtCo = Heavy Atom Count clangini
  double *AtWei; // list of atomic weights clangini
  double MolWei = 0.0; // Molecular Weight clangini
  int *CluIndex_sort; // Array Clu index number -> sorted Clu index number clangini
  std::string AlTySp;
  std::string FragNa_str; //C++ string equivalent to FragNa
#if  __cplusplus > 199711L
  std::unordered_map<std::string, int> FragNa_map;
#else
  std::map<std::string, int> FragNa_map;
#endif
  char buff[12]; // buffer for string manipulation
  bool align_flag; //align fragment
  double FrAlRef_m[3][3], *FrAlRef_rows[3], **FrAlRef,
        FrAlSet_m[3][3], *FrAlSet_rows[3], **FrAlSet;
  double AnglRo; //before it was a float. clangini
  double **FrCoor_NoAlign;
#ifdef USE_QUATERNION_ROTATION //clangini
  double SeFrAx[4]; //axis for fragment rotation
  Quaternion<double> q_SeFr; //calls default constructor
#endif
  /* CLANGINI 2016 END*/
#ifdef OMP
  char * numThreads;
#endif

  /*
     allocation in block size  -> to do realloc increase
   */

  int currentsize = _ALLOCSIZE; /* SEED starts counting at 1... */
  gettimeofday(&time_0,NULL);

  if(argc <= 1 )
  {
    printf("    ---------------------------------------------------------------------       \n");
    printf("               Solvation Energy for Exhaustive Docking (SEED)                   \n");
    printf("          N. Majeux, M. Scarsi, F. Dey, C. Langini and A. Caflisch              \n");
    printf("                        Department of Biochemistry                              \n");
    printf("                           University of Zurich                                 \n");
    printf("                         %s  (SEED 4.0.0)                                       \n",__DATE__);
    printf("\n    Copyright (C) 2017, Caflisch Lab, University of Zurich               \n");
    printf("                                                                           \n");
    printf("    SEED is free software: you can redistribute it and/or modify\n");
    printf("    it under the terms of the GNU General Public License as published by\n");
    printf("    the Free Software Foundation, either version 3 of the License, or\n");
    printf("    (at your option) any later version.\n");
    printf("\n");
    printf("    SEED is distributed in the hope that it will be useful,\n");
    printf("    but WITHOUT ANY WARRANTY; without even the implied warranty of\n");
    printf("    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n");
    printf("    GNU General Public License for more details.\n");
    printf("\n");
    printf("    You should have received a copy of the GNU General Public License\n");
    printf("    along with this program.  If not, see <http://www.gnu.org/licenses/>.\n");
    printf("    ---------------------------------------------------------------------       \n\n\n\n");

    printf("                WARNING No input file specified, exiting !\n\n");
    printf("                usage : seed [inputfile]\n\n");
    exit(12);
  }

  pi4 = M_PI*4; // clangini

  /* experimental ... OMP check :: */
#ifdef OMP
  numThreads = getenv ("SEED_NUM_THREADS");
  if(numThreads==NULL)
  {
    fprintf(stderr,"WARNING: \"SEED_NUM_THREADS\" is not set, resetting to 1\n");
    omp_set_num_threads(1);
  }
  else
  {
    if(atoi(numThreads)>omp_get_max_threads())
    {
      fprintf(stderr,"WARNING: \"SEED_NUM_THREADS\" %d is larger that maximal number of available threads %d, resetting to latter\n",atoi(numThreads),omp_get_max_threads());
      omp_set_num_threads(omp_get_max_threads());
      /* dbx check max num threads */
    }
    else
    {
      omp_set_num_threads(atoi(numThreads));
    }
  }
  /* t_omp_start = omp_get_wtime(); */
#endif
  /* CLANGINI 2016 */
  /* Make the directory outputs (if it does not exist)*/
  if ((stat("outputs",&DirExist) != 0)||(stat("outputs", &DirExist) == 0 && !S_ISDIR(DirExist.st_mode))){
    if (system("mkdir outputs") == -1)
      std::cout << "Cannot create outputs directory" << std::endl;
  }
  /* Make the directory scratch (if it does not exist)*/
  if ((stat("scratch",&DirExist) != 0)||(stat("scratch", &DirExist) == 0 && !S_ISDIR(DirExist.st_mode))){
    if(system("mkdir scratch") == -1)
      std::cout << "Cannot create outputs directory" << std::endl;
  }
  /* CLANGINI 2016 END */
  /* Read the input and parameters files */
  InpFil=argv[1];

  /* Check presence of input file */
  CheckFile(InpFil,'r');

  //SimWei=matrix(1,150,1,150); clangini
  SimWei=dmatrix(0,150,0,150); //Want to use atom element 0 for lone pair. clangini

  ReInFi(InpFil,RecFil,&BSResN,&BSReNu,FrFiNa,TREFiP,
      &SphAng,&SphPoN,&NuRoAx,&VdWFaB,&CoDieV,&CoDieP,&CoGrIn,&CoGrSi,
      OutFil,&BuEvFa,&FrMaEn,&PsSpRa,&GrSiCu_en,&FiNuMa,&GrInSo,
      &GrSiSo,&WaMoRa,&NPtSphere,&DielWa,&DielRe,ReDesoAlg,DesoMapAcc,
      DesoMapFile,FDexe,FDdir,&ReSurfDens_apol,&PtDensFr,&Sphere_apol,&NCutapolRatio,/*&NCutapol,*/
      &ScaleDeso,&ScaleVDW,SimWei,&SimExp,&SimCut,&FrMaEn_sd,
      &SimExp_sd,&SimCut_sd,&BSMeNu,&BSMeAN,&BSMeVE,CoGrAcc,CoGrFile,
      EvalEn,Solv_typ,SpPoCh_opt,SpPoCh_cent,&SpPoCh_rad,
      &SFDeso_fr,&SFDeso_re,&SFVWEn,&SFIntElec,&NuClusMem,&NuPosMem,&RedRPV_rp,
      &RedRPV_nkvRatio,&ScMaBump,&MuFaVdWCoff_ap,&NuLiEnClus,ApPoChoi,
      &VWGrIn,&VWGrSi,&BumpFaCut,VWGrAcc,VWGrFile,&MaxPosClus,&PrintLev,
      &NumbAT,&AtTyAr,&AtENAr,&VdWRad,&VdWEne,&BLAtTy,&distrPointBSNumb,
      &distrPointBS,&angle_rmin,&angle_rmax,&mult_fact_rmin,&mult_fact_rmax,
      EmpCorrB,gc_opt,&gc_reprke,&gc_cutclus,&gc_endifclus,&gc_weighneg,
      &gc_weighpos,&gc_maxwrite,write_pproc_opt,write_pproc_chm_opt,
      write_best_opt,write_sumtab_opt,write_best_sumtab_opt,&AtWei);/*clangini*/

  /* Check presence of parameter file */
  CheckFile(TREFiP,'r'); /* clangini This is superflous after already having read the parameters*/

  SpPoCh_rad_sq=SpPoCh_rad*SpPoCh_rad;
  /* SFWrNu_ar=ivector(1,FragNu); clangini 2016 ?? */
  /* ResN_fr=cmatrix(1,FragNu,1,10); clangini 2016 */

  /* Extract the output name for each fragment type from the fragments paths */
  /* FrFiNa_out=cmatrix(1,FragNu,1,_STRLENGTH); clangini*/
  /* ExtOutNam(FragNu,FrFiNa,FrFiNa_out); clangini */
  ExtOutNam(FrFiNa,FrFiNa_out);

  /* Open the input files, write informations and check the existence of some
     files */
  ChkExit=0;
  FPaOut=fopen(OutFil,"w"); /* OutFil is by default seed.out clangini */


  time(&runtime);
  fprintf (FPaOut,"Date and time: %s\n",ctime(&runtime)); //WARNING possible memory leak

  fprintf(FPaOut,"    ---------------------------------------------------------------------       \n");
  fprintf(FPaOut,"               Solvation Energy for Exhaustive Docking (SEED)                   \n");
  fprintf(FPaOut,"          N. Majeux, M. Scarsi, F. Dey, C. Langini and A. Caflisch              \n");
  fprintf(FPaOut,"                        Department of Biochemistry                              \n");
  fprintf(FPaOut,"                           University of Zurich                                 \n");
  fprintf(FPaOut,"                         %s  (SEED 4.0.0)                                       \n",__DATE__);
  fprintf(FPaOut,"\n    Copyright (C) 2017, Caflisch Lab, University of Zurich               \n");
  fprintf(FPaOut,"                                                                           \n");
  fprintf(FPaOut,"    SEED is free software: you can redistribute it and/or modify\n");
  fprintf(FPaOut,"    it under the terms of the GNU General Public License as published by\n");
  fprintf(FPaOut,"    the Free Software Foundation, either version 3 of the License, or\n");
  fprintf(FPaOut,"    (at your option) any later version.\n");
  fprintf(FPaOut,"\n");
  fprintf(FPaOut,"    SEED is distributed in the hope that it will be useful,\n");
  fprintf(FPaOut,"    but WITHOUT ANY WARRANTY; without even the implied warranty of\n");
  fprintf(FPaOut,"    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n");
  fprintf(FPaOut,"    GNU General Public License for more details.\n");
  fprintf(FPaOut,"\n");
  fprintf(FPaOut,"    You should have received a copy of the GNU General Public License\n");
  fprintf(FPaOut,"    along with this program.  If not, see <http://www.gnu.org/licenses/>.\n");
  fprintf(FPaOut,"    ---------------------------------------------------------------------       \n\n\n\n");


  fprintf(FPaOut,"Data of the input file :\n\n");
  FilePa=fopen(InpFil,"r"); /* seed.inp clangini */
  fgets_wrapper(StrLin,_STRLENGTH,FilePa);
  while (strncmp(StrLin,"end",3) && !feof(FilePa)) {
    fprintf(FPaOut,"%s",StrLin);
    fgets_wrapper(StrLin,_STRLENGTH,FilePa);
  }
  fclose(FilePa);
  fprintf(FPaOut,"\n");

  fprintf(FPaOut,"Data of the parameters file :\n\n");
  FilePa=fopen(TREFiP,"r"); /* seed.par clangini*/
  fgets_wrapper(StrLin,_STRLENGTH,FilePa);

  while (strncmp(StrLin,"end",3) && !feof(FilePa)   ) {
    fprintf(FPaOut,"%s",StrLin);
    fgets_wrapper(StrLin,_STRLENGTH,FilePa);
  }
  fclose(FilePa);
  fprintf(FPaOut,"\n");

  /* fprintf(FPaOut,"Numbering of the fragments :\n");
     for (i=1;i<=FragNu;i++)
     fprintf(FPaOut,"%3d  %s\n",i,&FrFiNa[i][1]);
     fprintf(FPaOut,"\n");

     if ((FilChk=fopen(RecFil,"r"))==NULL) {
     fprintf(FPaOut,"WARNING Was not able to open %s\n\n",RecFil);
     ChkExit=1;
     }
     else
     fclose(FilChk);

     for (i=1;i<=FragNu;i++) {
     if ((FilChk=fopen(FrFiNa[i]+1,"r"))==NULL) {
     fprintf(FPaOut,"WARNING Was not able to open %s\n\n",&FrFiNa[i][1]);
     ChkExit=1;
     }
     else
     fclose(FilChk);
     }

     if ((FilChk=fopen(TREFiP,"r"))==NULL) {
     fprintf(FPaOut,"WARNING Was not able to open %s\n\n",TREFiP);
     ChkExit=1;
     }
     else
     fclose(FilChk);

     if (!(strncmp(DesoMapAcc,"r",1))) {
     if ((FilChk=fopen(DesoMapFile,"r"))==NULL) {
     fprintf(FPaOut,"WARNING Was not able to open %s\n\n",DesoMapFile);
     ChkExit=1;
     }
     else
     fclose(FilChk);
     }

     if (!(strncmp(CoGrAcc,"r",1))) {
     if ((FilChk=fopen(CoGrFile,"r"))==NULL) {
     fprintf(FPaOut,"WARNING Was not able to open %s\n\n",CoGrFile);
     ChkExit=1;
     }
     else
     fclose(FilChk);
     }

     if (!(strncmp(VWGrAcc,"r",1))) {
     if ((FilChk=fopen(VWGrFile,"r"))==NULL) {
     fprintf(FPaOut,"WARNING Was not able to open %s\n\n",VWGrFile);
     ChkExit=1;
     }
     else
     fclose(FilChk);
     }

     for (i=1;i<=FragNu;i++)
     if (!((ApPoChoi[i][1]=='a')||(ApPoChoi[i][1]=='p')||(ApPoChoi[i][1]=='b'))) {
     fprintf(FPaOut,"WARNING The docking mode should be either a, p or b \n\n");
     ChkExit=1;
     }

     if (!((Solv_typ[0]=='f')||(Solv_typ[0]=='s')||(Solv_typ[0]=='p')||(Solv_typ[0]=='b'))) {
     fprintf(FPaOut,"WARNING The energy evaluation mode should be either f, s, p or b \n\n");
     ChkExit=1;
     }

     if (!((EvalEn[0]=='n')||(EvalEn[0]=='y'))) {
     fprintf(FPaOut,"WARNING The parameter for docking or energy evaluation mode \n");
     fprintf(FPaOut,"        should be either n or y \n\n");
     ChkExit=1;
     }

  if ((EvalEn[0]=='y')&&(Solv_typ[0]=='p'))
    Solv_typ[0]='b';

  fclose(FPaOut);
  FPaOut=fopen(OutFil,"a"); older version, modified by clangini */

    /* fprintf(FPaOut,"Numbering of the fragments :\n");
       for (i=1;i<=FragNu;i++)
       fprintf(FPaOut,"%3d  %s\n",i,&FrFiNa[i][1]);
       fprintf(FPaOut,"\n"); the numbering of the fragments is not known in advance */

    if ((FilChk=fopen(RecFil,"r"))==NULL) {
      fprintf(FPaOut,"WARNING Was not able to open receptor file %s\n\n",RecFil);
      ChkExit=1;
    }
    else
      fclose(FilChk);

  /* Fragment input file is now only one clangini*/
  if ((FilChk=fopen(FrFiNa,"r"))==NULL) {
    fprintf(FPaOut,"WARNING Was not able to open fragment file %s\n\n",FrFiNa);
    ChkExit=1;
  }
  else
    fclose(FilChk);


  if ((FilChk=fopen(TREFiP,"r"))==NULL) {
    fprintf(FPaOut,"WARNING Was not able to open parameter file %s\n\n",TREFiP);
    ChkExit=1;
  }
  else
    fclose(FilChk);

  if (!(strncmp(DesoMapAcc,"r",1))) {
    if ((FilChk=fopen(DesoMapFile,"r"))==NULL) {
      fprintf(FPaOut,"WARNING Was not able to open receptor desolvation grid %s\n\n",DesoMapFile);
      ChkExit=1;
    }
    else
      fclose(FilChk);
  }

  if (!(strncmp(CoGrAcc,"r",1))) {
    if ((FilChk=fopen(CoGrFile,"r"))==NULL) {
      fprintf(FPaOut,"WARNING Was not able to open Coulombic grid %s\n\n",CoGrFile);
      ChkExit=1;
    }
    else
      fclose(FilChk);
  }

  if (!(strncmp(VWGrAcc,"r",1))) {
    if ((FilChk=fopen(VWGrFile,"r"))==NULL) {
      fprintf(FPaOut,"WARNING Was not able to open Van der Waals grid %s\n\n",VWGrFile);
      ChkExit=1;
    }
    else
      fclose(FilChk);
  }

  /* Only one docking mode has to be chosen for all the fragments */
  if (!((ApPoChoi[0]=='a')||(ApPoChoi[0]=='p')||(ApPoChoi[0]=='b'))) {
    fprintf(FPaOut,"WARNING The docking mode should be either a, p or b \n\n");
    ChkExit=1;
  }

  if (!((Solv_typ[0]=='f')||(Solv_typ[0]=='s')||(Solv_typ[0]=='p')||(Solv_typ[0]=='b'))) {
    fprintf(FPaOut,"WARNING The energy evaluation mode should be either f, s, p or b \n\n");
    ChkExit=1;
  }

  if (!((EvalEn[0]=='d')||(EvalEn[0]=='e'))) {
    fprintf(FPaOut,"WARNING The parameter for docking or energy evaluation mode \n");
    fprintf(FPaOut,"        should be either d or e \n\n");
    ChkExit=1;
  }

  /*if ((EvalEn[0]=='e')&&(Solv_typ[0]=='p'))
    Solv_typ[0]='b';*/ //clangini April 2017 simplification of seed.par

  if((EvalEn[0] == 'e')&&(Solv_typ[0]!='s')) {
    std::cout << "Energy evalutaion (e) requested. Changing energy mode from ("
      << Solv_typ << ") to (s). Only slow energies will be computed."
      << std::endl;
    Solv_typ[0] = 's';
  } //clangini April 2017 simplification of seed.par

  fclose(FPaOut); /* Why do we close it in 'w' mode and re-open in 'a' mode? clangini */
  FPaOut=fopen(OutFil,"a");

  /* If any problem (WARNING), the program exits */
  if (ChkExit) {
    printf("Program exits\n");
    exit(0);
  }

  /* assign ResN_fr (name of fragment residue) */
  /* CheckRESN(FrFiNa,FragNu,ResN_fr); we do not use CheckRESN any more.
     ResN_fr will be read directly from the input file */

  fprintf(FPaOut,"Data for the receptor :\n\n");

  /* Read the receptor file */
  ReReFi_mol2(RecFil,&ReAtNu,&ReBdNu,&ReReNu,&ReAtEl,&ReCoor,&ReAtTy,
      &ReResN,&RePaCh,&ReBdAr); /* The receptor reading is not modified clangini*/

  /* KNEHANS 26 03 2012
     this function checks if there
     are residue indices in the input
     file which aren't corresponding
     to the indices in the receptor
     file; a message is issued and
     the program exits
   */


  for (i2=1;i2 <= BSResN; i2++) {
    int CheInpWroRe = 0;
    int CheInpDupRe = 0;
    /* for (i=1;i <= ReAtNu; i++)  orig  */
    for (i=2;i <= ReAtNu; i++){
      if(ReResN[i-1] != ReResN[i] || ReResN[i]==1){ /* modified to cover also first  residue, dey */
        if(ReResN[i] == BSReNu[i2]){
          CheInpWroRe++;
        }
      }
    }
    if (CheInpWroRe == 0){
      printf("WARNING! Input file residue index: %d (i10) is inconsistent with residue indices in the receptor file.\n", BSReNu[i2]);
      ChkExit=1;
    }
    for (i3=1; i3 <= BSResN; i3++){
      if(BSReNu[i3] == BSReNu[i2]){
        CheInpDupRe++;
      }

    }
    if (CheInpDupRe > 1){
      printf("WARNING! Input file residue index: %d (i10) is %d times in the input file.\n", BSReNu[i2], CheInpDupRe);
      ChkExit=1;
    }

    if (ChkExit) {
      printf("Program exits\n");
      exit(0);
    }
  }
  /*KNEHANS 2012 end*/

  /* Determine the first and last atoms and computes the total charge for each
     residue. Find the residue-representative atom which is the closest to the
     geometrical (for uncharged residue) or the charge (for charged residue)
     center */
  FiAtRes=ivector(1,ReReNu);
  LaAtRes=ivector(1,ReReNu);
  AtReprRes=ivector(1,ReReNu);
  TotChaRes=dvector(1,ReReNu);
  FeatRes(ReAtNu,ReCoor,ReReNu,ReResN,RePaCh,FiAtRes,LaAtRes,AtReprRes,
      TotChaRes,FPaOut);

  /* Create the receptor file in a kind of PDB format for UHBD.
     Only when UHBD is used */
  /*if (DesoMapAcc[0]!='r')
    Wr_pdb_uhbd(ReAtNu,ReCoor,RecFilPDB); clangini commented */

  /* Make the list of atoms which are in the binding site */
  MakBSAtList(ReAtNu,ReResN,BSResN,BSReNu,&BSAtNu,&BSAtLi);

  /* Assign the atom element numbers, van der Waals radii and energies
     for the receptor */
  ReVdWR=dvector(1,ReAtNu);
  ReVdWE=dvector(1,ReAtNu);
  ReAtEl_nu=ivector(1,ReAtNu);
  ReAtoTyp_nu=ivector(1,ReAtNu);
  if(! AssiRE(NumbAT,AtTyAr,VdWRad,VdWEne,ReAtNu,ReAtTy,ReVdWR,ReVdWE,FPaOut,
        AtENAr,ReAtEl_nu,ReAtoTyp_nu))
  {
    /* bail on missing types */
    printf("Program exits\n");
    exit(13);
  }

  /* Find the hybridization state for the carbon atoms of the receptor */
  HybReAt=ivector(1,ReAtNu);
  HybStAt(ReAtNu,ReAtEl_nu,ReCoor,ReBdNu,ReBdAr,HybReAt,FPaOut);

  /* Compute the square root of the van der Waals receptor energies
     and the power of 3 of the receptor van der Waals radii */
  ReVdWE_sr=dvector(1,ReAtNu);
  SqRoEn(ReAtNu,ReVdWE,ReVdWE_sr);
  ReVdWR_p3=dvector(1,ReAtNu);
  for (i=1;i<=ReAtNu;i++)
    ReVdWR_p3[i]=ReVdWR[i]*ReVdWR[i]*ReVdWR[i];

  /* Construct the list of the donor and receptor vectors for the binding site
     of the receptor */
  ReAcDo(BSAtNu,BSAtLi,ReAtEl_nu,ReCoor,HybReAt,ReBdNu,ReBdAr,SphAng,
      SphPoN,&ReAcNu,&ReDoNu,&ReDATy,&ReDAAt,&ReHydN,&ReVeCo,FPaOut,OutFil,
      BSMeNu,BSMeAN,BSMeVE);

  /* Reduce the number of polar vectors */
  PolVect_rec_01=ivector(1,ReAcNu+ReDoNu);

  Reduc_polvectre(ReAtNu,ReCoor,ReVdWR,ReVdWE_sr,ReAcNu,ReDoNu,ReDATy,
      ReDAAt,ReHydN,ReVeCo,RedRPV_rp,RedRPV_nkvRatio,PolVect_rec_01,
      FPaOut,distrPointBSNumb,distrPointBS,angle_rmin,
      angle_rmax,mult_fact_rmin,mult_fact_rmax,Sphere_apol);

  /* Find the largest van der Waals radius */
  VdWRaM(NumbAT,VdWRad,&LaVdWR);

  /* Find the maximal and minimal coordinates of the receptor */
  ReMMCo(ReAtNu,ReCoor,ReMaxC,ReMinC);

  /* Find the maximal and minimal coordinates of the binding site of the
     receptor */
  BSMMCo(ReCoor,BSAtNu,BSAtLi,BSMaxC,BSMinC);

  if ((Solv_typ[0]=='s')||(Solv_typ[0]=='b')||(Solv_typ[0]=='p')) {

    /* Construct the list of receptor atoms which are in the cubes of a grid
       for the checking of bumps in CoNuBu */
    ReLAIC(ReAtNu,ReCoor,LaVdWR,ReMaxC,ReMinC,CubNum,&CubFAI,&CubLAI,&CubLiA);

    /* Find which charged residues will always be taken into account for the
       evaluation of the energy */
    LiChResEn=ivector(1,ReReNu);
    ListChRes(BSResN,BSReNu,ReReNu,ReCoor,AtReprRes,TotChaRes,PsSpRa,FPaOut,
        &NuChResEn,LiChResEn);

    /* Construct the list of residue-representative atoms (of the receptor) which
       are in the cubes of a grid for the energy evaluation */
    for (i=1;i<=3;i++)
      CubNum_en[i]=ffloor((ReMaxC[i]-ReMinC[i])/GrSiCu_en)+1;
    CubFAI_en=i3tensor(1,CubNum_en[1],1,CubNum_en[2],1,CubNum_en[3]);
    CubLAI_en=i3tensor(1,CubNum_en[1],1,CubNum_en[2],1,CubNum_en[3]);
    CubLiA_en=ivector(1,ReReNu+10);
    ReLAIC_en(ReReNu,ReCoor,GrSiCu_en,ReMinC,CubNum_en,CubFAI_en,CubLAI_en,
        CubLiA_en,AtReprRes,FPaOut);

    /* Construct the pseudo-sphere for the evaluation of the energy */
    PsSpNC=ffloor(PsSpRa/GrSiCu_en)+1;
    PsSphe=i3tensor(1,2*PsSpNC+1,1,2*PsSpNC+1,1,2*PsSpNC+1);
    PsSpMa(PsSpRa,GrSiCu_en,PsSpNC,PsSphe);

  }

  if ((Solv_typ[0]=='f')||(Solv_typ[0]=='b')||(Solv_typ[0]=='p')) {

    gettimeofday(&time_1,NULL);
    /* times(&timevar); */
    /* time_1=timevar.tms_utime+timevar.tms_stime; */

    /* Compute or read the receptor part of the coulombic interaction on a grid */
    for (i=1;i<=3;i++)
      CoGPoN[i]=ffloor((BSMaxC[i]+CoGrIn-(BSMinC[i]-CoGrIn))/CoGrSi)+2;
    fprintf(FPaOut,"Number of points of the coulombic grid : %d\n\n",
        CoGPoN[1]*CoGPoN[2]*CoGPoN[3]);
    CoGrRP=d3tensor(1,CoGPoN[1],1,CoGPoN[2],1,CoGPoN[3]);
    if (CoGrAcc[0]=='r') {
      FilePa=fopen(CoGrFile,"r");
      fgets_wrapper(StrLin,_STRLENGTH,FilePa);
      sscanf(StrLin,"%d%d%d",&Nx,&Ny,&Nz);
      /* to do binary FilePa=fopen(CoGrFile,"rb");
         fread(&Nx,1,sizeof(int),FilePa);
         fread(&Ny,1,sizeof(int),FilePa);
         fread(&Nz,1,sizeof(int),FilePa);
       */
      if ((Nx!=CoGPoN[1])||(Ny!=CoGPoN[2])||(Nz!=CoGPoN[3])) {
        fprintf(FPaOut,"WARNING The Coulombic grid you want to read has not ");
        fprintf(FPaOut,"been created with the\n");
        fprintf(FPaOut,"        same input file parameters\n\n");
        exit(0);
      }
      for (i1=1;i1<=CoGPoN[1];i1++) {
        for (i2=1;i2<=CoGPoN[2];i2++) {
          for (i3=1;i3<=CoGPoN[3];i3++) {
            /* fread(&(CoGrRP[i1][i2][i3]),1,sizeof(double),FilePa); to do binary */

            fgets_wrapper(StrLin,_STRLENGTH,FilePa);
            sscanf(StrLin,"%lf",&(CoGrRP[i1][i2][i3]));

          }
        }
      }
      fclose(FilePa);
    }
    else {
      CoGReP(ReAtNu,ReCoor,RePaCh,CoDieV,CoDieP,CoGrIn,CoGrSi,BSMinC,CoGPoN,
          CoGrRP);
      printf("Receptor part of the coulombic interaction -> done\n");
      if (CoGrAcc[0]=='w') {
        FilePa=fopen(CoGrFile,"w"); /* to change to "wb" for binray grids*/
        /* FilePa=fopen(CoGrFile,"wb"); /\* to change to "wb" for binray grids*\/ */
        fprintf(FilePa,"%d %d %d\n",CoGPoN[1],CoGPoN[2],CoGPoN[3]);
        /*
           fwrite(&CoGPoN[1],sizeof(int),1,FilePa);
           fwrite(&CoGPoN[2],sizeof(int),1,FilePa);
           fwrite(&CoGPoN[3],sizeof(int),1,FilePa);
         */
        for (i1=1;i1<=CoGPoN[1];i1++) {
          for (i2=1;i2<=CoGPoN[2];i2++) {
            for (i3=1;i3<=CoGPoN[3];i3++) {
              fprintf(FilePa,"%.12f\n",CoGrRP[i1][i2][i3]);
              /* fwrite(&CoGrRP[i1][i2][i3],sizeof(double),1,FilePa); /\*  ->fixed size writing to do binary file writing*\/ */
            }
          }
        }
        fclose(FilePa);
      }
    }

    /* times(&timevar); */
    /* time_2=timevar.tms_utime+timevar.tms_stime; */

    /* times(&timevar); */
    /* time_3=timevar.tms_utime+timevar.tms_stime; */
    gettimeofday(&time_2,NULL);
    gettimeofday(&time_3,NULL);


    /* Compute or read the receptor part of the van der Waals interactions on a grid */
    for (i=1;i<=3;i++)
      VWGPoN[i]=ffloor((BSMaxC[i]+VWGrIn-(BSMinC[i]-VWGrIn))/VWGrSi)+2;
    //clangingi debug:
    //std::cout << "(BSMinC[i]-VWGrIn)  " << (BSMinC[1]-VWGrIn) << std::endl;
    //std::cout << "(BSMaxC[i]+VWGrIn)  " << (BSMaxC[1]+ VWGrIn) << std::endl;
    //std::cout << "VWGPoN[i]  " << VWGPoN[1] << std::endl;
    //end clangini debug

    fprintf(FPaOut,"Number of points of the van der Waals grid : %d\n\n",
        VWGPoN[1]*VWGPoN[2]*VWGPoN[3]);
    VWGrRP_at=d3tensor(1,VWGPoN[1],1,VWGPoN[2],1,VWGPoN[3]);
    VWGrRP_re=d3tensor(1,VWGPoN[1],1,VWGPoN[2],1,VWGPoN[3]);
    if (VWGrAcc[0]=='r') {
      FilePa=fopen(VWGrFile,"r");
      fgets_wrapper(StrLin,_STRLENGTH,FilePa);
      sscanf(StrLin,"%d%d%d",&Nx,&Ny,&Nz);
      if ((Nx!=VWGPoN[1])||(Ny!=VWGPoN[2])||(Nz!=VWGPoN[3])) {
        fprintf(FPaOut,"WARNING The van der Waals grid you want to read has not ");
        fprintf(FPaOut,"been created with the\n");
        fprintf(FPaOut,"        same input file parameters\n\n");
        exit(0);
      }
      for (i1=1;i1<=VWGPoN[1];i1++) {
        for (i2=1;i2<=VWGPoN[2];i2++) {
          for (i3=1;i3<=VWGPoN[3];i3++) {
            fgets_wrapper(StrLin,_STRLENGTH,FilePa);
            sscanf(StrLin,"%lf%lf",&(VWGrRP_at[i1][i2][i3]),&(VWGrRP_re[i1][i2][i3]));
          }
        }
      }
      fclose(FilePa);
    }
    else {
      VWGReP(ReAtNu,ReCoor,ReVdWR_p3,ReVdWE_sr,VWGrIn,VWGrSi,BSMinC,VWGPoN,
          ReMaxC,ReMinC,VWGrRP_at,VWGrRP_re);
      printf("Receptor part of the van der Waals interaction -> done\n");
      if (VWGrAcc[0]=='w') {
        FilePa=fopen(VWGrFile,"w");
        fprintf(FilePa,"%d %d %d\n",VWGPoN[1],VWGPoN[2],VWGPoN[3]);
        for (i1=1;i1<=VWGPoN[1];i1++) {
          for (i2=1;i2<=VWGPoN[2];i2++) {
            for (i3=1;i3<=VWGPoN[3];i3++) {
              fprintf(FilePa,"%.12f %.12f\n",VWGrRP_at[i1][i2][i3],VWGrRP_re[i1][i2][i3]);
            }
          }
        }
        fclose(FilePa);
      }
    }

    /* times(&timevar); */
    /* time_4=timevar.tms_utime+timevar.tms_stime; */
    gettimeofday(&time_4,NULL);

  }


  /*
     Hard-coded parameters (they all refer to solvation). They have been derived
     by comparison with PBEQ in CHARMM on several molecules as described in the
     following LaTeX file:


     \documentstyle[12pt,titlepage]{article}

     \voffset=-8mm
     \hoffset= 0mm
     \topmargin=0mm
     \headheight=0in
     \headsep=8mm
     \footheight=0in
     \footskip=0in
     \textheight=234mm
     \textwidth=148mm
     \oddsidemargin=0in
     \evensidemargin=0in
     \pagestyle{myheadings}
     \markboth{}{}

     \begin{document}

     \pagestyle{empty}

     \begin{tabular}{lrrrrrrrrr}
     \hline
     & & \multicolumn{2}{c} {  screened  } & & \multicolumn{2}{c} {  fragment  } & &
     \multicolumn{2}{c} {  receptor  } \\
     & & \multicolumn{2}{c} {  interaction  } & & \multicolumn{2}{c} {  desolvation  } & &
     \multicolumn{2}{c} {  desolvation  } \\
     \cline{3-4} \cline{6-7} \cline{9-10}
     model & & corr & slope & & corr & slope & & corr & slope \\
     \hline
     \\
     & & \multicolumn{8}{c} {{\bf thrombin} $\varepsilon_p$ = 1 (1025 data points)} \\
     SEED        & & 0.96 & 1.69 & & 0.91 & 1.19 & & 0.91 & 0.86 \\
     SEED\_corrB & & 0.95 & 0.95 & & 0.99 & 0.86 & &      &      \\
     \\
     & & \multicolumn{8}{c} {{\bf thrombin} $\varepsilon_p$ = 4 (1025)} \\
     SEED        & & 0.95 & 1.38 & & 0.93 & 1.15 & & 0.92 & 0.86 \\
     SEED\_corrB & & 0.88 & 0.81 & & 0.98 & 0.81 & &      &      \\
     \\
     & & \multicolumn{8}{c} {{\bf hiv} $\varepsilon_p$ = 1 (1490)} \\
     SEED        & & 0.97 & 1.44 & & 0.93 & 1.39 & & 0.96 & 0.89 \\
     SEED\_corrB & & 0.99 & 1.02 & & 0.98 & 0.88 & &      &      \\
     \\
     & & \multicolumn{8}{c} {{\bf hiv} $\varepsilon_p$ = 4 (1490)} \\
     SEED        & & 0.98 & 1.37 & & 0.95 & 1.31 & & 0.97 & 0.90 \\
     SEED\_corrB & & 0.99 & 0.98 & & 0.98 & 0.82 & &      &      \\
     \\
     & & \multicolumn{8}{c} {{\bf mdm2} $\varepsilon_p$ = 1 (613)} \\
     SEED        & & 0.97 & 1.41 & & 0.72 & 0.86 & & 0.85 & 0.55 \\
     SEED\_corrB & & 0.98 & 0.98 & & 0.97 & 0.90 & &      &      \\
     \\
     & & \multicolumn{8}{c} {{\bf mdm2} $\varepsilon_p$ = 4 (613)} \\
     SEED        & & 0.93 & 1.16 & & 0.78 & 0.92 & & 0.87 & 0.57 \\
     SEED\_corrB & & 0.95 & 0.83 & & 0.95 & 0.87 & &      &      \\
     \\
     & & \multicolumn{8}{c} {$\beta$-{\bf secretase} $\varepsilon_p$ = 1 (592)} \\
     SEED        & & 0.94 & 1.50 & & 0.92 & 1.57 & & 0.89 & 0.67 \\
     SEED\_corrB & & 0.99 & 1.02 & & 0.97 & 0.85 & &      &      \\
     \\
     & & \multicolumn{8}{c} {$\beta$-{\bf secretase} $\varepsilon_p$ = 4 (592)} \\
     SEED        & & 0.93 & 1.50 & & 0.94 & 1.46 & & 0.90 & 0.69 \\
     SEED\_corrB & & 0.96 & 1.04 & & 0.97 & 0.77 & &      &      \\
     \\
     & & \multicolumn{8}{c} {{\bf p38 MAP kinase} $\varepsilon_p$ = 1 (710)} \\
     SEED        & & 0.98 & 1.68 & & 0.90 & 1.65 & & 0.91 & 0.76 \\
     SEED\_corrB & & 0.99 & 0.98 & & 0.97 & 0.81 & &      &      \\
  \\
    & & \multicolumn{8}{c} {{\bf p38 MAP kinase} $\varepsilon_p$ = 4 (710)} \\
    SEED        & & 0.96 & 1.43 & & 0.93 & 1.51 & & 0.92 & 0.78 \\
    SEED\_corrB & & 0.94 & 0.85 & & 0.97 & 0.72 & &      &      \\
    \\
    \hline
    \end{tabular}

  \begin{tabular}{lrrrrrrrrr}
  \multicolumn{10}{l} {PB calculations were performed with the PBEQ module of CHARMM} \\
    \multicolumn{10}{l} {$\varepsilon_p = $ dielectric constant of the solute}
  \end{tabular}

  \newpage

    \begin{tabular}{lcccccc}
  \hline
    & &          & \multicolumn{1}{c}{total number}
  & \multicolumn{1}{c}{residues in the}
  & \multicolumn{1}{c}{number of}
  & \multicolumn{1}{c}{types of} \\
    \multicolumn{1}{c}{protein} & & \multicolumn{1}{c}{PDB code}
  & \multicolumn{1}{c}{of residues}
  & \multicolumn{1}{c}{binding site}
  & \multicolumn{1}{c}{frag. positions}
  & \multicolumn{1}{c}{fragments} \\
    \hline
    \\
    thrombin          & & 1hgt  & 251 & 251 & 1025 & 7$^a$ \\
    hiv               & & 5hvp? &  99 &  99 & 1490 & 7$^a$ \\
    mdm2              & & 1ycq  &  88 &  16 &  613 & 7$^a$ \\
    $\beta$-secretase & & 1fkn  & 391 &  21 &  592 & 15$^b$ \\
    p38 MAP kinase    & & 1a9u  & 351 &  32 &  710 & 15$^b$ \\
    \\
    \hline
    \end{tabular}

  \begin{tabular}{lcccccc}
  \multicolumn{7}{l} {$^a$ 1 apolar, 1 uncharged polar, 3 charged $-$, 2 charged $+$} \\
    \multicolumn{7}{l} {$^b$ 5 apolar, 5 uncharged polar, 3 charged $-$, 2 charged $+$}
  \end{tabular}

  \end{document}


  The average slopes are the following:

    screened         fragment         receptor
    interaction      desolvation      desolvation

    SEED          1.456             1.301           0.753
    SEED_corrB       0.946             0.829              "


    In the code all variables concerning the empirical correction term to the
    Coulomb field approximation for the accurate screened interaction and
    fragment desolvation energies have the suffix "corrB".

    */

    /* correction factor for receptor desolvation when using Coulomb approximation
       It affects the rec. desolvation of the slow and the fast algorithm.
       It is used also when calculating the desolvation resulting from a sphere
       on the surface of the frag according to the Coulomb approx. */
    corr_re_desoco = 1.33;

  if (EmpCorrB[0]=='y')
  {

    /* correction factor for fragment desolvation --> it affects only the fragment
       desolvation calculated with the slow algorithm (GB) */
    corr_fr_deso = 1.21;
    /* correction factor for screened interactions (only when using the GB formula -->
       --> it affects only slow algorithm) */
    corr_scrint = 1.06;

  }
  else
  {

    /* correction factor for fragment desolvation --> it affects only the fragment
       desolvation calculated with the slow algorithm (GB) */
    corr_fr_deso = 0.77;
    /* correction factor for screened interactions (only when using the GB formula -->
       --> it affects only slow algorithm) */
    corr_scrint = 0.69;

  }



  /* hard-coded parameters (they all refer to solvation). They have been derived
     by comparison with UHBD for thrombin with eps_w = 78.5 and eps_p = 1 ? */

  /* correction factor for receptor desolvation when using finite diff. approximation
     It affects the rec. desolvation of the slow and the fast algorithm. */
  corr_re_desofd = 0.39;
  /* surf point dens. for fast desolvation */
  SurfDens_deso = .5;
  /* The following is done so that the desolvations of the fast algorithm is
     indipendent of corr_re_desoco, corr_re_desofd and on the point density */
  sprintf(FD,"%s","fd\0");
  sprintf(CO,"%s","co\0");
  if ( strcmp(ReDesoAlg,CO) == 0 )
    corr_fast_deso = 0.24/(SurfDens_deso * corr_re_desoco);
  else if ( strcmp(ReDesoAlg,FD) == 0 ) /* this is no longer supported, ReDesoAlg is hard-coded to be "co" */
    corr_fast_deso = 0.074/(SurfDens_deso * corr_re_desofd);

  /* correction factor for fast desolvation
     (it depends only on the point density). Not used...
     corr_fast_deso = 0.19/SurfDens_deso;
   */



  /* Allocate some memory for variables that will be kept in use during all the program */
  /* all of this regards the receptor only clangini */
  ReRad=dvector(1,ReAtNu);
  ReRad2=dvector(1,ReAtNu);
  ReRadOut=dvector(1,ReAtNu);
  ReRadOut2=dvector(1,ReAtNu);

  surfpt_re=structpointvect(1,NPtSphere*ReAtNu);
  iatsrf_re=ivector(1,NPtSphere*ReAtNu);
  nsurf_re=ivector(1,ReAtNu);
  pointsrf_re=ivector(1,ReAtNu);
  for (ij=1;ij<=ReAtNu;ij++)
    pointsrf_re[ij]=0;
  for (iat=1;iat<=ReAtNu;iat++)
    nsurf_re[iat] = 0;

  ReRmax = MaxVector(ReVdWR,1,ReAtNu);
  NPtSphereMax = (int) (1. * pi4 * (ReRmax+WaMoRa) * (ReRmax+WaMoRa));
  surfpt_re_deso=structpointvect(1,NPtSphereMax*ReAtNu);
  iatsrf_re_deso=ivector(1,NPtSphereMax*ReAtNu);
  nsurf_re_deso=ivector(1,ReAtNu);
  pointsrf_re_deso=ivector(1,ReAtNu);
  for (ij=1;ij<=ReAtNu;ij++)
    pointsrf_re_deso[ij]=0;
  for (iat=1;iat<=ReAtNu;iat++)
    nsurf_re_deso[iat] = 0;

  /* Do all the precalculation necessary for continuum electrostatics */
  Solvation(ReAtNu,ReCoor,ReVdWE_sr,ReVdWR,ReRad,ReRad2,ReRadOut,ReRadOut2,
      ReMaxC,ReMinC,RePaCh,DielRe,DielWa,WaMoRa,GrSiSo,GrInSo,
      NPtSphere,ReResN,ReReNu,BSResN,BSReNu,ReDesoAlg,DesoMapAcc,
      DesoMapFile,RecFilPDB,FDexe,FDdir,&Min,&Max,&XGrid,
      &YGrid,&ZGrid,&GridMat,&NGridx,&NGridy,&NGridz,
      &Nsurfpt_re,surfpt_re,iatsrf_re,nsurf_re,pointsrf_re,
      &DeltaPrDeso,&nxminBS,&nyminBS,&nzminBS,&nxmaxBS,&nymaxBS,
      &nzmaxBS,&ReSurf_deso,&Nsurfpt_re_deso,surfpt_re_deso,
      iatsrf_re_deso,nsurf_re_deso,pointsrf_re_deso,SurfDens_deso,
      &nlist_re,&list_re,&MaxNeigh,&nstep_re,&cbmid_re,
      &midg_re,&rslop,&rscale_re,
      &ReSurf_apol,&Nsurfpt_re_apol,ReSurfDens_apol,
      Sphere_apol,NCutapolRatio,ScaleDeso,ScaleVDW, /*NCutapol,*/
      &apol_Vect_re,&ReApAt,&Nsurfpt_re_apol_BS,&Napol_Vect_re,
      &ReSelfVol,&Kelec,&Ksolv,&UnitVol,pi4,corr_re_desoco,
      corr_re_desofd,corr_fast_deso,
      distrPointBSNumb,distrPointBS,angle_rmin,angle_rmax,
      mult_fact_rmin,mult_fact_rmax,FPaOut,&ReSelfVol_corrB,EmpCorrB);



  FilePa=fopen("./outputs/apolar_rec_reduc.mol2","w");
  /*fprintf(FilePa,"# TRIPOS MOL2 file generated by WITNOTP\n"); clangini */
  fprintf(FilePa,"# TRIPOS MOL2 file generated by SEED\n"); /* clangini */
  fprintf(FilePa,"\n");
  fprintf(FilePa,"@<TRIPOS>MOLECULE\n");
  fprintf(FilePa,"apolar_rec_reduc\n");
  fprintf(FilePa,"%d 0 0 0 0\n",Napol_Vect_re);
  fprintf(FilePa,"****\n");
  fprintf(FilePa,"USER_CHARGES\n");
  fprintf(FilePa,"INVALID_CHARGES\n");
  fprintf(FilePa,"\n");
  fprintf(FilePa,"@<TRIPOS>ATOM\n");
  for (i1=1;i1<=Napol_Vect_re;i1++) {
    /*
          fprintf(FilePa,"%-5d %-3d%11.4f%11.4f%11.4f ***** 1    ****  0.000\n",
          i1*2-1,i1*2-1,apol_Vect_re[i1][1],apol_Vect_re[i1][2],
          apol_Vect_re[i1][3]);
          fprintf(FilePa,"%-5d %-3d%11.4f%11.4f%11.4f ***** 1    ****  0.000\n",
          i1*2,i1*2,apol_Vect_re[i1][4],apol_Vect_re[i1][5],
          apol_Vect_re[i1][6]);
     */
    fprintf(FilePa,"%-5d %-3d%11.4f%11.4f%11.4f ***** 1    ****  0.000\n",
        i1,i1,apol_Vect_re[i1][4],apol_Vect_re[i1][5],
        apol_Vect_re[i1][6]);
  }
  fclose(FilePa);


  /* Normalization of the vectors apol_Vect_re */
  for (i1=1;i1<=Napol_Vect_re;i1++)
    PoCoVe(apol_Vect_re[i1][1],apol_Vect_re[i1][2],apol_Vect_re[i1][3],
        apol_Vect_re[i1][4],apol_Vect_re[i1][5],apol_Vect_re[i1][6],1.0,
        &(apol_Vect_re[i1][4]),&(apol_Vect_re[i1][5]),
        &(apol_Vect_re[i1][6]));

  if(Nsurfpt_re_apol_BS==0)
    fprintf(stderr,"WARNING, number of kept apolar vectors == 0 !\n");/*dey*/

  fprintf(FPaOut,"Total number of apolar vectors generated on the ");
  fprintf(FPaOut,"receptor binding site : %d\n",Nsurfpt_re_apol_BS);
  fprintf(FPaOut,"Number of kept apolar vectors for seeding of ");
  fprintf(FPaOut,"apolar fragments : %d\n\n",Napol_Vect_re);

  fclose(FPaOut);

  /* times(&timevar); */
  /* time_5=timevar.tms_utime+timevar.tms_stime; */
  gettimeofday(&time_5,NULL);

  /* ***********************
     Loop over the fragments
   *********************** */
  /* The next part will be changed as we do not know in advance the number of fragments
     We loop until we get to the end of the fragment input file. clangini */

  /* dey new: memory allocation of arrays handled dynamically -> realloc(..) resizing */
  Index_ro=ivector(1,currentsize); /* up to now it is not clear to me what they are for. clangini */
  VW_f_ro=dvector(1,currentsize);
  VW_s_ro=dvector(1,currentsize);
  In_f_ro=dvector(1,currentsize);
  In_s_ro=dvector(1,currentsize);
  Dr_f_ro=dvector(1,currentsize);
  Dr_s_ro=dvector(1,currentsize);
  Df_f_ro=dvector(1,currentsize);
  Df_s_ro=dvector(1,currentsize);
  To_f_ro=dvector(1,currentsize);
  To_s_ro=dvector(1,currentsize);

  /* clangini 2016 */
  /* Setting up the table summary file */
  if (write_sumtab_opt[0]=='y'){ // Should introduce some check.
    strcpy(TabFil,"./outputs/seed_clus.dat");
    //char TabLin[_STRLENGTH];
    //std::ofstream TabOutStream;
    TabOutStream.open (TabFil);
    if(TabOutStream.is_open()){
      TabOutStream << std::left << std::setw(30) << "Name" << std::right
                   << std::setw(8)  << "Pose"
                   << std::setw(10) << "Cluster"
                   << std::setw(10) << "Fr_nu"
                   << std::setw(10) << "Tot"
                   << std::setw(10) << "ElinW"
                   << std::setw(10) << "rec_des"
                   << std::setw(10) << "frg_des"
                   << std::setw(10) << "vdW"
                   << std::setw(10) << "DElec"
                   << std::setw(10) << "DG_hydr"
                   << std::setw(10) << "Tot_eff"
                   << std::setw(10) << "vdW_eff"
                   << std::setw(10) << "Elec_eff"
                   << std::setw(10) << "HAC"
                   << std::setw(10) << "MW"
                   << std::endl;
    } else {
      std::cerr << "Unable to open file "<< TabFil
                << "for writing summary table." << std::endl;
    }
    TabOutStream.close();
  }
  if (write_best_sumtab_opt[0]=='y'){
    // Second summary table, with the best poses:
    strcpy(TabFil,"./outputs/seed_best.dat");
    TabOutStream.open (TabFil);
    if(TabOutStream.is_open()){
      TabOutStream << std::left << std::setw(30) << "Name" << std::right
                   << std::setw(8)  << "Pose"
                   << std::setw(10)  << "Cluster"
                   << std::setw(10) << "Fr_nu"
                   << std::setw(10) << "Tot"
                   << std::setw(10) << "ElinW"
                   << std::setw(10) << "rec_des"
                   << std::setw(10) << "frg_des"
                   << std::setw(10) << "vdW"
                   << std::setw(10) << "DElec"
                   << std::setw(10) << "DG_hydr"
                   << std::setw(10) << "Tot_eff"
                   << std::setw(10) << "vdW_eff"
                   << std::setw(10) << "Elec_eff"
                   << std::setw(10) << "HAC"
                   << std::setw(10) << "MW"
                   << std::endl;
    } else {
      std::cerr << "Unable to open file "<< TabFil
                << "for writing best poses summary table." << std::endl;
    }
    TabOutStream.close();
  }

  int CurFraTot = 0; /* Current fragment in the file (counting also skipped ones) clangini */
  int SkiFra = 0; /* Counter of skipped fragments clangini*/
  CurFra = 0; /* Current fragment clangini */
  int LstFra_f = 0; /* flag to signal that the last fragment was reached */

  std::ifstream FrInStream; /* declare input stream for the fragment */
  FrInStream.open(FrFiNa, std::ios_base::binary); /* open input stream for the fragment */
  if (FrInStream.fail()) {
    fprintf(stderr,"Cannot open specified input file %s\nProgram exits!\n",FrFiNa);
    exit(13);
  }
  std::streampos FrInPos = FrInStream.tellg(); /* Save the get position */

  if (write_pproc_opt[0]=='y'){
    /* FILE *FilePa;*/
    /* Create WriPat and open file to write a commented header */

    //sprintf(WriPat,"%s%s%s","./outputs/",&FrFiNa_out[1], /* should give the possibility to */
    //        "_clus_pproc.mol2\0");               /*use another output folder? */
    sprintf(WriPat,"%s%s%s","./outputs/",FrFiNa_out, /* should give the possibility to */
            "_clus.mol2\0");               /*use another output folder? */

    FilePa=fopen(WriPat,"w");
    fprintf(FilePa,"# TRIPOS MOL2 file generated by SEED\n\n");
    fclose(FilePa);
  }
  if (write_best_opt[0]=='y'){
    /* Create WriPat and open file to write a commented header */
    sprintf(WriPat,"%s%s%s","./outputs/",FrFiNa_out,
            "_best.mol2\0");
    FilePa=fopen(WriPat,"w");
    fprintf(FilePa,"# TRIPOS MOL2 file generated by SEED\n\n");
    fclose(FilePa);
  }
  while((FrInStream.eof() == 0)&&(LstFra_f == 0)) { /* loop until reach end of file -> loop over fragments clangini*/
    /*for (i=1;i<=FragNu;i++) */ /* open loop over fragments (has to be removed) clangini */
    /* control amount of printed "out of grid" error messages */

    //std::cout << "Entered loop over fragments!" << std::endl; /*clangini OK*/

    vdwErr = 0;
    clbErr = 0;

    /* CurFra=i; clangini */
    FPaOut=fopen(OutFil,"a");

    //fprintf(FPaOut,"-------------------------------------------------\n\n");//Moved after fragment reading. clangini
    //fprintf(FPaOut,"Data for the fragment type %d :\n\n",CurFraTot); //Moved after fragment reading. clangini


    /* times(&timevar); */
    /* time_6=timevar.tms_utime+timevar.tms_stime; */
    gettimeofday(&time_6,NULL);

    /* Read the fragment file of the current fragment (CurFra) without the
       coordinates */
    /* CLANGINI 2016 */

    /* ReFrFi_mol2 has been reimplemented in C++ clangini */
    /* Reads the next valid molecule. If it detects any problems it skips the molecule. clangini */
    LstFra_f = ReFrFi_mol2(&FrInStream,&FrInPos,&SkiFra,&CurFraTot,FragNa,
                           FragNa_str,
                           &FrAtNu,&FrBdNu,&FrAtEl,&FrCoor,&FrAtTy,&FrSyAtTy,
                           &FrPaCh,&FrBdAr,&FrBdTy,FrSubN,FrSubC,&FrCoNu,
                           &SubNa,AlTySp);
    if (LstFra_f){
      std::cerr << "\tNo more fragments in input file " << FrFiNa << std::endl;
      fclose(FPaOut);
      break;
    }

    //CurFraTot++; //This is updated in ReFrFi_mol2
    if ( ++FragNa_map[FragNa_str] > 1){
      sprintf(buff, "$%d",FragNa_map[FragNa_str]);
      strcat(FragNa,buff);
    }

    CurFra++; /* After succesfully reading, we update the current fragment index */
    i = CurFra; /* Not to change i to CurFra everywhere in the following clangini*/
    fprintf(FPaOut,"-------------------------------------------------\n\n");//Moved here clangini
    fprintf(FPaOut,"Data for the fragment type %d :\n\n",CurFraTot); //Moved here clangini
    // Here we rotate the fragment to a common framework:
    FrAlRef = convert2pp_return(1,3,1,3,FrAlRef_m[0],FrAlRef_rows);
    FrAlSet = convert2pp_return(1,3,1,3,FrAlSet_m[0],FrAlSet_rows);
    align_flag = set_align_ref(FrCoor,FrAtNu,FrAlRef,FrAlSet); //setting the reference
    //print_pp(1,3,1,3,FrAlRef);
    //align_flag = false;
    if(align_flag && (EvalEn[0] == 'd')){
      // Find rotation quaternion
      //std::cout<<"Old: "<<FrCoor[1][1]<<" "<<FrCoor[1][2]<<" "<<FrCoor[1][3]<<" "<<std::endl;
#ifdef DEBUG_CL
      sprintf(WriPat,"%s","rotate_test.mol2\0"); // clangini
      FilePa=fopen(WriPat,"a");
      append_pose_to_mol2(FilePa,FragNa,FrAtNu,FrBdNu,1,FrAtEl,FrCoor,
                          1,FrSyAtTy,FrAtTy,CurFra,FrBdAr,
                          FrBdTy,1,0.0,
                          FrPaCh,SubNa,AlTySp);
      fclose(FilePa);
#endif
      struct_align(FrCoor,FrAtNu,FrAlRef,FrAlSet,3);
#ifdef DEBUG_CL
      sprintf(WriPat,"%s","rotate_test.mol2\0"); // clangini
      FilePa=fopen(WriPat,"a");
      append_pose_to_mol2(FilePa,FragNa,FrAtNu,FrBdNu,1,FrAtEl,FrCoor,
                          1,FrSyAtTy,FrAtTy,CurFra,FrBdAr,
                          FrBdTy,1,0.0,
                          FrPaCh,SubNa,AlTySp);
      fclose(FilePa);
      //std::cout<<"New: "<<FrCoor[1][1]<<" "<<FrCoor[1][2]<<" "<<FrCoor[1][3]<<" "<<std::endl;
      //std::cout<<"New: "<<FrCoor[2][1]<<" "<<FrCoor[2][2]<<" "<<FrCoor[2][3]<<" "<<std::endl;
      //std::cout<<"New: "<<FrCoor[3][1]<<" "<<FrCoor[3][2]<<" "<<FrCoor[3][3]<<" "<<std::endl;
#endif
    }
    else if(align_flag && (EvalEn[0] == 'e')){
      FrCoor_NoAlign=dmatrix(1,FrAtNu,1,3);//original coordinates before
        //alignment. clangini

      for (int ii=1;ii <= FrAtNu;ii++){
        for(int jj=1;jj <=3; jj++){
          FrCoor_NoAlign[ii][jj] = FrCoor[ii][jj];
        }
      }
      struct_align(FrCoor,FrAtNu,FrAlRef,FrAlSet,3);
    }
    else if (!align_flag && (EvalEn[0] == 'e')){
      FrCoor_NoAlign=dmatrix(1,FrAtNu,1,3);//original coordinates before
        //alignment. clangini

      for (int ii=1;ii <= FrAtNu;ii++){
        for(int jj=1;jj <=3; jj++){
          FrCoor_NoAlign[ii][jj] = FrCoor[ii][jj];
        }
      }
      std::cout << "Fragment " << FragNa
      << " for energy evaluation run has not been pre-aligned." << std::endl;
    }
    else {
      std::cout << "Fragment " << FragNa << " has not been pre-aligned." << std::endl;
    }
    /* CLANGINI 2016 end */

    /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
    /* Lo1 : Loop over the conformations     !! No further indentation */
    /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
    /* Now we will have only one conformation per fragment. Different conformations
       are considered different fragments. clangini */

    for (Ind_num_cn=1;Ind_num_cn<=FrCoNu;Ind_num_cn++) { /*This loop is no longer needed clangini*/
      /*Ind_num_cn = 1; clangini */

      /* Read the coordinates of the current conformation Ind_num_cn of the current
         fragment type CurFra */
      /* ReFrFi_coo_mol2(CurFra,Ind_num_cn,FrFiNa,FrAtNu,FrCoor); clangini*/

      /*fprintf(FPaOut,"\nResidue name of the fragment   %s\n",&FrFiNa_out[CurFra][1]); clangini */
      fprintf(FPaOut,"\nFragment name is: %s\n", FragNa);

      /* Compute the maximal number of tolerated bumps */
      BumpMa=ffloor(ScMaBump*((double) FrAtNu));
      if ((EvalEn[0]=='d')&&((Solv_typ[0]=='s')||(Solv_typ[0]=='b')))
        fprintf(FPaOut,"\nMaximal number of tolerated bumps : %d\n\n",BumpMa);

      /* Assign the atom element numbers, van der Waals radii and energies
         for the current fragment */
      FrVdWR=dvector(1,FrAtNu);
      FrVdWE=dvector(1,FrAtNu);
      if (Ind_num_cn==1)
        FrAtEl_nu=ivector(1,FrAtNu);
      FrAtoTyp_nu=ivector(1,FrAtNu);
      if(!AssiRE(NumbAT,AtTyAr,VdWRad,VdWEne,FrAtNu,FrAtTy,FrVdWR,FrVdWE,FPaOut,
            AtENAr,FrAtEl_nu,FrAtoTyp_nu))
      {
        /* bail on missing types */
        printf("Program exits\n");
        exit(13);
      }

      /*---------------------------------------------------------
        PART NEEDED FOR CONTINUUM ELECTROSTATICS
        ---------------------------------------------------------*/

      // clangini 30-01-17 START
      /*Rotate the fragment to a reference frame in order to ensure
      reproducibility when seeding the same fragment prepared with
      with different absolute coordinates.*/

      // clangini 30-01-17 END

      /* Get the charge radii and the internal distances for the current fragment */
      FrRad=dvector(1,FrAtNu);
      FrRad2=dvector(1,FrAtNu);
      FrRadOut=dvector(1,FrAtNu);
      FrRadOut2=dvector(1,FrAtNu);
      Frdist2=dmatrix(1,FrAtNu,1,FrAtNu);
      nn = get_Ch_Rad_Fr(FrAtNu,FrCoor,FrVdWR,WaMoRa,FrRad,FrRad2,FrRadOut,
          FrRadOut2,Frdist2,&Rmin,&FrRmax);
          //Frdist2 is passed as **double and not as ***double as it is
          //allocated in the main and not in get_Ch_Rad_Fr. clangini

      /* Get the extremes of frag coor */
      vect=dvector(1,FrAtNu);
      getColumnFrom2DArray(FrCoor, 1, 1, FrAtNu, vect);
      FrMinC[1] = MinVector(vect, 1, FrAtNu);
      FrMaxC[1] = MaxDVector(vect, 1, FrAtNu);
      getColumnFrom2DArray(FrCoor, 2, 1, FrAtNu, vect);
      FrMinC[2] = MinVector(vect, 1, FrAtNu);
      FrMaxC[2] = MaxDVector(vect, 1, FrAtNu);
      getColumnFrom2DArray(FrCoor, 3, 1, FrAtNu, vect);
      FrMinC[3] = MinVector(vect, 1, FrAtNu);
      FrMaxC[3] = MaxDVector(vect, 1, FrAtNu);
      free_dvector(vect,1,FrAtNu);

      if ((Solv_typ[0]=='s')||(Solv_typ[0]=='b')) {
        /* Make all the preparation necessary for the slow and accurate continuum
           elec algorithm */

        /* Place points over the frag surface to describe its SAS1 (to get
           the volume occupied by the fragment - slow algorithm) */
        surfpt_fr=structpointvect(1,NPtSphere*FrAtNu);
        iatsrf_fr=ivector(1,NPtSphere*FrAtNu);
        pointsrf_fr=ivector(1,FrAtNu);
        for (ij=1;ij<=FrAtNu;ij++)
          pointsrf_fr[ij]=0;
        nsurf_fr=ivector(1,FrAtNu);
        for (iat=1;iat<=FrAtNu;iat++)
          nsurf_fr[iat] = 0;
        nn = Surf_Grid(FrAtNu,FrCoor,FrMaxC,FrMinC,FrRad,WaMoRa,NPtSphere,
            &Nsurfpt_fr,surfpt_fr,iatsrf_fr,nsurf_fr,pointsrf_fr);
        /* Now allocate space for the SAS1 surface resulting from the fragment
           together with the receptor neighbour atoms ("extended fragment")
           It will be used in ElecFrag */
        surfpt_ex=structpointvect(1,NPtSphere * (ReAtNu+FrAtNu));

        /* Calculate the solvation energy of the fragment without the receptor
           Use a very fine grid spacing (0.1) */
        nn = FragSolvEn(FrAtNu,FrCoor,FrPaCh,FrVdWR,FrRadOut,
            FrRadOut2,Frdist2,Nsurfpt_fr,surfpt_fr,WaMoRa,0.1,
            Ksolv,pi4,&FrSolvEn,EmpCorrB,FPaOut);
        /* fprintf(FPaOut,"Dielectric value = %4.1f -> %4.1f transfer energy for the fragment (%s) : ",
            DielRe,DielWa,&ResN_fr[i][1]); */
        fprintf(FPaOut,"Dielectric value = %4.1f -> %4.1f transfer energy for the fragment (%s) : ",
            DielRe,DielWa,FragNa); /*clangini*/
        fprintf(FPaOut,"%10.4f\n",FrSolvEn);
        /*    printf("Solv En. = %f\n",FrSolvEn); */
      }

      if ((Solv_typ[0]=='f')||(Solv_typ[0]=='b')||(Solv_typ[0]=='p')) {
        /* Make all the preparation necessary for the fast and approximate continuum
           elec algorithm */

        /* Place points over the frag surface to describe its SAS2 (for a
           fast evaluation of the frag desolvation) */
        NPtSphereMax_Fr = (int) (SurfDens_deso * pi4 * (FrRmax+WaMoRa) * 2);
        //Why we do not use 4*pi*radius^2 ??? clangini
        /* * 2 necessary otherwise seg.fault with e.g. tribrommethane */

        surfpt_fr_deso_orig=structpointvect(1,NPtSphereMax_Fr*FrAtNu);
        /* surfpt_fr_deso will be used later in ElecFragFast */
        surfpt_fr_deso=structpointvect(1,NPtSphereMax_Fr*FrAtNu);
        iatsrf_fr_deso=ivector(1,NPtSphereMax_Fr*FrAtNu);
        pointsrf_fr_deso=ivector(1,FrAtNu);
        for (ij=1;ij<=FrAtNu;ij++)
          pointsrf_fr_deso[ij]=0;
        nsurf_fr_deso=ivector(1,FrAtNu);
        for (iat=1;iat<=FrAtNu;iat++)
          nsurf_fr_deso[iat] = 0;


        nn = Surf_Grid_unif_dens(FrAtNu,FrCoor,FrMaxC,FrMinC,FrRad,
            FrRadOut2,WaMoRa,SurfDens_deso,
            pi4,&Nsurfpt_fr_deso,surfpt_fr_deso_orig,
            iatsrf_fr_deso,nsurf_fr_deso,
            pointsrf_fr_deso);

        //std::cout<<"Number of points on SAS: "<<Nsurfpt_fr_deso<<std::endl;

        /* Calculate the frag desolvation resulting from the removal of a sphere of
           water having the center on all the surfpt_fr_deso_orig SAS2 points */
        nn = FragDesoSurf(FrAtNu,FrCoor,FrPaCh,FrVdWR,FrRad2,FrRadOut,
            FrRadOut2,Nsurfpt_fr_deso,surfpt_fr_deso_orig,
            WaMoRa,0.3,DielRe,DielWa,pi4,corr_fast_deso,
            corr_re_desoco,&FrSurf_deso);
      }
      /*---------------------------------------------------------
        END OF THE PART NEEDED FOR CONTINUUM ELECTROSTATICS
        ---------------------------------------------------------*/

      /* Find the aliphatic hydrogens */
      if (Ind_num_cn==1)
        AliHyd=ivector(1,FrAtNu);
      FiAlHy(FrAtNu,FrBdNu,FrAtEl_nu,FrBdAr,AliHyd,&NonAliHy);

      /* Find the hybridization state for the carbon atoms of the current fragment */
      HybFrAt=ivector(1,FrAtNu);
      HybStAt(FrAtNu,FrAtEl_nu,FrCoor,FrBdNu,FrBdAr,HybFrAt,FPaOut);

      /* Compute the square root of the van der Waals current fragment energies
         and the power of 3 of the fragment van der Waals radii */
      FrVdWE_sr=dvector(1,FrAtNu);
      SqRoEn(FrAtNu,FrVdWE,FrVdWE_sr);
      FrVdWR_p3=dvector(1,FrAtNu);
      for (j=1;j<=FrAtNu;j++)
        FrVdWR_p3[j]=FrVdWR[j]*FrVdWR[j]*FrVdWR[j];

      /* Construct the list of the donor and receptor vectors for the current
         fragment */
      FrAcDo(FrAtNu,FrAtEl_nu,FrCoor,HybFrAt,FrBdNu,FrBdAr,&FrAcNu,&FrDoNu,
          &FrDATy,&FrDAAt,&FrHydN,&FrVeCo,FPaOut,OutFil); //Couldn't this be done
          //only if polar or polar/apolar seeding is requested? clangini

      /* Find the symmetries in the fragment i.e. find the undistinguishable fragment
         atoms */
      UndisAt_fr=ivector(1,FrAtNu);
      FindFrSym(FrAtNu,FrCoor,FrAtTy,UndisAt_fr);
      /*fprintf(FPaOut,"Vectors considered or not (1,0) on the frag. atoms");
      fprintf(FPaOut," because of rotat. symmetry :\n");
      for (j=1;j<=FrAtNu;j++)
        fprintf(FPaOut,"%7d %7d\n",j,UndisAt_fr[j]);
      fprintf(FPaOut,"\n");*/ //clangini commented

      fclose(FPaOut);

      //std::cout << "We evaluated the vectors" << std::endl; //clangini OK

      FPaOut=fopen(OutFil,"a");

      SeFrCo=dmatrix(1,FrAtNu,1,3);
      RoSFCo=dmatrix(1,FrAtNu,1,3);
      SDFrRe_ps=dmatrix(1,FrAtNu,1,ReAtNu);
      SDFrRe_ps_elec=dmatrix(1,FrAtNu,1,ReAtNu);
      ChFrRe_ps_elec=dmatrix(1,FrAtNu,1,ReAtNu);

      if (Ind_num_cn==1)
        FrNuVdW_ap=-1;
      //std::cout << "Before selecting ApPoChoi mode" << std::endl; // clangini OK
      //std::cout << "FrAcNu+FrDoNu: "<< FrAcNu+FrDoNu << std::endl; // clangini OK
      //std::cout << "EvalEn[0]: "<< EvalEn[0] << std::endl; // clangini OK
      //std::cout << "ApPoChoi[0]: "<< ApPoChoi[0] << std::endl; // clangini OK
      if (((FrAcNu+FrDoNu)!=0)&&(EvalEn[0]!='e')&&
          (ApPoChoi[0]=='p')) {
        /* --------------------------- */
        /* Seeding of a polar fragment */
        /* --------------------------- */
        fprintf(FPaOut,"This fragment is considered as a polar fragment\n\n");

        /* Compute the total number of fragment positions */
        ToNFrP=0;
        for (j=1;j<=ReAcNu+ReDoNu;j++) {
          for (k=1;k<=FrAcNu+FrDoNu;k++) {
            if ((ReDATy[j]!=FrDATy[k])&&(UndisAt_fr[FrDAAt[k]])&&
                (PolVect_rec_01[j]))
              ToNFrP=ToNFrP+1;
          }
        }
        ToNFrP=ToNFrP*NuRoAx*FrCoNu;

        if (Ind_num_cn==1) {
          FrCoPo=ppdvector(1,ToNFrP);
          ConfArr=ivector(1,ToNFrP);
        }

        if (Ind_num_cn==1) {
          SFWrNu=0;
          FrNuTo=0;
          FrNuPB=0;
          FrNuSp=0;
        }
        for (j=1;j<=ReAcNu+ReDoNu;j++) { //j loops over receptor vectors. clangini
          for (k=1;k<=FrAcNu+FrDoNu;k++) {//k loops over fragment vectors. clangini
            if ((ReDATy[j]!=FrDATy[k])&&     /* vector types must be different */
                (UndisAt_fr[FrDAAt[k]])&&    /* atom accepted (symmetry) */
                (PolVect_rec_01[j])) {       /* rec. vector accept. (probe vdW) */

              RecHyN=ReHydN[j];
              FraHyN=FrHydN[k];

              /* Seed fragment */
              SeedFr(j,ReVeCo,k,FrVeCo,FrAtNu,FrCoor,ReDATy,SeFrCo,FrAtoTyp_nu,
                  FrDAAt,ReAtoTyp_nu,ReDAAt,BLAtTy,FPaOut);

#ifdef USE_QUATERNION_ROTATION //clangini
              //AnglRo = 6.283185307179586476925286766559/NuRoAx;
              AnglRo = M_PI*2/NuRoAx;
              //AnglRo=6.2831854/NuRoAx;
              //AnglRo=boost::math::constants::pi<double>()*2/NuRoAx;

              SeFrAx[1]=ReCoor[ReDAAt[j]][1]-SeFrCo[FrDAAt[k]][1];
              SeFrAx[2]=ReCoor[ReDAAt[j]][2]-SeFrCo[FrDAAt[k]][2];
              SeFrAx[3]=ReCoor[ReDAAt[j]][3]-SeFrCo[FrDAAt[k]][3];
              NormVe(&SeFrAx[1],&SeFrAx[2],&SeFrAx[3]);
              q_SeFr.fromAngleAxis(AnglRo,SeFrAx);

              // Test for precision problems
              //Quaternion<double> aaa;
              //aaa.fromAngleAxis(0.1,SeFrAx);

              for (int ii=1;ii<=FrAtNu;ii++) {
                RoSFCo[ii][1] = SeFrCo[ii][1];
                RoSFCo[ii][2] = SeFrCo[ii][2];
                RoSFCo[ii][3] = SeFrCo[ii][3];
              }

              // Test for precision problems
              /*for (int ii=1;ii<=FrAtNu;ii++){
                aaa.quatConjugateVecRef(RoSFCo[ii],SeFrCo[FrDAAt[k]]);
              }*/

#ifdef DEBUG_CL
              sprintf(WriPat,"%s","quaternion_start.mol2\0"); // clangini
              FilePa=fopen(WriPat,"a");
              append_pose_to_mol2(FilePa,FragNa,FrAtNu,FrBdNu,1,FrAtEl,RoSFCo,
                                  1,FrSyAtTy,FrAtTy,CurFra,FrBdAr,
                                  FrBdTy,1,0.0,
                                  FrPaCh,SubNa,AlTySp);
              fclose(FilePa);
              double **RoSFCo_dbg;
              RoSFCo_dbg=dmatrix(1,FrAtNu,1,3);
              for (int ii=1;ii<=FrAtNu;ii++) {
                //Need to add back the
                RoSFCo_dbg[ii][1] = SeFrCo[ii][1];
                RoSFCo_dbg[ii][2] = SeFrCo[ii][2];
                RoSFCo_dbg[ii][3] = SeFrCo[ii][3];
              }
              //double AnglRo_double = 6.283185307179586476925286766559/NuRoAx;
              int lll = NuRoAx;
              double AnglRo_double = M_PI*2/NuRoAx;
              //double AnglRo_double = M_PI*2/NuRoAx*lll;
              q_SeFr.fromAngleAxis(AnglRo_double,SeFrAx);
              //std::cout << q_SeFr <<std::endl;
              for (int ii=1; ii <= NuRoAx; ii++){
                for (int jj=1;jj<=FrAtNu;jj++){
                  q_SeFr.quatConjugateVecRef<double>(RoSFCo_dbg[jj],
                                      SeFrCo[FrDAAt[k]][1],
                                      SeFrCo[FrDAAt[k]][2],
                                      SeFrCo[FrDAAt[k]][3]);
                }
                //std::cout << q_SeFr <<std::endl;
                //q_SeFr.norm_inplace();
              }
              sprintf(WriPat,"%s","quaternion_end.mol2\0"); // clangini
              FilePa=fopen(WriPat,"a");
              append_pose_to_mol2_double(FilePa,FragNa,FrAtNu,FrBdNu,1,FrAtEl,
                                  RoSFCo_dbg,
                                  1,FrSyAtTy,FrAtTy,CurFra,FrBdAr,
                                  FrBdTy,1,0.0,
                                  FrPaCh,SubNa,AlTySp);
              fclose(FilePa);
              free_dmatrix(RoSFCo_dbg,1,FrAtNu,1,3);
              q_SeFr.fromAngleAxis(AnglRo,SeFrAx);
#endif

#else
#ifdef DEBUG_CL
              sprintf(WriPat,"%s","euler_start.mol2\0"); // clangini
              FilePa=fopen(WriPat,"a");
              append_pose_to_mol2(FilePa,FragNa,FrAtNu,FrBdNu,1,FrAtEl,SeFrCo,
                                  1,FrSyAtTy,FrAtTy,CurFra,FrBdAr,
                                  FrBdTy,1,0.0,
                                  FrPaCh,SubNa,AlTySp);
              fclose(FilePa);

              double **RoSFCo_dbg;
              RoSFCo_dbg=dmatrix(1,FrAtNu,1,3);
              for (int ii=1;ii<=FrAtNu;ii++) {
                //Need to add back the
                RoSFCo_dbg[ii][1] = SeFrCo[ii][1];
                RoSFCo_dbg[ii][2] = SeFrCo[ii][2];
                RoSFCo_dbg[ii][3] = SeFrCo[ii][3];
              }

              int lll = NuRoAx;
              double AnglRo_double = M_PI*2/NuRoAx*lll;
              RoSeFr(j,ReDAAt,ReCoor,k,FrDAAt,FrAtNu,SeFrCo,
                AnglRo_double,RoSFCo_dbg);
              sprintf(WriPat,"%s","euler_end.mol2\0"); // clangini
              FilePa=fopen(WriPat,"a");
              append_pose_to_mol2_double(FilePa,FragNa,FrAtNu,FrBdNu,1,FrAtEl,
                                  RoSFCo_dbg,
                                  1,FrSyAtTy,FrAtTy,CurFra,FrBdAr,
                                  FrBdTy,1,0.0,
                                  FrPaCh,SubNa,AlTySp);
              fclose(FilePa);
              free_dmatrix(RoSFCo_dbg,1,FrAtNu,1,3);
#endif
#endif

              /* Make rotations */
              for (l=1;l<=NuRoAx;l++) {
#ifdef USE_QUATERNION_ROTATION //clangini
                FrNuTo = FrNuTo+1;
                for (int ii=1;ii<=FrAtNu;ii++){
                  q_SeFr.quatConjugateVecRef(RoSFCo[ii],SeFrCo[FrDAAt[k]]);
                }
#else
                //AnglRo=6.2831854/NuRoAx*l; clangini
                //AnglRo = 6.283185307179586476925286766559/NuRoAx*l; //clangini
                AnglRo = M_PI*2/NuRoAx*l; //clangini
                FrNuTo=FrNuTo+1;

                /* RoSFCo is the position of the rotated fragment */
                RoSeFr(j,ReDAAt,ReCoor,k,FrDAAt,FrAtNu,SeFrCo,AnglRo,RoSFCo);
#endif

                VW_f=0.0;
                VW_s=0.0;
                In_f=0.0;
                In_s=0.0;
                Dr_f=0.0;
                Dr_s=0.0;
                Df_f=0.0;
                Df_s=0.0;
                To_f=0.0;
                To_s=0.0;


                /* ----- SLOW or BOTH methods ----- */
                if ((Solv_typ[0]=='s')||(Solv_typ[0]=='b')) {

                  SpPoCh_bool=1;

                  /* Only select a fragment if it is in the specified sphere. In order to do that,
                     compute the geometric center of the fragment and its distance to the sphere
                     center */
                  if (SpPoCh_opt[0]=='y') {

                    SpPoCh_gc[1]=0.0;
                    SpPoCh_gc[2]=0.0;
                    SpPoCh_gc[3]=0.0;
                    for (i1=1;i1<=FrAtNu;i1++) {
                      SpPoCh_gc[1]=SpPoCh_gc[1]+RoSFCo[i1][1];
                      SpPoCh_gc[2]=SpPoCh_gc[2]+RoSFCo[i1][2];
                      SpPoCh_gc[3]=SpPoCh_gc[3]+RoSFCo[i1][3];
                    }
                    SpPoCh_gc[1]=SpPoCh_gc[1]/FrAtNu;
                    SpPoCh_gc[2]=SpPoCh_gc[2]/FrAtNu;
                    SpPoCh_gc[3]=SpPoCh_gc[3]/FrAtNu;

                    SpPoCh_dist=DistSq(SpPoCh_gc[1],SpPoCh_gc[2],SpPoCh_gc[3],
                        SpPoCh_cent[1],SpPoCh_cent[2],SpPoCh_cent[3]);

                    if (SpPoCh_dist<=SpPoCh_rad_sq)
                      FrNuSp=FrNuSp+1;
                    else
                      SpPoCh_bool=0;

                  }

                  if (SpPoCh_bool) {

                    /* Compute the # bumps */
                    FrAcce=CoNuBu(CubNum,CubFAI,CubLAI,CubLiA,LaVdWR,ReMaxC,ReMinC,
                        ReCoor,ReVdWR,FrAtNu,RoSFCo,FrVdWR,VdWFaB,BumpMa,
                        RecHyN,FraHyN,BuEvFa);

                    if (FrAcce) {

                      FrNuPB=FrNuPB+1;

                      Rot_Tran(FrAtNu,FrCoor,RoSFCo,Tr,U1,U2);

                      /* Compute the squared distances between the fragment atoms and the
                         corresponding receptor atoms of the pseudo-sphere procedure, one list
                         for the vdW energy and another for the electrostatic interaction energy.
                         In the case of the electrostatic interaction energy, it also provides the
                         corresponding partial or total charges */
                      SqDisFrRe_ps(FrAtNu,RoSFCo,ReCoor,ReMinC,GrSiCu_en,
                          CubNum_en,CubFAI_en,CubLAI_en,CubLiA_en,
                          PsSpNC,PsSphe,SDFrRe_ps,ReAtNu,PsSpRa,
                          RePaCh,ReReNu,AtReprRes,FiAtRes,LaAtRes,
                          TotChaRes,NuChResEn,LiChResEn,
                          SDFrRe_ps_elec,ChFrRe_ps_elec);

                      /* Compute vdW energy */
                      PsSpEE(FrAtNu,ReAtNu,ReVdWE_sr,FrVdWE_sr,
                          ReVdWR,FrVdWR,&VWEnEv_ps,SDFrRe_ps);

                      /* Compute receptor desolvation, fragment desolvation and receptor-fragment
                         interaction (with screening effect) energies (slow method) */
                      ElecFrag(ReAtNu,ReCoor,RePaCh,ChFrRe_ps_elec,
                          ReRad,ReRad2,ReRadOut,
                          ReRadOut2,surfpt_re,nsurf_re,
                          pointsrf_re,ReSelfVol,FrAtNu,RoSFCo,FrCoor,
                          FrPaCh,FrRad,FrRad2,FrRadOut,FrRadOut2,
                          Frdist2,SDFrRe_ps_elec,FrMinC,FrMaxC,&FrSolvEn,
                          Nsurfpt_fr,surfpt_fr,
                          nsurf_fr,pointsrf_fr,surfpt_ex,Tr,U1,U2,WaMoRa,
                          GrSiSo,NPtSphere,Min,Max,XGrid,YGrid,ZGrid,
                          NGridx,NGridy,NGridz,GridMat,
                          DeltaPrDeso,Kelec,Ksolv,UnitVol,
                          pi4,nxminBS,nyminBS,nzminBS,nxmaxBS,nymaxBS,
                          nzmaxBS,corr_scrint,corr_fr_deso,&ReDesoElec,
                          &ReFrIntElec,&FrDesoElec,ReSelfVol_corrB,EmpCorrB,FPaOut);

                      VW_s=SFVWEn*VWEnEv_ps;
                      In_s=SFIntElec*ReFrIntElec;
                      Dr_s=SFDeso_re*ReDesoElec;
                      Df_s=SFDeso_fr*FrDesoElec;
                      To_s=VW_s+In_s+Dr_s+Df_s;


                      /* Check energy cutoff */
                      if (To_s<=FrMaEn) {

                        SFWrNu=SFWrNu+1;

                        FrCoPo[SFWrNu]=dmatrix(1,FrAtNu,1,3);
                        for (i1=1;i1<=FrAtNu;i1++) {
                          for (i2=1;i2<=3;i2++) {
                            FrCoPo[SFWrNu][i1][i2]=RoSFCo[i1][i2];
                          }
                        }

                        ConfArr[SFWrNu]=Ind_num_cn;

                        if (Solv_typ[0]=='b') {

                          /* Van der Waals energy on a grid using the geometric mean approximation */
                          VW_f=VWFrEn(FrAtNu,RoSFCo,BSMinC,VWGrIn,VWGrSi,FrVdWR_p3,
                              FrVdWE_sr,VWGPoN,VWGrRP_at,VWGrRP_re,FPaOut,&vdwErr);

                          /* Coulombic interaction energy on a grid */
                          In_f=CoFrEn(FrAtNu,FrPaCh,RoSFCo,BSMinC,CoGrIn,CoGrSi,
                              CoGPoN,CoGrRP,FPaOut,&clbErr);

                          /* Compute receptor desolvation, fragment desolvation energies (fast method) */
                          ElecFragFast(ReAtNu,ReCoor,ReRadOut,ReRadOut2,
                              surfpt_re_deso,nsurf_re_deso,ReSurf_deso,
                              pointsrf_re_deso,nlist_re,list_re,
                              nstep_re,cbmid_re,midg_re,rscale_re,
                              FrAtNu,RoSFCo,FrCoor,FrRad,FrRadOut,
                              FrRadOut2,FrMinC,FrMaxC,FrRmax,
                              Nsurfpt_fr_deso,surfpt_fr_deso_orig,
                              surfpt_fr_deso,nsurf_fr_deso,
                              FrSurf_deso,pointsrf_fr_deso,
                              Tr,U1,U2,WaMoRa,rslop,
                              MaxNeigh,&ReDesoElec,&FrDesoElec);

                          VW_f=SFVWEn*VW_f;
                          In_f=SFIntElec*In_f;
                          Dr_f=SFDeso_re*ReDesoElec;
                          Df_f=SFDeso_fr*FrDesoElec;
                          To_f=VW_f+In_f+Dr_f+Df_f;


                        }

                        if(SFWrNu + 1>= currentsize)
                        {
                          /*  to do  */
                          /* #ifdef OMP */
                          /* #pragma omp critical */
                          /* #endif */
                          currentsize += _ALLOCSIZE;
                          VW_f_ro = dvecresize(VW_f_ro,currentsize+1); /* seed starts counting at 1 */
                          VW_s_ro = dvecresize(VW_s_ro,currentsize+1);
                          In_f_ro = dvecresize(In_f_ro,currentsize+1);
                          In_s_ro = dvecresize(In_s_ro,currentsize+1);
                          Dr_f_ro = dvecresize(Dr_f_ro,currentsize+1);
                          Dr_s_ro = dvecresize(Dr_s_ro,currentsize+1);
                          Df_f_ro = dvecresize(Df_f_ro,currentsize+1);
                          Df_s_ro = dvecresize(Df_s_ro,currentsize+1);
                          To_f_ro = dvecresize(To_f_ro,currentsize+1);
                          To_s_ro = dvecresize(To_s_ro,currentsize+1);
                          Index_ro = ivecresize(Index_ro,currentsize+1);

                        }

                        Index_ro[SFWrNu] = SFWrNu;
                        VW_f_ro[SFWrNu] = VW_f;
                        VW_s_ro[SFWrNu] = VW_s;
                        In_f_ro[SFWrNu] = In_f;
                        In_s_ro[SFWrNu] = In_s;
                        Dr_f_ro[SFWrNu] = Dr_f;
                        Dr_s_ro[SFWrNu] = Dr_s;
                        Df_f_ro[SFWrNu] = Df_f;
                        Df_s_ro[SFWrNu] = Df_s;
                        To_f_ro[SFWrNu] = To_f;
                        To_s_ro[SFWrNu] = To_s;
                      }

                    }

                  }

                }


                /* ----- FAST method ----- */
                if ((Solv_typ[0]=='f')||(Solv_typ[0]=='p')) {

                  SpPoCh_bool=1;

                  /* Only select a fragment if it is in the specified sphere. In order to do that,
                     compute the geometric center of the fragment and its distance to the sphere
                     center */
                  if (SpPoCh_opt[0]=='y') {

                    SpPoCh_gc[1]=0.0;
                    SpPoCh_gc[2]=0.0;
                    SpPoCh_gc[3]=0.0;
                    for (i1=1;i1<=FrAtNu;i1++) {
                      SpPoCh_gc[1]=SpPoCh_gc[1]+RoSFCo[i1][1];
                      SpPoCh_gc[2]=SpPoCh_gc[2]+RoSFCo[i1][2];
                      SpPoCh_gc[3]=SpPoCh_gc[3]+RoSFCo[i1][3];
                    }
                    SpPoCh_gc[1]=SpPoCh_gc[1]/FrAtNu;
                    SpPoCh_gc[2]=SpPoCh_gc[2]/FrAtNu;
                    SpPoCh_gc[3]=SpPoCh_gc[3]/FrAtNu;

                    SpPoCh_dist=DistSq(SpPoCh_gc[1],SpPoCh_gc[2],SpPoCh_gc[3],
                        SpPoCh_cent[1],SpPoCh_cent[2],SpPoCh_cent[3]);

                    if (SpPoCh_dist<=SpPoCh_rad_sq)
                      FrNuSp=FrNuSp+1;
                    else
                      SpPoCh_bool=0;

                  }

                  if (SpPoCh_bool) {

                    /* Van der Waals energy on a grid using the geometric mean approximation */
                    VW_f=VWFrEn(FrAtNu,RoSFCo,BSMinC,VWGrIn,VWGrSi,FrVdWR_p3,
                        FrVdWE_sr,VWGPoN,VWGrRP_at,VWGrRP_re,FPaOut,&vdwErr);

                    VW_f=SFVWEn*VW_f;

                    /* The van der Waals energy acts as a bump checking */
                    if (VW_f<=BumpFaCut) {

                      FrNuPB=FrNuPB+1;

                      Rot_Tran(FrAtNu,FrCoor,RoSFCo,Tr,U1,U2);

                      /* Coulombic interaction energy on a grid */
                      In_f=CoFrEn(FrAtNu,FrPaCh,RoSFCo,BSMinC,CoGrIn,CoGrSi,
                          CoGPoN,CoGrRP,FPaOut,&clbErr);

                      /* Compute receptor desolvation, fragment desolvation energies (fast method) */
#ifdef XFAST
                      if(SFIntElec*In_f +  VW_f <= FrMaEn) /*dbx check remove*/
                      {
#endif
                        ElecFragFast(ReAtNu,ReCoor,ReRadOut,ReRadOut2,
                            surfpt_re_deso,nsurf_re_deso,ReSurf_deso,
                            pointsrf_re_deso,nlist_re,list_re,
                            nstep_re,cbmid_re,midg_re,rscale_re,
                            FrAtNu,RoSFCo,FrCoor,FrRad,FrRadOut,
                            FrRadOut2,FrMinC,FrMaxC,FrRmax,
                            Nsurfpt_fr_deso,surfpt_fr_deso_orig,
                            surfpt_fr_deso,nsurf_fr_deso,FrSurf_deso,
                            pointsrf_fr_deso,
                            Tr,U1,U2,WaMoRa,rslop,
                            MaxNeigh,&ReDesoElec,&FrDesoElec);


                        Dr_f=SFDeso_re*ReDesoElec;
                        Df_f=SFDeso_fr*FrDesoElec;
#ifdef XFAST
                      }
                      else{ Df_f=1000; Dr_f=1000;} /*dbx check remove*/
#endif
                      In_f=SFIntElec*In_f;
                      To_f=VW_f+In_f+Dr_f+Df_f;

                      /* Check energy cutoff */
                      if (To_f<=FrMaEn) {

                        SFWrNu=SFWrNu+1;

                        FrCoPo[SFWrNu]=dmatrix(1,FrAtNu,1,3);
                        for (i1=1;i1<=FrAtNu;i1++) {
                          for (i2=1;i2<=3;i2++) {
                            FrCoPo[SFWrNu][i1][i2]=RoSFCo[i1][i2];
                          }
                        }

                        ConfArr[SFWrNu]=Ind_num_cn;


                        if(SFWrNu + 1>= currentsize)
                        {
                          /*  to do  */
                          /* #ifdef OMP */
                          /* #pragma omp critical */
                          /* #endif */
                          currentsize += _ALLOCSIZE;
                          VW_f_ro = dvecresize(VW_f_ro,currentsize+1);
                          VW_s_ro = dvecresize(VW_s_ro,currentsize+1);
                          In_f_ro = dvecresize(In_f_ro,currentsize+1);
                          In_s_ro = dvecresize(In_s_ro,currentsize+1);
                          Dr_f_ro = dvecresize(Dr_f_ro,currentsize+1);
                          Dr_s_ro = dvecresize(Dr_s_ro,currentsize+1);
                          Df_f_ro = dvecresize(Df_f_ro,currentsize+1);
                          Df_s_ro = dvecresize(Df_s_ro,currentsize+1);
                          To_f_ro = dvecresize(To_f_ro,currentsize+1);
                          To_s_ro = dvecresize(To_s_ro,currentsize+1);
                          Index_ro = ivecresize(Index_ro,currentsize+1);
                        }

                        Index_ro[SFWrNu] = SFWrNu;
                        VW_f_ro[SFWrNu] = VW_f;
                        VW_s_ro[SFWrNu] = VW_s;
                        In_f_ro[SFWrNu] = In_f;
                        In_s_ro[SFWrNu] = In_s;
                        Dr_f_ro[SFWrNu] = Dr_f;
                        Dr_s_ro[SFWrNu] = Dr_s;
                        Df_f_ro[SFWrNu] = Df_f;
                        Df_s_ro[SFWrNu] = Df_s;
                        To_f_ro[SFWrNu] = To_f;
                        To_s_ro[SFWrNu] = To_s;
                        //debug clangini start
                        /*if (SFWrNu == 39 || SFWrNu == 43|| SFWrNu == 8
                        || SFWrNu == 7|| SFWrNu == 14){
                          std::cout << SFWrNu << ":  "<< l
                          << "  To_f:  " << To_f << std::endl;
                        }*/
                        //debug clangini end

                      }

                    }

                  }

                }


              }


            }
          }
        }

      }
      else if ((EvalEn[0]!='e')&&(ApPoChoi[0]=='a')) {
        /* ----------------------------- */
        /* Seeding of an apolar fragment */
        /* ----------------------------- */
        fprintf(FPaOut,"This fragment is considered as an apolar fragment\n\n");

        /* ---------------------------------------------------------------*/
        /* This part is needed to generate VECTOR FOR APOLAR FRAGMENTS !! */

        /* Make another grid on the surface of the fragment to describe its SAS
           (this one will be used for the vectors of apolar fragments) */

        /*      printf("\n\tGeneration of the SAS for the vectors of apolar
                fragments\n"); */

        /*
           dey mod 26072006 :
           NPtSphereMax_Fr TOO SMALL with :
           NPtSphereMax_Fr = (int) (PtDensFr * pi4 * (FrRmax+WaMoRa));
           -->> seg. fault with e.g. dibrommethane

thus :
NPtSphereMax_Fr = (int) (PtDensFr * pi4 * (FrRmax+WaMoRa) * 2);
eventually also adapt :
NPtSphereMax = (int) (1. * pi4 * (ReRmax+WaMoRa) * (ReRmax+WaMoRa));
NPtSphereMax_Fr = (int) (SurfDens_deso * pi4 * (FrRmax+WaMoRa));
         */
        NPtSphereMax_Fr = (int) (PtDensFr * pi4 * (FrRmax+WaMoRa) * 2);
        /* * 2 necessary otherwise seg.fault with e.g. tribrommethane */
        //std::cout<<"NPtSphereMax_Fr: "<<NPtSphereMax_Fr<<std::endl; clangini
        //std::cout<<"dens*4pi*r^2: "<< //clangini
        //   (int) (PtDensFr * pi4 * (FrRmax+WaMoRa) * (FrRmax+WaMoRa))<<std::endl;

        surfpt_fr_apol=structpointvect(1,NPtSphereMax_Fr*FrAtNu);
        iatsrf_fr_apol=ivector(1,NPtSphereMax_Fr*FrAtNu);
        pointsrf_fr_apol=ivector(1,FrAtNu);
        for (ij=1;ij<=FrAtNu;ij++)
          pointsrf_fr_apol[ij]=0;
        nsurf_fr_apol=ivector(1,FrAtNu);
        for (iat=1;iat<=FrAtNu;iat++)
          nsurf_fr_apol[iat] = 0;


        nn = Surf_Grid_unif_dens(FrAtNu,FrCoor,FrMaxC,FrMinC,FrRad,FrRadOut2,
            Sphere_apol,PtDensFr,pi4,&Nsurfpt_fr_apol,
            surfpt_fr_apol,iatsrf_fr_apol,
            nsurf_fr_apol,pointsrf_fr_apol);

        apol_Vect_fr = dmatrix(1,Nsurfpt_fr_apol,1,6);
        FrApAt=ivector(1,Nsurfpt_fr_apol);
        for (isph=1;isph<=Nsurfpt_fr_apol;isph++) {
          for (j=1;j<=3;j++)
            apol_Vect_fr[isph][j+3]=FrCoor[iatsrf_fr_apol[isph]][j];//one end of
            //the apolar vector is the position of the atom generating the vector. clangini
          apol_Vect_fr[isph][1] = (double) surfpt_fr_apol[isph].x;
          apol_Vect_fr[isph][2] = (double) surfpt_fr_apol[isph].y;
          apol_Vect_fr[isph][3] = (double) surfpt_fr_apol[isph].z;
          FrApAt[isph]=iatsrf_fr_apol[isph];
        }

        free_ivector(pointsrf_fr_apol,1,FrAtNu);
        free_ivector(nsurf_fr_apol,1,FrAtNu);
        free_ivector(iatsrf_fr_apol,1,NPtSphereMax_Fr*FrAtNu);
        free_structpointvect(surfpt_fr_apol,1,NPtSphereMax_Fr*FrAtNu);

        /* Normalization of the vectors apol_Vect_fr */
        for (i1=1;i1<=Nsurfpt_fr_apol;i1++)
          PoCoVe(apol_Vect_fr[i1][4],apol_Vect_fr[i1][5],apol_Vect_fr[i1][6],
              apol_Vect_fr[i1][1],apol_Vect_fr[i1][2],apol_Vect_fr[i1][3],
              1.0,&(apol_Vect_fr[i1][1]),&(apol_Vect_fr[i1][2]),
              &(apol_Vect_fr[i1][3]));

        fprintf(FPaOut,"Numb. of frag. vectors (surf. points) used to seed ");
        fprintf(FPaOut,"this apolar fragment : %d\n",Nsurfpt_fr_apol);
        fprintf(FPaOut,"This number is then reduced if undistinguishable");
        fprintf(FPaOut," fragment atoms\n\n");

        /*
            FilePa=fopen("./outputs/apolar_frag.mol2","w");
            fprintf(FilePa,"# TRIPOS MOL2 file generated by WITNOTP\n");
            fprintf(FilePa,"\n");
            fprintf(FilePa,"@<TRIPOS>MOLECULE\n");
            fprintf(FilePa,"apolar_frag\n");
            fprintf(FilePa,"%d 0 0 0 0\n",2*Nsurfpt_fr_apol);
            fprintf(FilePa,"****\n");
            fprintf(FilePa,"USER_CHARGES\n");
            fprintf(FilePa,"INVALID_CHARGES\n");
            fprintf(FilePa,"\n");
            fprintf(FilePa,"@<TRIPOS>ATOM\n");

            for (i1=1;i1<=Nsurfpt_fr_apol;i1++) {
            fprintf(FilePa,"%-5d %-3d%11.4f%11.4f%11.4f ***** 1    ****  0.000\n",
            i1*2-1,i1*2-1,apol_Vect_fr[i1][1],apol_Vect_fr[i1][2],
            apol_Vect_fr[i1][3]);
            fprintf(FilePa,"%-5d %-3d%11.4f%11.4f%11.4f ***** 1    ****  0.000\n",
            i1*2,i1*2,apol_Vect_fr[i1][4],apol_Vect_fr[i1][5],
            apol_Vect_fr[i1][6]);
            }
            fclose(FilePa);
         */

        /* End of the part needed to generate VECTOR FOR APOLAR FRAGMENTS !! */
        /* ------------------------------------------------------------------*/


        VdWCoff_ap=MuFaVdWCoff_ap*((double) FrAtNu);

        /* Compute the total number of fragment positions */
        ToNFrP=0;
        for (j=1;j<=Napol_Vect_re;j++) {
          for (k=1;k<=Nsurfpt_fr_apol;k++) {
            if (UndisAt_fr[FrApAt[k]])
              ToNFrP=ToNFrP+1;//all possible combinations of receptor and
              //fragment vectors. clangini
          }
        }
        ToNFrP=ToNFrP*NuRoAx*FrCoNu;

        if (Ind_num_cn==1) {
          FrCoPo=ppdvector(1,ToNFrP);
          ConfArr=ivector(1,ToNFrP);
        }

        if (Ind_num_cn==1) {
          SFWrNu=0;
          FrNuTo=0;
          FrNuPB=0;
          FrNuSp=0;
          FrNuVdW_ap=0;
        }
        for (j=1;j<=Napol_Vect_re;j++) {
          for (k=1;k<=Nsurfpt_fr_apol;k++) {
            if (UndisAt_fr[FrApAt[k]]) {        /* accepted atom (symmetry) */

              /* Seed fragment */
              SeedFr_ap(j,apol_Vect_re,ReApAt,k,apol_Vect_fr,FrApAt,FrCoor,
                  ReVdWR,FrVdWR,FrAtNu,FPaOut,SeFrCo);
#ifdef USE_QUATERNION_ROTATION //clangini
              //AnglRo = 6.283185307179586476925286766559/NuRoAx;
              AnglRo = M_PI*2/NuRoAx;//+0.0000000000000001;
              //AnglRo = 6.2831853/NuRoAx; //clangini debug
              SeFrAx[1]=ReCoor[ReApAt[j]][1]-SeFrCo[FrApAt[k]][1];
              SeFrAx[2]=ReCoor[ReApAt[j]][2]-SeFrCo[FrApAt[k]][2];
              SeFrAx[3]=ReCoor[ReApAt[j]][3]-SeFrCo[FrApAt[k]][3];
              NormVe(&SeFrAx[1],&SeFrAx[2],&SeFrAx[3]);
              q_SeFr.fromAngleAxis(AnglRo,SeFrAx);

              for (int ii=1;ii<=FrAtNu;ii++) {
                //Need to add back the
                RoSFCo[ii][1] = SeFrCo[ii][1];
                RoSFCo[ii][2] = SeFrCo[ii][2];
                RoSFCo[ii][3] = SeFrCo[ii][3];
              }
#endif

              /* Make rotations */
              for (l=1;l<=NuRoAx;l++) {
#ifdef USE_QUATERNION_ROTATION //clangini
                FrNuTo = FrNuTo+1;
                for (int ii=1;ii<=FrAtNu;ii++){
                  q_SeFr.quatConjugateVecRef(RoSFCo[ii],SeFrCo[FrApAt[k]]);
                }
#else
                AnglRo = M_PI*2/NuRoAx*l;

                FrNuTo=FrNuTo+1;

                /* RoSFCo is the position of the rotated fragment */
                RoSeFr(j,ReApAt,ReCoor,k,FrApAt,FrAtNu,SeFrCo,AnglRo,RoSFCo);
#endif

                VW_f=0.0;
                VW_s=0.0;
                In_f=0.0;
                In_s=0.0;
                Dr_f=0.0;
                Dr_s=0.0;
                Df_f=0.0;
                Df_s=0.0;
                To_f=0.0;
                To_s=0.0;


                /* ----- SLOW or BOTH methods ----- */
                if ((Solv_typ[0]=='s')||(Solv_typ[0]=='b')) {

                  SpPoCh_bool=1;

                  /* Only select a fragment if it is in the specified sphere. In order to do that,
                     compute the geometric center of the fragment and its distance to the sphere
                     center */
                  if (SpPoCh_opt[0]=='y') {

                    SpPoCh_gc[1]=0.0;
                    SpPoCh_gc[2]=0.0;
                    SpPoCh_gc[3]=0.0;
                    for (i1=1;i1<=FrAtNu;i1++) {
                      SpPoCh_gc[1]=SpPoCh_gc[1]+RoSFCo[i1][1];
                      SpPoCh_gc[2]=SpPoCh_gc[2]+RoSFCo[i1][2];
                      SpPoCh_gc[3]=SpPoCh_gc[3]+RoSFCo[i1][3];
                    }
                    SpPoCh_gc[1]=SpPoCh_gc[1]/FrAtNu;
                    SpPoCh_gc[2]=SpPoCh_gc[2]/FrAtNu;
                    SpPoCh_gc[3]=SpPoCh_gc[3]/FrAtNu;

                    SpPoCh_dist=DistSq(SpPoCh_gc[1],SpPoCh_gc[2],SpPoCh_gc[3],
                        SpPoCh_cent[1],SpPoCh_cent[2],SpPoCh_cent[3]);

                    if (SpPoCh_dist<=SpPoCh_rad_sq)
                      FrNuSp=FrNuSp+1;
                    else
                      SpPoCh_bool=0;

                  }

                  if (SpPoCh_bool) {

                    /* Compute the # bumps */
                    FrAcce=CoNuBu(CubNum,CubFAI,CubLAI,CubLiA,LaVdWR,ReMaxC,ReMinC,
                        ReCoor,ReVdWR,FrAtNu,RoSFCo,FrVdWR,VdWFaB,BumpMa,
                        0,0,BuEvFa);

                    if (FrAcce) {

                      FrNuPB=FrNuPB+1;

                      /* Compute the squared distances between the fragment atoms and the
                         corresponding receptor atoms of the pseudo-sphere procedure, one list
                         for the vdW energy and another for the electrostatic interaction energy.
                         In the case of the electrostatic interaction energy, it also provides the
                         corresponding partial or total charges */
                      SqDisFrRe_ps(FrAtNu,RoSFCo,ReCoor,ReMinC,GrSiCu_en,
                          CubNum_en,CubFAI_en,CubLAI_en,CubLiA_en,
                          PsSpNC,PsSphe,SDFrRe_ps,ReAtNu,PsSpRa,
                          RePaCh,ReReNu,AtReprRes,FiAtRes,LaAtRes,
                          TotChaRes,NuChResEn,LiChResEn,
                          SDFrRe_ps_elec,ChFrRe_ps_elec);

                      /* Compute vdW energy */
                      PsSpEE(FrAtNu,ReAtNu,ReVdWE_sr,FrVdWE_sr,
                          ReVdWR,FrVdWR,&VWEnEv_ps,SDFrRe_ps);

                      VW_s=SFVWEn*VWEnEv_ps;

                      /* The vdW energy has to be smaller than a cutoff */
                      if (VW_s<=VdWCoff_ap) {

                        FrNuVdW_ap=FrNuVdW_ap+1;

                        Rot_Tran(FrAtNu,FrCoor,RoSFCo,Tr,U1,U2);

                        /* Compute receptor desolvation, fragment desolvation and receptor-fragment
                           interaction (with screening effect) energies (slow method) */
                        ElecFrag(ReAtNu,ReCoor,RePaCh,ChFrRe_ps_elec,
                            ReRad,ReRad2,ReRadOut,
                            ReRadOut2,surfpt_re,nsurf_re,
                            pointsrf_re,ReSelfVol,FrAtNu,RoSFCo,FrCoor,
                            FrPaCh,FrRad,FrRad2,FrRadOut,FrRadOut2,
                            Frdist2,SDFrRe_ps_elec,FrMinC,FrMaxC,
                            &FrSolvEn,Nsurfpt_fr,surfpt_fr,
                            nsurf_fr,pointsrf_fr,surfpt_ex,Tr,U1,U2,WaMoRa,
                            GrSiSo,NPtSphere,Min,Max,XGrid,YGrid,ZGrid,
                            NGridx,NGridy,NGridz,GridMat,
                            DeltaPrDeso,Kelec,Ksolv,UnitVol,
                            pi4,nxminBS,nyminBS,nzminBS,nxmaxBS,nymaxBS,
                            nzmaxBS,corr_scrint,corr_fr_deso,&ReDesoElec,
                            &ReFrIntElec,&FrDesoElec,ReSelfVol_corrB,EmpCorrB,FPaOut);

                        In_s=SFIntElec*ReFrIntElec;
                        Dr_s=SFDeso_re*ReDesoElec;
                        Df_s=SFDeso_fr*FrDesoElec;
                        To_s=VW_s+In_s+Dr_s+Df_s;


                        /* Check energy cutoff */
                        if (To_s<=FrMaEn) {

                          SFWrNu=SFWrNu+1;

                          FrCoPo[SFWrNu]=dmatrix(1,FrAtNu,1,3);
                          for (i1=1;i1<=FrAtNu;i1++) {
                            for (i2=1;i2<=3;i2++) {
                              FrCoPo[SFWrNu][i1][i2]=RoSFCo[i1][i2];
                            }
                          }

                          ConfArr[SFWrNu]=Ind_num_cn;

                          if (Solv_typ[0]=='b') {

                            /* Van der Waals energy on a grid using the geometric mean approximation */
                            VW_f=VWFrEn(FrAtNu,RoSFCo,BSMinC,VWGrIn,VWGrSi,FrVdWR_p3,
                                FrVdWE_sr,VWGPoN,VWGrRP_at,VWGrRP_re,FPaOut,&vdwErr);

                            /* Coulombic interaction energy on a grid */
                            In_f=CoFrEn(FrAtNu,FrPaCh,RoSFCo,BSMinC,CoGrIn,CoGrSi,
                                CoGPoN,CoGrRP,FPaOut,&clbErr);

                            /* Compute receptor desolvation, fragment desolvation energies (fast method) */
                            ElecFragFast(ReAtNu,ReCoor,ReRadOut,
                                ReRadOut2,surfpt_re_deso,nsurf_re_deso,
                                ReSurf_deso,pointsrf_re_deso,nlist_re,
                                list_re,nstep_re,cbmid_re,midg_re,
                                rscale_re,FrAtNu,RoSFCo,FrCoor,FrRad,
                                FrRadOut,FrRadOut2,FrMinC,FrMaxC,FrRmax,
                                Nsurfpt_fr_deso,surfpt_fr_deso_orig,
                                surfpt_fr_deso,nsurf_fr_deso,
                                FrSurf_deso,pointsrf_fr_deso,Tr,U1,U2,
                                WaMoRa,rslop,MaxNeigh,
                                &ReDesoElec,&FrDesoElec);

                            VW_f=SFVWEn*VW_f;
                            In_f=SFIntElec*In_f;
                            Dr_f=SFDeso_re*ReDesoElec;
                            Df_f=SFDeso_fr*FrDesoElec;
                            To_f=VW_f+In_f+Dr_f+Df_f;

                          }


                          if(SFWrNu + 1>= currentsize)
                          {
                            /*  to do  */
                            /* #ifdef OMP */
                            /* #pragma omp critical */
                            /* #endif */
                            currentsize += _ALLOCSIZE;
                            VW_f_ro = dvecresize(VW_f_ro,currentsize+1);
                            VW_s_ro = dvecresize(VW_s_ro,currentsize+1);
                            In_f_ro = dvecresize(In_f_ro,currentsize+1);
                            In_s_ro = dvecresize(In_s_ro,currentsize+1);
                            Dr_f_ro = dvecresize(Dr_f_ro,currentsize+1);
                            Dr_s_ro = dvecresize(Dr_s_ro,currentsize+1);
                            Df_f_ro = dvecresize(Df_f_ro,currentsize+1);
                            Df_s_ro = dvecresize(Df_s_ro,currentsize+1);
                            To_f_ro = dvecresize(To_f_ro,currentsize+1);
                            To_s_ro = dvecresize(To_s_ro,currentsize+1);
                            Index_ro = ivecresize(Index_ro,currentsize+1);
                          }

                          Index_ro[SFWrNu] = SFWrNu;
                          VW_f_ro[SFWrNu] = VW_f;
                          VW_s_ro[SFWrNu] = VW_s;
                          In_f_ro[SFWrNu] = In_f;
                          In_s_ro[SFWrNu] = In_s;
                          Dr_f_ro[SFWrNu] = Dr_f;
                          Dr_s_ro[SFWrNu] = Dr_s;
                          Df_f_ro[SFWrNu] = Df_f;
                          Df_s_ro[SFWrNu] = Df_s;
                          To_f_ro[SFWrNu] = To_f;
                          To_s_ro[SFWrNu] = To_s;

                        }

                      }

                    }

                  }

                }


                /* ----- FAST method ----- */
                if ((Solv_typ[0]=='f')||(Solv_typ[0]=='p')) {

                  SpPoCh_bool=1;

                  /* Only select a fragment if it is in the specified sphere. In order to do that,
                     compute the geometric center of the fragment and its distance to the sphere
                     center */
                  if (SpPoCh_opt[0]=='y') {

                    SpPoCh_gc[1]=0.0;
                    SpPoCh_gc[2]=0.0;
                    SpPoCh_gc[3]=0.0;
                    for (i1=1;i1<=FrAtNu;i1++) {
                      SpPoCh_gc[1]=SpPoCh_gc[1]+RoSFCo[i1][1];
                      SpPoCh_gc[2]=SpPoCh_gc[2]+RoSFCo[i1][2];
                      SpPoCh_gc[3]=SpPoCh_gc[3]+RoSFCo[i1][3];
                    }
                    SpPoCh_gc[1]=SpPoCh_gc[1]/FrAtNu;
                    SpPoCh_gc[2]=SpPoCh_gc[2]/FrAtNu;
                    SpPoCh_gc[3]=SpPoCh_gc[3]/FrAtNu;

                    SpPoCh_dist=DistSq(SpPoCh_gc[1],SpPoCh_gc[2],SpPoCh_gc[3],
                        SpPoCh_cent[1],SpPoCh_cent[2],SpPoCh_cent[3]);

                    if (SpPoCh_dist<=SpPoCh_rad_sq)
                      FrNuSp=FrNuSp+1;
                    else
                      SpPoCh_bool=0;

                  }

                  if (SpPoCh_bool) {


                    /* Van der Waals energy on a grid using the geometric mean approximation */
                    VW_f=VWFrEn(FrAtNu,RoSFCo,BSMinC,VWGrIn,VWGrSi,FrVdWR_p3,
                        FrVdWE_sr,VWGPoN,VWGrRP_at,VWGrRP_re,FPaOut,&vdwErr);

                    VW_f=SFVWEn*VW_f;

                    /* The van der Waals energy acts as a bump checking */
                    if (VW_f<=BumpFaCut) {

                      FrNuPB=FrNuPB+1;

                      /* The vdW energy has to be smaller than a cutoff */
                      if (VW_f<=VdWCoff_ap) {

                        FrNuVdW_ap=FrNuVdW_ap+1;

                        Rot_Tran(FrAtNu,FrCoor,RoSFCo,Tr,U1,U2);

                        /* Coulombic interaction energy on a grid */
                        In_f=CoFrEn(FrAtNu,FrPaCh,RoSFCo,BSMinC,CoGrIn,CoGrSi,
                            CoGPoN,CoGrRP,FPaOut,&clbErr);

                        /* Compute receptor desolvation, fragment desolvation energies (fast method) */
#ifdef XFAST
                        if(SFIntElec*In_f +  VW_f <= FrMaEn) /*dbx check remove*/
                        {
#endif
                          ElecFragFast(ReAtNu,ReCoor,ReRadOut,
                              ReRadOut2,surfpt_re_deso,nsurf_re_deso,
                              ReSurf_deso,pointsrf_re_deso,nlist_re,
                              list_re,nstep_re,cbmid_re,midg_re,
                              rscale_re,FrAtNu,RoSFCo,FrCoor,FrRad,
                              FrRadOut,FrRadOut2,FrMinC,FrMaxC,FrRmax,
                              Nsurfpt_fr_deso,surfpt_fr_deso_orig,
                              surfpt_fr_deso,nsurf_fr_deso,FrSurf_deso,
                              pointsrf_fr_deso,Tr,U1,U2,
                              WaMoRa,rslop,MaxNeigh,&ReDesoElec,&FrDesoElec);

                          Dr_f=SFDeso_re*ReDesoElec;
                          Df_f=SFDeso_fr*FrDesoElec;
#ifdef XFAST
                        }
                        else /* dbx check remove */
                        {
                          Dr_f=1000;
                          Df_f=1000;
                        }
#endif
                        In_f=SFIntElec*In_f;
                        To_f=VW_f+In_f+Dr_f+Df_f;

                        /* Check energy cutoff */
                        if (To_f<=FrMaEn) {

                          SFWrNu=SFWrNu+1;

                          FrCoPo[SFWrNu]=dmatrix(1,FrAtNu,1,3);
                          for (i1=1;i1<=FrAtNu;i1++) {
                            for (i2=1;i2<=3;i2++) {
                              FrCoPo[SFWrNu][i1][i2]=RoSFCo[i1][i2];
                            }
                          }

                          ConfArr[SFWrNu]=Ind_num_cn;


                          if(SFWrNu + 1>= currentsize)
                          {
                            /*  to do  */
                            /* #ifdef OMP */
                            /* #pragma omp critical */
                            /* #endif */
                            currentsize += _ALLOCSIZE;
                            VW_f_ro = dvecresize(VW_f_ro,currentsize+1);
                            VW_s_ro = dvecresize(VW_s_ro,currentsize+1);
                            In_f_ro = dvecresize(In_f_ro,currentsize+1);
                            In_s_ro = dvecresize(In_s_ro,currentsize+1);
                            Dr_f_ro = dvecresize(Dr_f_ro,currentsize+1);
                            Dr_s_ro = dvecresize(Dr_s_ro,currentsize+1);
                            Df_f_ro = dvecresize(Df_f_ro,currentsize+1);
                            Df_s_ro = dvecresize(Df_s_ro,currentsize+1);
                            To_f_ro = dvecresize(To_f_ro,currentsize+1);
                            To_s_ro = dvecresize(To_s_ro,currentsize+1);
                            Index_ro = ivecresize(Index_ro,currentsize+1);
                          }

                          Index_ro[SFWrNu] = SFWrNu;
                          VW_f_ro[SFWrNu] = VW_f;
                          VW_s_ro[SFWrNu] = VW_s;
                          In_f_ro[SFWrNu] = In_f;
                          In_s_ro[SFWrNu] = In_s;
                          Dr_f_ro[SFWrNu] = Dr_f;
                          Dr_s_ro[SFWrNu] = Dr_s;
                          Df_f_ro[SFWrNu] = Df_f;
                          Df_s_ro[SFWrNu] = Df_s;
                          To_f_ro[SFWrNu] = To_f;
                          To_s_ro[SFWrNu] = To_s;
                          //debug clangini start
                          //if (SFWrNu == 1552 || SFWrNu == 3436|| SFWrNu == 6295
                          //|| SFWrNu == 13581|| SFWrNu == 5982||SFWrNu == 14650){
                          //  std::cout << SFWrNu << ":  "<< l
                          //  << "  To_s:  " << To_s << std::endl;
                          //}
                          //debug clangini end

                        }

                      }

                    }

                  }

                }


              }


            }
          }
        }

        free_dmatrix(apol_Vect_fr,1,Nsurfpt_fr_apol,1,6);
        free_ivector(FrApAt,1,Nsurfpt_fr_apol);

      }
      else if ((EvalEn[0]!='e')&&
          (ApPoChoi[0]=='b')) { // SEEDING both polar and apolar START. clangini
        /* ------------------------------------- */
        /* Seeding in both polar and apolar ways */
        /* ------------------------------------- */
        fprintf(FPaOut,"This fragment is considered as a polar and apolar fragment\n\n");


        /* Compute the total number of fragment positions for the polar seeding */
        ToNFrP=0;
        for (j=1;j<=ReAcNu+ReDoNu;j++) {
          for (k=1;k<=FrAcNu+FrDoNu;k++) {
            if ((ReDATy[j]!=FrDATy[k])&&(UndisAt_fr[FrDAAt[k]])&&
                (PolVect_rec_01[j]))
              ToNFrP=ToNFrP+1;
          }
        }
        ToNFrP=ToNFrP*NuRoAx*FrCoNu;

        /* ---------------------------------------------------------------*/
        /* This part is needed to generate VECTOR FOR APOLAR FRAGMENTS !! */

        /* Make another grid on the surface of the fragment to describe its SAS
           (this one will be used for the vectors of apolar fragments) */

        /*      printf("\n\tGeneration of the SAS for the vectors of apolar
                fragments\n"); */
        NPtSphereMax_Fr = (int) (PtDensFr * pi4 * (FrRmax+WaMoRa) * 2);
        /* * 2 necessary otherwise seg.fault with e.g. tribrommethane */

        surfpt_fr_apol=structpointvect(1,NPtSphereMax_Fr*FrAtNu);
        iatsrf_fr_apol=ivector(1,NPtSphereMax_Fr*FrAtNu);
        pointsrf_fr_apol=ivector(1,FrAtNu);
        for (ij=1;ij<=FrAtNu;ij++)
          pointsrf_fr_apol[ij]=0;
        nsurf_fr_apol=ivector(1,FrAtNu);
        for (iat=1;iat<=FrAtNu;iat++)
          nsurf_fr_apol[iat] = 0;

        nn = Surf_Grid_unif_dens(FrAtNu,FrCoor,FrMaxC,FrMinC,FrRad,FrRadOut2,
            Sphere_apol,PtDensFr,pi4,&Nsurfpt_fr_apol,
            surfpt_fr_apol,iatsrf_fr_apol,
            nsurf_fr_apol,pointsrf_fr_apol);

        apol_Vect_fr = dmatrix(1,Nsurfpt_fr_apol,1,6);
        FrApAt=ivector(1,Nsurfpt_fr_apol);
        for (isph=1;isph<=Nsurfpt_fr_apol;isph++) {
          for (j=1;j<=3;j++)
            apol_Vect_fr[isph][j+3]=FrCoor[iatsrf_fr_apol[isph]][j];
          apol_Vect_fr[isph][1] = (double) surfpt_fr_apol[isph].x;
          apol_Vect_fr[isph][2] = (double) surfpt_fr_apol[isph].y;
          apol_Vect_fr[isph][3] = (double) surfpt_fr_apol[isph].z;
          FrApAt[isph]=iatsrf_fr_apol[isph];
        }

        free_ivector(pointsrf_fr_apol,1,FrAtNu);
        free_ivector(nsurf_fr_apol,1,FrAtNu);
        free_ivector(iatsrf_fr_apol,1,NPtSphereMax_Fr*FrAtNu);
        free_structpointvect(surfpt_fr_apol,1,NPtSphereMax_Fr*FrAtNu);

        /* Normalization of the vectors apol_Vect_fr */
        for (i1=1;i1<=Nsurfpt_fr_apol;i1++)
          PoCoVe(apol_Vect_fr[i1][4],apol_Vect_fr[i1][5],apol_Vect_fr[i1][6],
              apol_Vect_fr[i1][1],apol_Vect_fr[i1][2],apol_Vect_fr[i1][3],
              1.0,&(apol_Vect_fr[i1][1]),&(apol_Vect_fr[i1][2]),
              &(apol_Vect_fr[i1][3]));

        fprintf(FPaOut,"Numb. of frag. vectors (surf. points) used to seed ");
        fprintf(FPaOut,"this apolar fragment : %d\n",Nsurfpt_fr_apol);
        fprintf(FPaOut,"This number is then reduced if undistinguishable");
        fprintf(FPaOut," fragment atoms\n\n");

        /*
            FilePa=fopen("./outputs/apolar_frag.mol2","w");
            fprintf(FilePa,"# TRIPOS MOL2 file generated by WITNOTP\n");
            fprintf(FilePa,"\n");
            fprintf(FilePa,"@<TRIPOS>MOLECULE\n");
            fprintf(FilePa,"apolar_frag\n");
            fprintf(FilePa,"%d 0 0 0 0\n",2*Nsurfpt_fr_apol);
            fprintf(FilePa,"****\n");
            fprintf(FilePa,"USER_CHARGES\n");
            fprintf(FilePa,"INVALID_CHARGES\n");
            fprintf(FilePa,"\n");
            fprintf(FilePa,"@<TRIPOS>ATOM\n");

            for (i1=1;i1<=Nsurfpt_fr_apol;i1++) {
            fprintf(FilePa,"%-5d %-3d%11.4f%11.4f%11.4f ***** 1    ****  0.000\n",
            i1*2-1,i1*2-1,apol_Vect_fr[i1][1],apol_Vect_fr[i1][2],
            apol_Vect_fr[i1][3]);
            fprintf(FilePa,"%-5d %-3d%11.4f%11.4f%11.4f ***** 1    ****  0.000\n",
            i1*2,i1*2,apol_Vect_fr[i1][4],apol_Vect_fr[i1][5],
            apol_Vect_fr[i1][6]);
            }
            fclose(FilePa);
         */

        /* End of the part needed to generate VECTOR FOR APOLAR FRAGMENTS !! */
        /* ------------------------------------------------------------------*/

        VdWCoff_ap=MuFaVdWCoff_ap*((double) FrAtNu);

        /* Compute the total number of fragment positions for the apolar seeding */
        ToNFrP_ap=0;
        for (j=1;j<=Napol_Vect_re;j++) {
          for (k=1;k<=Nsurfpt_fr_apol;k++) {
            if (UndisAt_fr[FrApAt[k]])
              ToNFrP_ap=ToNFrP_ap+1;
          }
        }
        ToNFrP_ap=ToNFrP_ap*NuRoAx*FrCoNu;

        ToNFrP=ToNFrP+ToNFrP_ap;

        if (Ind_num_cn==1) {
          FrCoPo=ppdvector(1,ToNFrP);
          ConfArr=ivector(1,ToNFrP);
        }

        if (Ind_num_cn==1) {
          SFWrNu=0;
          FrNuTo=0;
          FrNuPB=0;
          FrNuSp=0;
          FrNuVdW_ap=0;
        }

        //std::cout << "Right before polar seeding" << std::endl; //clangini OK
        /* - - - - - - - */
        /* POLAR SEEDING */
        /* - - - - - - - */
        for (j=1;j<=ReAcNu+ReDoNu;j++) {
          for (k=1;k<=FrAcNu+FrDoNu;k++) {
            if ((ReDATy[j]!=FrDATy[k])&&     /* vector types must be different */
                (UndisAt_fr[FrDAAt[k]])&&    /* atom accepted (symmetry) */
                (PolVect_rec_01[j])) {       /* rec. vector accept. (probe vdW) */

              RecHyN=ReHydN[j];
              FraHyN=FrHydN[k];

              /* Seed fragment */
              SeedFr(j,ReVeCo,k,FrVeCo,FrAtNu,FrCoor,ReDATy,SeFrCo,FrAtoTyp_nu,
                  FrDAAt,ReAtoTyp_nu,ReDAAt,BLAtTy,FPaOut);

#ifdef USE_QUATERNION_ROTATION //clangini
              AnglRo = M_PI*2/NuRoAx;

              SeFrAx[1]=ReCoor[ReDAAt[j]][1]-SeFrCo[FrDAAt[k]][1];
              SeFrAx[2]=ReCoor[ReDAAt[j]][2]-SeFrCo[FrDAAt[k]][2];
              SeFrAx[3]=ReCoor[ReDAAt[j]][3]-SeFrCo[FrDAAt[k]][3];
              NormVe(&SeFrAx[1],&SeFrAx[2],&SeFrAx[3]);
              q_SeFr.fromAngleAxis(AnglRo,SeFrAx);

              for (int ii=1;ii<=FrAtNu;ii++) {
                //Need to add back the
                RoSFCo[ii][1] = SeFrCo[ii][1];
                RoSFCo[ii][2] = SeFrCo[ii][2];
                RoSFCo[ii][3] = SeFrCo[ii][3];
              }
#endif

              /* Make rotations */
              for (l=1;l<=NuRoAx;l++) {
#ifdef USE_QUATERNION_ROTATION //clangini
                FrNuTo = FrNuTo+1;
                for (int ii=1;ii<=FrAtNu;ii++){
                  q_SeFr.quatConjugateVecRef(RoSFCo[ii],SeFrCo[FrDAAt[k]]);
                }
#else
                AnglRo=M_PI*2/NuRoAx*l;

                FrNuTo=FrNuTo+1;

                /* RoSFCo is the position of the rotated fragment */
                RoSeFr(j,ReDAAt,ReCoor,k,FrDAAt,FrAtNu,SeFrCo,AnglRo,RoSFCo);
#endif

                VW_f=0.0;
                VW_s=0.0;
                In_f=0.0;
                In_s=0.0;
                Dr_f=0.0;
                Dr_s=0.0;
                Df_f=0.0;
                Df_s=0.0;
                To_f=0.0;
                To_s=0.0;


                /* ----- SLOW or BOTH methods ----- */
                if ((Solv_typ[0]=='s')||(Solv_typ[0]=='b')) {

                  SpPoCh_bool=1;

                  /* Only select a fragment if it is in the specified sphere. In order to do that,
                     compute the geometric center of the fragment and its distance to the sphere
                     center */
                  if (SpPoCh_opt[0]=='y') {

                    SpPoCh_gc[1]=0.0;
                    SpPoCh_gc[2]=0.0;
                    SpPoCh_gc[3]=0.0;
                    for (i1=1;i1<=FrAtNu;i1++) {
                      SpPoCh_gc[1]=SpPoCh_gc[1]+RoSFCo[i1][1];
                      SpPoCh_gc[2]=SpPoCh_gc[2]+RoSFCo[i1][2];
                      SpPoCh_gc[3]=SpPoCh_gc[3]+RoSFCo[i1][3];
                    }
                    SpPoCh_gc[1]=SpPoCh_gc[1]/FrAtNu;
                    SpPoCh_gc[2]=SpPoCh_gc[2]/FrAtNu;
                    SpPoCh_gc[3]=SpPoCh_gc[3]/FrAtNu;

                    SpPoCh_dist=DistSq(SpPoCh_gc[1],SpPoCh_gc[2],SpPoCh_gc[3],
                        SpPoCh_cent[1],SpPoCh_cent[2],SpPoCh_cent[3]);

                    if (SpPoCh_dist<=SpPoCh_rad_sq)
                      FrNuSp=FrNuSp+1;
                    else
                      SpPoCh_bool=0;

                  }

                  if (SpPoCh_bool) {

                    /* Compute the # bumps */
                    FrAcce=CoNuBu(CubNum,CubFAI,CubLAI,CubLiA,LaVdWR,ReMaxC,ReMinC,
                        ReCoor,ReVdWR,FrAtNu,RoSFCo,FrVdWR,VdWFaB,BumpMa,
                        RecHyN,FraHyN,BuEvFa);

                    if (FrAcce) {

                      FrNuPB=FrNuPB+1;

                      Rot_Tran(FrAtNu,FrCoor,RoSFCo,Tr,U1,U2);

                      /* Compute the squared distances between the fragment atoms and the
                         corresponding receptor atoms of the pseudo-sphere procedure, one list
                         for the vdW energy and another for the electrostatic interaction energy.
                         In the case of the electrostatic interaction energy, it also provides the
                         corresponding partial or total charges */
                      SqDisFrRe_ps(FrAtNu,RoSFCo,ReCoor,ReMinC,GrSiCu_en,
                          CubNum_en,CubFAI_en,CubLAI_en,CubLiA_en,
                          PsSpNC,PsSphe,SDFrRe_ps,ReAtNu,PsSpRa,
                          RePaCh,ReReNu,AtReprRes,FiAtRes,LaAtRes,
                          TotChaRes,NuChResEn,LiChResEn,
                          SDFrRe_ps_elec,ChFrRe_ps_elec);

                      /* Compute vdW energy */
                      PsSpEE(FrAtNu,ReAtNu,ReVdWE_sr,FrVdWE_sr,
                          ReVdWR,FrVdWR,&VWEnEv_ps,SDFrRe_ps);

                      /* Compute receptor desolvation, fragment desolvation and receptor-fragment
                         interaction (with screening effect) energies (slow method) */
                      ElecFrag(ReAtNu,ReCoor,RePaCh,ChFrRe_ps_elec,
                          ReRad,ReRad2,ReRadOut,
                          ReRadOut2,surfpt_re,nsurf_re,
                          pointsrf_re,ReSelfVol,FrAtNu,RoSFCo,FrCoor,
                          FrPaCh,FrRad,FrRad2,FrRadOut,FrRadOut2,
                          Frdist2,SDFrRe_ps_elec,FrMinC,FrMaxC,&FrSolvEn,
                          Nsurfpt_fr,surfpt_fr,
                          nsurf_fr,pointsrf_fr,surfpt_ex,Tr,U1,U2,WaMoRa,
                          GrSiSo,NPtSphere,Min,Max,XGrid,YGrid,ZGrid,
                          NGridx,NGridy,NGridz,GridMat,
                          DeltaPrDeso,Kelec,Ksolv,UnitVol,
                          pi4,nxminBS,nyminBS,nzminBS,nxmaxBS,nymaxBS,
                          nzmaxBS,corr_scrint,corr_fr_deso,&ReDesoElec,
                          &ReFrIntElec,&FrDesoElec,ReSelfVol_corrB,EmpCorrB,FPaOut);

                      VW_s=SFVWEn*VWEnEv_ps;
                      In_s=SFIntElec*ReFrIntElec;
                      Dr_s=SFDeso_re*ReDesoElec;
                      Df_s=SFDeso_fr*FrDesoElec;
                      To_s=VW_s+In_s+Dr_s+Df_s;



                      /* Check energy cutoff */
                      if (To_s<=FrMaEn) {

                        SFWrNu=SFWrNu+1;

                        FrCoPo[SFWrNu]=dmatrix(1,FrAtNu,1,3);
                        for (i1=1;i1<=FrAtNu;i1++) {
                          for (i2=1;i2<=3;i2++) {
                            FrCoPo[SFWrNu][i1][i2]=RoSFCo[i1][i2];
                          }
                        }

                        ConfArr[SFWrNu]=Ind_num_cn;

                        if (Solv_typ[0]=='b') {

                          /* Van der Waals energy on a grid using the geometric mean approximation */
                          VW_f=VWFrEn(FrAtNu,RoSFCo,BSMinC,VWGrIn,VWGrSi,FrVdWR_p3,
                              FrVdWE_sr,VWGPoN,VWGrRP_at,VWGrRP_re,FPaOut,&vdwErr);

                          /* Coulombic interaction energy on a grid */
                          In_f=CoFrEn(FrAtNu,FrPaCh,RoSFCo,BSMinC,CoGrIn,CoGrSi,
                              CoGPoN,CoGrRP,FPaOut,&clbErr);

                          /* Compute receptor desolvation, fragment desolvation energies (fast method) */
                          ElecFragFast(ReAtNu,ReCoor,ReRadOut,ReRadOut2,
                              surfpt_re_deso,nsurf_re_deso,ReSurf_deso,
                              pointsrf_re_deso,nlist_re,list_re,
                              nstep_re,cbmid_re,midg_re,rscale_re,
                              FrAtNu,RoSFCo,FrCoor,FrRad,FrRadOut,
                              FrRadOut2,FrMinC,FrMaxC,FrRmax,
                              Nsurfpt_fr_deso,surfpt_fr_deso_orig,
                              surfpt_fr_deso,nsurf_fr_deso,
                              FrSurf_deso,pointsrf_fr_deso,
                              Tr,U1,U2,WaMoRa,rslop,
                              MaxNeigh,&ReDesoElec,&FrDesoElec);

                          VW_f=SFVWEn*VW_f;
                          In_f=SFIntElec*In_f;
                          Dr_f=SFDeso_re*ReDesoElec;
                          Df_f=SFDeso_fr*FrDesoElec;
                          To_f=VW_f+In_f+Dr_f+Df_f;

                        }


                        if(SFWrNu + 1>= currentsize)
                        {
                          /*  to do  */
                          /* #ifdef OMP */
                          /* #pragma omp critical */
                          /* #endif */
                          currentsize += _ALLOCSIZE;
                          VW_f_ro = dvecresize(VW_f_ro,currentsize+1);
                          VW_s_ro = dvecresize(VW_s_ro,currentsize+1);
                          In_f_ro = dvecresize(In_f_ro,currentsize+1);
                          In_s_ro = dvecresize(In_s_ro,currentsize+1);
                          Dr_f_ro = dvecresize(Dr_f_ro,currentsize+1);
                          Dr_s_ro = dvecresize(Dr_s_ro,currentsize+1);
                          Df_f_ro = dvecresize(Df_f_ro,currentsize+1);
                          Df_s_ro = dvecresize(Df_s_ro,currentsize+1);
                          To_f_ro = dvecresize(To_f_ro,currentsize+1);
                          To_s_ro = dvecresize(To_s_ro,currentsize+1);
                          Index_ro = ivecresize(Index_ro,currentsize+1);
                        }

                        Index_ro[SFWrNu] = SFWrNu;
                        VW_f_ro[SFWrNu] = VW_f;
                        VW_s_ro[SFWrNu] = VW_s;
                        In_f_ro[SFWrNu] = In_f;
                        In_s_ro[SFWrNu] = In_s;
                        Dr_f_ro[SFWrNu] = Dr_f;
                        Dr_s_ro[SFWrNu] = Dr_s;
                        Df_f_ro[SFWrNu] = Df_f;
                        Df_s_ro[SFWrNu] = Df_s;
                        To_f_ro[SFWrNu] = To_f;
                        To_s_ro[SFWrNu] = To_s;

                      }

                    }

                  }

                }


                /* ----- FAST method ----- */
                if ((Solv_typ[0]=='f')||(Solv_typ[0]=='p')) {

                  SpPoCh_bool=1;

                  /* Only select a fragment if it is in the specified sphere. In order to do that,
                     compute the geometric center of the fragment and its distance to the sphere
                     center */
                  if (SpPoCh_opt[0]=='y') {

                    SpPoCh_gc[1]=0.0;
                    SpPoCh_gc[2]=0.0;
                    SpPoCh_gc[3]=0.0;
                    for (i1=1;i1<=FrAtNu;i1++) {
                      SpPoCh_gc[1]=SpPoCh_gc[1]+RoSFCo[i1][1];
                      SpPoCh_gc[2]=SpPoCh_gc[2]+RoSFCo[i1][2];
                      SpPoCh_gc[3]=SpPoCh_gc[3]+RoSFCo[i1][3];
                    }
                    SpPoCh_gc[1]=SpPoCh_gc[1]/FrAtNu;
                    SpPoCh_gc[2]=SpPoCh_gc[2]/FrAtNu;
                    SpPoCh_gc[3]=SpPoCh_gc[3]/FrAtNu;

                    SpPoCh_dist=DistSq(SpPoCh_gc[1],SpPoCh_gc[2],SpPoCh_gc[3],
                        SpPoCh_cent[1],SpPoCh_cent[2],SpPoCh_cent[3]);

                    if (SpPoCh_dist<=SpPoCh_rad_sq)
                      FrNuSp=FrNuSp+1;
                    else
                      SpPoCh_bool=0;

                  }

                  if (SpPoCh_bool) {

                    /* Van der Waals energy on a grid using the geometric mean approximation */
                    VW_f=VWFrEn(FrAtNu,RoSFCo,BSMinC,VWGrIn,VWGrSi,FrVdWR_p3,
                        FrVdWE_sr,VWGPoN,VWGrRP_at,VWGrRP_re,FPaOut,&vdwErr);

                    VW_f=SFVWEn*VW_f;

                    /* The van der Waals energy acts as a bump checking */
                    if (VW_f<=BumpFaCut) {

                      FrNuPB=FrNuPB+1;

                      Rot_Tran(FrAtNu,FrCoor,RoSFCo,Tr,U1,U2);

                      /* Coulombic interaction energy on a grid */
                      In_f=CoFrEn(FrAtNu,FrPaCh,RoSFCo,BSMinC,CoGrIn,CoGrSi,
                          CoGPoN,CoGrRP,FPaOut,&clbErr);

                      /* Compute receptor desolvation, fragment desolvation energies (fast method) */
                      ElecFragFast(ReAtNu,ReCoor,ReRadOut,ReRadOut2,
                          surfpt_re_deso,nsurf_re_deso,ReSurf_deso,
                          pointsrf_re_deso,nlist_re,list_re,
                          nstep_re,cbmid_re,midg_re,rscale_re,
                          FrAtNu,RoSFCo,FrCoor,FrRad,FrRadOut,
                          FrRadOut2,FrMinC,FrMaxC,FrRmax,
                          Nsurfpt_fr_deso,surfpt_fr_deso_orig,
                          surfpt_fr_deso,nsurf_fr_deso,FrSurf_deso,
                          pointsrf_fr_deso,
                          Tr,U1,U2,WaMoRa,rslop,
                          MaxNeigh,&ReDesoElec,&FrDesoElec);

                      In_f=SFIntElec*In_f;
                      Dr_f=SFDeso_re*ReDesoElec;
                      Df_f=SFDeso_fr*FrDesoElec;
                      To_f=VW_f+In_f+Dr_f+Df_f;

                      /* Check energy cutoff */
                      if (To_f<=FrMaEn) {

                        SFWrNu=SFWrNu+1;

                        FrCoPo[SFWrNu]=dmatrix(1,FrAtNu,1,3);
                        for (i1=1;i1<=FrAtNu;i1++) {
                          for (i2=1;i2<=3;i2++) {
                            FrCoPo[SFWrNu][i1][i2]=RoSFCo[i1][i2];
                          }
                        }

                        ConfArr[SFWrNu]=Ind_num_cn;


                        if(SFWrNu + 1>= currentsize)
                        {
                          /*  to do  */
                          /* #ifdef OMP */
                          /* #pragma omp critical */
                          /* #endif */
                          currentsize += _ALLOCSIZE;
                          VW_f_ro = dvecresize(VW_f_ro,currentsize+1);
                          VW_s_ro = dvecresize(VW_s_ro,currentsize+1);
                          In_f_ro = dvecresize(In_f_ro,currentsize+1);
                          In_s_ro = dvecresize(In_s_ro,currentsize+1);
                          Dr_f_ro = dvecresize(Dr_f_ro,currentsize+1);
                          Dr_s_ro = dvecresize(Dr_s_ro,currentsize+1);
                          Df_f_ro = dvecresize(Df_f_ro,currentsize+1);
                          Df_s_ro = dvecresize(Df_s_ro,currentsize+1);
                          To_f_ro = dvecresize(To_f_ro,currentsize+1);
                          To_s_ro = dvecresize(To_s_ro,currentsize+1);
                          Index_ro = ivecresize(Index_ro,currentsize+1);
                        }

                        Index_ro[SFWrNu] = SFWrNu;
                        VW_f_ro[SFWrNu] = VW_f;
                        VW_s_ro[SFWrNu] = VW_s;
                        In_f_ro[SFWrNu] = In_f;
                        In_s_ro[SFWrNu] = In_s;
                        Dr_f_ro[SFWrNu] = Dr_f;
                        Dr_s_ro[SFWrNu] = Dr_s;
                        Df_f_ro[SFWrNu] = Df_f;
                        Df_s_ro[SFWrNu] = Df_s;
                        To_f_ro[SFWrNu] = To_f;
                        To_s_ro[SFWrNu] = To_s;

                      }

                    }

                  }

                }


              }


            }
          }
        }

        //std::cout << "Right before apolar seeding" << std::endl; //clangini OK
        /* - - - - - - - - */
        /*  APOLAR SEEDING */
        /* - - - - - - - - */
        for (j=1;j<=Napol_Vect_re;j++) {
          for (k=1;k<=Nsurfpt_fr_apol;k++) {
            if (UndisAt_fr[FrApAt[k]]) {        /* accepted atom (symmetry) */

              /* Seed fragment */
              SeedFr_ap(j,apol_Vect_re,ReApAt,k,apol_Vect_fr,FrApAt,FrCoor,
                  ReVdWR,FrVdWR,FrAtNu,FPaOut,SeFrCo);
#ifdef USE_QUATERNION_ROTATION //clangini
              AnglRo = M_PI*2/NuRoAx;

              SeFrAx[1]=ReCoor[ReApAt[j]][1]-SeFrCo[FrApAt[k]][1];
              SeFrAx[2]=ReCoor[ReApAt[j]][2]-SeFrCo[FrApAt[k]][2];
              SeFrAx[3]=ReCoor[ReApAt[j]][3]-SeFrCo[FrApAt[k]][3];
              NormVe(&SeFrAx[1],&SeFrAx[2],&SeFrAx[3]);
              q_SeFr.fromAngleAxis(AnglRo,SeFrAx);

              for (int ii=1;ii<=FrAtNu;ii++) {
                //Need to add back the
                RoSFCo[ii][1] = SeFrCo[ii][1];
                RoSFCo[ii][2] = SeFrCo[ii][2];
                RoSFCo[ii][3] = SeFrCo[ii][3];
              }
#endif

              /* Make rotations */
              for (l=1;l<=NuRoAx;l++) {
#ifdef USE_QUATERNION_ROTATION //clangini
                FrNuTo = FrNuTo+1;
                for (int ii=1;ii<=FrAtNu;ii++){
                  q_SeFr.quatConjugateVecRef(RoSFCo[ii],SeFrCo[FrApAt[k]]);
                }
#else
                AnglRo=M_PI*2/NuRoAx*l;

                FrNuTo=FrNuTo+1;

                /* RoSFCo is the position of the rotated fragment */
                RoSeFr(j,ReApAt,ReCoor,k,FrApAt,FrAtNu,SeFrCo,AnglRo,RoSFCo);
#endif

                VW_f=0.0;
                VW_s=0.0;
                In_f=0.0;
                In_s=0.0;
                Dr_f=0.0;
                Dr_s=0.0;
                Df_f=0.0;
                Df_s=0.0;
                To_f=0.0;
                To_s=0.0;


                /* ----- SLOW or BOTH methods ----- */
                if ((Solv_typ[0]=='s')||(Solv_typ[0]=='b')) {

                  SpPoCh_bool=1;

                  /* Only select a fragment if it is in the specified sphere. In order to do that,
                     compute the geometric center of the fragment and its distance to the sphere
                     center */
                  if (SpPoCh_opt[0]=='y') {

                    SpPoCh_gc[1]=0.0;
                    SpPoCh_gc[2]=0.0;
                    SpPoCh_gc[3]=0.0;
                    for (i1=1;i1<=FrAtNu;i1++) {
                      SpPoCh_gc[1]=SpPoCh_gc[1]+RoSFCo[i1][1];
                      SpPoCh_gc[2]=SpPoCh_gc[2]+RoSFCo[i1][2];
                      SpPoCh_gc[3]=SpPoCh_gc[3]+RoSFCo[i1][3];
                    }
                    SpPoCh_gc[1]=SpPoCh_gc[1]/FrAtNu;
                    SpPoCh_gc[2]=SpPoCh_gc[2]/FrAtNu;
                    SpPoCh_gc[3]=SpPoCh_gc[3]/FrAtNu;

                    SpPoCh_dist=DistSq(SpPoCh_gc[1],SpPoCh_gc[2],SpPoCh_gc[3],
                        SpPoCh_cent[1],SpPoCh_cent[2],SpPoCh_cent[3]);

                    if (SpPoCh_dist<=SpPoCh_rad_sq)
                      FrNuSp=FrNuSp+1;
                    else
                      SpPoCh_bool=0;

                  }

                  if (SpPoCh_bool) {

                    /* Compute the # bumps */
                    FrAcce=CoNuBu(CubNum,CubFAI,CubLAI,CubLiA,LaVdWR,ReMaxC,ReMinC,
                        ReCoor,ReVdWR,FrAtNu,RoSFCo,FrVdWR,VdWFaB,BumpMa,
                        0,0,BuEvFa);

                    if (FrAcce) {

                      FrNuPB=FrNuPB+1;

                      /* Compute the squared distances between the fragment atoms and the
                         corresponding receptor atoms of the pseudo-sphere procedure, one list
                         for the vdW energy and another for the electrostatic interaction energy.
                         In the case of the electrostatic interaction energy, it also provides the
                         corresponding partial or total charges */
                      SqDisFrRe_ps(FrAtNu,RoSFCo,ReCoor,ReMinC,GrSiCu_en,
                          CubNum_en,CubFAI_en,CubLAI_en,CubLiA_en,
                          PsSpNC,PsSphe,SDFrRe_ps,ReAtNu,PsSpRa,
                          RePaCh,ReReNu,AtReprRes,FiAtRes,LaAtRes,
                          TotChaRes,NuChResEn,LiChResEn,
                          SDFrRe_ps_elec,ChFrRe_ps_elec);

                      /* Compute vdW energy */
                      PsSpEE(FrAtNu,ReAtNu,ReVdWE_sr,FrVdWE_sr,
                          ReVdWR,FrVdWR,&VWEnEv_ps,SDFrRe_ps);

                      VW_s=SFVWEn*VWEnEv_ps;

                      /* The vdW energy has to be smaller than a cutoff */
                      if (VW_s<=VdWCoff_ap) {

                        FrNuVdW_ap=FrNuVdW_ap+1;

                        Rot_Tran(FrAtNu,FrCoor,RoSFCo,Tr,U1,U2);

                        /* Compute receptor desolvation, fragment desolvation and receptor-fragment
                           interaction (with screening effect) energies (slow method) */
                        ElecFrag(ReAtNu,ReCoor,RePaCh,ChFrRe_ps_elec,
                            ReRad,ReRad2,ReRadOut,
                            ReRadOut2,surfpt_re,nsurf_re,
                            pointsrf_re,ReSelfVol,FrAtNu,RoSFCo,FrCoor,
                            FrPaCh,FrRad,FrRad2,FrRadOut,FrRadOut2,
                            Frdist2,SDFrRe_ps_elec,FrMinC,FrMaxC,
                            &FrSolvEn,Nsurfpt_fr,surfpt_fr,
                            nsurf_fr,pointsrf_fr,surfpt_ex,Tr,U1,U2,WaMoRa,
                            GrSiSo,NPtSphere,Min,Max,XGrid,YGrid,ZGrid,
                            NGridx,NGridy,NGridz,GridMat,
                            DeltaPrDeso,Kelec,Ksolv,UnitVol,
                            pi4,nxminBS,nyminBS,nzminBS,nxmaxBS,nymaxBS,
                            nzmaxBS,corr_scrint,corr_fr_deso,&ReDesoElec,
                            &ReFrIntElec,&FrDesoElec,ReSelfVol_corrB,EmpCorrB,FPaOut);

                        In_s=SFIntElec*ReFrIntElec;
                        Dr_s=SFDeso_re*ReDesoElec;
                        Df_s=SFDeso_fr*FrDesoElec;
                        To_s=VW_s+In_s+Dr_s+Df_s;


                        /* Check energy cutoff */
                        if (To_s<=FrMaEn) {

                          SFWrNu=SFWrNu+1;

                          FrCoPo[SFWrNu]=dmatrix(1,FrAtNu,1,3);
                          for (i1=1;i1<=FrAtNu;i1++) {
                            for (i2=1;i2<=3;i2++) {
                              FrCoPo[SFWrNu][i1][i2]=RoSFCo[i1][i2];
                            }
                          }

                          ConfArr[SFWrNu]=Ind_num_cn;

                          if (Solv_typ[0]=='b') {

                            /* Van der Waals energy on a grid using the geometric mean approximation */
                            VW_f=VWFrEn(FrAtNu,RoSFCo,BSMinC,VWGrIn,VWGrSi,FrVdWR_p3,
                                FrVdWE_sr,VWGPoN,VWGrRP_at,VWGrRP_re,FPaOut,&vdwErr);

                            /* Coulombic interaction energy on a grid */
                            In_f=CoFrEn(FrAtNu,FrPaCh,RoSFCo,BSMinC,CoGrIn,CoGrSi,
                                CoGPoN,CoGrRP,FPaOut,&clbErr);

                            /* Compute receptor desolvation, fragment desolvation energies (fast method) */
                            ElecFragFast(ReAtNu,ReCoor,ReRadOut,
                                ReRadOut2,surfpt_re_deso,nsurf_re_deso,
                                ReSurf_deso,pointsrf_re_deso,nlist_re,
                                list_re,nstep_re,cbmid_re,midg_re,
                                rscale_re,FrAtNu,RoSFCo,FrCoor,FrRad,
                                FrRadOut,FrRadOut2,FrMinC,FrMaxC,FrRmax,
                                Nsurfpt_fr_deso,surfpt_fr_deso_orig,
                                surfpt_fr_deso,nsurf_fr_deso,
                                FrSurf_deso,pointsrf_fr_deso,Tr,U1,U2,
                                WaMoRa,rslop,MaxNeigh,
                                &ReDesoElec,&FrDesoElec);

                            VW_f=SFVWEn*VW_f;
                            In_f=SFIntElec*In_f;
                            Dr_f=SFDeso_re*ReDesoElec;
                            Df_f=SFDeso_fr*FrDesoElec;
                            To_f=VW_f+In_f+Dr_f+Df_f;

                          }


                          if(SFWrNu + 1>= currentsize)
                          {
                            /*  to do  */
                            /* #ifdef OMP */
                            /* #pragma omp critical */
                            /* #endif */
                            currentsize += _ALLOCSIZE;
                            VW_f_ro = dvecresize(VW_f_ro,currentsize+1);
                            VW_s_ro = dvecresize(VW_s_ro,currentsize+1);
                            In_f_ro = dvecresize(In_f_ro,currentsize+1);
                            In_s_ro = dvecresize(In_s_ro,currentsize+1);
                            Dr_f_ro = dvecresize(Dr_f_ro,currentsize+1);
                            Dr_s_ro = dvecresize(Dr_s_ro,currentsize+1);
                            Df_f_ro = dvecresize(Df_f_ro,currentsize+1);
                            Df_s_ro = dvecresize(Df_s_ro,currentsize+1);
                            To_f_ro = dvecresize(To_f_ro,currentsize+1);
                            To_s_ro = dvecresize(To_s_ro,currentsize+1);
                            Index_ro = ivecresize(Index_ro,currentsize+1);
                          }

                          Index_ro[SFWrNu] = SFWrNu;
                          VW_f_ro[SFWrNu] = VW_f;
                          VW_s_ro[SFWrNu] = VW_s;
                          In_f_ro[SFWrNu] = In_f;
                          In_s_ro[SFWrNu] = In_s;
                          Dr_f_ro[SFWrNu] = Dr_f;
                          Dr_s_ro[SFWrNu] = Dr_s;
                          Df_f_ro[SFWrNu] = Df_f;
                          Df_s_ro[SFWrNu] = Df_s;
                          To_f_ro[SFWrNu] = To_f;
                          To_s_ro[SFWrNu] = To_s;

                        }

                      }

                    }

                  }

                }


                /* ----- FAST method ----- */
                if ((Solv_typ[0]=='f')||(Solv_typ[0]=='p')) {

                  SpPoCh_bool=1;

                  /* Only select a fragment if it is in the specified sphere. In order to do that,
                     compute the geometric center of the fragment and its distance to the sphere
                     center */
                  if (SpPoCh_opt[0]=='y') {

                    SpPoCh_gc[1]=0.0;
                    SpPoCh_gc[2]=0.0;
                    SpPoCh_gc[3]=0.0;
                    for (i1=1;i1<=FrAtNu;i1++) {
                      SpPoCh_gc[1]=SpPoCh_gc[1]+RoSFCo[i1][1];
                      SpPoCh_gc[2]=SpPoCh_gc[2]+RoSFCo[i1][2];
                      SpPoCh_gc[3]=SpPoCh_gc[3]+RoSFCo[i1][3];
                    }
                    SpPoCh_gc[1]=SpPoCh_gc[1]/FrAtNu;
                    SpPoCh_gc[2]=SpPoCh_gc[2]/FrAtNu;
                    SpPoCh_gc[3]=SpPoCh_gc[3]/FrAtNu;

                    SpPoCh_dist=DistSq(SpPoCh_gc[1],SpPoCh_gc[2],SpPoCh_gc[3],
                        SpPoCh_cent[1],SpPoCh_cent[2],SpPoCh_cent[3]);

                    if (SpPoCh_dist<=SpPoCh_rad_sq)
                      FrNuSp=FrNuSp+1;
                    else
                      SpPoCh_bool=0;

                  }

                  if (SpPoCh_bool) {

                    /* Van der Waals energy on a grid using the geometric mean approximation */
                    VW_f=VWFrEn(FrAtNu,RoSFCo,BSMinC,VWGrIn,VWGrSi,FrVdWR_p3,
                        FrVdWE_sr,VWGPoN,VWGrRP_at,VWGrRP_re,FPaOut,&vdwErr);

                    VW_f=SFVWEn*VW_f;

                    /* The van der Waals energy acts as a bump checking */
                    if (VW_f<=BumpFaCut) {

                      FrNuPB=FrNuPB+1;

                      /* The vdW energy has to be smaller than a cutoff */
                      if (VW_f<=VdWCoff_ap) {

                        FrNuVdW_ap=FrNuVdW_ap+1;

                        Rot_Tran(FrAtNu,FrCoor,RoSFCo,Tr,U1,U2);

                        /* Coulombic interaction energy on a grid */
                        In_f=CoFrEn(FrAtNu,FrPaCh,RoSFCo,BSMinC,CoGrIn,CoGrSi,
                            CoGPoN,CoGrRP,FPaOut,&clbErr);

                        /* Compute receptor desolvation, fragment desolvation energies (fast method) */
#ifdef XFAST
                        if(SFIntElec*In_f +  VW_f <= FrMaEn) /*dbx check remove*/
                        {
#endif
                          ElecFragFast(ReAtNu,ReCoor,ReRadOut,
                              ReRadOut2,surfpt_re_deso,nsurf_re_deso,
                              ReSurf_deso,pointsrf_re_deso,nlist_re,
                              list_re,nstep_re,cbmid_re,midg_re,
                              rscale_re,FrAtNu,RoSFCo,FrCoor,FrRad,
                              FrRadOut,FrRadOut2,FrMinC,FrMaxC,FrRmax,
                              Nsurfpt_fr_deso,surfpt_fr_deso_orig,
                              surfpt_fr_deso,nsurf_fr_deso,FrSurf_deso,
                              pointsrf_fr_deso,Tr,U1,U2,
                              WaMoRa,rslop,MaxNeigh,&ReDesoElec,
                              &FrDesoElec);

                          Dr_f=SFDeso_re*ReDesoElec;
                          Df_f=SFDeso_fr*FrDesoElec;
#ifdef XFAST
                        }
                        else
                        {
                          Dr_f=1000;
                          Df_f=1000;
                        }
#endif
                        In_f=SFIntElec*In_f;
                        To_f=VW_f+In_f+Dr_f+Df_f;

                        /* Check energy cutoff */
                        if (To_f<=FrMaEn) {

                          SFWrNu=SFWrNu+1;

                          FrCoPo[SFWrNu]=dmatrix(1,FrAtNu,1,3);
                          for (i1=1;i1<=FrAtNu;i1++) {
                            for (i2=1;i2<=3;i2++) {
                              FrCoPo[SFWrNu][i1][i2]=RoSFCo[i1][i2];
                            }
                          }

                          ConfArr[SFWrNu]=Ind_num_cn;

                          if(SFWrNu + 1>= currentsize)
                          {
                            /*  to do  */
                            /* #ifdef OMP */
                            /* #pragma omp critical */
                            /* #endif */
                            currentsize += _ALLOCSIZE;
                            VW_f_ro = dvecresize(VW_f_ro,currentsize+1);
                            VW_s_ro = dvecresize(VW_s_ro,currentsize+1);
                            In_f_ro = dvecresize(In_f_ro,currentsize+1);
                            In_s_ro = dvecresize(In_s_ro,currentsize+1);
                            Dr_f_ro = dvecresize(Dr_f_ro,currentsize+1);
                            Dr_s_ro = dvecresize(Dr_s_ro,currentsize+1);
                            Df_f_ro = dvecresize(Df_f_ro,currentsize+1);
                            Df_s_ro = dvecresize(Df_s_ro,currentsize+1);
                            To_f_ro = dvecresize(To_f_ro,currentsize+1);
                            To_s_ro = dvecresize(To_s_ro,currentsize+1);
                            Index_ro = ivecresize(Index_ro,currentsize+1);
                          }

                          Index_ro[SFWrNu] = SFWrNu;
                          VW_f_ro[SFWrNu] = VW_f;
                          VW_s_ro[SFWrNu] = VW_s;
                          In_f_ro[SFWrNu] = In_f;
                          In_s_ro[SFWrNu] = In_s;
                          Dr_f_ro[SFWrNu] = Dr_f;
                          Dr_s_ro[SFWrNu] = Dr_s;
                          Df_f_ro[SFWrNu] = Df_f;
                          Df_s_ro[SFWrNu] = Df_s;
                          To_f_ro[SFWrNu] = To_f;
                          To_s_ro[SFWrNu] = To_s;

                        }

                      }

                    }

                  }

                }


              }


            }
          }
        }

        free_dmatrix(apol_Vect_fr,1,Nsurfpt_fr_apol,1,6);
        free_ivector(FrApAt,1,Nsurfpt_fr_apol);

        //std::cout << "End of seeding" << std::endl; //clangini OK
      }
      else if (EvalEn[0]=='e') {
        /* ------------------------------------------------------------------------ */
        /* Evaluating the energy for the fragment positions given in the input file */
        /* ------------------------------------------------------------------------ */

        /* Artificial memory allocation */
        ToNFrP=0;
        if (Ind_num_cn==1) {
          FrCoPo=ppdvector(1,ToNFrP);
          ConfArr=ivector(1,ToNFrP);
        }

        SFWrNu=0;

        /*for (j=1;j<=FrAtNu;j++) {
          RoSFCo[j][1]=FrCoor[j][1];
          RoSFCo[j][2]=FrCoor[j][2];
          RoSFCo[j][3]=FrCoor[j][3];
        }*/

        for (j=1;j<=FrAtNu;j++) {
          RoSFCo[j][1]=FrCoor_NoAlign[j][1];
          RoSFCo[j][2]=FrCoor_NoAlign[j][2];
          RoSFCo[j][3]=FrCoor_NoAlign[j][3];
        } //clangini. The seeded frag coordinates correspond to the original
        //coordinates without pre-alignment.

        /* Make the translation vector and the rotation matrices for ElecFrag */
        // Now FrCoor is no longer identical to RoSFCo
        /*for (j=1;j<=3;j++)
          Tr[j]=0.0;
        for (j=1;j<=3;j++) {
          for (k=1;k<=3;k++) {
            if (j==k) {
              U1[j][k]=1.0;
              U2[j][k]=1.0;
            }
            else {
              U1[j][k]=0.0;
              U2[j][k]=0.0;
            }
          }
        }*/

        Rot_Tran(FrAtNu,FrCoor,RoSFCo,Tr,U1,U2);//clangini

        VW_f=0.0;
        VW_s=0.0;
        In_f=0.0;
        In_s=0.0;
        Dr_f=0.0;
        Dr_s=0.0;
        Df_f=0.0;
        Df_s=0.0;
        To_f=0.0;
        To_s=0.0;

        /* Check whether the whole ligand lies in the grid(s) */
        /* Modification 09.08.2004 NM */
        MaxGrIn=GrInSo;
        if (Solv_typ[0]!='s')
        {
          if (VWGrIn>MaxGrIn) {MaxGrIn=VWGrIn;}
          if (CoGrIn>MaxGrIn) {MaxGrIn=CoGrIn;}
        }
        ChkInGrid=1;
        for (j=1;j<=FrAtNu;j++)
        {
          for (k=1;k<=3;k++)
          {
            if (RoSFCo[j][k]>(BSMaxC[k]+MaxGrIn)) {ChkInGrid=0;}
            if (RoSFCo[j][k]<(BSMinC[k]-MaxGrIn)) {ChkInGrid=0;}
          }
        }
        if (!ChkInGrid)
          fprintf(FPaOut,"WARNING Fragment %s is outside the grid(s)\n\n",FragNa); /*clangini*/
          /*fprintf(FPaOut,"WARNING Fragment %s is outside the grid(s)\n\n",
              &FrFiNa_out[CurFra][1]); clangini*/

        /* ----- SLOW or BOTH methods ----- */
        if (((Solv_typ[0]=='s')||(Solv_typ[0]=='b'))&&(ChkInGrid)) {

          /* Compute the squared distances between the fragment atoms and the
             corresponding receptor atoms of the pseudo-sphere procedure, one list
             for the vdW energy and another for the electrostatic interaction energy.
             In the case of the electrostatic interaction energy, it also provides the
             corresponding partial or total charges */
          SqDisFrRe_ps(FrAtNu,RoSFCo,ReCoor,ReMinC,GrSiCu_en,
              CubNum_en,CubFAI_en,CubLAI_en,CubLiA_en,
              PsSpNC,PsSphe,SDFrRe_ps,ReAtNu,PsSpRa,
              RePaCh,ReReNu,AtReprRes,FiAtRes,LaAtRes,
              TotChaRes,NuChResEn,LiChResEn,
              SDFrRe_ps_elec,ChFrRe_ps_elec);

          /* Compute vdW energy */
          PsSpEE(FrAtNu,ReAtNu,ReVdWE_sr,FrVdWE_sr,
              ReVdWR,FrVdWR,&VWEnEv_ps,SDFrRe_ps);

          /* Compute receptor desolvation, fragment desolvation and receptor-fragment
             interaction (with screening effect) energies (slow method) */
          ElecFrag(ReAtNu,ReCoor,RePaCh,ChFrRe_ps_elec,ReRad,ReRad2,
              ReRadOut,ReRadOut2,surfpt_re,nsurf_re,
              pointsrf_re,ReSelfVol,FrAtNu,RoSFCo,FrCoor,
              FrPaCh,FrRad,FrRad2,FrRadOut,FrRadOut2,Frdist2,
              SDFrRe_ps_elec,FrMinC,FrMaxC,&FrSolvEn,Nsurfpt_fr,
              surfpt_fr,nsurf_fr,pointsrf_fr,surfpt_ex,Tr,U1,U2,
              WaMoRa,GrSiSo,NPtSphere,Min,Max,XGrid,YGrid,ZGrid,
              NGridx,NGridy,NGridz,GridMat,DeltaPrDeso,Kelec,Ksolv,
              UnitVol,pi4,nxminBS,nyminBS,nzminBS,nxmaxBS,nymaxBS,
              nzmaxBS,corr_scrint,corr_fr_deso,&ReDesoElec,
              &ReFrIntElec,&FrDesoElec,ReSelfVol_corrB,EmpCorrB,FPaOut);

          VW_s=SFVWEn*VWEnEv_ps;
          In_s=SFIntElec*ReFrIntElec;
          Dr_s=SFDeso_re*ReDesoElec;
          Df_s=SFDeso_fr*FrDesoElec;
          To_s=VW_s+In_s+Dr_s+Df_s;


        }

        /* ----- FAST or BOTH methods ----- */
        if (((Solv_typ[0]=='f')||(Solv_typ[0]=='b'))&&(ChkInGrid)) {

          /* Van der Waals energy on a grid using the geometric mean approximation */
          VW_f=VWFrEn(FrAtNu,RoSFCo,BSMinC,VWGrIn,VWGrSi,FrVdWR_p3,
              FrVdWE_sr,VWGPoN,VWGrRP_at,VWGrRP_re,FPaOut,&vdwErr);

          /* Coulombic interaction energy on a grid */
          In_f=CoFrEn(FrAtNu,FrPaCh,RoSFCo,BSMinC,CoGrIn,CoGrSi,
              CoGPoN,CoGrRP,FPaOut,&clbErr);

          /* Compute receptor desolvation, fragment desolvation energies (fast method) */
          ElecFragFast(ReAtNu,ReCoor,
              ReRadOut,ReRadOut2,surfpt_re_deso,nsurf_re_deso,
              ReSurf_deso,pointsrf_re_deso,nlist_re,list_re,
              nstep_re,cbmid_re,midg_re,rscale_re,FrAtNu,
              RoSFCo,FrCoor,FrRad,FrRadOut,FrRadOut2,FrMinC,FrMaxC,
              FrRmax,Nsurfpt_fr_deso,surfpt_fr_deso_orig,
              surfpt_fr_deso,nsurf_fr_deso,
              FrSurf_deso,pointsrf_fr_deso,Tr,U1,U2,
              WaMoRa,rslop,MaxNeigh,&ReDesoElec,&FrDesoElec);

          VW_f=SFVWEn*VW_f;
          In_f=SFIntElec*In_f;
          Dr_f=SFDeso_re*ReDesoElec;
          Df_f=SFDeso_fr*FrDesoElec;
          To_f=VW_f+In_f+Dr_f+Df_f;

        }
        /* OLD OUTPUT FULL
           fprintf(FPaOut,"Position of Frag_num %d :\n\n",i);

           fprintf(FPaOut,"                          intermolecular                    ");
           fprintf(FPaOut,"electrostatic_desolvation            total_energy\n");
           fprintf(FPaOut,"                 van_der_Waals       electrostatic         ");
           fprintf(FPaOut,"receptor            fragment\n");
           fprintf(FPaOut,"                fast      acc.      fast      acc.      ");
           fprintf(FPaOut,"fast      acc.      fast      acc.      fast      acc.\n\n");

           fprintf(FPaOut,"      Conf      VW_f      VW_s      In_f      In_s");
           fprintf(FPaOut,"      Dr_f      Dr_s      Df_f      Df_s      To_f      To_s\n");

           fprintf(FPaOut,
           "     %5d%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n",
           Ind_num_cn,VW_f,VW_s,In_f,In_s,Dr_f,Dr_s,Df_f,Df_s,
           To_f,To_s);

         */

        /* NEW OUTPUT, fast energies not printed */

        fprintf(FPaOut,"                 intermolecular        electrostat_desolv.");
        fprintf(FPaOut,"        Total\n");
        fprintf(FPaOut,"  Conform.     vdWaals  electrost");
        fprintf(FPaOut,"     receptor    fragment       energy\n\n");

        fprintf(FPaOut,
            "%10d%12.2f%11.2f%13.2f%12.2f%13.2f\n",
            Ind_num_cn,VW_s,In_s,Dr_s,Df_s,
            To_s);
        // if energy evaluation run is requested, the energies are saved in
        // the best output summary table. clangini
        HeAtCo = count_heavy_atom(FrAtEl_nu, FrAtNu);
        MolWei = molecular_weight(FrAtEl_nu, FrAtNu, AtWei);
        if (write_best_sumtab_opt[0]=='y'){
          // append to _best_pproc summary table
          strcpy(TabFil,"./outputs/seed_best.dat");
          TabOutStream.open (TabFil, std::ios::out | std::ios::app); // append mode
          if(TabOutStream.is_open()){
            //for (j=1;j<=((NuPosMem<NuPosSdCl)?NuPosMem:NuPosSdCl);j++) {
            sprintf(TabLin,
                  "%-30s%8d%10d%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10d%10.4f",
                  FragNa,1,
                  1,
                  To_s,
                  In_s,
                  Dr_s,
                  Df_s,
                  VW_s,
                  (In_s-FrSolvEn),
                  FrSolvEn,(To_s/HeAtCo),
                  (VW_s/HeAtCo),
                  (In_s/HeAtCo),
                  HeAtCo,MolWei);
            TabOutStream << TabLin << std::endl;
            //}
          } else {
            std::cerr << "Unable to write to file "<< TabFil << std::endl;
          }
          TabOutStream.close();
        }
        free_dmatrix(FrCoor_NoAlign,1,FrAtNu,1,3);//clangini
      }
      else {
        /* --------------------- */
        /* Situation of conflict */
        /* --------------------- */
        /* Conflict between the fact that the fragment cannot be seeded as polar by
           SEED (because of the implemented rules) and the user choice of seeding it
           as a polar fragment */
        fprintf(FPaOut,"WARNING The fragment type %d cannot be ",CurFra);
        fprintf(FPaOut,"seeded as a polar fragment\n\n");

        /* Artificial memory allocation */
        ToNFrP=0;
        if (Ind_num_cn==1) {
          FrCoPo=ppdvector(1,ToNFrP);
          ConfArr=ivector(1,ToNFrP);
        }

        SFWrNu=0;

        FrNuTo=0;
        FrNuPB=0;
        FrNuSp=0;

      }


      free_dmatrix(SeFrCo,1,FrAtNu,1,3);
      free_dmatrix(RoSFCo,1,FrAtNu,1,3);
      free_dmatrix(SDFrRe_ps,1,FrAtNu,1,ReAtNu);
      free_dmatrix(SDFrRe_ps_elec,1,FrAtNu,1,ReAtNu);
      free_dmatrix(ChFrRe_ps_elec,1,FrAtNu,1,ReAtNu);

      free_dmatrix(Frdist2,1,FrAtNu,1,FrAtNu);
      free_dvector(FrRadOut2,1,FrAtNu);
      free_dvector(FrRadOut,1,FrAtNu);
      free_dvector(FrRad2,1,FrAtNu);
      free_dvector(FrRad,1,FrAtNu);
      if ((Solv_typ[0]=='s')||(Solv_typ[0]=='b')) {
        free_structpointvect(surfpt_fr,1,NPtSphere*FrAtNu);
        free_ivector(iatsrf_fr,1,NPtSphere*FrAtNu);
        free_ivector(pointsrf_fr,1,FrAtNu);
        free_ivector(nsurf_fr,1,FrAtNu);
        free_structpointvect(surfpt_ex,1,NPtSphere * (ReAtNu+FrAtNu));
      }
      if ((Solv_typ[0]=='f')||(Solv_typ[0]=='b')||(Solv_typ[0]=='p')) {
        NPtSphereMax_Fr = (int) (SurfDens_deso * pi4 * (FrRmax+WaMoRa) * 2);
        /* * 2 necessary otherwise seg.fault with e.g. tribrommethane */

        free_ivector(iatsrf_fr_deso,1,NPtSphereMax_Fr*FrAtNu);
        free_ivector(nsurf_fr_deso,1,FrAtNu);
        free_structpointvect(surfpt_fr_deso,1,NPtSphereMax_Fr*FrAtNu);
        free_ivector(pointsrf_fr_deso,1,FrAtNu);
        free_dvector(FrSurf_deso,1,Nsurfpt_fr);
        free_structpointvect(surfpt_fr_deso_orig,1,NPtSphereMax_Fr*FrAtNu);
      }

      free_ivector(FrDATy,1,FrAcNu+FrDoNu);
      free_ivector(FrDAAt,1,FrAcNu+FrDoNu);
      free_ivector(FrHydN,1,FrAcNu+FrDoNu);
      free_dmatrix(FrVeCo,1,FrAcNu+FrDoNu,1,6);
      free_dvector(FrVdWR,1,FrAtNu);
      free_dvector(FrVdWE,1,FrAtNu);
      free_dvector(FrVdWE_sr,1,FrAtNu);
      free_dvector(FrVdWR_p3,1,FrAtNu);
      free_ivector(HybFrAt,1,FrAtNu);
      free_ivector(UndisAt_fr,1,FrAtNu);
      free_ivector(FrAtoTyp_nu,1,FrAtNu);
      /* Lo1 : End */
    } /* end of for (Ind_num_cn=1;Ind_num_cn<=FrCoNu;Ind_num_cn++) clangini*/
    //std::cout << "End of for cycle on conformations" << std::endl; //clangini OK
    /* times(&timevar); */
    /* time_7=timevar.tms_utime+timevar.tms_stime; */
    gettimeofday(&time_7,NULL);
    //std::cout << "i: "<< i << std::endl; //clangini
    //std::cout << "SFWrNu: "<< SFWrNu << std::endl; //clangini
    /*SFWrNu_ar[i]=SFWrNu; clangini */
    SFWrNu_ar=SFWrNu; //clangini

    fprintf(FPaOut,"\n");
    if (EvalEn[0]!='e') {
      fprintf(FPaOut,"Total number of generated fragments of type %d (%s) : %d\n",
          i,FragNa,FrNuTo); /* clangini */
      //std::cout << FrNuTo << " " << ToNFrP << std::endl; //clangini
      //std::cout << "HERE" << std::endl; //clangini
      /*fprintf(FPaOut,"Total number of generated fragments of type %d (%s) : %d\n",
          i,&FrFiNa_out[CurFra][1],FrNuTo); clangini*/
      if (SpPoCh_opt[0]=='y')
        fprintf(FPaOut,"Fragments that passed the sphere checking : %d\n",FrNuSp);
      fprintf(FPaOut,"Fragments that passed the bump checking : %d\n",FrNuPB);
      if (FrNuVdW_ap>=0)
        fprintf(FPaOut,"Fragments that passed the vdW energy cutoff : %d\n",
            FrNuVdW_ap);
      fprintf(FPaOut,"Fragments that passed the total energy cutoff : %d\n",
          SFWrNu);
      fprintf(FPaOut,"CPU time in sec. for the seeding of %s : %.2f\n",FragNa,
          ((time_7.tv_sec  - time_6.tv_sec) * 1000000u +
           time_7.tv_usec - time_6.tv_usec) / 1.e6);/* clangini */
      /*fprintf(FPaOut,"CPU time in sec. for the seeding of %s : %.2f\n",
          &FrFiNa_out[CurFra][1],
          ((time_7.tv_sec  - time_6.tv_sec) * 1000000u +
           time_7.tv_usec - time_6.tv_usec) / 1.e6); clangini */
      /* &ResN_fr[i][1],(time_7-time_6)*0.01); */
      fprintf(FPaOut,"\n");
    }
    /* !!! From here, the value of SFWrNu might change */
    SFWrNu_init=SFWrNu;
    SFWrNu=((SFWrNu<MaxPosClus)?SFWrNu:MaxPosClus);

    if (EvalEn[0]!='e')
      fprintf(FPaOut,"Number of positions conserved for clustering : %d\n\n",SFWrNu);

    fclose(FPaOut);

    /* Co2 : Continue with the current fragment only if seeded positions of the
       fragment were found */
    if (SFWrNu_init) {

      /* Read the various energies in the output file, sort the total energies
         and cluster the fragments */

      /* times(&timevar); */
      /* time_9=timevar.tms_utime+timevar.tms_stime; */
      gettimeofday(&time_9,NULL);

      Coulo_ro=dvector(1,SFWrNu_init);
      Vande_ro=dvector(1,SFWrNu_init);
      TotEn_ro=dvector(1,SFWrNu_init);
      TotEn_ro_cp=dvector(1,SFWrNu_init);

      /* REMOVED  -> reading of energies from "do_not_touch_ener" */
      /* ReaOut(CurFra,SFWrNu_init,Index_ro,VW_f_ro,VW_s_ro,In_f_ro,In_s_ro,  */
      /* 	     Dr_f_ro,Dr_s_ro,Df_f_ro,Df_s_ro,To_f_ro,To_s_ro);  */

      /* Fill some arrays which are used in the next part (writing files
         and clustering) */
      if ((Solv_typ[0]=='f')||(Solv_typ[0]=='p')) {
        for (j=1;j<=SFWrNu_init;j++) {
          Vande_ro[j]=VW_f_ro[j];
          Coulo_ro[j]=In_f_ro[j];
          TotEn_ro[j]=To_f_ro[j];
        }
      }
      else if ((Solv_typ[0]=='s')||(Solv_typ[0]=='b')) {
        for (j=1;j<=SFWrNu_init;j++) {
          Vande_ro[j]=VW_s_ro[j];
          Coulo_ro[j]=In_s_ro[j];
          TotEn_ro[j]=To_s_ro[j];
        }
      }

      /* Sorting the fragments with respect to the total energy */
      for (j=1;j<=SFWrNu_init;j++)
        TotEn_ro_cp[j]=TotEn_ro[j];
      Sort(SFWrNu_init,Index_ro,TotEn_ro_cp); // This is the sorting according to fast energy. clangini

      FPaOut=fopen(OutFil,"a");
      /*OLD OUTPUT FULL
        fprintf(FPaOut,"Sorted frag. of type %d (%s) :\n",CurFra,&ResN_fr[i][1]);
        if (Solv_typ[0]=='b')

        fprintf(FPaOut,"   Num Fr_nu  Conf    VW_f    VW_s    In_f    In_s");
        fprintf(FPaOut,"    Dr_f    Dr_s    Df_f    Df_s    To_f    To_s\n");
        for (j=1;j<=((NuLiEnClus<SFWrNu)?NuLiEnClus:SFWrNu);j++)
        fprintf(FPaOut,
        "%6d%6d%6d%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n",
        j,Index_ro[j],ConfArr[Index_ro[j]],
        VW_f_ro[Index_ro[j]],VW_s_ro[Index_ro[j]],In_f_ro[Index_ro[j]],
        In_s_ro[Index_ro[j]],Dr_f_ro[Index_ro[j]],Dr_s_ro[Index_ro[j]],
        Df_f_ro[Index_ro[j]],Df_s_ro[Index_ro[j]],To_f_ro[Index_ro[j]],
        To_s_ro[Index_ro[j]]);
       */

      /*NEW OUTPUT, shortened */

      fprintf(FPaOut,"Sorted frag. of type %d (%s) :\n",CurFra,FragNa); /* clangini */
      /*fprintf(FPaOut,"Sorted frag. of type %d (%s) :\n",CurFra,&FrFiNa_out[CurFra][1]); clangini */
      if (Solv_typ[0]=='b') {

        fprintf(FPaOut,"   Num    Fr_nu  Conf    VW_s    In_s");
        fprintf(FPaOut,"    Dr_s    Df_s    To_s\n");
        for (j=1;j<=((NuLiEnClus<SFWrNu)?NuLiEnClus:SFWrNu);j++)
          fprintf(FPaOut,
              "%6d%9d%6d%8.2f%8.2f%8.2f%8.2f%8.2f\n", /* There was a bug here. added %8.2f clangini*/
              j,Index_ro[j],ConfArr[Index_ro[j]],
              VW_s_ro[Index_ro[j]],
              In_s_ro[Index_ro[j]],Dr_s_ro[Index_ro[j]],
              Df_s_ro[Index_ro[j]],
              To_s_ro[Index_ro[j]]);

      }
      else if (Solv_typ[0]=='f') {

        fprintf(FPaOut,"   Num    Fr_nu  Conf    VW_f    In_f");
        fprintf(FPaOut,"    Dr_f    Df_f    To_f\n");
        for (j=1;j<=((NuLiEnClus<SFWrNu)?NuLiEnClus:SFWrNu);j++)
          fprintf(FPaOut,
              "%6d%9d%6d%8.2f%8.2f%8.2f%8.2f%8.2f\n",
              j,Index_ro[j],ConfArr[Index_ro[j]],
              VW_f_ro[Index_ro[j]],In_f_ro[Index_ro[j]],
              Dr_f_ro[Index_ro[j]],
              Df_f_ro[Index_ro[j]],To_f_ro[Index_ro[j]]);

      }
      else if ((Solv_typ[0]=='p')&&(PrintLev>0)) {

        fprintf(FPaOut,"   Num    Fr_nu  Conf    VW_f    In_f");
        fprintf(FPaOut,"    Dr_f    Df_f    To_f\n");
        for (j=1;j<=((NuLiEnClus<SFWrNu)?NuLiEnClus:SFWrNu);j++)
          fprintf(FPaOut,
              "%6d%9d%6d%8.2f%8.2f%8.2f%8.2f%8.2f\n",
              j,Index_ro[j],ConfArr[Index_ro[j]],
              VW_f_ro[Index_ro[j]],In_f_ro[Index_ro[j]],
              Dr_f_ro[Index_ro[j]],
              Df_f_ro[Index_ro[j]],To_f_ro[Index_ro[j]]);

      }
      else if (Solv_typ[0]=='s') {

        fprintf(FPaOut,"   Num    Fr_nu  Conf    VW_s    In_s");
        fprintf(FPaOut,"    Dr_s    Df_s    To_s\n");
        for (j=1;j<=((NuLiEnClus<SFWrNu)?NuLiEnClus:SFWrNu);j++)
          fprintf(FPaOut,
              "%6d%9d%6d%8.2f%8.2f%8.2f%8.2f%8.2f\n",
              j,Index_ro[j],ConfArr[Index_ro[j]],
              VW_s_ro[Index_ro[j]],
              In_s_ro[Index_ro[j]],Dr_s_ro[Index_ro[j]],
              Df_s_ro[Index_ro[j]],To_s_ro[Index_ro[j]]);

      }
      fprintf(FPaOut,"\n");
      fclose(FPaOut);

      if (Solv_typ[0]!='p')     /* Modification NM 26-05-2004 */
      {
        if (FrAtNu>=3){ }
          /* write_charmm(CurFra,SFWrNu,FrAtNu,Index_ro,Coulo_ro,
              Vande_ro,TotEn_ro,FrAtTy,FrCoPo,ResN_fr,
              FrFiNa_out); clangini function no longer needed*/
              /*To be changed. ResN_fr does not exist any more. clangini*/
      }

      /* If both (SLOW and FAST) methods are used, show the two scoring lists */
      /*
         if (Solv_typ[0]=='b') {
         Index_both=ivector(1,SFWrNu);
         TotEn_both=vector(1,SFWrNu);
         for (j=1;j<=SFWrNu;j++) {
         Index_both[j]=j;
         TotEn_both[j]=To_f_ro[j];
         }
         Sort(SFWrNu,Index_both,TotEn_both);
         FPaOut=fopen(OutFil,"a");
         fprintf(FPaOut,"Two scoring lists for both (SLOW and FAST) methods\n");
         for (j=1;j<=((NuLiEnClus<SFWrNu)?NuLiEnClus:SFWrNu);j++) {
         LoopTe=1;
         for (i1=1;((i1<=SFWrNu)&&LoopTe);i1++) {
         if (Index_both[i1]==Index_ro[j]) {
         fprintf(FPaOut,"%7d%7d\n",j,i1);
         LoopTe=0;
         }
         }
         }
         fprintf(FPaOut,"\n");
         fclose(FPaOut);
         free_ivector(Index_both,1,SFWrNu);
         free_vector(TotEn_both,1,SFWrNu);
         }
       */

      /* Clustering with the GSEAL method in two steps */

      FrCoor_clus=dmatrix(1,FrAtNu,1,3);

      /* First clustering */

      /* Compute the normalization for the similarity number, which is the same
         for all the fragments of the same type and of the same conformation,
         because the translation and rotation are done with a fragment which is
         kept rigid */ //Normalization is max(S_AA,S_BB) clangini
      /* clangini START */
      FraSim_no=dvector(1,FrCoNu); // Different conformer (not poses!) of the same molecule can have different similarities! clangini
                                  // In the case of the new seed this is always 1 as different conformers are treated as
                                  // different fragments. clangini
      for (j=1;j<=FrCoNu;j++) { /*Now FrCoNu is always 1*/
        /*ReFrFi_coo_mol2(CurFra,j,FrFiNa,FrAtNu,FrCoor_clus); want to avoid to reread it*/
        for (i1 = 1; i1 <= FrAtNu; i1++){ /*copy matrix (is it really necessary??)*/
          FrCoor_clus[i1][1] = FrCoor[i1][1];
          FrCoor_clus[i1][2] = FrCoor[i1][2];
          FrCoor_clus[i1][3] = FrCoor[i1][3];
        }
        FraSim_no[j]=0.0;
        for (i1=1;i1<=FrAtNu;i1++) {
          for (i2=1;i2<=FrAtNu;i2++) {
             if(SimWei[FrAtEl_nu[i1]][FrAtEl_nu[i2]]!=0)
             {
               Dist_sq=DistSq(FrCoor_clus[i1][1],FrCoor_clus[i1][2],
                   FrCoor_clus[i1][3],FrCoor_clus[i2][1],
                   FrCoor_clus[i2][2],FrCoor_clus[i2][3]);
               FraSim_no[j]=FraSim_no[j]+SimWei[FrAtEl_nu[i1]][FrAtEl_nu[i2]]*
                 expf(-SimExp*Dist_sq);
             }
           }
         }
      }
      /*FraSim_no=vector(1,FrCoNu);
      for (j=1;j<=FrCoNu;j++) {
        ReFrFi_coo_mol2(CurFra,j,FrFiNa,FrAtNu,FrCoor_clus);
        FraSim_no[j]=0.0;
        for (i1=1;i1<=FrAtNu;i1++) {
          for (i2=1;i2<=FrAtNu;i2++) {
            if(SimWei[FrAtEl_nu[i1]][FrAtEl_nu[i2]]!=0)
            {
              Dist_sq=DistSq(FrCoor_clus[i1][1],FrCoor_clus[i1][2],
                  FrCoor_clus[i1][3],FrCoor_clus[i2][1],
                  FrCoor_clus[i2][2],FrCoor_clus[i2][3]);
              FraSim_no[j]=FraSim_no[j]+SimWei[FrAtEl_nu[i1]][FrAtEl_nu[i2]]*
                expf(-SimExp*Dist_sq);
            }
          }
        }
      } I think it is unnecessary to reread FrCoor clangini*/
      /* clangini END */
      // Here first clustering starts. clangini
      FraEqu=ivector(1,SFWrNu_init);
      for (i1=1;i1<=SFWrNu_init;i1++)
        FraEqu[i1]=i1;
      ClusNu=SFWrNu; // start with number of clusters = number of poses to consider
      for (i1=1;i1<=SFWrNu;i1++) {
        if (FraEqu[Index_ro[i1]]==Index_ro[i1]) {
          for (i2=i1+1;i2<=SFWrNu;i2++) {
            if (FraEqu[Index_ro[i2]]==Index_ro[i2]) {
              Simila(Index_ro[i1],Index_ro[i2],FrAtNu,FrAtEl_nu,FrCoPo,SimWei,
                     SimExp,&FraSim);
              FraSim_no_max=(FraSim_no[ConfArr[Index_ro[i1]]]>
                             FraSim_no[ConfArr[Index_ro[i2]]])?
                             FraSim_no[ConfArr[Index_ro[i1]]]:
                             FraSim_no[ConfArr[Index_ro[i2]]];
              if ((FraSim/FraSim_no_max)>SimCut) {
                FraEqu[Index_ro[i2]]=Index_ro[i1];
                ClusNu=ClusNu-1;
              }
            }
          }
        }
      }
      // actual clustering ends, but still need to rework indexes. clangini
      /* Making of the indexes list and of the cluster list done with GSEAL */
      ClusIn=ivector(1,ClusNu);
      ClusLi=ivector(1,SFWrNu);
      ClusVa_1=0;
      ClusVa_2=0;
      for (i1=1;i1<=SFWrNu;i1++) {
        if (FraEqu[Index_ro[i1]]==Index_ro[i1]) { // if (representative of cluster). clangini
          ClusVa_1=ClusVa_1+1;
          ClusLi[ClusVa_1]=Index_ro[i1];
          for (i2=i1+1;i2<=SFWrNu;i2++) { // loop looking for the other members of the cluster
            if (FraEqu[Index_ro[i2]]==Index_ro[i1]) {
              ClusVa_1=ClusVa_1+1;
              ClusLi[ClusVa_1]=Index_ro[i2];
            }
          }
          ClusVa_2=ClusVa_2+1;
          ClusIn[ClusVa_2]=ClusVa_1;
        }
      }

      FPaOut=fopen(OutFil,"a");
      fprintf(FPaOut,"Number of clusters after first clustering ");
      fprintf(FPaOut,"for the fragment type %d (%s) : %d\n",
          CurFra,FragNa,ClusNu);
      /*fprintf(FPaOut,"for the fragment type %d (%s) : %d\n",
          CurFra,&FrFiNa_out[CurFra][1],ClusNu); clangini*/
      /*
         ClusVa_1=0;
         ClusVa_2=0;
         ClusVa_wr=0;
         for (i1=1;i1<=ClusNu;i1++) {
         ClusVa_1=ClusVa_1+1;
         ClusVa_2=ClusVa_2+1;
         ClusVa_wr=ClusVa_wr+1;
         if (ClusVa_wr<=NuLiEnClus)
         fprintf(FPaOut,"Representative : %d\n",ClusLi[ClusVa_1]);
         for (i2=ClusVa_1+1;i2<=ClusIn[ClusVa_2];i2++) {
         ClusVa_1=ClusVa_1+1;
         ClusVa_wr=ClusVa_wr+1;
         if (ClusVa_wr<=NuLiEnClus)
         fprintf(FPaOut,"                      %d\n",ClusLi[ClusVa_1]);
         }
         }
       */
      fprintf(FPaOut,"\n");
      fclose(FPaOut);

      /* Second clustering */

      /* Compute the normalization for the similarity number */
      /* clangini START */
      FraSim_no_sd=dvector(1,FrCoNu);
      for (j=1;j<=FrCoNu;j++) {
        /*ReFrFi_coo_mol2(CurFra,j,FrFiNa,FrAtNu,FrCoor_clus);*/
        for (i1 = 1; i1 <= FrAtNu; i1++){ /*copy matrix (is it really necessary??)*/
          FrCoor_clus[i1][1] = FrCoor[i1][1];
          FrCoor_clus[i1][2] = FrCoor[i1][2];
          FrCoor_clus[i1][3] = FrCoor[i1][3];
        }
        FraSim_no_sd[j]=0.0;
        for (i1=1;i1<=FrAtNu;i1++) {
          for (i2=1;i2<=FrAtNu;i2++) {
            if( SimWei[FrAtEl_nu[i1]][FrAtEl_nu[i2]] != 0 )
            {
              Dist_sq=DistSq(FrCoor_clus[i1][1],FrCoor_clus[i1][2],
                  FrCoor_clus[i1][3],FrCoor_clus[i2][1],
                  FrCoor_clus[i2][2],FrCoor_clus[i2][3]);
              FraSim_no_sd[j]=FraSim_no_sd[j]+SimWei[FrAtEl_nu[i1]][FrAtEl_nu[i2]]*
                expf(-SimExp_sd*Dist_sq);
            }
          }
        }
      }

      /*FraSim_no_sd=vector(1,FrCoNu);
      for (j=1;j<=FrCoNu;j++) {
        ReFrFi_coo_mol2(CurFra,j,FrFiNa,FrAtNu,FrCoor_clus);
        FraSim_no_sd[j]=0.0;
        for (i1=1;i1<=FrAtNu;i1++) {
          for (i2=1;i2<=FrAtNu;i2++) {
            if( SimWei[FrAtEl_nu[i1]][FrAtEl_nu[i2]] != 0 )
            {
              Dist_sq=DistSq(FrCoor_clus[i1][1],FrCoor_clus[i1][2],
                  FrCoor_clus[i1][3],FrCoor_clus[i2][1],
                  FrCoor_clus[i2][2],FrCoor_clus[i2][3]);
              FraSim_no_sd[j]=FraSim_no_sd[j]+SimWei[FrAtEl_nu[i1]][FrAtEl_nu[i2]]*
                expf(-SimExp_sd*Dist_sq);
            }
          }
        }
      } clangini*/
      /* clangini END*/
      // Second clustering starts. clangini
      FraEqu_sd=ivector(1,SFWrNu_init);
      for (i1=1;i1<=SFWrNu_init;i1++)
        FraEqu_sd[i1]=i1;
      ClusVa_1=1;
      ClusVa_2=0;
      for (i1=1;i1<=ClusNu;i1++) { // seems i1 = ClusVa_2 always! // loops over first clusters
        ClusVa_2=ClusVa_2+1;
        if (TotEn_ro[ClusLi[ClusVa_1]]<=FrMaEn_sd) { // I do second clustering only if
          //representative (first pose of the cluster -
          //which is also the most energetically favorable) has energy below cutoff. clangini
          for (i2=ClusVa_1;i2<=ClusIn[ClusVa_2];i2++) {
            if (FraEqu_sd[ClusLi[i2]]==ClusLi[i2]) { // check if not already assigned to another (second) cluster
              for (i3=i2+1;i3<=ClusIn[ClusVa_2];i3++) {
                if (FraEqu_sd[ClusLi[i3]]==ClusLi[i3]) { // check if not already assigned to another (second) cluster
                  Simila(ClusLi[i2],ClusLi[i3],FrAtNu,FrAtEl_nu,FrCoPo,SimWei,
                         SimExp_sd,&FraSim);
                  FraSim_no_max=(FraSim_no_sd[ConfArr[ClusLi[i2]]]>
                                 FraSim_no_sd[ConfArr[ClusLi[i3]]])?
                                 FraSim_no_sd[ConfArr[ClusLi[i2]]]:
                                 FraSim_no_sd[ConfArr[ClusLi[i3]]]; // The presence of ConfArr is no longer needed in my case
                                                                   // as differnet conformers are considered as different fragments
                  if ((FraSim/FraSim_no_max)>SimCut_sd)
                    FraEqu_sd[ClusLi[i3]]=ClusLi[i2];
                }
              }
            }
          }
        }
        ClusVa_1=ClusIn[ClusVa_2]+1;
      }
      // Second clustering ends.
      FPaOut=fopen(OutFil,"a");
      if (PrintLev>1) {
        fprintf(FPaOut,"Clusters after second clustering (with energy cutoff) ");
        fprintf(FPaOut,"for fragment type %d (%s) :\n",CurFra,FragNa); /* clangini */
        /*fprintf(FPaOut,"for fragment type %d (%s) :\n",CurFra,&FrFiNa_out[CurFra][1]); clangini */
      }
      ClusLi_sd=ivector(1,SFWrNu);
      ClusLi_sd_01=ivector(1,SFWrNu);
      ClusLi_sd_01_reduc=ivector(1,SFWrNu);
      ClusLi_sd_pproc=ivector(1,SFWrNu);
      for (i1=1;i1<=SFWrNu;i1++) {
        ClusLi_sd_01[i1]=0;
        ClusLi_sd_01_reduc[i1]=0;
        ClusLi_sd_pproc[i1]=0;
      }
      ClusVa_1=1;
      ClusVa_2=0;
      ClusVa_3=0;
      ClusVa_wr=0;
      for (i1=1;i1<=ClusNu;i1++) {//ClusNu=#of clusters after first clustering. clangini
        ClusVa_2=ClusVa_2+1;
        if (TotEn_ro[ClusLi[ClusVa_1]]<=FrMaEn_sd) {
          PrintVa=1;
          for (i2=ClusVa_1;i2<=ClusIn[ClusVa_2];i2++) {//loop within poses of first cluster
            if (FraEqu_sd[ClusLi[i2]]==ClusLi[i2]) { // check if representative of second cluster
              if (PrintVa) {
                ClusVa_wr=ClusVa_wr+1;
                if ((ClusVa_wr<=NuLiEnClus)&&(PrintLev>1))
                  fprintf(FPaOut,"                 %d\n",ClusLi[i2]);
                PrintVa=0;
                ClusVa_3=ClusVa_3+1;
                ClusLi_sd[ClusVa_3]=ClusLi[i2];
                ClusLi_sd_01[ClusVa_3]=1; // Representative of second GSEAL
                ClusLi_sd_01_reduc[ClusVa_3]=1;
                ClusLi_sd_pproc[ClusVa_3]=2; // Representative of first GSEAL
                ClusVa_4=1;
              } else { // else = if(PrintVa == 0)
                ClusVa_wr=ClusVa_wr+1;
                if ((ClusVa_wr<=NuLiEnClus)&&(PrintLev>1))
                  fprintf(FPaOut,"                      %d\n",ClusLi[i2]);
                ClusVa_3=ClusVa_3+1;
                ClusLi_sd[ClusVa_3]=ClusLi[i2];
                ClusLi_sd_01[ClusVa_3]=1; // Representative of second GSEAL
                if (ClusVa_4<NuClusMem) { // if not exceeding the maximum number
                                          //of cluster members
                  ClusLi_sd_01_reduc[ClusVa_3]=1;
                  ClusLi_sd_pproc[ClusVa_3]=1;
                  ClusVa_4=ClusVa_4+1;
                }
              }
              for (i3=i2+1;i3<=ClusIn[ClusVa_2];i3++) { // loop over the next poses of same first cluster
                if (FraEqu_sd[ClusLi[i3]]==ClusLi[i2]) {// check if belongs to the same second cluster
                  ClusVa_wr=ClusVa_wr+1;
                  if ((ClusVa_wr<=NuLiEnClus)&&(PrintLev>1))
                    fprintf(FPaOut,"                           %d\n",
                            ClusLi[i3]);
                  ClusVa_3=ClusVa_3+1;
                  ClusLi_sd[ClusVa_3]=ClusLi[i3];
                }
              }
            }
          }
        }
        ClusVa_1=ClusIn[ClusVa_2]+1;
      }
      fprintf(FPaOut,"\n");
      fclose(FPaOut);

      if (Solv_typ[0]!='p')     /* Modification NM 26-05-2004 */
      {
        if (NonAliHy>=3 && write_pproc_chm_opt[0]=='y')
        {

          /*write_chm_clus(CurFra,SFWrNu,FrAtNu,AliHyd,Coulo_ro,Vande_ro,TotEn_ro,
              FrAtTy,FrCoPo,ClusLi_sd,ClusLi_sd_01,NonAliHy,ResN_fr,
              FrFiNa_out); clangini *//*Need to be removed. clangini*/
          /*write_chm_clus_reduc(CurFra,SFWrNu,FrAtNu,AliHyd,Coulo_ro,Vande_ro,
              TotEn_ro,FrAtTy,FrCoPo,ClusLi_sd,
              ClusLi_sd_01_reduc,NonAliHy,ResN_fr,FrFiNa_out); clangini*/
              /*Need to be removed. clangini*/
        }
      }

      /* times(&timevar); */
      /* time_10=timevar.tms_utime+timevar.tms_stime; */
      gettimeofday(&time_10,NULL);

      FPaOut=fopen(OutFil,"a");
      fprintf(FPaOut,"CPU time in sec. for the sorting and clusterings of ");
      /*fprintf(FPaOut,"%s : %.2f\n\n",&FrFiNa_out[CurFra][1],
          ((time_10.tv_sec  - time_9.tv_sec) * 1000000u +
           time_10.tv_usec - time_9.tv_usec) / 1.e6); clangini */
      fprintf(FPaOut,"%s : %.2f\n\n",FragNa,
          ((time_10.tv_sec  - time_9.tv_sec) * 1000000u +
           time_10.tv_usec - time_9.tv_usec) / 1.e6); /* clangini */

      /* (time_10-time_9)*0.01); */
      fclose(FPaOut);

      if (Solv_typ[0]!='p')     /* Modification NM 26-05-2004 */ /* This has to be modified. clangini */
      {
        /*for (i1=1;i1<=((SFWrNu<FiNuMa)?SFWrNu:FiNuMa);i1++) {
          sprintf(WriPat,"%s%s%s%d%s","./outputs/",&FrFiNa_out[CurFra][1],"_",
              Index_ro[i1],".mol2\0");
          sprintf(WriNam,"%s%d%s%d%s","seeded_frag_",CurFra,"_",
              Index_ro[i1],"\0");
          write_mol2(WriPat,WriNam,FrAtNu,FrBdNu,FrAtEl,FrCoPo[Index_ro[i1]],
              FrAtTy,FrPaCh,FrBdAr,FrBdTy,FrSubN,FrSubC,ResN_fr,CurFra);
        } clangini write_mol2 has to be modified */
      }



      /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
      /* Postprocessing : Compute the energy with the slow method for some members   */
      /*                  of the second clustering step created with the fast method */
      /* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

      /*
         FPaOut=fopen(OutFil,"a");
         fprintf(FPaOut,"\nPostProcess\n");
         fprintf(FPaOut,"   i1     ClusLi_sd[i1]      ClusLi_sd_pproc[i1]\n");
         for (i1=1;i1<=SFWrNu;i1++) {
         fprintf(FPaOut,"    %6d  %6d  %6d\n",i1,ClusLi_sd[i1],ClusLi_sd_pproc[i1]);
         }
         fprintf(FPaOut,"\n\n\n");
         fclose(FPaOut);
       */

      if (Solv_typ[0]=='p') {

        FPaOut=fopen(OutFil,"a");

        /* times(&timevar); */
        /* time_11=timevar.tms_utime+timevar.tms_stime; */
        gettimeofday(&time_11,NULL);

        /* Loop over the conformations */
        for (Ind_num_cn=1;Ind_num_cn<=FrCoNu;Ind_num_cn++) { /* This loop is no longer needed. clangini */
        /* Ind_num_cn should be now always be 1 clangini */
          /* Read the coordinates of the current conformation Ind_num_cn of the current
             fragment type CurFra */

          /*ReFrFi_coo_mol2(CurFra,Ind_num_cn,FrFiNa,FrAtNu,FrCoor); Think FrCoor is still in memory*/

          /* Assign the atom element numbers, van der Waals radii and energies
             for the current fragment */
          FrVdWR=dvector(1,FrAtNu);
          FrVdWE=dvector(1,FrAtNu);
          FrAtoTyp_nu=ivector(1,FrAtNu);
          if(! AssiRE(NumbAT,AtTyAr,VdWRad,VdWEne,FrAtNu,FrAtTy,FrVdWR,FrVdWE,FPaOut,
                AtENAr,FrAtEl_nu,FrAtoTyp_nu) )
          {
            /* bail on missing types */
            printf("Program exits\n");
            exit(13);
          }

          /*---------------------------------------------------------
            PART NEEDED FOR CONTINUUM ELECTROSTATICS
            ---------------------------------------------------------*/

          /* Get the charge radii and the internal distances for the current fragment */
          FrRad=dvector(1,FrAtNu);
          FrRad2=dvector(1,FrAtNu);
          FrRadOut=dvector(1,FrAtNu);
          FrRadOut2=dvector(1,FrAtNu);
          Frdist2=dmatrix(1,FrAtNu,1,FrAtNu);
          nn = get_Ch_Rad_Fr(FrAtNu,FrCoor,FrVdWR,WaMoRa,FrRad,FrRad2,FrRadOut,
              FrRadOut2,Frdist2,&Rmin,&FrRmax);

          /* Get the extremes of frag coor */
          vect=dvector(1,FrAtNu);
          getColumnFrom2DArray(FrCoor, 1, 1, FrAtNu, vect);
          FrMinC[1] = MinVector(vect, 1, FrAtNu);
          FrMaxC[1] = MaxDVector(vect, 1, FrAtNu);
          getColumnFrom2DArray(FrCoor, 2, 1, FrAtNu, vect);
          FrMinC[2] = MinVector(vect, 1, FrAtNu);
          FrMaxC[2] = MaxDVector(vect, 1, FrAtNu);
          getColumnFrom2DArray(FrCoor, 3, 1, FrAtNu, vect);
          FrMinC[3] = MinVector(vect, 1, FrAtNu);
          FrMaxC[3] = MaxDVector(vect, 1, FrAtNu);
          free_dvector(vect,1,FrAtNu);

          /* Preparation for the slow method */
          /* Place points over the frag surface to describe its SAS1 (to get
             the volume occupied by the fragment - slow algorithm) */
          surfpt_fr=structpointvect(1,NPtSphere*FrAtNu);
          iatsrf_fr=ivector(1,NPtSphere*FrAtNu);
          pointsrf_fr=ivector(1,FrAtNu);
          for (ij=1;ij<=FrAtNu;ij++)
            pointsrf_fr[ij]=0;
          nsurf_fr=ivector(1,FrAtNu);
          for (iat=1;iat<=FrAtNu;iat++)
            nsurf_fr[iat] = 0;
          nn = Surf_Grid(FrAtNu,FrCoor,FrMaxC,FrMinC,FrRad,WaMoRa,NPtSphere,
              &Nsurfpt_fr,surfpt_fr,iatsrf_fr,nsurf_fr,pointsrf_fr);

          /* Now allocate space for the S&R surface resulting from the fragment
             together with the receptor ("extended fragment") */
          surfpt_ex=structpointvect(1,NPtSphere * (ReAtNu+FrAtNu));

          /* Calculate the solvation energy of the fragment without the receptor
             Use a very fine grid spacing (0.1) */
          nn = FragSolvEn(FrAtNu,FrCoor,FrPaCh,FrVdWR,FrRadOut,
              FrRadOut2,Frdist2,Nsurfpt_fr,surfpt_fr,WaMoRa,0.1,
              Ksolv,pi4,&FrSolvEn,EmpCorrB,FPaOut);
          /*fprintf(FPaOut,"Dielectric value = %4.1f -> %4.1f transfer energy for the fragment (%s) : ",
              DielRe,DielWa,&FrFiNa_out[CurFra][1]);*/
          fprintf(FPaOut,"Dielectric value = %4.1f -> %4.1f transfer energy for the fragment (%s) : ",
              DielRe,DielWa,FragNa); /* clangini */
          fprintf(FPaOut,"%10.4f\n\n",FrSolvEn);
          /*---------------------------------------------------------
            END OF THE PART NEEDED FOR CONTINUUM ELECTROSTATICS
            ---------------------------------------------------------*/

          /* Find the aliphatic hydrogens */
          FiAlHy(FrAtNu,FrBdNu,FrAtEl_nu,FrBdAr,AliHyd,&NonAliHy);

          /* Compute the square root of the van der Waals current fragment energies
             and the power of 3 of the fragment van der Waals radii */
          FrVdWE_sr=dvector(1,FrAtNu);
          SqRoEn(FrAtNu,FrVdWE,FrVdWE_sr);
          FrVdWR_p3=dvector(1,FrAtNu);
          for (j=1;j<=FrAtNu;j++)
            FrVdWR_p3[j]=FrVdWR[j]*FrVdWR[j]*FrVdWR[j];


          RoSFCo=dmatrix(1,FrAtNu,1,3);
          SDFrRe_ps=dmatrix(1,FrAtNu,1,ReAtNu);
          SDFrRe_ps_elec=dmatrix(1,FrAtNu,1,ReAtNu);
          ChFrRe_ps_elec=dmatrix(1,FrAtNu,1,ReAtNu);


          /* Loop over the list of the second clustering step and compute the energy
             with the slow method for the kept members of the current conformation */

          /* dey test -> shared FrAtNu,ReCoor,CubNum_en,GrSiCu_en,CubLAI_en,CubLiA_en */
          /* does not work yet: */
          /* none) shared(VW_s_ro,In_s_ro,Dr_f_ro,Df_s_ro,To_s_ro) */
          /* dbx check first private */
          /* #ifdef OMP  */
          /* #pragma omp parallel for default(shared) private(i1,i2,Tr,U1,U2,VWEnEv_ps,FrSolvEn,ReDesoElec,ReFrIntElec,FrDesoElec) firstprivate(RoSFCo) */
          /* #endif */
          for (i1=1;i1<=SFWrNu;i1++) {

            if (((ClusLi_sd_pproc[i1]==2)||(ClusLi_sd_pproc[i1]==1))&&
                (ConfArr[ClusLi_sd[i1]]==Ind_num_cn)) {

              for (i2=1;i2<=FrAtNu;i2++) {
                RoSFCo[i2][1]=FrCoPo[ClusLi_sd[i1]][i2][1];
                RoSFCo[i2][2]=FrCoPo[ClusLi_sd[i1]][i2][2];
                RoSFCo[i2][3]=FrCoPo[ClusLi_sd[i1]][i2][3];
              }

              //std::cout << "ClusLi_sd[i1] =========== : " << ClusLi_sd[i1] << std::endl;
              Rot_Tran(FrAtNu,FrCoor,RoSFCo,Tr,U1,U2);

              /* Compute the squared distances between the fragment atoms and the
                 corresponding receptor atoms of the pseudo-sphere procedure, one list
                 for the vdW energy and another for the electrostatic interaction energy.
                 In the case of the electrostatic interaction energy, it also provides the
                 corresponding partial or total charges */
              SqDisFrRe_ps(FrAtNu,RoSFCo,ReCoor,ReMinC,GrSiCu_en,
                  CubNum_en,CubFAI_en,CubLAI_en,CubLiA_en,
                  PsSpNC,PsSphe,SDFrRe_ps,ReAtNu,PsSpRa,
                  RePaCh,ReReNu,AtReprRes,FiAtRes,LaAtRes,
                  TotChaRes,NuChResEn,LiChResEn,
                  SDFrRe_ps_elec,ChFrRe_ps_elec);

              /* Compute vdW energy */
              PsSpEE(FrAtNu,ReAtNu,ReVdWE_sr,FrVdWE_sr,
                  ReVdWR,FrVdWR,&VWEnEv_ps,SDFrRe_ps);

              /* Compute receptor desolvation, fragment desolvation and receptor-fragment
                 interaction (with screening effect) energies (slow method) */
              ElecFrag(ReAtNu,ReCoor,RePaCh,ChFrRe_ps_elec,
                  ReRad,ReRad2,ReRadOut,
                  ReRadOut2,surfpt_re,nsurf_re,
                  pointsrf_re,ReSelfVol,FrAtNu,RoSFCo,FrCoor,
                  FrPaCh,FrRad,FrRad2,FrRadOut,FrRadOut2,
                  Frdist2,SDFrRe_ps_elec,FrMinC,FrMaxC,&FrSolvEn,
                  Nsurfpt_fr,surfpt_fr,
                  nsurf_fr,pointsrf_fr,surfpt_ex,Tr,U1,U2,WaMoRa,
                  GrSiSo,NPtSphere,Min,Max,XGrid,YGrid,ZGrid,
                  NGridx,NGridy,NGridz,GridMat,
                  DeltaPrDeso,Kelec,Ksolv,UnitVol,
                  pi4,nxminBS,nyminBS,nzminBS,nxmaxBS,nymaxBS,
                  nzmaxBS,corr_scrint,corr_fr_deso,&ReDesoElec,
                  &ReFrIntElec,&FrDesoElec,ReSelfVol_corrB,EmpCorrB,FPaOut);

              VW_s_ro[ClusLi_sd[i1]]=SFVWEn*VWEnEv_ps;
              In_s_ro[ClusLi_sd[i1]]=SFIntElec*ReFrIntElec;
              Dr_s_ro[ClusLi_sd[i1]]=SFDeso_re*ReDesoElec;
              Df_s_ro[ClusLi_sd[i1]]=SFDeso_fr*FrDesoElec;
              To_s_ro[ClusLi_sd[i1]]=VW_s_ro[ClusLi_sd[i1]]+In_s_ro[ClusLi_sd[i1]]+
                                     Dr_s_ro[ClusLi_sd[i1]]+Df_s_ro[ClusLi_sd[i1]];

            }

          }

          free_dmatrix(RoSFCo,1,FrAtNu,1,3);
          free_dmatrix(SDFrRe_ps,1,FrAtNu,1,ReAtNu);
          free_dmatrix(SDFrRe_ps_elec,1,FrAtNu,1,ReAtNu);
          free_dmatrix(ChFrRe_ps_elec,1,FrAtNu,1,ReAtNu);
          free_structpointvect(surfpt_ex,1,NPtSphere * (ReAtNu+FrAtNu));
          free_dmatrix(Frdist2,1,FrAtNu,1,FrAtNu);
          free_dvector(FrRadOut2,1,FrAtNu);
          free_dvector(FrRadOut,1,FrAtNu);
          free_dvector(FrRad2,1,FrAtNu);
          free_dvector(FrRad,1,FrAtNu);
          free_structpointvect(surfpt_fr,1,NPtSphere*FrAtNu);
          free_ivector(iatsrf_fr,1,NPtSphere*FrAtNu);
          free_ivector(pointsrf_fr,1,FrAtNu);
          free_ivector(nsurf_fr,1,FrAtNu);
          free_dvector(FrVdWR,1,FrAtNu);
          free_dvector(FrVdWE,1,FrAtNu);
          free_dvector(FrVdWE_sr,1,FrAtNu);
          free_dvector(FrVdWR_p3,1,FrAtNu);
          free_ivector(FrAtoTyp_nu,1,FrAtNu);

        } /* End of for (Ind_num_cn=1;Ind_num_cn<=FrCoNu;Ind_num_cn++) clangini*/

        /* times(&timevar); */
        /* time_12=timevar.tms_utime+timevar.tms_stime; */
        gettimeofday(&time_12,NULL);

        /* Compute NuSdClKe and NuPosSdCl */
        NuSdClKe=0;
        NuPosSdCl=0;
        //std::cout << "SFWrNu = " << SFWrNu << std::endl; //clangini
        for (j=1;j<=SFWrNu;j++) {
          //std::cout << "j = " << j << " ClusLi_sd_pproc[j] = " <<ClusLi_sd_pproc[j]<< std::endl; //clangini
          if (ClusLi_sd_pproc[j]==2) { // if it is a first cluster representative
            NuSdClKe=NuSdClKe+1;
            NuPosSdCl=NuPosSdCl+1;
          }
          if (ClusLi_sd_pproc[j]==1) // if it is a second cluster representative (reduced number)
            NuPosSdCl=NuPosSdCl+1;
        }

        /* Prepare lists of kept positions in second clusters for sorting
           IntVar1 numbering on kept positions
           IntVar2 numbering on clusters */
        FrPosAr_pproc=ivector(1,NuPosSdCl);
        SdClusAr_pproc=ivector(1,NuPosSdCl);
        TotEnSdClus_pproc=dvector(1,NuPosSdCl);

        IntVar1=0;
        IntVar2=0;
        for (j=1;j<=SFWrNu;j++) {
          if (ClusLi_sd_pproc[j]==2) { // if representative of first clustering
                                       //(and hence of second as well) clangini
            IntVar1=IntVar1+1;
            IntVar2=IntVar2+1;
            FrPosAr_pproc[IntVar1]=ClusLi_sd[j];
            SdClusAr_pproc[IntVar1]=IntVar2;
            TotEnSdClus_pproc[IntVar1]=To_s_ro[ClusLi_sd[j]];
            if(To_s_ro[ClusLi_sd[j]] > FrMaEn){
              PosToRem.push_back(IntVar1);
            }
          }
          if (ClusLi_sd_pproc[j]==1) {// if representative of second clustering
                                      //only but not of first. clangini
            IntVar1=IntVar1+1;
            FrPosAr_pproc[IntVar1]=ClusLi_sd[j];
            SdClusAr_pproc[IntVar1]=IntVar2;
            TotEnSdClus_pproc[IntVar1]=To_s_ro[ClusLi_sd[j]];
            if(To_s_ro[ClusLi_sd[j]] > FrMaEn){
              PosToRem.push_back(IntVar1);
            }
          }
        }

        /* Sort the kept positions of the kept second clusters with respect to
           the total energy computed with the slow method */
        Index_pproc=ivector(1,NuPosSdCl); // new index after sorting
        for (j=1;j<=NuPosSdCl;j++)
          Index_pproc[j]=j;
        Sort(NuPosSdCl,Index_pproc,TotEnSdClus_pproc);//This is the sorting on the slow energy. clangini

        /* Obtain new order of clusters with ordered members with respect to
           the total energy (slow method) */
        // Sort is sorting all the kept positions, here we split the ordered poses
        // into the clusters they belong to. clangini
        CluIndex_sort=ivector(1,NuSdClKe);// clangini
        for (j=1;j<= NuSdClKe;j++){
          CluIndex_sort[j] = j; // clangini
        }

        FrPosAr_sort=ivector(1,NuPosSdCl);
        SdClusAr_sort=ivector(1,NuPosSdCl);
        FlagAr=ivector(1,NuPosSdCl);
        for (j=1;j<=NuPosSdCl;j++){
          FlagAr[j]=1;
        }

        IntVar1=0; // index of cluster (energy ordered)
        IntVar2=0; // index
        for (i1=1;i1<=NuPosSdCl;i1++) {
          if (FlagAr[Index_pproc[i1]]) {
            IntVar2=IntVar2+1;
            // std::cout << "IntVar2 " << IntVar2 << "\n";
            CluIndex_sort[SdClusAr_pproc[Index_pproc[i1]]]=IntVar2; //clangini
            // std::cout<<"CluIndex_sort["<<SdClusAr_pproc[Index_pproc[i1]]<<"] = "
                    //  <<IntVar2<<std::endl;//clangini
            for (i2=1;i2<=NuPosSdCl;i2++) {
              if ((FlagAr[Index_pproc[i2]])&&
                  (SdClusAr_pproc[Index_pproc[i2]]==SdClusAr_pproc[Index_pproc[i1]])) {
                IntVar1=IntVar1+1;
                FrPosAr_sort[IntVar1]=FrPosAr_pproc[Index_pproc[i2]];
                SdClusAr_sort[IntVar1]=IntVar2;
                FlagAr[Index_pproc[i2]]=0;
                // std::cout << "IntVar1 " << IntVar1 << "\n";
                // std::cout << "FrPosAr_sort[IntVar1] " << FrPosAr_sort[IntVar1] << "\n";
                // std::cout << "SdClusAr_sort[IntVar1] " << SdClusAr_sort[IntVar1] << "\n";
              }
            }
          }
        }
        /* OLD OUTPUT, FULL LENGTH
           fprintf(FPaOut,"Total number of conserved second clusters : %d\n",NuSdClKe);
           fprintf(FPaOut,"Total number of conserved positions : %d\n\n",NuPosSdCl);

           fprintf(FPaOut,"Postprocessed clusters for fragment type %d (%s) :\n\n",
           CurFra,&ResN_fr[i][1]);

           fprintf(FPaOut,"                            intermolecular             ");
           fprintf(FPaOut,"electrostatic_desolvation      total_energy\n");
           fprintf(FPaOut,"                     van_der_Waals   electrostatic      ");
           fprintf(FPaOut,"receptor        fragment\n");
           fprintf(FPaOut,"                      fast    acc.    fast    acc.    ");
           fprintf(FPaOut,"fast    acc.    fast    acc.    fast    acc.\n\n");

           fprintf(FPaOut,"   Num Fr_nu  Conf    VW_f    VW_s    In_f    In_s");
           fprintf(FPaOut,"    Dr_f    Dr_s    Df_f    Df_s    To_f    To_s\n\n");

           IntVar1=1;
           fprintf(FPaOut,"%5d\n",IntVar1);
           for (j=1;j<=((NuLiEnClus<NuPosSdCl)?NuLiEnClus:NuPosSdCl);j++) {
           if (SdClusAr_sort[j]!=IntVar1) {
           fprintf(FPaOut,"\n");
           IntVar1=SdClusAr_sort[j];
           fprintf(FPaOut,"%5d\n",IntVar1);
           }
           fprintf(FPaOut,
           "%6d%6d%6d%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n",
           j,FrPosAr_sort[j],ConfArr[FrPosAr_sort[j]],
           VW_f_ro[FrPosAr_sort[j]],VW_s_ro[FrPosAr_sort[j]],In_f_ro[FrPosAr_sort[j]],
           In_s_ro[FrPosAr_sort[j]],Dr_f_ro[FrPosAr_sort[j]],Dr_s_ro[FrPosAr_sort[j]],
           Df_f_ro[FrPosAr_sort[j]],Df_s_ro[FrPosAr_sort[j]],To_f_ro[FrPosAr_sort[j]],
           To_s_ro[FrPosAr_sort[j]]);
           }
           fprintf(FPaOut,"\n");
         */
        /* NEW OUTPUT, short*/

        fprintf(FPaOut,"Total number of conserved second clusters : %d\n",NuSdClKe);
        fprintf(FPaOut,"Total number of conserved positions : %d\n\n",NuPosSdCl);

        /*fprintf(FPaOut,"Postprocessed clusters for fragment type %d (%s) :\n\n",
            CurFra,&FrFiNa_out[CurFra][1]); clangini*/
        fprintf(FPaOut,"Postprocessed clusters for fragment type %d (%s) :\n\n",
            CurFra,FragNa); /* clangini */
        /* This part of seed.out should also go in the summary table. clangini */
        fprintf(FPaOut,"                             intermolecular       ");
        fprintf(FPaOut,"electrostat_desolv.       Total\n");

        fprintf(FPaOut,"   Num    Fr_nu  Conform.  vdWaals  electrost     receptor   fragment");
        fprintf(FPaOut,"       energy\n\n");

        IntVar1=1;
        fprintf(FPaOut,"%5d\n",IntVar1);
        for (j=1;j<=((NuLiEnClus<NuPosSdCl)?NuLiEnClus:NuPosSdCl);j++) {
          if (SdClusAr_sort[j]!=IntVar1) {
            fprintf(FPaOut,"\n");
            IntVar1=SdClusAr_sort[j];
            fprintf(FPaOut,"%5d\n",IntVar1);
          }
          fprintf(FPaOut,
              "%6d%9d%6d%13.2f%11.2f%13.2f%11.2f%13.2f\n",
              j,FrPosAr_sort[j],ConfArr[FrPosAr_sort[j]],
              VW_s_ro[FrPosAr_sort[j]],In_s_ro[FrPosAr_sort[j]],Dr_s_ro[FrPosAr_sort[j]],
              Df_s_ro[FrPosAr_sort[j]],To_s_ro[FrPosAr_sort[j]]);
        }
        fprintf(FPaOut,"\n");

        /* START modification clangini 17/01/17 */
        /* Ordered postprocessed output to be written in pproc.mol2 and seed_summary.dat */
        fprintf(FPaOut,"Best postprocessed poses for fragment type %d (%s) :\n\n",
                CurFra,FragNa);
        fprintf(FPaOut,"                                   intermolecular       ");
        fprintf(FPaOut,"electrostat_desolv.       Total\n");

        fprintf(FPaOut,"   Num   Clu    Fr_nu  Conform.  vdWaals  electrost     ");
        fprintf(FPaOut,"receptor   fragment       energy\n\n");

        //IntVar1=1;
        //fprintf(FPaOut,"%5d\n",IntVar1);
        for (j=1;j<=((NuPosMem<NuPosSdCl)?NuPosMem:NuPosSdCl);j++) {
        // if NuPosMem is too large -> reduce it to NuPosSdCl
          //if (SdClusAr_sort[j]!=IntVar1) {
          //  fprintf(FPaOut,"\n");
          //  IntVar1=SdClusAr_sort[j];
          //  fprintf(FPaOut,"%5d\n",IntVar1);
          //}

          fprintf(FPaOut,
              "%6d%6d%9d%6d%13.2f%11.2f%13.2f%11.2f%13.2f\n",
              j, CluIndex_sort[SdClusAr_pproc[Index_pproc[j]]],
              FrPosAr_pproc[Index_pproc[j]],
              ConfArr[FrPosAr_pproc[Index_pproc[j]]],VW_s_ro[FrPosAr_pproc[Index_pproc[j]]],
              In_s_ro[FrPosAr_pproc[Index_pproc[j]]],Dr_s_ro[FrPosAr_pproc[Index_pproc[j]]],
              Df_s_ro[FrPosAr_pproc[Index_pproc[j]]],To_s_ro[FrPosAr_pproc[Index_pproc[j]]]);
        }
        fprintf(FPaOut,"\n");
        /* END modification clangini 17/01/17 */

        if (PosToRem.size() > 0){
          fprintf(FPaOut,"The poses ");
          for(std::vector<int>::iterator it = PosToRem.begin(); it != PosToRem.end(); ++it) {
            fprintf(FPaOut, "%d, ", FrPosAr_pproc[*it]);
          }
          fprintf(FPaOut, "will not be written to the output files as they have a total slow energy > %lf\n", FrMaEn);
        }


        /* Append lines to summary table. START clangini */
        /* Calculate HAC (Heavy Atom Count) and MW (Molecular Weight)*/
        HeAtCo = count_heavy_atom(FrAtEl_nu, FrAtNu);
        MolWei = molecular_weight(FrAtEl_nu, FrAtNu, AtWei);
        if (write_sumtab_opt[0]=='y'){
          // append to _pproc summary table
          //std::ofstream TabOutStream;
          strcpy(TabFil,"./outputs/seed_clus.dat");
          TabOutStream.open (TabFil, std::ios::out | std::ios::app); // append mode
          if(TabOutStream.is_open()){
            for (j=1;j<=((NuLiEnClus<NuPosSdCl)?NuLiEnClus:NuPosSdCl);j++) {
              if(To_s_ro[FrPosAr_sort[j]] <= FrMaEn){
                sprintf(TabLin,
                    "%-30s%8d%10d%10d%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10d%10.4f",
                    FragNa,j,SdClusAr_sort[j],FrPosAr_sort[j],To_s_ro[FrPosAr_sort[j]],In_s_ro[FrPosAr_sort[j]],
                    Dr_s_ro[FrPosAr_sort[j]],Df_s_ro[FrPosAr_sort[j]],
                    VW_s_ro[FrPosAr_sort[j]],(In_s_ro[FrPosAr_sort[j]]-FrSolvEn),
                    FrSolvEn,(To_s_ro[FrPosAr_sort[j]]/HeAtCo),
                    (VW_s_ro[FrPosAr_sort[j]]/HeAtCo),(In_s_ro[FrPosAr_sort[j]]/HeAtCo),
                    HeAtCo,MolWei);
                TabOutStream << TabLin << std::endl;
              }
            }
          } else {
            std::cerr << "Unable to write to file "<< TabFil << std::endl;
          }
          TabOutStream.close();
        }
        if (write_best_sumtab_opt[0]=='y'){
          // append to _best_pproc summary table
          strcpy(TabFil,"./outputs/seed_best.dat");
          TabOutStream.open (TabFil, std::ios::out | std::ios::app); // append mode
          if(TabOutStream.is_open()){
            for (j=1;j<=((NuPosMem<NuPosSdCl)?NuPosMem:NuPosSdCl);j++) {
              if(To_s_ro[FrPosAr_pproc[Index_pproc[j]]] < FrMaEn){
                sprintf(TabLin,
                    "%-30s%8d%10d%10d%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10d%10.4f",
                    FragNa,j,
                    CluIndex_sort[SdClusAr_pproc[Index_pproc[j]]],
                    FrPosAr_pproc[Index_pproc[j]],
                    To_s_ro[FrPosAr_pproc[Index_pproc[j]]],
                    In_s_ro[FrPosAr_pproc[Index_pproc[j]]],
                    Dr_s_ro[FrPosAr_pproc[Index_pproc[j]]],
                    Df_s_ro[FrPosAr_pproc[Index_pproc[j]]],
                    VW_s_ro[FrPosAr_pproc[Index_pproc[j]]],
                    (In_s_ro[FrPosAr_pproc[Index_pproc[j]]]-FrSolvEn),
                    FrSolvEn,(To_s_ro[FrPosAr_pproc[Index_pproc[j]]]/HeAtCo),
                    (VW_s_ro[FrPosAr_pproc[Index_pproc[j]]]/HeAtCo),
                    (In_s_ro[FrPosAr_pproc[Index_pproc[j]]]/HeAtCo),
                    HeAtCo,MolWei);
                TabOutStream << TabLin << std::endl;
              }
            }
          } else {
            std::cerr << "Unable to write to file "<< TabFil << std::endl;
          }
          TabOutStream.close();
        }
        /* Append lines summary table. END clangini */

        /* Write charmm format files */
        if (write_pproc_opt[0]=='y') {

          if(write_pproc_chm_opt[0]=='y' && NonAliHy>=3)
          {
            /*write_chm_clus_pproc(CurFra,NuPosSdCl,FrAtNu,AliHyd,In_s_ro,VW_s_ro,
                To_s_ro,FrAtTy,FrCoPo,FrPosAr_sort,
                NonAliHy,ResN_fr,FrFiNa_out); clangini to be fixed*/ /*Need to be removed. clangini*/
          }
          /*CLANGINI 2016 START*/
          /*CL write_mol2_clus_pproc(CurFra,NuPosSdCl,FrAtNu,FrAtTy,FrCoPo,FrPosAr_sort,ResN_fr,
            FrFiNa_out,FrBdNu,FrAtEl,FrPaCh,FrBdAr,FrBdTy,FragNa,
            SdClusAr_sort,To_s_ro);*/
          /*write_mol2_clus_pproc_separate(CurFra,NuPosSdCl,FrAtNu,FrAtTy,FrSyAtTy,FrCoPo,
              FrPosAr_sort,ResN_fr,FrFiNa_out,FrBdNu,FrAtEl,FrPaCh,
              FrBdAr,FrBdTy,FragNa,SdClusAr_sort,To_s_ro);*/
          sprintf(WriPat,"%s%s%s","./outputs/",FrFiNa_out,
                  "_clus.mol2\0"); // clangini
          FilePa=fopen(WriPat,"a");
          for (j=1;j<=((NuLiEnClus<NuPosSdCl)?NuLiEnClus:NuPosSdCl);j++){
#ifdef DEBUG_CL
            append_pose_to_mol2_double(FilePa,FragNa,FrAtNu,
                                FrBdNu,j,FrAtEl,FrCoPo[FrPosAr_sort[j]],
                                FrPosAr_sort[j],FrSyAtTy,FrAtTy,CurFra,FrBdAr,
                                FrBdTy,SdClusAr_sort[j],To_s_ro[FrPosAr_sort[j]],
                                FrPaCh,SubNa,AlTySp);
#else
            if(To_s_ro[FrPosAr_sort[j]] <= FrMaEn){
              append_pose_to_mol2(FilePa,FragNa,/*FragNa_map[FragNa_str],*/FrAtNu,
                                  FrBdNu,j,FrAtEl,FrCoPo,
                                  FrPosAr_sort[j],FrSyAtTy,FrAtTy,CurFra,FrBdAr,
                                  FrBdTy,SdClusAr_sort[j],To_s_ro[FrPosAr_sort[j]],
                                  FrPaCh,SubNa,AlTySp);
            }
#endif
          }
          fclose(FilePa);
          //append_to_mol2(CurFra,NuPosSdCl,FrAtNu,FrAtTy,FrSyAtTy,FrCoPo,
          //               FrPosAr_sort,SubNa,FrFiNa_out,FrBdNu,FrAtEl,FrPaCh,
          //               FrBdAr,FrBdTy,FragNa,SdClusAr_sort,To_s_ro,WriPat);
          /*CLANGINI 2016 END*/
        }

        // clangini 18/01/17 START
        if (write_best_opt[0]=='y') {// Write *best_pproc* mol2 file
          sprintf(WriPat,"%s%s%s","./outputs/",FrFiNa_out,
                  "_best.mol2\0");
          FilePa=fopen(WriPat,"a");
          for (j=1;j<=((NuPosMem<NuPosSdCl)?NuPosMem:NuPosSdCl);j++){
#ifdef DEBUG_CL
            append_pose_to_mol2_double(FilePa,FragNa,FrAtNu,
                                FrBdNu,j,FrAtEl,
                                FrCoPo[FrPosAr_pproc[Index_pproc[j]]],
                                FrPosAr_pproc[Index_pproc[j]],FrSyAtTy,
                                FrAtTy,CurFra,FrBdAr,FrBdTy,
                                CluIndex_sort[SdClusAr_pproc[Index_pproc[j]]],
                                To_s_ro[FrPosAr_pproc[Index_pproc[j]]],
                                FrPaCh,SubNa,AlTySp);
#else
            if(To_s_ro[FrPosAr_pproc[Index_pproc[j]]] < FrMaEn){
              append_pose_to_mol2(FilePa,FragNa,/*FragNa_map[FragNa_str],*/FrAtNu,
                                  FrBdNu,j,FrAtEl,FrCoPo,
                                  FrPosAr_pproc[Index_pproc[j]],FrSyAtTy,
                                  FrAtTy,CurFra,FrBdAr,FrBdTy,
                                  CluIndex_sort[SdClusAr_pproc[Index_pproc[j]]],
                                  To_s_ro[FrPosAr_pproc[Index_pproc[j]]],
                                  FrPaCh,SubNa,AlTySp);
            }
#endif
          }
          fclose(FilePa);
        }
        // clangini 18/01/17 END

        ReprSdClAr=ivector(1,NuSdClKe);
        IntVar1=1;
        ReprSdClAr[IntVar1]=FrPosAr_sort[1];
        for (j=1;j<=NuPosSdCl;j++) {
          if (SdClusAr_sort[j]!=IntVar1) {
            IntVar1=SdClusAr_sort[j];
            ReprSdClAr[IntVar1]=FrPosAr_sort[j];
          }
        }
        if (write_pproc_opt[0]=='y'&& write_pproc_chm_opt[0]=='y') {

          /*write_chm_clus_pprocr(CurFra,NuSdClKe,FrAtNu,FrAtTy,FrCoPo,ResN_fr,
              FrFiNa_out,ReprSdClAr); no longer used clangini*/ /*Needs to be removed clangini*/
        }
        if (gc_opt[0]=='y')
        {
          /*GeomCenter_FFLD(CurFra,FrFiNa_out,FrAtNu,FrAtEl_nu,NuSdClKe,ReprSdClAr,
              To_s_ro,FrCoPo,ResN_fr,FragNa,gc_reprke,gc_cutclus,
              gc_endifclus,gc_weighneg,gc_weighpos,gc_maxwrite);
              Do we want to keep support for FFLD? clangini */
        }
        free_ivector(ReprSdClAr,1,NuSdClKe);


        /* times(&timevar); */
        /* time_13=timevar.tms_utime+timevar.tms_stime; */
        gettimeofday(&time_13,NULL);

        fprintf(FPaOut,"CPU time in sec. for the postprocessing ");
        fprintf(FPaOut,"(energy with slow method) : %.2f\n",
            ((time_12.tv_sec  - time_11.tv_sec) * 1000000u +
             time_12.tv_usec - time_11.tv_usec) / 1.e6);

        /* (time_12-time_11)*0.01); */
        fprintf(FPaOut,"CPU time in sec. for the sorting of the ");
        fprintf(FPaOut,"postprocessed positions : %.2f\n\n",
            ((time_13.tv_sec  - time_12.tv_sec) * 1000000u +
             time_13.tv_usec - time_12.tv_usec) / 1.e6);

        /* (time_13-time_12)*0.01); */

        free_ivector(FrPosAr_pproc,1,NuPosSdCl);
        free_ivector(SdClusAr_pproc,1,NuPosSdCl);
        free_dvector(TotEnSdClus_pproc,1,NuPosSdCl);
        free_ivector(Index_pproc,1,NuPosSdCl);
        free_ivector(FrPosAr_sort,1,NuPosSdCl);
        free_ivector(SdClusAr_sort,1,NuPosSdCl);
        free_ivector(FlagAr,1,NuPosSdCl);

        fclose(FPaOut);

      } /* End of if (Solv_typ[0]=='p') clangini */

      free_dvector(Coulo_ro,1,SFWrNu_init);
      free_dvector(Vande_ro,1,SFWrNu_init);
      free_dvector(TotEn_ro,1,SFWrNu_init);
      free_dvector(TotEn_ro_cp,1,SFWrNu_init);

      free_ivector(FraEqu,1,SFWrNu_init);
      free_ivector(ClusIn,1,ClusNu);
      free_ivector(ClusLi,1,SFWrNu);
      free_ivector(FraEqu_sd,1,SFWrNu_init);
      free_ivector(ClusLi_sd,1,SFWrNu);
      free_ivector(ClusLi_sd_01,1,SFWrNu);
      free_ivector(ClusLi_sd_01_reduc,1,SFWrNu);
      free_ivector(ClusLi_sd_pproc,1,SFWrNu);
      free_dmatrix(FrCoor_clus,1,FrAtNu,1,3);
      free_dvector(FraSim_no,1,FrCoNu);
      free_dvector(FraSim_no_sd,1,FrCoNu);

      /* Co2 : End */
    } /*End of if (SFWrNu_init) ... clangini*/
    else { /* else: no seeded positions found clangini*/
      FPaOut=fopen(OutFil,"a");
      if (EvalEn[0]!='e') {
        fprintf(FPaOut,"WARNING No seeded positions were found for ");
        fprintf(FPaOut,"the fragment type %d\n\n",CurFra);
      }
      fclose(FPaOut);
    }

    /* times(&timevar); */
    /* time_14=timevar.tms_utime+timevar.tms_stime; */
    gettimeofday(&time_14,NULL);

    FPaOut=fopen(OutFil,"a");
    /*fprintf(FPaOut,"Total CPU time in sec. for the seeding of fragment %d (%s) : %.2f\n",
        CurFra,&FrFiNa_out[CurFra][1],
        ((time_14.tv_sec  - time_6.tv_sec) * 1000000u +
         time_14.tv_usec - time_6.tv_usec) / 1.e6); clangini */
    fprintf(FPaOut,"Total CPU time in sec. for the seeding of fragment %d (%s) : %.2f\n",
        CurFra,FragNa,
        ((time_14.tv_sec  - time_6.tv_sec) * 1000000u +
         time_14.tv_usec - time_6.tv_usec) / 1.e6); /* clangini */

    /* (time_14-time_6)*0.01); */
    fclose(FPaOut);

    for (i1=1;i1<=SFWrNu_init;i1++)
      free_dmatrix(FrCoPo[i1],1,FrAtNu,1,3);
    free_ppdvector(FrCoPo,1,ToNFrP);
    free_ivector(ConfArr,1,ToNFrP);

    free_dmatrix(FrCoor,1,FrAtNu,1,3);
    free_cmatrix(FrAtEl,1,FrAtNu,1,5);
    free_cmatrix(FrAtTy,1,FrAtNu,1,7);
    free_cmatrix(FrSyAtTy,1,FrAtNu,1,7); /*clangini 2016*/
    free_cmatrix(SubNa,1,FrAtNu,1,10); /*clangini 2016*/
    free_dvector(FrPaCh,1,FrAtNu);
    free_imatrix(FrBdAr,1,FrBdNu,1,2);
    free_cmatrix(FrBdTy,1,FrBdNu,1,4);
    free_ivector(AliHyd,1,FrAtNu);
    free_ivector(FrAtEl_nu,1,FrAtNu);
  } /* End of while(!FrInStream.eof()&&!LstFra_f)  clangini */


    /* Collect the energies of all the fragment types, sort them with respect to
       the total energy and write them in a file */
    /*
       Sort_all(FragNu,SFWrNu_ar,ResN_fr,Solv_typ);  -> REMOVED previously requires do_not_touch_ener
     */

    /* times(&timevar); */
    /* time_8=timevar.tms_utime+timevar.tms_stime; */
    gettimeofday(&time_8,NULL);

    FPaOut=fopen(OutFil,"a");
    fprintf(FPaOut,"-------------------------------------------------\n\n");
    fprintf(FPaOut,"CPU time in sec. for creating/reading Coulombic grid : %.2f\n",
        /* (time_2-time_1)*0.01); */
      ((time_2.tv_sec  - time_1.tv_sec) * 1000000u +
       time_2.tv_usec - time_1.tv_usec) / 1.e6);
    fprintf(FPaOut,"CPU time in sec. for creating/reading van der Waals grid : %.2f\n",
        ((time_4.tv_sec  - time_3.tv_sec) * 1000000u +
         time_4.tv_usec - time_3.tv_usec) / 1.e6);
    /* (time_4-time_3)*0.01); */
    fprintf(FPaOut,"CPU time in sec. up to the loop over the fragments : %.2f\n",
        ((time_5.tv_sec  - time_0.tv_sec) * 1000000u +
         time_5.tv_usec - time_0.tv_usec) / 1.e6);
    /* time_5*0.01);/xxxxx */
    fprintf(FPaOut,"Total CPU time in sec. : %.2f\n",
        ((time_8.tv_sec  - time_0.tv_sec) * 1000000u +
         time_8.tv_usec - time_0.tv_usec) / 1.e6);
    /* time_8*0.01); */
    fclose(FPaOut);


    /* clean-up */
    free_ivector(Index_ro,1,currentsize);
    free_dvector(VW_f_ro,1,currentsize);
    free_dvector(VW_s_ro,1,currentsize);
    free_dvector(In_f_ro,1,currentsize);
    free_dvector(In_s_ro,1,currentsize);
    free_dvector(Dr_f_ro,1,currentsize);
    free_dvector(Dr_s_ro,1,currentsize);
    free_dvector(Df_f_ro,1,currentsize);
    free_dvector(Df_s_ro,1,currentsize);
    free_dvector(To_f_ro,1,currentsize);
    free_dvector(To_s_ro,1,currentsize);


    /* Free the dynamically allocated memory */
    free_dvector(ReSelfVol,1,ReAtNu);
    free_c3tensor(GridMat,1,NGridx+1,1,NGridy+1,1,NGridz+1);
    free_ivector(ReApAt,1,Napol_Vect_re);
    free_dmatrix(apol_Vect_re,1,Napol_Vect_re,1,6);
    free_dvector(ReSurf_apol,1,Nsurfpt_re_apol);
    free_dvector(ReRadOut2,1,ReAtNu);
    free_dvector(ReRadOut,1,ReAtNu);
    free_dvector(ReRad2,1,ReAtNu);
    free_dvector(ReRad,1,ReAtNu);
    free_ivector(nsurf_re,1,ReAtNu);
    free_ivector(pointsrf_re,1,ReAtNu);
    free_ivector(iatsrf_re,1,NPtSphere*ReAtNu);
    free_structpointvect(surfpt_re,1,NPtSphere*ReAtNu);
    free_structpointvect(surfpt_re_deso,1,NPtSphereMax*ReAtNu);
    free_ivector(iatsrf_re_deso,1,NPtSphereMax*ReAtNu);
    free_ivector(nsurf_re_deso,1,ReAtNu);
    free_ivector(pointsrf_re_deso,1,ReAtNu);
    free_d3tensor(DeltaPrDeso,1,NGridx,1,NGridy,1,NGridz);
    free_dvector(ZGrid,1,NGridz);
    free_dvector(YGrid,1,NGridy);
    free_dvector(XGrid,1,NGridx);
    free_dvector(ReSurf_deso,1,Nsurfpt_re);
    free_imatrix(list_re,1,MaxNeigh,1,nstep_re*nstep_re*nstep_re);
    free_i3tensor(nlist_re,1,nstep_re,1,nstep_re,1,nstep_re);
    free_dvector(ReSelfVol_corrB,1,ReAtNu);

    free_ivector(BSReNu,1,BSResN);
    /*free_cmatrix(FrFiNa,1,FragNu+CorrFiNumb,1,_STRLENGTH); Now FrFiNa is a char vector */
    free_cmatrix(ReAtEl,1,ReAtNu,1,5);
    free_dmatrix(ReCoor,1,ReAtNu,1,3);
    free_cmatrix(ReAtTy,1,ReAtNu,1,7);
    free_ivector(ReResN,1,ReAtNu);
    free_dvector(RePaCh,1,ReAtNu);
    free_ivector(ReDATy,1,ReAcNu+ReDoNu);
    free_ivector(ReDAAt,1,ReAcNu+ReDoNu);
    free_ivector(ReHydN,1,ReAcNu+ReDoNu);
    free_dmatrix(ReVeCo,1,ReAcNu+ReDoNu,1,6);
    free_cmatrix(AtTyAr,1,NumbAT,1,7);
    free_ivector(AtENAr,1,NumbAT);
    free_dvector(VdWRad,1,NumbAT);
    free_dvector(VdWEne,1,NumbAT);
    free_dvector(ReVdWR,1,ReAtNu);
    free_dvector(ReVdWE,1,ReAtNu);
    if ((Solv_typ[0]=='s')||(Solv_typ[0]=='b')||(Solv_typ[0]=='p')) {
      free_i3tensor(CubFAI,1,CubNum[1],1,CubNum[2],1,CubNum[3]);
      free_i3tensor(CubLAI,1,CubNum[1],1,CubNum[2],1,CubNum[3]);
      free_ivector(CubLiA,1,ReAtNu+10);
      free_i3tensor(CubFAI_en,1,CubNum_en[1],1,CubNum_en[2],1,CubNum_en[3]);
      free_i3tensor(CubLAI_en,1,CubNum_en[1],1,CubNum_en[2],1,CubNum_en[3]);
      free_ivector(CubLiA_en,1,ReReNu+10);
      free_i3tensor(PsSphe,1,2*PsSpNC+1,1,2*PsSpNC+1,1,2*PsSpNC+1);
      free_ivector(LiChResEn,1,ReReNu);
    }
    if ((Solv_typ[0]=='f')||(Solv_typ[0]=='b')||(Solv_typ[0]=='p')) {
      free_d3tensor(CoGrRP,1,CoGPoN[1],1,CoGPoN[2],1,CoGPoN[3]);
      free_d3tensor(VWGrRP_at,1,VWGPoN[1],1,VWGPoN[2],1,VWGPoN[3]);
      free_d3tensor(VWGrRP_re,1,VWGPoN[1],1,VWGPoN[2],1,VWGPoN[3]);
    }
    free_dvector(ReVdWE_sr,1,ReAtNu);
    /*free_vector(FrMaEn,1,FragNu+CorrFiNumb); clangini CHECK*/
    //free_matrix(SimWei,1,150,1,150); clangini
    free_dmatrix(SimWei,0,150,0,150); //Want to use atom element 0 for lone pair. clangini
    /*free_vector(FrMaEn_sd,1,FragNu+CorrFiNumb); clangini CHECK*/
    free_ivector(ReAtEl_nu,1,ReAtNu);
    free_ivector(BSAtLi,1,BSAtNu);
    free_imatrix(ReBdAr,1,ReBdNu,1,2);
    free_ivector(HybReAt,1,ReAtNu);
    free_ivector(BSMeAN,1,BSMeNu);
    free_dmatrix(BSMeVE,1,BSMeNu,1,3);
    /*free_ivector(SFWrNu_ar,1,FragNu); clangini*//* -CorrFiNumb not needed !*/
    /*free_cmatrix(ResN_fr,1,FragNu,1,10); clangini*//* -CorrFiNumb not needed !*/
    /* free_cmatrix(FrFiNa_out,1,FragNu,1,_STRLENGTH); allocated statically clangini */
    free_ivector(PolVect_rec_01,1,ReAcNu+ReDoNu);
    free_dmatrix(BLAtTy,1,NumbAT,1,NumbAT);
    free_ivector(ReAtoTyp_nu,1,ReAtNu);
    free_ivector(FiAtRes,1,ReReNu);
    free_ivector(LaAtRes,1,ReReNu);
    free_ivector(AtReprRes,1,ReReNu);
    free_dvector(TotChaRes,1,ReReNu);
    /*free_cmatrix(ApPoChoi,1,FragNu+CorrFiNumb,1,2); clangini CHECK*/
    free_dvector(ReVdWR_p3,1,ReAtNu);
    free_dmatrix(distrPointBS,1,distrPointBSNumb,1,3);


    return 0;

  }
