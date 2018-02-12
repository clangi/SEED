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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nrutil.h"
#include "funct.h"

#ifdef ENABLE_MPI
  #include <mpi.h>
  #ifndef MASTERRANK
    #define MASTERRANK 0
  #endif
#endif

#ifndef _STRLENGTH
#define _STRLENGTH 500
#endif

void Solvation(int ReAtNu,double **ReCoor,double *ReVdWE_sr,double *ReVdWR,
               double *ReRad,double *ReRad2,double *ReRadOut,double *ReRadOut2,
               double *ReMaxC,double *ReMinC,double *RePaCh,double DielRe,
               double DielWa,double WaMoRa,double GrSiSo,double GrInSo,
               int NPtSphere,int *ReResN,int ReReNu,int BSResN,
               int *BSReNu,char *ReDesoAlg,char *DesoMapAcc,
               char *DesoMapFile,char *RecFilPDB,
               char *FDexe,char *FDdir,struct point *PMin,struct point *PMax,
               double **XGrid,double **YGrid,double **ZGrid,char ****GridMat,
               int *PNGridx,int *PNGridy,int *PNGridz,int *PNsurfpt_re,
               struct point *surfpt_re,int *iatsrf_re,int *nsurf_re,
               int *pointsrf_re,double ****DeltaPrDeso,
               int *PnxminBS,int *PnyminBS,int *PnzminBS,
               int *PnxmaxBS,int *PnymaxBS,int *PnzmaxBS,
               double **ReSurf_deso,int *PNsurfpt_re_deso,
               struct point *surfpt_re_deso,int *iatsrf_re_deso,
               int *nsurf_re_deso,int *pointsrf_re_deso,double SurfDens_deso,
               int ****nlist_re,int ***list_re,
               int *PMaxNeigh,int *Pnstep_re,struct point *Pcbmid_re,
               int *Pmidg_re,double *Prslop,double *Prscale,
               double **ReSurf_apol,int *PNsurfpt_re_apol,
               double ReSurfDens_apol,double Sphere_apol,double NCutapolRatio,/*int NCutapol,*/
               double ScaleDeso,double ScaleVDW,
               double ***apol_Vect_re,int **ReApAt,int *Nsurfpt_re_apol_BS,
               int *PNapol_Vect_re,double **SelfVol,
               double *PKelec,double *PKsolv,double *PUnitVol,double pi4,
               double corr_re_desoco,double corr_re_desofd,double corr_fast_deso,
	             int distrPointBSNumb,double **distrPointBS,double angle_rmin,
	             double angle_rmax,double mult_fact_rmin,double mult_fact_rmax,
	             FILE *FPaOut,double **SelfVol_corrB,char *EmpCorrB)

/*#####################################################################################
Continuum Electrostatics: It sets up grids and precalculates many quantities
that will be used in the loops over the fragments
######################################################################################*/

/*####################################################################################
int ReAtNu -------------- Tot # rec atoms
double **ReCoor ---------- Rec Coordinates
double *ReVdWE_sr -------- Rec vdW energies
double *ReVdWR ----------- Rec vdW radii
double *ReRad ----------- Rec charge radii (=vdW radii apart "enclosed" H)
double *ReRad2 ---------- (Rec charge radii)^2
double *ReRadOut -------- Rec charge radii + WaMoRa
double *ReRadOut2 ------- (Rec charge radii + WaMoRa)^2
double *ReMaxC ----------- Max (along x,y,z) of ReCoor
double *ReMinC ----------- Min (along x,y,z) of ReCoor
double *RePaCh ----------- Rec partial charges
double DielRe ----------- Dielectric constant of the rec (and of the fragments)
double DielWa ----------- Dielectric constant of the water
double WaMoRa ----------- Radius of the water molecule
double GrSiSo ----------- Size of the 3D grid used for cont. electrostatics
double GrInSo ----------- Margin left along each dimension (positive and
                          negative) of the rec for building the 3D grid box
int NPtSphere ---------- Amount of points placed over each atom to generate
                         the SAS (built in order to obtain the volume
                         enclosed by the MS)
int *ReResN ------------ # rec residue
int ReReNu ------------- Tot number of rec residues
int BSResN ------------- Tot number of rec residues in the Binding Site
int *BSReNu ------------ # rec residue in the Binding Site
char *ReDesoAlg -------- Type of algorithm to calculate rec desolvation
                         (fd = finite diff. / co = coulomb approx.)
char *DesoMapAcc ------- Access to deolvation map file
                         (w = calculate and write / n = calculate / r = read)
char *DesoMapFile ------ Name of the desolvation map file
char *RecFilPDB -------- name of the PDB file of the rec
                         (by default is ./outputs/receptor_uhbd.pdb)
char *FDexe ------------ Name of the UHBD executable file
char *FDdir ------------ Directory where to put the outputs of UHBD
struct point *PMin ----- Min coor of the grid box
struct point *PMax ----- Max coor of the grid box
double **XGrid --------- X coor of the grid points
double **YGrid --------- Y coor of the grid points
double **ZGrid --------- Z coor of the grid points
char ****GridMat ------- Matrix telling if a grid point is occupied by the
                         rec (o), empty (e), or if it belongs to the interface
                         between SAS and MS (s)
int *PNGridx ----------- Tot # of rec grid points along x
int *PNGridy ----------- Tot # of rec grid points along y
int *PNGridz ----------- Tot # of rec grid points along z

==============================================================================
                  SAS1 is a SAS built to obtain the volume enclosed by the MS
                  SAS2 is a SAS built for a fast evaluation of the desolvation
                  SAS3 is a SAS built to estimate the surface hydrophobicity
==============================================================================

int *PNsurfpt_re ------- Tot # of points over the rec SAS1
struct point *surfpt_re - Coor of points over SAS1
int *iatsrf_re --------- iatsrf_re[n] = atom number generating the nth SAS1 point
int *nsurf_re ---------- nsurf_re[n] = amount of SAS1 surface points
                         generated by rec tom n
int *pointsrf_re ------- pointsrf_re[n] = first rec SAS1 point (in the
                         list surfpt_re) that is generated by atom n
double ****DeltaPrDeso - Elec rec desolvation due to the occupation of a grid point
int *PnxminBS ---------- grid point (along x) where the BS starts
int *PnyminBS ---------- grid point (along y) where the BS starts
int *PnzminBS ---------- grid point (along z) where the BS starts
int *PnxmaxBS ---------- grid point (along x) where the BS ends
int *PnymaxBS ---------- grid point (along y) where the BS ends
int *PnzmaxBS ---------- grid point (along z) where the BS ends
double **ReSurf_deso --- ReSurf_deso[n] = Desolvation of a sphere placed
                         on the nth point of the rec SAS2
int *PNsurfpt_re_deso -- Tot # of points over the rec SAS2
struct point *surfpt_re_deso - Coor of points over SAS2
int *iatsrf_re_deso ---- iatsrf_re_deso[n] = atom number generating the
                         nth SAS2 point
int *nsurf_re_deso ----- nsurf_re_deso[n] = amount of SAS2 surface points
                         generated by atom n
int *pointsrf_re_deso -- pointsrf_re_deso[n] = first rec SAS2 point (in the list
                         surfpt_re_deso) generated by atom n
double SurfDens_deso --- Surface point density ( in A^(-2) ) to generate SAS2
int ****nlist_re ------- nlist_re[ix][iy][iz]: # of rec atoms associated with cube
                         grid point ix,iy,iz. It is a special grid created to
                         have a fast neighbour list for the fast desolvation
int ***list_re --------- list[m][n(ix,iy,iz)]: given the grid
                         point ix,iy,iz "list" tells which rec atom numbers are
                         associated with this point (m goes from
                         1 to nlist_re[ix][iy][iz] )
int *PMaxNeigh --------- Max # of atoms that can be associated to the same
                        grid point. It is intriduced to facilitate
                        memory allocation
int *Pnstep_re --------- Dimension of the grid box for the rec neighbour list. It
                         is a cubic box of even elements on each side
struct point *Pcbmid_re - Coor of the center of the grid box for the
                          neighbour list
int *Pmidg_re ---------- *Pnstep_re/2. + 1
double *Prslop --------- Parameter introduced to define the grid box for the
                         neighbour list
double *Prscale -------- 2.*(MaximumRecChargeRadius + WaMoRa) + *Prslop
double **ReSurf_apol --- ReSurf_apol[n] = vdW + Desolvation of a sphere placed
                         on the nth point of the SAS3
int *PNsurfpt_re_apol -- Tot # of points over the rec SAS3
double ReSurfDens_apol - Surface point density ( in A^(-2) ) to generate SAS3
                         for the receptor
double Sphere_apol ----- Dimension of the probe sphere used to generate SAS3
int NCutapol ----------- Amount of most hydrophobic rec points used to seed
                         apolar fragments -->> deprecated
int NCutapolRatio ------ Ratio of most hydrophobic rec points used to seed
                         apolar fragments
double ScaleDeso ------- Scaling factor for vdW interactions in the calculation
                         of ReSurf_apol[n]
double ScaleVDW -------- Scaling factor for elec desolvation in the calculation
                         of ReSurf_apol[n]
double ***apol_Vect_re -- most hydrophobic rec points used to seed
                         apolar fragments
int **ReApAt ----------- ReApAt[n] = number of the atom that has generated
                         the nth most hydrophobic point
int *Nsurfpt_re_apol_BS - Tot # of points over the rec SAS3
int *PNapol_Vect_re ---- # of most hydrophobic points used to seed the apolar
                         vectors. It is an input value, but if the value chosen
                         is smaller than the tot # of hydrophobic points, it is
                         fixed to the tot #
double **SelfVol ------- SelfVol[iat] = integral of 1/r^4 over the receptor
                         volume for atom iat
double *PKelec --------- Constant
double *PKsolv --------- Constant
double *PUnitVol ------- Volume of the grid element for cont. elec.
double pi4 ------------- 4 * greekpi
double corr_re_desoco -- Correction factor for slow rec elec desolvation
                         (coulomb approx.)
double corr_re_desofd -- Correction factor for slow rec elec desolvation
                         (fd approx.)
double corr_fast_deso -- Correction factor for fast elec desolvation
###########################################################################*/
{
  double **apol_Vect_re_loc;
  double *EffRad,*SelfEn,*ReSurf_apol_loc,*vdWSurfPt,SelfEnTot,IntEnTot,
      Rmin,Rmax,rmax,slen;
  int i,j,k,ix,iy,iz,iat,NGrid,NGridx,NGridy,NGridz,*ResFirstAtom,
      *NAtom_per_Res,*nsurf_re_apol,*pointsrf_re_apol,*iatsrf_re_apol,
      *isurf_More_apol,NPtSphereMax,nn; /*,ntime;*/
  struct point Min,Max,*surfpt_re_apol,len;
  char FD[3],CO[3];

/*##########################################################################
double **apol_Vect_re_loc --- most hydrophobic rec points used to seed
                             apolar fragments (local)
double *EffRad ------------- EffRad[n] = effective radius of rec atom n
double *SelfEn ------------- SelfEn[n] = self energy of rec atom n
double *ReSurf_apol_loc ---- ReSurf_apol_loc[n] = vdW + Desolvation of a sphere
                             placed on the nth point of the SAS3 (local)
double *vdWSurfPt ---------- vdWSurfPt[n] = vdW of a sphere placed on the
                             nth point of the SAS3
double SelfEnTot ----------- Total self energy of the isolated receptor
double IntEnTot ------------ Total interaction energy of the isolated receptor
double Rmin ---------------- Smallest rec charge radius
double Rmax ---------------- Largest rec charge radius
double rmax ---------------- Largest rec charge radius
double slen ---------------- Max dimension of the the grid box for the
                             rec neighbour list
int i,j,k,ix,iy,iz,iat ----- Multipurpose indices
int NGrid ------------------ NGridx*NGridy*NGridz
int NGridx(y,z) ------------ Tot # of grid points along x(y,z)
int *ResFirstAtom ---------- ResFirstAtom[n] = atom number of the first atom
                             belonging to residue n
int *NAtom_per_Res --------- NAtom_per_Res[n] = # of atoms of residue n
int *nsurf_re_apol --------- nsurf_re_apol[n] = amount of SAS3 surface
                             points generated by atom n
int *pointsrf_re_apol ------ pointsrf_re_apol[n] = first rec SAS3 point (in
                             the list surfpt_re_apol) generated by atom n
int *iatsrf_re_apol --------
int *iatsrf_re_apol -------- iatsrf_re_apol[n] = atom number generating the
                             nth SAS3 point
int *isurf_More_apol ------- isurf_More_apol[n] = rec SAS3 point with the nth
                             highest hydrophobicity
int NPtSphereMax ----------- Max # of points per atom given a certain
                             surface point density (ReSurfDens_apol)
int nn,ntime --------------- Multipurpose variables
struct point Min ----------- Min coor of the grid box
struct point Max ----------- Max coor of the grid box
struct point *surfpt_re_apol -  Coor of points over SAS3
struct point len ----------- ReMaxC - ReMinC
##############################################################################*/

  /* Get the # grid pts needed to cover the receptor molecule along the
     three axis */

  printf("\n\tEvaluation of grid dimension...\n");
  nn = get_Grid_Dim(ReAtNu,ReCoor,ReVdWR,WaMoRa,GrSiSo,GrInSo,&NGridx,&NGridy,
                    &NGridz,&NGrid,&Min,&Max,PUnitVol);

  *PNGridx = NGridx;
  *PNGridy = NGridy;
  *PNGridz = NGridz;
  PMin->x = Min.x;
  PMin->y = Min.y;
  PMin->z = Min.z;
  PMax->x = Max.x;
  PMax->y = Max.y;
  PMax->z = Max.z;

  /* Get the cartesian coords of the grid pts (after having
     allocated space) */

  printf("\n\tCartesian coords of the grid...\n");
  *XGrid=dvector(1,NGridx);
  *YGrid=dvector(1,NGridy);
  *ZGrid=dvector(1,NGridz);
  nn = Coor_Grid_Pts(GrSiSo,NGridx,NGridy,NGridz,Min,
                    *XGrid,*YGrid,*ZGrid);
/*
  printf("\n\tXGrid(1) %f YGrid(1) %f ZGrid(1) %f\n",(*XGrid)[1],
         (*YGrid)[1],(*ZGrid)[1]);
  printf("\n\tXGrid(NGridx) %f YGrid(NGridy) %f ZGrid(NGridz) %f\n",
         (*XGrid)[NGridx],(*YGrid)[NGridy],(*ZGrid)[NGridz]);
  printf("\n\tNGridx %d NGridy %d NGridz %d\n",
         NGridx,NGridy,NGridz);
  printf("\n\tRec. Min x %f Rec. Max x %f Extent %f\n",Min.x+GrInSo,
                             Max.x-GrInSo,Max.x-GrInSo-Min.x-GrInSo);
  printf("\n\tRec. Min y %f Rec. Max y %f Extent %f\n",Min.y+GrInSo,
                             Max.y-GrInSo,Max.y-GrInSo-Min.y-GrInSo);
  printf("\n\tRec. Min z %f Rec. Max z %f Extent %f\n",Min.z+GrInSo,
                             Max.z-GrInSo,Max.z-GrInSo-Min.z-GrInSo);
*/

  /* Create array ReRad (array of charge radii: they are the vdW
     radii unless an atom (usually H) is completely inclosed into
     another atom: in this case its radius is enlarged) and also
     some other useful arrays: ReRad2 (ReRad^2), ReRadOut (ReRad+WaMoRa),
     ReRadOut2 ([ReRad+WaMoRa]^2). Before, allocate space... */

  printf("\n\tCharge radii...\n");
  nn = get_Ch_Rad(ReAtNu,ReCoor,ReVdWR,WaMoRa,ReRad,ReRad2,
                  ReRadOut,ReRadOut2,&Rmin,&Rmax);
  /* GridMat is the matrix that tells us whether a gr pt is occupied */

  *GridMat = c3tensor(1,NGridx+1,1,NGridy+1,1,NGridz+1);
/* Initialize GridMat */
  for (i=1;i<=NGridx+1;i++)
    for (j=1;j<=NGridy+1;j++)
      for (k=1;k<=NGridz+1;k++)
        (*GridMat)[i][j][k] = 'e';

/* Make the map (GridMat) of the 3D grid points occupied by the volume
   enclosed by the rec SAS  */
  printf("\n\tMap of volume enclosed by the SAS...\n");
  nn = SAS_Volume(ReAtNu,ReCoor,ReRadOut,ReRadOut2,Min,GrSiSo,
                 1,1,1,NGridx,NGridy,NGridz,*GridMat);

/* Place points on the surface of the receptor to describe its SAS */
/* FIRST SAS (SAS1): needed to define the volume enclosed by the MS of the receptor */
  printf("\n\tGeneration of the SAS of the receptor...\n");
  nn = Surf_Grid(ReAtNu,ReCoor,ReMaxC,ReMinC,ReRad,WaMoRa,NPtSphere,
                PNsurfpt_re,surfpt_re,iatsrf_re,nsurf_re,pointsrf_re);

/* Place a sphere of radius WaMoRa on every surface grid point and
   mark as empty all the volume grid points falling inside the sphere */
  printf("\n\tGeneration of volume enclosed by the MS...\n");
  *SelfVol=dvector(1,ReAtNu);
  *SelfVol_corrB=dvector(1,ReAtNu);
  for (iat=1;iat<=ReAtNu;iat++)
  {
    (*SelfVol)[iat] = 0.;
    (*SelfVol_corrB)[iat] = 0.;
  }
  nn = Excl_Grid(ReAtNu,ReCoor,Min,ReRadOut,ReRadOut2,WaMoRa,GrSiSo,
                 1,1,1,NGridx,NGridy,NGridz,*PUnitVol,*GridMat,
                 *PNsurfpt_re,surfpt_re,*SelfVol,*SelfVol_corrB,EmpCorrB);

/* Calculation of the integral 1 over r4 (GB formula) for each atom */
  printf("\n\tCalculation of the integral 1 over r4 for each atom\n");
  printf("\tover the volume enclosed by the MS...\n");
  nn = Get_Self_Vol(ReAtNu,ReCoor,ReRadOut2,GrSiSo,*XGrid,*YGrid,*ZGrid,1,1,1,
                    NGridx,NGridy,NGridz,*PUnitVol,*GridMat,*SelfVol,
		                *SelfVol_corrB,EmpCorrB);

/* Evaluation of the self energy */
  EffRad=dvector(1,ReAtNu);
  SelfEn=dvector(1,ReAtNu);
  *PKsolv = 332.0716 * ( 1./DielWa - 1./DielRe );
  *PKelec = 332.0716 / DielRe;
  SelfEnTot = 0.;
  for (iat=1;iat<=ReAtNu;iat++){


    if (EmpCorrB[0]!='y')
      EffRad[iat] = 1. / ( 1./ReRadOut[iat] - (*SelfVol)[iat]/pi4 ); // eq. (14) from [Scarsi et al. 1997]. clangini
    else
    {

	     EffRad[iat] = 1./( (-1.*(1./ReRadOut[iat] - (*SelfVol)[iat]/pi4))
			    + 3.0*sqrt( (1./(2.*ReRadOut[iat]*ReRadOut[iat])) -
				       ((*SelfVol_corrB)[iat]/pi4) ) )
	        + 0.215;


	/*
	  Dey exception handling :
	  in rare cases the expression :
	  (1./(2.*ReRadOut[iat]*ReRadOut[iat])) - ((*SelfVol_corrB)[iat]/pi4)
	  can become < 0 -> the sqrt() function cannot be evaluated, which leads
	  to "nan" values or the effective born radius is smaller than 0

	*/
	     if(EffRad[iat]<=0 || isnan(EffRad[iat])) {
#ifndef NOWARN
	       fprintf(FPaOut,"WARNING could not calculate effective born radius for atom %d, using standard approach\n",iat);
#endif
	       EffRad[iat] = 1. / ( 1./ReRadOut[iat] - (*SelfVol)[iat]/pi4 );
	     }



    }


    SelfEn[iat] = *PKsolv * RePaCh[iat] * RePaCh[iat] / (2. * EffRad[iat]); // eq. (12) in [Scarsi et al. 1997]. clangini
    SelfEnTot += SelfEn[iat];

  }




  printf("\n\tSelf Energy %f\n",SelfEnTot);

/* Evaluation of the interaction energy */
  printf("\n\tElectrostatic interactions...\n");
  nn = GB_int(ReAtNu,ReCoor,RePaCh,EffRad,*PKsolv,&IntEnTot);


  printf("\n\tSolvation Interaction Energy %f\n",IntEnTot);

/* Get some infos about the dimension and the location of the binding site...\n") */
  printf("\n\tGet Binding Site Infos...\n");
  NAtom_per_Res=ivector(1,ReReNu);
  ResFirstAtom=ivector(1,ReReNu);
  nn = Get_BS_Info(ReAtNu,ReCoor,ReRadOut,Min,GrSiSo,GrInSo,ReResN,
                   BSResN,BSReNu,NGridx,NGridy,NGridz,PnxminBS,
                   PnyminBS,PnzminBS,PnxmaxBS,PnymaxBS,PnzmaxBS,
                   ResFirstAtom,NAtom_per_Res);

/* Evaluation of the desolvation of the receptor due to the occupation
   of a grid cube in the binding site */
  printf("\n\tElectric Displacement and receptor desolvation...\n");
  *DeltaPrDeso = d3tensor(1,NGridx,1,NGridy,1,NGridz);
  for (ix=1;ix<=NGridx;ix++)
    for (iy=1;iy<=NGridy;iy++)
      for (iz=1;iz<=NGridz;iz++)
        (*DeltaPrDeso)[ix][iy][iz] = 0.0000000000000000;
  sprintf(FD,"%s","fd"); /* prev "fd\0" */
  sprintf(CO,"%s","co"); /* prev "co\0" */
  if (DesoMapAcc[0] != 'r') {
    if ( strcmp(ReDesoAlg,FD) == 0 )
/* Evaluation of the desolvation in the BS by finite difference calculations */
      nn = Calc_D_uhbd(ReAtNu,ReRad,Min,Max,RePaCh,WaMoRa,NPtSphere,
                       DielRe,DielWa,GrSiSo,GrInSo,pi4,NGridx,NGridy,NGridz,
                       *PnxminBS,*PnyminBS,*PnzminBS,*PnxmaxBS,*PnymaxBS,
                       *PnzmaxBS,*XGrid,*YGrid,*ZGrid,*GridMat,*PUnitVol,
                       corr_re_desofd,RecFilPDB,FDexe,FDdir,DesoMapAcc,
                       DesoMapFile,*DeltaPrDeso);
    else if ( strcmp(ReDesoAlg,CO) == 0 ) {
/* Evaluation of the desolvation in the BS by Coulmomb approximation */
      nn = Calc_D_Coul(ReAtNu,ReCoor,ReRad2,Min,Max,RePaCh,DielRe,DielWa,
                       GrSiSo,pi4,*PnxminBS,*PnyminBS,*PnzminBS,
                       *PnxmaxBS,*PnymaxBS,*PnzmaxBS,*XGrid,*YGrid,*ZGrid,
                       *GridMat,*PUnitVol,corr_re_desoco,*DeltaPrDeso);
    }
  }
  if (DesoMapAcc[0] != 'n')
/* If requested, read or write the map for the desolvation of the receptor */
    nn = Read_Write_Desol_Map(WaMoRa,NPtSphere,DielRe,DielWa,GrSiSo,GrInSo,
                              NGridx,NGridy,NGridz,*PnxminBS,*PnyminBS,
                              *PnzminBS,*PnxmaxBS,*PnymaxBS,*PnzmaxBS,
                              *XGrid,*YGrid,*ZGrid,*GridMat,*PUnitVol,
                              ReDesoAlg,DesoMapAcc,DesoMapFile,*DeltaPrDeso);

/* Fast method to evaluate receptor desolvation */
  printf("\n\tFast method to evaluate receptor desolvation\n");
/* get some values for the cubing grid that is used for the SAS nighbour list
   All this is done to allocate the memory for nlist_re and list_re*/
  len.x = ReMaxC[1] - ReMinC[1];
  len.y = ReMaxC[2] - ReMinC[2];
  len.z = ReMaxC[3] - ReMinC[3];
  slen = (len.x>len.y) ? len.x : len.y;
  slen = (slen>len.z) ? slen : len.z;
  rmax = MaxDVector(ReRad, 1, ReAtNu);
  *Prslop = .2;
  *Prscale = 2.*(rmax + WaMoRa) + *Prslop;
  *Pnstep_re = (int) (( slen + 2.*(rmax + WaMoRa) ) / (*Prscale) + 0.5) + 2;
  *Pnstep_re = ( ((*Pnstep_re >> 1) << 1) == *Pnstep_re) ?
               *Pnstep_re : *Pnstep_re + 1;
  *nlist_re=i3tensor(1,*Pnstep_re,1,*Pnstep_re,1,*Pnstep_re);
  *PMaxNeigh = 500;
  *list_re=imatrix(1,*PMaxNeigh,1,(*Pnstep_re)*(*Pnstep_re)*(*Pnstep_re));

/* Place point over the receptor SAS with uniform density and get also the neighbour
   list (nlist_re,list_re) over a grid */
/* SECOND SAS (SAS2): needed to estimate the approximate receptor desolvation
   at the surface */
  nn = Surf_Grid_unif_dens_neighlist(ReAtNu,ReCoor,ReMinC,ReRad,
                                     ReRadOut2,WaMoRa,SurfDens_deso,pi4,
                                     PNsurfpt_re_deso,surfpt_re_deso,
                                     iatsrf_re_deso,nsurf_re_deso,
                                     pointsrf_re_deso,len,slen,
                                     *nlist_re,*list_re,*Pnstep_re,Pcbmid_re,
                                     Pmidg_re,*Prscale);
  *ReSurf_deso = dvector(1,*PNsurfpt_re_deso);
  nn = Fast_Desol_Surf(Min,WaMoRa,GrSiSo,1,1,1,NGridx,NGridy,NGridz,
                      *PNsurfpt_re_deso,surfpt_re_deso,*DeltaPrDeso,
                      corr_fast_deso,*ReSurf_deso);

/* THIRD SAS (SAS3): needed to identify the more hydrophobic zones
   on the receptor surface */
  printf("\n\tGeneration of the SAS for the Desolvation Surface\n");
  NPtSphereMax = (int) (ReSurfDens_apol * pi4 * (Rmax+WaMoRa) * (Rmax+WaMoRa));
  surfpt_re_apol=structpointvect(1,NPtSphereMax*ReAtNu);
  iatsrf_re_apol=ivector(1,NPtSphereMax*ReAtNu);
  nsurf_re_apol=ivector(1,ReAtNu);
  pointsrf_re_apol=ivector(1,ReAtNu);
  for (iat=1;iat<=ReAtNu;iat++)
    pointsrf_re_apol[iat]=0;
  for (iat=1;iat<=ReAtNu;iat++)
    nsurf_re_apol[iat] = 0;
/* Generation of a SAS around the receptor */
  nn = Surf_Grid_unif_dens(ReAtNu,ReCoor,ReMaxC,ReMinC,ReRad,ReRadOut2,
                           Sphere_apol,ReSurfDens_apol,pi4,PNsurfpt_re_apol,
                           surfpt_re_apol,iatsrf_re_apol,nsurf_re_apol,
                           pointsrf_re_apol);

/* Evaluation of  the vdw interaction of the probe sphere with the receptor */
  printf("\n\tvdW Interaction Surface\n");
  vdWSurfPt = dvector(1,*PNsurfpt_re_apol);
  for (i=1;i<=*PNsurfpt_re_apol;i++)
    (vdWSurfPt)[i] = 0.;
  nn = vdW_Surf(ReAtNu,ReCoor,ReVdWE_sr,ReVdWR,Sphere_apol,BSResN,
                BSReNu,NAtom_per_Res,ResFirstAtom,surfpt_re_apol,
                nsurf_re_apol,pointsrf_re_apol,vdWSurfPt);

/* Evaluation of  the desolvation of the receptor operated by the removal
   of a sphere placed with the center over a SAS3 surface grid point */
  printf("\n\tDesolvation Surface\n");
  ReSurf_apol_loc = dvector(1,*PNsurfpt_re_apol);
  *ReSurf_apol = ReSurf_apol_loc;
  for (i=1;i<=*PNsurfpt_re_apol;i++)
/*    ReSurf_apol_loc[i] = 0.; */
    ReSurf_apol_loc[i] = ScaleVDW * vdWSurfPt[i];
  isurf_More_apol=ivector(1,NPtSphereMax*ReAtNu);
  nn = Desol_Surf(ReAtNu,Min,Sphere_apol,GrSiSo,NGridx,NGridy,NGridz,BSResN,
                  BSReNu,NAtom_per_Res,ResFirstAtom,surfpt_re_apol,
                  nsurf_re_apol,pointsrf_re_apol,*DeltaPrDeso,
                  ReSurf_apol_loc,/*NCutapol,*/NCutapolRatio,ScaleDeso,
                  isurf_More_apol,PNapol_Vect_re,Nsurfpt_re_apol_BS,
		              iatsrf_re_apol,distrPointBSNumb,distrPointBS,angle_rmin,
		              angle_rmax,mult_fact_rmin,mult_fact_rmax,ReCoor,
		              FPaOut);

/* Creation of vectors for the placement of apolar fragments */
  printf("\n\tCreate vectors for the placement of apolar fragments\n");
  apol_Vect_re_loc = dmatrix(1,*PNapol_Vect_re,1,6);
  *apol_Vect_re = apol_Vect_re_loc;
  *ReApAt=ivector(1,*PNapol_Vect_re);
  for (i=1;i<=*PNapol_Vect_re;i++) {
    for (j=1;j<=3;j++)
      apol_Vect_re_loc[i][j]=ReCoor[iatsrf_re_apol[isurf_More_apol[i]]][j];
    apol_Vect_re_loc[i][4] = (double) surfpt_re_apol[isurf_More_apol[i]].x;
    apol_Vect_re_loc[i][5] = (double) surfpt_re_apol[isurf_More_apol[i]].y;
    apol_Vect_re_loc[i][6] = (double) surfpt_re_apol[isurf_More_apol[i]].z;
    (*ReApAt)[i]=iatsrf_re_apol[isurf_More_apol[i]];
  }

/* Free memory */
  free_dvector(EffRad,1,ReAtNu);
  free_ivector(isurf_More_apol,1,NPtSphere*ReAtNu);
  free_ivector(nsurf_re_apol,1,ReAtNu);
  free_ivector(pointsrf_re_apol,1,ReAtNu);
  free_ivector(iatsrf_re_apol,1,NPtSphereMax*ReAtNu);
  free_structpointvect(surfpt_re_apol,1,NPtSphereMax*ReAtNu);
  free_ivector(ResFirstAtom,1,ReReNu);
  free_ivector(NAtom_per_Res,1,ReReNu);
  free_dvector(SelfEn,1,ReAtNu);

  free_dvector(vdWSurfPt,1,*PNsurfpt_re_apol); /* prev. mem. leak*/

  printf("\n\tEnd Solv\n");
}

int get_Grid_Dim(int ReAtNu,double **ReCoor,double *ReVdWR,double WaMoRa,
                 double GrSiSo,double GrInSo,int *PNGridx,int *PNGridy,
                 int *PNGridz,int *PNGrid,struct point *PMin,
                 struct point *PMax,double *PUnitVol)
/*########################################################################
Prepare the grid over the receptor: find max and min over the 3 coords,
then calculate the # grid pts that are necessary (NGridx,NGridy,NGridz),
then calculate the volume occupied by a grid pt (UnitVol)
##########################################################################*/

/*#####################################################################
int ReAtNu -------------- Tot # rec (frag) atoms
double **ReCoor ---------- Rec (frag) Coordinates
double *ReVdWR ----------- Rec (frag) vdW radii
double WaMoRa ----------- Radius of the water molecule
double GrSiSo ----------- Size of the 3D grid used for cont. electrostatics
double GrInSo ----------- Margin left along each dimension (positive and
                          negative) of the rec (frag) for building the 3D grid box
int *PNGridx ------------ Tot # of grid points along x
int *PNGridy ------------ Tot # of grid points along y
int *PNGridz ------------ Tot # of grid points along z
int *PNGrid ------------- NGridx*NGridy*NGridz
struct point *PMin ------ Min coor of the rec (frag) grid box
struct point *PMax ------ Max coor of the rec (frag) grid box
double *PUnitVol -------- Volume of the rec (frag) grid element for cont. elec.
#####################################################################*/

{
  double *vect1,*vect2;
  int i;

  vect1=dvector(1,ReAtNu);
  vect2=dvector(1,ReAtNu);
  for (i=1; i<=3; i++) {
    getColumnFrom2DArray(ReCoor, i, 1, ReAtNu, vect1);
    SumVectors(vect1, ReVdWR, 1, ReAtNu, vect2);
    if (i == 1) {
      PMax->x = MaxDVector(vect2, 1, ReAtNu);
      PMax->x += (WaMoRa + GrInSo);
    }
    else if (i == 2) {
      PMax->y=MaxDVector(vect2, 1, ReAtNu);
      PMax->y += (WaMoRa + GrInSo);
    }
    else if (i == 3) {
      PMax->z=MaxDVector(vect2, 1, ReAtNu);
      PMax->z += (WaMoRa + GrInSo);
    }
    SubstractVectors(vect1, ReVdWR, 1, ReAtNu, vect2);
    if (i == 1) {
      PMin->x=MinVector(vect2, 1, ReAtNu);
      PMin->x -= (WaMoRa + GrInSo);
    }
    else if (i == 2) {
      PMin->y=MinVector(vect2, 1, ReAtNu);
      PMin->y -= (WaMoRa + GrInSo);
    }
    else if (i == 3) {
      PMin->z=MinVector(vect2, 1, ReAtNu);
      PMin->z -= (WaMoRa + GrInSo);
    }
  }
  free_dvector(vect2,1,ReAtNu);
  free_dvector(vect1,1,ReAtNu);
  *PNGridx = ( (PMax->x-PMin->x)/GrSiSo + 1);
  //clangini debug:
  //std::cout << "PMin->x  " << PMin->x << std::endl;
  //std::cout << "PMax->x  " << PMax->x << std::endl;
  //std::cout << "*PNGridx  " << *PNGridx<< std::endl;
  //clangini debug end
  *PNGridy = ( (PMax->y-PMin->y)/GrSiSo + 1);
  *PNGridz = ( (PMax->z-PMin->z)/GrSiSo + 1);
  *PNGrid = (*PNGridx) * (*PNGridy) * (*PNGridz);
  *PUnitVol = GrSiSo*GrSiSo*GrSiSo;
  if (i == 4 )
    return 1;
  else
    return 0;
}

int Coor_Grid_Pts(double GrSiSo,int NGridx,int NGridy,int NGridz,
                  struct point Min,double *XGrid,double *YGrid,double *ZGrid)
/*##########################################
Get the coor of the grid points
###########################################*/

/*##########################################
double GrSiSo ----------- Size of the 3D grid used for cont. electrostatics
int NGridx -------------- Tot # of grid points along x
int NGridy -------------- Tot # of grid points along y
int NGridz -------------- Tot # of grid points along z
struct point Min -------- Min coor of the rec (frag) grid box
double *XGrid ----------- X coor of the grid points
double *YGrid ----------- Y coor of the grid points
double *ZGriz ----------- Z coor of the grid points
###########################################*/
{
  int i;

  for (i=1;i<=NGridx;i++)
    XGrid[i] = GrSiSo*(i-0.5) + Min.x;
  for (i=1;i<=NGridy;i++)
    YGrid[i] = GrSiSo*(i-0.5) + Min.y;
  for (i=1;i<=NGridz;i++)
    ZGrid[i] = GrSiSo*(i-0.5) + Min.z;
  if (i == NGridz+1 )
    return 1;
  else
    return 0;
}

int get_Ch_Rad(int ReAtNu,double **ReCoor,double *ReVdWR,double WaMoRa,
               double *ReRad,double *ReRad2,double *ReRadOut,
               double *ReRadOut2,double *PRmin,double *PRmax)
/*##########################################
Calculate the charge radii for the receptor
###########################################*/

/*##########################################
int ReAtNu -------------- Tot # rec atoms
double **ReCoor ---------- Rec Coordinates
double *ReVdWR ----------- Rec vdW radii
double WaMoRa ----------- Radius of the water molecule
double *ReRad ----------- Rec charge radii (=vdW radii apart "enclosed" H)
double *ReRad2 ---------- (Rec charge radii)^2
double *ReRadOut -------- Rec charge radii + WaMoRa
double *ReRadOut2 ------- (Rec charge radii + WaMoRa)^2
double *PRmin ----------- Smallest rec charge radius
double *PRmax ----------- Largest rec charge radius
###########################################*/
{
  int i,j;
  double d,d2;

  for (i=1;i<=ReAtNu;i++)
    ReRad[i] =  ReVdWR[i];

  *PRmin = 9999.;
  *PRmax = -9999.;
  for (i=1;i<=ReAtNu;i++) {
    for (j=i+1;j<=ReAtNu;j++) {
      d2 = (ReCoor[i][1]-ReCoor[j][1])*(ReCoor[i][1]-ReCoor[j][1]) +
           (ReCoor[i][2]-ReCoor[j][2])*(ReCoor[i][2]-ReCoor[j][2]) +
           (ReCoor[i][3]-ReCoor[j][3])*(ReCoor[i][3]-ReCoor[j][3]);
      d = sqrtf (d2);
      if ( d + ReVdWR[i] < ReVdWR[j] )
        ReRad[i] = ReVdWR[j] - d;
      else if ( d + ReVdWR[j] < ReVdWR[i] )
        ReRad[j] = ReVdWR[i] - d;
    }
    *PRmin = (*PRmin < ReRad[i]) ? *PRmin : ReRad[i];
    *PRmax = (*PRmax > ReRad[i]) ? *PRmax : ReRad[i];
    ReRad2[i] = ReRad[i]*ReRad[i];
    ReRadOut[i] = ReRad[i] + WaMoRa;
    ReRadOut2[i] = ReRadOut[i]*ReRadOut[i];
  }
  if (i == ReAtNu+1 )
    return 1;
  else
    return 0;
}

int get_Ch_Rad_Fr(int ReAtNu,double **ReCoor,double *ReVdWR,double WaMoRa,
                  double *ReRad,double *ReRad2,double *ReRadOut,
                  double *ReRadOut2,double **dist2,double *PRmin,double *PRmax)
/*##########################################
Calculate the charge radii for the fragment (the diff. with get_Ch_Rad
is that the interatomic distances are transmitted out).
###########################################*/

// The only radii different from vdW are those of atoms fully enclosed in
//the vdW sphere of the another atom

/*##########################################
int ReAtNu -------------- Tot # atoms
double **ReCoor ---------- Coordinates
double *ReVdWR ----------- vdW radii
double WaMoRa ----------- Radius of the water molecule
double *ReRad ----------- Charge radii (=vdW radii apart "enclosed" H)
double *ReRad2 ---------- (Charge radii)^2
double *ReRadOut -------- Charge radii + WaMoRa
double *ReRadOut2 ------- (Charge radii + WaMoRa)^2
double **dist2 ----------- squared interatomic distances
double *PRmin ----------- Smallest charge radius
double *PRmax ----------- Largest charge radius
###########################################*/
{
  int i,j;
  double d;

  for (i=1;i<=ReAtNu;i++)
    ReRad[i] =  ReVdWR[i];

  *PRmin = 9999.;
  *PRmax = -9999.;
  for (i=1;i<=ReAtNu;i++) {
    for (j=i+1;j<=ReAtNu;j++) {
      dist2[i][j]=(ReCoor[i][1]-ReCoor[j][1])*(ReCoor[i][1]-ReCoor[j][1]) +
                  (ReCoor[i][2]-ReCoor[j][2])*(ReCoor[i][2]-ReCoor[j][2]) +
                  (ReCoor[i][3]-ReCoor[j][3])*(ReCoor[i][3]-ReCoor[j][3]);
      d = sqrtf (dist2[i][j]);
      if ( d + ReVdWR[i] < ReVdWR[j] ){
        ReRad[i] = ReVdWR[j] - d;
        std::cout << "Correction for atom "<< i <<std::endl;
      }
      else if ( d + ReVdWR[j] < ReVdWR[i] ){
        ReRad[j] = ReVdWR[i] - d;
        std::cout << "Correction for atom "<< j <<std::endl;
      }
    }
    *PRmin = (*PRmin < ReRad[i]) ? *PRmin : ReRad[i];
    *PRmax = (*PRmax > ReRad[i]) ? *PRmax : ReRad[i];
    ReRad2[i] = ReRad[i]*ReRad[i];
    ReRadOut[i] = ReRad[i] + WaMoRa;
    ReRadOut2[i] = ReRadOut[i]*ReRadOut[i];

  }

  // clangini debug:
  /*std::cout << "ReRad list: " << std::endl;
  for (i=1;i<=ReAtNu;i++){
    std::cout << "ReRad[" << i << "] = " << ReRad[i] << std::endl;
  }*/
  // clangini debug end

  if (i == ReAtNu+1 )
    return 1;
  else
    return 0;
}

int SAS_Volume(int ReAtNu,double **ReCoor,double *ReRadOut,
               double *ReRadOut2,struct point Min,
               double GrSiSo,int NStartGridx,int NStartGridy,
               int NStartGridz,int NGridx,int NGridy,int NGridz,
               char ***GridMat)
/*##########################################
Get the volume enclosed by the SAS
###########################################*/

/*##########################################
int ReAtNu -------------- Tot # rec (frag) atoms
double **ReCoor ---------- Rec (frag) coordinates
double *ReRadOut -------- Rec (frag) charge radii + WaMoRa
double *ReRadOut2 ------- (Rec (frag) charge radii + WaMoRa)^2
struct point Min -------- Min coor of the grid box
double GrSiSo ----------- Size of the 3D grid used for cont. electrostatics
int NStartGridx --------- First grid point along x (it is 1)
int NStartGridy --------- First grid point along y (it is 1)
int NStartGridz --------- First grid point along z (it is 1)
int  NGridx ------------- Tot # of grid points along x
int  NGridy ------------- Tot # of grid points along y
int  NGridz ------------- Tot # of grid points along z
char ****GridMat -------- Matrix telling if a grid point is occupied (o),
                          empty (e), or if it belongs to the interface
                          between SAS and MS (s)
###########################################*/
{
  int iat,ix,iy,iz,ixmin,iymin,izmin,ixmax,iymax,izmax;
  double xtemp,x2temp,ytemp,xy2temp,ztemp,r2;

/* Calculate grid points occupied by volume enclosed in SAS */
  for (iat=1;iat<=ReAtNu;iat++) { // for each atom, mark the occupied locations. clangini
    ixmin = ( (ReCoor[iat][1] - ReRadOut[iat] - Min.x) / GrSiSo + 1 );
    iymin = ( (ReCoor[iat][2] - ReRadOut[iat] - Min.y) / GrSiSo + 1 );
    izmin = ( (ReCoor[iat][3] - ReRadOut[iat] - Min.z) / GrSiSo + 1 );
    ixmax = ( (ReCoor[iat][1] + ReRadOut[iat] - Min.x) / GrSiSo + 1 );
    iymax = ( (ReCoor[iat][2] + ReRadOut[iat] - Min.y) / GrSiSo + 1 );
    izmax = ( (ReCoor[iat][3] + ReRadOut[iat] - Min.z) / GrSiSo + 1 );
    ixmin = (ixmin > NStartGridx) ? ixmin : NStartGridx;
    iymin = (iymin > NStartGridy) ? iymin : NStartGridy;
    izmin = (izmin > NStartGridz) ? izmin : NStartGridz;
    ixmax = (ixmax < NGridx) ? ixmax : NGridx;
    iymax = (iymax < NGridy) ? iymax : NGridy;
    izmax = (izmax < NGridz) ? izmax : NGridz;

    xtemp = GrSiSo * (ixmin - 1.5) + Min.x - ReCoor[iat][1];
    for (ix=ixmin;ix<=ixmax;ix++) {
      xtemp += GrSiSo;
      x2temp = xtemp*xtemp;
      ytemp = GrSiSo * (iymin - 1.5) + Min.y - ReCoor[iat][2];
      for (iy=iymin;iy<=iymax;iy++) {
        ytemp += GrSiSo;
        xy2temp = ytemp*ytemp + x2temp;
        if (ReRadOut2[iat] > xy2temp) {
          ztemp = GrSiSo * (izmin - 1.5) + Min.z - ReCoor[iat][3];
          for (iz=izmin;iz<=izmax;iz++) {
            ztemp += GrSiSo;
            if (GridMat[ix][iy][iz] != 'o') {
              r2 = ztemp * ztemp + xy2temp;
              if (ReRadOut2[iat] > r2)
                GridMat[ix][iy][iz] = 'o';
            }
          }
        }
      }
    }
  }

  if (iat == ReAtNu+1 )
    return 1;
  else
    return 0;
}

int Surf_Grid_unif_dens(int ReAtNu,double **ReCoor,double *ReMaxC,double *ReMinC,
                        double *ReRad,double *ReRadOut2,double WaMoRa,
                        double ReSurfDens_apol,double pi4,int *PNsurfpt,
                        struct point *surfpt,int *iatsrf,int *nsurf,
                        int *pointsrf)
/*##########################################
Place points over the SAS assuming a uniform surface density
###########################################*/

/*##########################################
int ReAtNu -------------- Tot # atoms
double **ReCoor ---------- Coordinates
double *ReMaxC ----------- Max (along x,y,z) of ReCoor
double *ReMinC ----------- Min (along x,y,z) of ReCoor
double *ReRad ----------- Charge radii (=vdW radii apart "enclosed" H)
double *ReRadOut2 ------- (Charge radii + WaMoRa)^2
double WaMoRa ----------- Radius of the water molecule
double ReSurfDens_apol -- Surface point density ( in A^(-2) ) to generate SAS
                         for the receptor
double pi4 -------------- 4 * greekpi
int *PNsurfpt ----------- Tot # of points over SAS
struct point *surfpt ---- Coor of points over SAS
int *iatsrf ------------- iatsrf[n] = atom number generating the nth SAS point
int *nsurf -------------- nsurf_re[n] = amount of SAS surface points
                          generated by atom n
int *pointsrf ----------- pointsrf_re[n] = first rec SAS point (in the
                          list surfpt) that is generated by atom n
###########################################*/
{
    double rmax,slen,rslop,rscale,cbscal ,*cbrdp,*cbrdp2,d2;
  int iat,ix,iy,iz,nstep,midg,nnnmax,NPtSphere_old,isph,*NPtSphereiat,
      nx,ny,nz,inayb,listmp,***nlist,**list,first;
  struct point len,cbmid,*ccrd,*SphCoor,sphpt,d,midgpt;

/*#########################################################################
double rmax ---------- Largest charge radius
double slen ---------- Max dimension of the the grid box for the
                      neighbour list
double rslop --------- Parameter introduced to define the grid box for the
                      neighbour list
double rscale -------- 2.*(MaximumChargeRadius + WaMoRa) + rslop
double cbscal -------- 1./rscale
double *cbrdp -------- cbrdp[iat] = cbscal * ( ReRad[iat] + WaMoRa )
double *cbrdp2 ------- cbrdp^2
double d2
int iat,ix,iy,iz ----- multipurpose indices
int nstep ------------ Dimension of the grid box for the neighbour list. It
                       is a cubic box of even elements on each side
int midg ------------- nstep/2. + 1
int nnnmax ----------- max # of neighbours belonging to a grid point
int NPtSphere_old ---- tmp variable
int isph ------------- index
int *NPtSphereiat ---- NPtSphereiat[n] = # of points over atom n
                      to generate SAS
int nx,ny,nz --------- indices
int inayb ------------ index
int listmp ----------- tmp variable
int ***nlist --------- nlist_re[ix][iy][iz]: # of atoms associated with cube
                       grid point ix,iy,iz
int **list ----------- list[m][n(ix,iy,iz)]: given the grid
                       point ix,iy,iz "list" tells which atom numbers are
                       associated with this point (m goes from
                       1 to nlist_re[ix][iy][iz] )
int first ------------ tmp variable
struct point len ----- ReMaxC - ReMinC
struct point cbmid --- Coor of the center of the grid box for the
                       neighbour list
struct point *ccrd --- cubing grid coords
struct point *SphCoor - Coor of the points placed around one atom
struct point sphpt --- Tmp variable
struct point d ------- Tmp variable
struct point midgpt -- nstep/2. + 1
##############################################################*/

  rmax = MaxDVector(ReRad, 1, ReAtNu);

/* Get some useful values */

  len.x = ReMaxC[1] - ReMinC[1];
  len.y = ReMaxC[2] - ReMinC[2];
  len.z = ReMaxC[3] - ReMinC[3];
  slen = (len.x>len.y) ? len.x : len.y;
  slen = (slen>len.z) ? slen : len.z;

/* Begin Shrake and Rupley surface generation
   Determine neighbour list-grid position and dimension */

  rslop = .2;
  rscale = 2.*(rmax + WaMoRa) + rslop;
  cbscal = 1./rscale;
  cbmid.x = ReMinC[1] + len.x/2.;
  cbmid.y = ReMinC[2] + len.y/2.;
  cbmid.z = ReMinC[3] + len.z/2.;

  nstep = (int) (( slen + 2.*(rmax + WaMoRa) ) * cbscal + 0.5) + 2;
  nstep = ( ((nstep >> 1) << 1) == nstep) ? nstep : nstep + 1;

  midg = ( nstep/2. + 1 );
  midgpt.x = nstep/2. + 1;
  midgpt.y = nstep/2. + 1;
  midgpt.z = nstep/2. + 1;

 /* convert all atom coords to cubing grid coords (ccrd): */

  cbrdp=dvector(1,ReAtNu);
  cbrdp2=dvector(1,ReAtNu);

  ccrd=structpointvect(1,ReAtNu);

  for (iat=1;iat<=ReAtNu;iat++) {
    cbrdp[iat] = cbscal * ( ReRad[iat] + WaMoRa );
    cbrdp2[iat] = cbrdp[iat] * cbrdp[iat];

    ccrd[iat].x = (ReCoor[iat][1] - cbmid.x)*cbscal + midg;
    ccrd[iat].y = (ReCoor[iat][2] - cbmid.y)*cbscal + midg;
    ccrd[iat].z = (ReCoor[iat][3] - cbmid.z)*cbscal + midg;

  }

 /* NOW GENERATE LIST OF ATOMS ASSOCIATED WITH EACH GRID POINT,
    BY TRUNCATION:
    keep track of max number of atoms associated with any grid point */

  nnnmax = -1000;
 /* initialize nlist: */
  nlist=i3tensor(1,nstep,1,nstep,1,nstep);
  for (ix=1;ix<=nstep;ix++) {
    for (iy=1;iy<=nstep;iy++) {
      for (iz=1;iz<=nstep;iz++) {
        nlist[ix][iy][iz] = 0;
      }
    }
  }

  list=imatrix(1,500,1,nstep*nstep*nstep);
  for (iat=1;iat<=ReAtNu;iat++) {
/* generate truncated cube coords of each atom, and use as subscript */
    ix = (int) ccrd[iat].x;
    iy = (int) ccrd[iat].y;
    iz = (int) ccrd[iat].z;

/* number of atoms associated with cube grid point ix,iy,iz: */
    ++nlist[ix][iy][iz];
    nnnmax = ( nlist[ix][iy][iz]>nnnmax ) ? nlist[ix][iy][iz] : nnnmax;

/* store index of atom associated with cube grid point ix,iy,iz: */
    list[nlist[ix][iy][iz]][nstep*nstep*(ix-1) + nstep*(iy-1) + iz]
    = iat;
  }

/* LOOP OVER ATOMS; FOR EACH ATOM, GENERATE CORRESPONDING SPHERE;
   DETERMINE GRID POINT ASSOCIATED WITH ATOM; LOOP OVER SPHERE POINTS
   OF ATOM; FOR EACH SPHERE POINT, LOOP OVER NEIGHBORING ATOMS
   (AS DETERMINED FROM LIST ARRAY); IF CURRENT SURFACE POINT IS
   EXCLUDED BY A NEIGHBOR, MOVE ON TO NEXT POINT ON CURRENT ATOM;
   IF NOT, MOVE ON TO NEXT NEIGHBOR ATOM; IF NO ATOMS EXCLUDE
   SURFACE POINT, ADD COORDS OF SURFACE POINT TO SRF ARRAY. */

/* Initialize to zero the number of S&R surface points */
  *PNsurfpt=0;

/* Begin loop over atoms */
  NPtSphereiat=ivector(1,ReAtNu);
  for (iat=1;iat<=ReAtNu;iat++) {
    first = 1;
    NPtSphereiat[iat] = (int) (ReSurfDens_apol * pi4 * ReRadOut2[iat]);

/* The following assignement is because NPtSphere can be a bit changed
   by SpherePoints(&NPtSphere) */
    NPtSphere_old = NPtSphereiat[iat];

    SphCoor=SpherePoints(&(NPtSphereiat[iat]));

/* Begin loop over surface point of atom iat: get the cubing coords
   of every surface point */
    for (isph=1;isph<=NPtSphereiat[iat];isph++) {
      sphpt = pt_plus_pt(c_per_pt(cbrdp[iat],SphCoor[isph]),
                         ccrd[iat]);

/* determine grid point associated with current sphere surface point: */
      nx = (int) (sphpt.x + 0.5);
      ny = (int) (sphpt.y + 0.5);
      nz = (int) (sphpt.z + 0.5);

/* loop over neighbor grid points.  The cubing grid is set up so that the
   only atoms which could "knock out" this possible protein surface point
   are those associated with the atom's own grid point, or that grid
   point's neighbors. */
      for (iz=nz-1;iz<=nz;iz++) {
        if (iz>=1 || iz <= nstep) {//to avoid to use an invalid index. clangini

          for (iy=ny-1;iy<=ny;iy++) {
            if (iy>=1 || iy <= nstep) {

              for (ix=nx-1;ix<=nx;ix++) {
                if (ix>=1 || ix <= nstep) {




/* loop over atoms associated with current grid point: */
                  for (inayb=1;inayb<=nlist[ix][iy][iz];inayb++) {

                    if ( (listmp = list[inayb][nstep*nstep*(ix-1) + nstep*(iy-1) + iz])
                      != iat ) {
                      d=pt_minus_pt(ccrd[listmp],sphpt);
                      d2 = pt_scal_pt(d,d);
                      if (d2<cbrdp2[listmp])
                        goto newpoint;
                    }
/* end loop over atoms of current neighbor grid point (inayb): */
                  }
                }
/* end loops over neighbor grid points for this possible surface point:
   ix loop: */
              }
            }
/* iy loop: */
          }
        }
/* iz loop: */
      }

/* if get here, sphere point has survived the neighbors, and is accepted
   as a Shrake and Rupley surface point */
      iatsrf[++(*PNsurfpt)] = iat;
      if (first) {
        pointsrf[iat] = *PNsurfpt;
        first = 0;
      }

/* convert sphpt coords back to real space and save in surface array */
      surfpt[*PNsurfpt] =
       pt_plus_pt(c_per_pt(rscale,pt_minus_pt(sphpt,midgpt)),cbmid);

/* end loop over sphere points */
      ++nsurf[iat];
newpoint:
      {
	continue;
      }
    }
/* end loop over atoms: */

    free_structpointvect(SphCoor,1,NPtSphere_old);
  }


/* Free memory */

  free_imatrix(list,1,500,1,nstep*nstep*nstep);
  free_i3tensor(nlist,1,nstep,1,nstep,1,nstep);
  free_structpointvect(ccrd,1,ReAtNu);
  free_dvector(cbrdp2,1,ReAtNu);

  free_ivector(NPtSphereiat,1,ReAtNu); /* dey check memory leak */

/*
  dey 26072006
  seg. fault with e.g. dibrommethane
  when
  int *PNsurfpt and   struct point *surfpt
  are allocated too few memory in main()
  -> modified "NPtSphereMax_XXX" in main()
*/
  free_dvector(cbrdp,1,ReAtNu);

  if (iat == ReAtNu+1 )
    return 1;
  else
    return 0;
}

int Surf_Grid_unif_dens_neighlist(int ReAtNu,double **ReCoor,double *ReMinC,
                                  double *ReRad,double *ReRadOut2,
                                  double WaMoRa,double ReSurfDens_deso,
                                  double pi4,int *PNsurfpt,
                                  struct point *surfpt,int *iatsrf,int *nsurf,
                                  int *pointsrf,struct point len,double slen,
                                  int ***nlist,int **list,
                                  int nstep,struct point *Pcbmid,int *Pmidg,
                                  double rscale)
/*##########################################
Place points over the SAS assuming a uniform surface density and transmit the
atom-to-grid neighbour list
###########################################*/

/*##########################################
int ReAtNu -------------- Tot # atoms
double **ReCoor ---------- Coordinates
double *ReMinC ----------- Min (along x,y,z) of ReCoor
double *ReRad ----------- Charge radii (=vdW radii apart "enclosed" H)
double *ReRadOut2 ------- (Charge radii + WaMoRa)^2
double WaMoRa ----------- Radius of the water molecule
double ReSurfDens_deso -- Surface point density ( in A^(-2) ) to generate SAS
                         for the receptor
double pi4 -------------- 4 * greekpi
int *PNsurfpt ----------- Tot # of points over SAS
struct point *surfpt ---- Coor of points over SAS
int *iatsrf ------------- iatsrf[n] = atom number generating the nth SAS point
int *nsurf -------------- nsurf_re[n] = amount of SAS surface points
                          generated by atom n
int *pointsrf ----------- pointsrf_re[n] = first rec SAS point (in the
                          list surfpt) that is generated by atom n
struct point len -------- ReMaxC - ReMinC
double slen ------------- Max dimension of the the grid box for the
                          neighbour list
int ***nlist ------------ nlist_re[ix][iy][iz]: # of atoms associated with cube
                          grid point ix,iy,iz
int **list -------------- list[m][n(ix,iy,iz)]: given the grid
                          point ix,iy,iz "list" tells which atom numbers are
                          associated with this point (m goes from
                          1 to nlist_re[ix][iy][iz] )
int nstep --------------- Dimension of the grid box for the neighbour list. It
                          is a cubic box of even elements on each side
struct point cbmid ------ Coor of the center of the grid box for the
                          neighbour list
int midg ---------------- nstep/2. + 1
double rscale ----------- 2.*(MaximumChargeRadius + WaMoRa) + rslop
###########################################*/

{
  double *cbrdp,*cbrdp2,d2;
  int iat,ix,iy,iz,nnnmax,NPtSphere_old,isph,*NPtSphereiat,
      nx,ny,nz,inayb,listmp,first;
  struct point *ccrd,*SphCoor,sphpt,d,midgpt;

/*#########################################################################
double *cbrdp -------- cbrdp[iat] = cbscal * ( ReRad[iat] + WaMoRa )
double *cbrdp2 ------- cbrdp^2
double d2
int iat,ix,iy,iz ----- multipurpose indices
int nnnmax ----------- max # of neighbours belonging to a grid point
int NPtSphere_old ---- tmp variable
int isph ------------- index
int *NPtSphereiat ---- NPtSphereiat[n] = # of points over atom n
                      to generate SAS
int nx,ny,nz --------- indices
int inayb ------------ index
int listmp ----------- tmp variable
int first ------------ tmp variable
struct point *ccrd --- cubing grid coords
struct point *SphCoor - Coor of the points placed around one atom
struct point sphpt --- Tmp variable
struct point d ------- Tmp variable
struct point midgpt -- nstep/2. + 1
##############################################################*/

/* get some values for the cubing grid */
  Pcbmid->x = ReMinC[1] + len.x/2.;
  Pcbmid->y = ReMinC[2] + len.y/2.;
  Pcbmid->z = ReMinC[3] + len.z/2.;

  *Pmidg = (  nstep/2. + 1 );
  midgpt.x = nstep/2. + 1;
  midgpt.y = nstep/2. + 1;
  midgpt.z = nstep/2. + 1;

 /* convert all atom coords to cubing grid coords (ccrd): */

  cbrdp=dvector(1,ReAtNu);
  cbrdp2=dvector(1,ReAtNu);
  ccrd=structpointvect(1,ReAtNu);

  for (iat=1;iat<=ReAtNu;iat++) {
    cbrdp[iat] = ( ReRad[iat] + WaMoRa ) / (rscale);
    cbrdp2[iat] = cbrdp[iat] * cbrdp[iat];

    ccrd[iat].x = (ReCoor[iat][1] - Pcbmid->x) / (rscale) + *Pmidg;
    ccrd[iat].y = (ReCoor[iat][2] - Pcbmid->y) / (rscale) + *Pmidg;
    ccrd[iat].z = (ReCoor[iat][3] - Pcbmid->z) / (rscale) + *Pmidg;
  }

 /* NOW GENERATE LIST OF ATOMS ASSOCIATED WITH EACH GRID POINT,
    BY TRUNCATION:
    keep track of max number of atoms associated with any grid point */

  nnnmax = -1000;
 /* initialize nlist: */
  for (ix=1;ix<=nstep;ix++) {
    for (iy=1;iy<=nstep;iy++) {
      for (iz=1;iz<=nstep;iz++) {
        nlist[ix][iy][iz] = 0;
      }
    }
  }

  for (iat=1;iat<=ReAtNu;iat++) {
/* generate truncated cube coords of each atom, and use as subscript */
    ix = (int) ccrd[iat].x;
    iy = (int) ccrd[iat].y;
    iz = (int) ccrd[iat].z;

/* number of atoms associated with cube grid point ix,iy,iz: */
    ++nlist[ix][iy][iz];
    nnnmax = ( nlist[ix][iy][iz]>nnnmax ) ? nlist[ix][iy][iz] : nnnmax;

/* store index of atom associated with cube grid point ix,iy,iz: */
    list[nlist[ix][iy][iz]][nstep*nstep*(ix-1) + nstep*(iy-1) + iz]
    = iat;
  }

/* LOOP OVER ATOMS; FOR EACH ATOM, GENERATE CORRESPONDING SPHERE;
   DETERMINE GRID POINT ASSOCIATED WITH ATOM; LOOP OVER SPHERE POINTS
   OF ATOM; FOR EACH SPHERE POINT, LOOP OVER NEIGHBORING ATOMS
   (AS DETERMINED FROM LIST ARRAY); IF CURRENT SURFACE POINT IS
   EXCLUDED BY A NEIGHBOR, MOVE ON TO NEXT POINT ON CURRENT ATOM;
   IF NOT, MOVE ON TO NEXT NEIGHBOR ATOM; IF NO ATOMS EXCLUDE
   SURFACE POINT, ADD COORDS OF SURFACE POINT TO SRF ARRAY. */

/* Initialize to zero the number of S&R surface points */
  *PNsurfpt=0;

/* Begin loop over atoms */
  NPtSphereiat=ivector(1,ReAtNu);
  for (iat=1;iat<=ReAtNu;iat++) {
    first = 1;
    NPtSphereiat[iat] = (int) (ReSurfDens_deso * pi4 * ReRadOut2[iat]);

/* The following assignement is because NPtSphere can be a bit changed
   by SpherePoints(&NPtSphere) */
    NPtSphere_old = NPtSphereiat[iat];

    SphCoor=SpherePoints(&(NPtSphereiat[iat]));

/* Begin loop over surface point of atom iat: get the cubing coords
   of every surface point */
    for (isph=1;isph<=NPtSphereiat[iat];isph++) {
      sphpt = pt_plus_pt(c_per_pt(cbrdp[iat],SphCoor[isph]),
                         ccrd[iat]);

/* determine grid point associated with current sphere surface point: */
      nx = (int) (sphpt.x + 0.5);
      ny = (int) (sphpt.y + 0.5);
      nz = (int) (sphpt.z + 0.5);

/* loop over neighbor grid points.  The cubing grid is set up so that the
   only atoms which could "knock out" this possible protein surface point
   are those associated with the atom's own grid point, or that grid
   point's neighbors. */
      for (iz=nz-1;iz<=nz;iz++) {
        if (iz>=1 || iz <= nstep) {

          for (iy=ny-1;iy<=ny;iy++) {
            if (iy>=1 || iy <= nstep) {

              for (ix=nx-1;ix<=nx;ix++) {
                if (ix>=1 || ix <= nstep) {

/* loop over atoms associated with current grid point: */
                  for (inayb=1;inayb<=nlist[ix][iy][iz];inayb++) {
                    if ( (listmp =
                      list[inayb][nstep*nstep*(ix-1) + nstep*(iy-1) + iz])
                      != iat ) {
                      d=pt_minus_pt(ccrd[listmp],sphpt);
                      d2 = pt_scal_pt(d,d);
                      if (d2<cbrdp2[listmp])
                        goto newpoint;
                    }
/* end loop over atoms of current neighbor grid point (inayb): */
                  }
                }
/* end loops over neighbor grid points for this possible surface point:
   ix loop: */
              }
            }
/* iy loop: */
          }
        }
/* iz loop: */
      }

/* if get here, sphere point has survived the neighbors, and is accepted
   as a Shrake and Rupley surface point */
      iatsrf[++(*PNsurfpt)] = iat;
      if (first) {
        pointsrf[iat] = *PNsurfpt;
        first = 0;
      }

/* convert sphpt coords back to real space and save in surface array */
      surfpt[*PNsurfpt] =
       pt_plus_pt(c_per_pt(rscale,pt_minus_pt(sphpt,midgpt)),*Pcbmid);

/* end loop over sphere points */
      ++nsurf[iat];
newpoint:
      continue;
    }
/* end loop over atoms: */
    free_structpointvect(SphCoor,1,NPtSphere_old);
  }

/* Free memory */
  free_ivector(NPtSphereiat,1,ReAtNu); /* prev. mem. leak.*/
  free_structpointvect(ccrd,1,ReAtNu);
  free_dvector(cbrdp2,1,ReAtNu);
  free_dvector(cbrdp,1,ReAtNu);

  if (iat == ReAtNu+1 )
    return 1;
  else
    return 0;
}

int Surf_Grid(int ReAtNu,double **ReCoor,double *ReMaxC,double *ReMinC,
               double *ReRad,double WaMoRa,int NPtSphere,int *PNsurfpt,
               struct point *surfpt,int *iatsrf,int *nsurf,int *pointsrf)
/*##########################################
Place points over the SAS assuming a given # of points per atom (NPtSphere)
###########################################*/

/*##########################################
int ReAtNu -------------- Tot # atoms
double **ReCoor ---------- Coordinates
double *ReMaxC ----------- Max (along x,y,z) of ReCoor
double *ReMinC ----------- Min (along x,y,z) of ReCoor
double *ReRad ----------- Charge radii (=vdW radii apart "enclosed" H)
double WaMoRa ----------- Radius of the water molecule
int NPtSphere ----------- # of points per atom
int *PNsurfpt ----------- Tot # of points over SAS
struct point *surfpt ---- Coor of points over SAS
int *iatsrf ------------- iatsrf[n] = atom number generating the nth SAS point
int *nsurf -------------- nsurf_re[n] = amount of SAS surface points
                          generated by atom n
int *pointsrf ----------- pointsrf_re[n] = first rec SAS point (in the
                          list surfpt) that is generated by atom n
###########################################*/
{
  double rmax,slen,rslop,rscale,cbscal,*cbrdp,*cbrdp2,d2;
  int iat,ix,iy,iz,nstep,midg,nnnmax,NPtSphere_old,isph,
      nx,ny,nz,inayb,listmp,***nlist,**list,first;
  struct point len,cbmid,*ccrd,*SphCoor,sphpt,d,midgpt;

/*#########################################################################
double rmax ---------- Largest charge radius
double slen ---------- Max dimension of the the grid box for the
                      neighbour list
double rslop --------- Parameter introduced to define the grid box for the
                      neighbour list
double rscale -------- 2.*(MaximumChargeRadius + WaMoRa) + rslop
double cbscal -------- 1./rscale
double *cbrdp -------- cbrdp[iat] = cbscal * ( ReRad[iat] + WaMoRa )
double *cbrdp2 ------- cbrdp^2
double d2
int iat,ix,iy,iz ----- multipurpose indices
int nstep ------------ Dimension of the grid box for the neighbour list. It
                       is a cubic box of even elements on each side
int midg ------------- nstep/2. + 1
int nnnmax ----------- max # of neighbours belonging to a grid point
int NPtSphere_old ---- tmp variable
int isph ------------- index
int *NPtSphereiat ---- NPtSphereiat[n] = # of points over atom n
                      to generate SAS
int nx,ny,nz --------- indices
int inayb ------------ index
int listmp ----------- tmp variable
int ***nlist --------- nlist_re[ix][iy][iz]: # of atoms associated with cube
                       grid point ix,iy,iz
int **list ----------- list[m][n(ix,iy,iz)]: given the grid
                       point ix,iy,iz "list" tells which atom numbers are
                       associated with this point (m goes from
                       1 to nlist_re[ix][iy][iz] )
int first ------------ tmp variable
struct point len ----- ReMaxC - ReMinC
struct point cbmid --- Coor of the center of the grid box for the
                       neighbour list
struct point *ccrd --- cubing grid coords
struct point *SphCoor - Coor of the points placed around one atom
struct point sphpt --- Tmp variable
struct point d ------- Tmp variable
struct point midgpt -- nstep/2. + 1
##############################################################*/

  rmax = MaxDVector(ReRad, 1, ReAtNu);

/* Get some useful values */

  len.x = ReMaxC[1] - ReMinC[1];
  len.y = ReMaxC[2] - ReMinC[2];
  len.z = ReMaxC[3] - ReMinC[3];
  slen = (len.x>len.y) ? len.x : len.y;
  slen = (slen>len.z) ? slen : len.z;

/* Begin Shrake and Rupley surface generation
   Determine neighbour list-grid position and dimension */

  rslop = .2;
  rscale = 2.*(rmax + WaMoRa) + rslop;
  cbscal = 1./rscale;
  cbmid.x = ReMinC[1] + len.x/2.;
  cbmid.y = ReMinC[2] + len.y/2.;
  cbmid.z = ReMinC[3] + len.z/2.;

  nstep = (int) (( slen + 2.*(rmax + WaMoRa) ) * cbscal + 0.5) + 2;
  nstep = ( ((nstep >> 1) << 1) == nstep) ? nstep : nstep + 1;

  midg = nstep/2. + 1;
  midgpt.x = nstep/2. + 1;
  midgpt.y = nstep/2. + 1;
  midgpt.z = nstep/2. + 1;

 /* convert all atom coords to cubing grid coords (ccrd): */

  cbrdp=dvector(1,ReAtNu);
  cbrdp2=dvector(1,ReAtNu);
  ccrd=structpointvect(1,ReAtNu);

  for (iat=1;iat<=ReAtNu;iat++) {
    cbrdp[iat] = cbscal * ( ReRad[iat] + WaMoRa );
    cbrdp2[iat] = cbrdp[iat] * cbrdp[iat];

    ccrd[iat].x = (ReCoor[iat][1] - cbmid.x)*cbscal + midg;
    ccrd[iat].y = (ReCoor[iat][2] - cbmid.y)*cbscal + midg;
    ccrd[iat].z = (ReCoor[iat][3] - cbmid.z)*cbscal + midg;
  }

 /* NOW GENERATE LIST OF ATOMS ASSOCIATED WITH EACH GRID POINT,
    BY TRUNCATION:
    keep track of max number of atoms associated with any grid point */

  nnnmax = -1000;
 /* initialize nlist: */
  nlist=i3tensor(1,nstep,1,nstep,1,nstep);
  for (ix=1;ix<=nstep;ix++) {
    for (iy=1;iy<=nstep;iy++) {
      for (iz=1;iz<=nstep;iz++) {
        nlist[ix][iy][iz] = 0;
      }
    }
  }

  list=imatrix(1,500,1,nstep*nstep*nstep);
  for (iat=1;iat<=ReAtNu;iat++) {
/* generate truncated cube coords of each atom, and use as subscript */
    ix = (int) ccrd[iat].x;
    iy = (int) ccrd[iat].y;
    iz = (int) ccrd[iat].z;

/* number of atoms associated with cube grid point ix,iy,iz: */
    ++nlist[ix][iy][iz];
    nnnmax = ( nlist[ix][iy][iz]>nnnmax ) ? nlist[ix][iy][iz] : nnnmax;

/* store index of atom associated with cube grid point ix,iy,iz: */
    list[nlist[ix][iy][iz]][nstep*nstep*(ix-1) + nstep*(iy-1) + iz]
    = iat;
  }

/* The following assignement is because NPtSphere can be a bit changed
   by SpherePoints(&NPtSphere) */
  NPtSphere_old = NPtSphere;

/* Get (in SphCoor) the coordinates of about NPtSphere unifromly
   distributed over a unit sphere */
  SphCoor=SpherePoints(&NPtSphere);

/* LOOP OVER ATOMS; FOR EACH ATOM, GENERATE CORRESPONDING SPHERE;
   DETERMINE GRID POINT ASSOCIATED WITH ATOM; LOOP OVER SPHERE POINTS
   OF ATOM; FOR EACH SPHERE POINT, LOOP OVER NEIGHBORING ATOMS
   (AS DETERMINED FROM LIST ARRAY); IF CURRENT SURFACE POINT IS
   EXCLUDED BY A NEIGHBOR, MOVE ON TO NEXT POINT ON CURRENT ATOM;
   IF NOT, MOVE ON TO NEXT NEIGHBOR ATOM; IF NO ATOMS EXCLUDE
   SURFACE POINT, ADD COORDS OF SURFACE POINT TO SRF ARRAY. */

/* Initialize to zero the number of S&R surface points */
  *PNsurfpt=0;

/* Begin loop over atoms */
  for (iat=1;iat<=ReAtNu;iat++) {
    first = 1;

/* Begin loop over surface point of atom iat: get the cubing coords
   of every surface point */
    for (isph=1;isph<=NPtSphere;isph++) {
      sphpt = pt_plus_pt(c_per_pt(cbrdp[iat],SphCoor[isph]),
                         ccrd[iat]);

/* determine grid point associated with current sphere surface point: */
      nx = (int) (sphpt.x + 0.5);
      ny = (int) (sphpt.y + 0.5);
      nz = (int) (sphpt.z + 0.5);

/* loop over neighbor grid points.  The cubing grid is set up so that the
   only atoms which could "knock out" this possible protein surface point
   are those associated with the atom's own grid point, or that grid
   point's neighbors. */
      for (iz=nz-1;iz<=nz;iz++) {
        if (iz>=1 || iz <= nstep) {

          for (iy=ny-1;iy<=ny;iy++) {
            if (iy>=1 || iy <= nstep) {

              for (ix=nx-1;ix<=nx;ix++) {
                if (ix>=1 || ix <= nstep) {

/* loop over atoms associated with current grid point: */
                  for (inayb=1;inayb<=nlist[ix][iy][iz];inayb++) {
                    if ( (listmp =
                      list[inayb][nstep*nstep*(ix-1) + nstep*(iy-1) + iz])
                      != iat ) {
                      d=pt_minus_pt(ccrd[listmp],sphpt);
                      d2 = pt_scal_pt(d,d);
                      if (d2<cbrdp2[listmp])
                        goto newpoint;
                    }
/* end loop over atoms of current neighbor grid point (inayb): */
                  }
                }
/* end loops over neighbor grid points for this possible surface point:
   ix loop: */
              }
            }
/* iy loop: */
          }
        }
/* iz loop: */
      }

/* if get here, sphere point has survived the neighbors, and is accepted
   as a Shrake and Rupley surface point */
      iatsrf[++(*PNsurfpt)] = iat;
      if (first) {
        pointsrf[iat] = *PNsurfpt;
        first = 0;
      }

/* convert sphpt coords back to real space and save in surface array */
      surfpt[*PNsurfpt] =
       pt_plus_pt(c_per_pt(rscale,pt_minus_pt(sphpt,midgpt)),cbmid);

/* end loop over sphere points */
      ++nsurf[iat];
newpoint:
      continue;
    }
/* end loop over atoms: */
  }

/* Free memory */
  free_structpointvect(SphCoor,1,NPtSphere_old);
  free_imatrix(list,1,500,1,nstep*nstep*nstep);
  free_i3tensor(nlist,1,nstep,1,nstep,1,nstep);
  free_structpointvect(ccrd,1,ReAtNu);
  free_dvector(cbrdp2,1,ReAtNu);
  free_dvector(cbrdp,1,ReAtNu);

  if (iat == ReAtNu+1 )
    return 1;
  else
    return 0;
}

struct point *SpherePoints(int *PNPtSphere)
/*##########################################
Place approximately *PNPtSphere points over a sphere (and return
the new *PNPtSphere value)
############################################*/
{
  double theta,stheta,ctheta,psi,cpsi,psistp,spsi,frac,thtstp;
  int ntheta, npsmax,i, npsi,nt,np;
  struct point *SphCoor;

  frac = *PNPtSphere / 4.;
  ntheta = (sqrt(3.1415926535897931 * frac));
  npsmax = (2. * ntheta);
  thtstp = 3.1415926535897931 / ntheta;
  i=1;

  SphCoor=structpointvect(1,*PNPtSphere);
  for (nt=1;nt<=ntheta;nt++) {
    theta = thtstp * nt;
    stheta = sinf(theta);
    ctheta = cosf(theta);
    npsi = (int) (stheta*npsmax + 0.5);
    if (npsi != 0) {
      psistp = 6.2831853071795862 / npsi;
      for (np=1;np<=npsi;np++) {
        psi = np * psistp;
        cpsi = cosf(psi);
        spsi = sinf(psi);
        SphCoor[i].x = cpsi * stheta;
        SphCoor[i].y = spsi * stheta;
        SphCoor[i].z = ctheta;
        i = i + 1;
      }
    }
  }

  *PNPtSphere = i - 1;

  return SphCoor;
}

int Excl_Grid(int ReAtNu,double **ReCoor,struct point Min,double *ReRadOut,
              double *ReRadOut2,double WaMoRa,double GrSiSo,
              int NStartGridx,int NStartGridy,int NStartGridz,
              int NGridx,int NGridy,int NGridz,
              double UnitVol,char ***GridMat,int Nsurfpt,
              struct point *surfpt,double *SelfVol,double *SelfVol_corrB,
	            char *EmpCorrB)
/*########################################################
Place a sphere (WaMoRa) over each SAS point (surfpt) and set to empty
(GridMat = 'e') all the grid points occupied by the sphere
##########################################################*/

/*##########################################
int ReAtNu -------------- Tot # rec (frag) atoms
double **ReCoor ---------- Rec (frag) coordinates
struct point Min -------- Min coor of the grid box
double *ReRadOut -------- Rec (frag)  charge radii + WaMoRa
double *ReRadOut2 ------- (Rec (frag)  charge radii + WaMoRa)^2
double WaMoRa ----------- Radius of the water molecule
double GrSiSo ----------- Size of the 3D grid used for cont. electrostatics
int NStartGridx --------- First grid point along x (1)
int NStartGridy --------- First grid point along y (1)
int NStartGridz --------- First grid point along z (1)
int  NGridx ------------- Tot # of grid points along x
int  NGridy ------------- Tot # of grid points along y
int  NGridz ------------- Tot # of grid points along z
double UnitVol ---------- Volume of the grid element for cont. elec.
char ***GridMat --------- Matrix telling if a grid point is occupied (o),
                          empty (e), or if it belongs to the interface
                          between SAS and MS (s)
int Nsurfpt ------------- Tot # of points over SAS
struct point *surfpt ---- Coor of points over SAS
double *SelfVol --------- SelfVol[iat] = integral of 1/r^4 over the solute
                          volume for rec (frag) atom iat
###########################################*/
{
  int i,iat,ix,iy,iz,ixmin,iymin,izmin,ixmax,iymax,izmax;
  double WaMoRa2,xtemp,x2temp,ytemp,xy2temp,ztemp,r2,r4;

  WaMoRa2 = WaMoRa * WaMoRa;

/* Loop over S&R surface points, place a probe sphere over each of them
   and set to 's' all the 3D grid point covered by it */

  for (i=1;i<=Nsurfpt;i++) {

/* Calculate the extremes, in the grid frame, of the cube containing the
   probe sphere. They are integer! */

    ixmin =( (surfpt[i].x - WaMoRa - Min.x) / GrSiSo + 1 );
    iymin =( (surfpt[i].y - WaMoRa - Min.y) / GrSiSo + 1 );
    izmin =( (surfpt[i].z - WaMoRa - Min.z) / GrSiSo + 1 );
    ixmax =( (surfpt[i].x + WaMoRa - Min.x) / GrSiSo + 1 );
    iymax =( (surfpt[i].y + WaMoRa - Min.y) / GrSiSo + 1 );
    izmax =( (surfpt[i].z + WaMoRa - Min.z) / GrSiSo + 1 );
    ixmin = (ixmin > NStartGridx) ? ixmin : NStartGridx;
    iymin = (iymin > NStartGridy) ? iymin : NStartGridy;
    izmin = (izmin > NStartGridz) ? izmin : NStartGridz;
    ixmax = (ixmax < NGridx) ? ixmax : NGridx;
    iymax = (iymax < NGridy) ? iymax : NGridy;
    izmax = (izmax < NGridz) ? izmax : NGridz;

/* xtemp,ytemp,ztemp are cartesian coordinates of the grid pt ix,iy,iz,
   in a frame with the origin in the suface point surfpt[i] */

    xtemp = GrSiSo * (ixmin - 1.5) + Min.x - surfpt[i].x;
    for (ix=ixmin;ix<=ixmax;ix++) {
      xtemp += GrSiSo;
      x2temp = xtemp*xtemp;
      ytemp = GrSiSo * (iymin - 1.5) + Min.y - surfpt[i].y;
      for (iy=iymin;iy<=iymax;iy++) {
        ytemp += GrSiSo;
        xy2temp = ytemp*ytemp + x2temp;

/*  Check if we are already out of the sphere */

        if (WaMoRa2 > xy2temp) {
          ztemp = GrSiSo * (izmin - 1.5) + Min.z - surfpt[i].z;
          for (iz=izmin;iz<=izmax;iz++) {
            ztemp += GrSiSo;
            if (GridMat[ix][iy][iz] == 'o') {
              r2 = ztemp * ztemp + xy2temp;

/* Check if our grid point is inside the sphere. If yes change GridMat */

              if (WaMoRa2 > r2){
                GridMat[ix][iy][iz] = 's';//'s';
              }
            }
          }
        }
      }
    }
  }

/* Loop over the atom and when they include a 3D grid point marked 's'
   assign a negative contribution to the integral Selfvol. This is
   done to balance a successive assignment of the same contribution
   with a sign +. It is done to make things faster and simpler. */

  //int jumppoint; //debug
  for (iat=1;iat<=ReAtNu;iat++) {
    //jumppoint = 0;
    ixmin = ( (ReCoor[iat][1] - ReRadOut[iat] - Min.x) / GrSiSo + 1 );
    iymin = ( (ReCoor[iat][2] - ReRadOut[iat] - Min.y) / GrSiSo + 1 );
    izmin = ( (ReCoor[iat][3] - ReRadOut[iat] - Min.z) / GrSiSo + 1 );
    ixmax = ( (ReCoor[iat][1] + ReRadOut[iat] - Min.x) / GrSiSo + 1 );
    iymax = ( (ReCoor[iat][2] + ReRadOut[iat] - Min.y) / GrSiSo + 1 );
    izmax = ( (ReCoor[iat][3] + ReRadOut[iat] - Min.z) / GrSiSo + 1 );
    ixmin = (ixmin > NStartGridx) ? ixmin : NStartGridx;
    iymin = (iymin > NStartGridy) ? iymin : NStartGridy;
    izmin = (izmin > NStartGridz) ? izmin : NStartGridz;
    ixmax = (ixmax < NGridx) ? ixmax : NGridx;
    iymax = (iymax < NGridy) ? iymax : NGridy;
    izmax = (izmax < NGridz) ? izmax : NGridz;

    xtemp = GrSiSo * (ixmin - 1.5) + Min.x - ReCoor[iat][1];
    for (ix=ixmin;ix<=ixmax;ix++) {
      xtemp += GrSiSo;
      x2temp = xtemp*xtemp;
      ytemp = GrSiSo * (iymin - 1.5) + Min.y - ReCoor[iat][2];
      for (iy=iymin;iy<=iymax;iy++) {
        ytemp += GrSiSo;
        xy2temp = ytemp*ytemp + x2temp;
        if (ReRadOut2[iat] > xy2temp) {
          ztemp = GrSiSo * (izmin - 1.5) + Min.z - ReCoor[iat][3];
          for (iz=izmin;iz<=izmax;iz++) {
            ztemp += GrSiSo;
            r2 = ztemp * ztemp + xy2temp;
            if (ReRadOut2[iat] > r2) {
              if (GridMat[ix][iy][iz] == 's') { // 's'
                r4 = r2 * r2;
                SelfVol[iat] -= UnitVol/r4;
		            if (EmpCorrB[0]=='y')
		              SelfVol_corrB[iat] -= UnitVol/(r4*sqrt(r2));
                /*std::cout << "here subtracting volume for atom " << iat << "\n";*/
                //jumppoint++;
              }
            }
          }
        }
      }
    }
    //std::cout << "Subtracted surf point for atom " << iat << " = " << jumppoint <<"\n";
  }
  if (iat == ReAtNu+1 )
    return 1;
  else
    return 0;
}

int Get_Self_Vol(int ReAtNu,double **ReCoor,double *ReRadOut2,
                 double GrSiSo,double *XGrid,double *YGrid,double *ZGrid,
                 int NStartGridx,int NStartGridy,int NStartGridz,
                 int NGridx,int NGridy,int NGridz,
                 double UnitVol,char ***GridMat,double *SelfVol,
		             double *SelfVol_corrB,char *EmpCorrB)
/*##########################################
Calculate the integral of 1/r^4 for each atom over the
solute volume defined by GridMat
###########################################*/

/*##########################################
int ReAtNu -------------- Tot # atoms
double **ReCoor ---------- Coordinates
double *ReRadOut2 ------- (Rec charge radii + WaMoRa)^2
double GrSiSo ----------- Size of the 3D grid used for cont. electrostatics
double *XGrid ----------- X coor of the grid points
double *YGrid ----------- Y coor of the grid points
double *ZGrid ----------- Z coor of the grid points
int NStartGridx --------- First grid point along x (1)
int NStartGridy --------- First grid point along y (1)
int NStartGridz --------- First grid point along z (1)
int  NGridx ------------- Tot # of grid points along x
int  NGridy ------------- Tot # of grid points along y
int  NGridz ------------- Tot # of grid points along z
double UnitVol ---------- Volume of the grid element for cont. elec.
char ***GridMat --------- Matrix telling if a grid point is occupied by the
                         rec (o), empty (e), or if it belongs to the interface
                         between SAS and MS (s)
double *SelfVol --------- SelfVol[iat] = integral of 1/r^4 over the solute
                          volume for atom iat
###########################################*/
{
  int ix,iy,iz,ixf,iyf,izf,Neigh,*NeighList,iat,jat,filled,dolist;
  double UnitVol27,cutoff,cutoff2,cutoff2_grid,r2,*SelfVolTmp,*SelfVolTmp_corrB;

  cutoff = 5.;      /* cutoff between fine and coarse 3D grid */
  cutoff2 = cutoff * cutoff;
  cutoff2_grid = (cutoff+sqrt(3.)*GrSiSo)*(cutoff+sqrt(3.)*GrSiSo);
  UnitVol27 = UnitVol * 27.;

  NeighList=ivector(1,ReAtNu);
  SelfVolTmp=dvector(1,ReAtNu);
  SelfVolTmp_corrB=dvector(1,ReAtNu);
  for (iat=1;iat<=ReAtNu;iat++)
  {
    SelfVolTmp[iat]=0.0000000000000000;
    SelfVolTmp_corrB[iat]=0.0000000000000000;
  }

/* The grid points falling at a distance smaller than cutoff are accounted
as they are. The ones falling at a larger distance are grouped in cubes of
27 elements and are accounted according to the fact that the central grid
point or not.
 NeighList gives the list af atoms that are close to a 27-element cube */

/* Loop over the centers of the cubes of 27 elements */
  for (ix=NStartGridx+1; ix<=NGridx; ix+=3) {
    for (iy=NStartGridy+1; iy<=NGridy; iy+=3) {
      for (iz=NStartGridz+1; iz<=NGridz; iz+=3) {
        filled = 0;
        dolist = 0;
        if (GridMat[ix][iy][iz] == 'o') {
          filled = 1;
          dolist = 1;
        }
        else {
/* Loop inside the cubes of 27 elements to check if at least one is occupied */
          for (ixf=ix-1;ixf<=ix+1;ixf++) {
            for (iyf=iy-1;iyf<=iy+1;iyf++) {
              for (izf=iz-1;izf<=iz+1;izf++) {
                if (GridMat[ixf][iyf][izf] == 'o') {
                  dolist = 1;
                  goto do_neigh_list; // I think this goto is simply a "multiple break to get out of the nested loops. clangini
                }
              }
            }
          }
        }
do_neigh_list:
/* If come here it means that at least one point of the
27-element cube is occupied */
/* Now loop inside the cube and over the atoms and assign the contribution
to the integral of 1/r^4 */
        if (dolist) {
          Neigh = 0;
          for (iat=1;iat<=ReAtNu;iat++) {
            r2 = (XGrid[ix]-ReCoor[iat][1])*(XGrid[ix]-ReCoor[iat][1])+
                 (YGrid[iy]-ReCoor[iat][2])*(YGrid[iy]-ReCoor[iat][2])+
                 (ZGrid[iz]-ReCoor[iat][3])*(ZGrid[iz]-ReCoor[iat][3]);
            if ( r2 > cutoff2_grid ) {
              if (filled) {
                SelfVolTmp[iat] += UnitVol27 / (r2 * r2);
		            if (EmpCorrB[0]=='y')
		              SelfVolTmp_corrB[iat] += UnitVol27 / (r2 * r2 * sqrt(r2));
              }
            }
            else if (r2 < cutoff2_grid && r2 > cutoff2) {
              if (filled) {
                SelfVolTmp[iat] += UnitVol27 / (r2 * r2);
		          if (EmpCorrB[0]=='y')
		            SelfVolTmp_corrB[iat] += UnitVol27 / (r2 * r2 * sqrt(r2));
	            }
              ++Neigh;
              NeighList[Neigh] = iat;
            }
            else {
              ++Neigh;
              NeighList[Neigh] = iat;
            }
          }
          for (ixf=ix-1;ixf<=ix+1;ixf++) {
            for (iyf=iy-1;iyf<=iy+1;iyf++) {
              for (izf=iz-1;izf<=iz+1;izf++) {
                if (GridMat[ixf][iyf][izf] == 'o') {
                  for (jat=1;jat<=Neigh;jat++) {
                    r2 = (XGrid[ixf]-ReCoor[NeighList[jat]][1]) *
                         (XGrid[ixf]-ReCoor[NeighList[jat]][1]) +
                         (YGrid[iyf]-ReCoor[NeighList[jat]][2]) *
                         (YGrid[iyf]-ReCoor[NeighList[jat]][2]) +
                         (ZGrid[izf]-ReCoor[NeighList[jat]][3]) *
                         (ZGrid[izf]-ReCoor[NeighList[jat]][3]);
                    if ( r2 < cutoff2 && r2 > ReRadOut2[NeighList[jat]] ) {
                      SelfVolTmp[NeighList[jat]] += UnitVol / (r2 * r2);
		                  if (EmpCorrB[0]=='y')
		                    SelfVolTmp_corrB[NeighList[jat]] += UnitVol / (r2 * r2 * sqrt(r2));
		                }
                  }
                }
              }
            }
          }
        } // end of if(dolist)
      }
    }
  }

  for (iat=1;iat<=ReAtNu;iat++)
  {
    SelfVol[iat] += SelfVolTmp[iat];
    SelfVol_corrB[iat] += SelfVolTmp_corrB[iat];
  }

  free_dvector(SelfVolTmp,1,ReAtNu);
  free_dvector(SelfVolTmp_corrB,1,ReAtNu);
  free_ivector(NeighList,1,ReAtNu);

  if (ix > NGridx )
    return 1;
  else
    return 0;
}

int GB_int(int ReAtNu,double **ReCoor,double *RePaCh,double *EffRad,
           double Ksolv,double *PIntEnTot)
/*##########################################
Calculate the interaction energy according to the GB formula
###########################################*/

/*##########################################
int ReAtNu -------------- Tot # atoms
double **ReCoor ---------- Coordinates
double *RePaCh ---------- partial charges
double *EffRad ---------- Effective radii
double Ksolv ------------ Constant
double *PIntEnTot ------- Tot interaction energy
###########################################*/
{
  int iat,jat;
  double Rij,Chiat,Inte,d2;

  *PIntEnTot = 0.;
  for (iat=1;iat<=ReAtNu;iat++) {
    if (RePaCh[iat] != 0. ) {
      Chiat = Ksolv * RePaCh[iat];
      for (jat=iat+1;jat<=ReAtNu;jat++) {
        if (RePaCh[jat] != 0. ) {  /* MODIFICATION 13/11/99 */
          d2 = (ReCoor[iat][1]-ReCoor[jat][1])*(ReCoor[iat][1]-ReCoor[jat][1]) +
               (ReCoor[iat][2]-ReCoor[jat][2])*(ReCoor[iat][2]-ReCoor[jat][2]) +
               (ReCoor[iat][3]-ReCoor[jat][3])*(ReCoor[iat][3]-ReCoor[jat][3]);
          Rij = EffRad[iat] * EffRad[jat];


          Inte = Chiat*RePaCh[jat] / sqrtf (d2 + Rij * exp(-d2/(4.*Rij))); // eq. (15) of [Scarsi et. al 1997]. clangini
          *PIntEnTot += Inte;
        }
      }
    }
  }

  if (iat == ReAtNu+1 )
    return 1;
  else
    return 0;
}

int Get_BS_Info(int ReAtNu,double **ReCoor,double *ReRadOut,
                struct point Min,double GrSiSo,double GrInSo,int *ReResN,
                int BSResN,int *BSReNu,int NGridx,
                int NGridy,int NGridz,int *PnxminBS,int *PnyminBS,
                int *PnzminBS,int *PnxmaxBS,int *PnymaxBS,
                int *PnzmaxBS,int *ResFirstAtom,int *NAtom_per_Res)
/*##########################################
Get informations about the size and the location of the Binding Site
###########################################*/

/*##########################################
int ReAtNu -------------- Tot # rec atoms
double **ReCoor ---------- Rec coordinates
double *ReRadOut -------- Rec charge radii + WaMoRa
struct point Min -------- Min coor of the grid box
double GrSiSo ----------- Size of the 3D grid used for cont. electrostatics
double GrInSo ----------- Margin left along each dimension (positive and
                         negative) of the rec for building the 3D grid box
int *ReResN ------------ # rec residue
int BSResN ------------- Tot number of rec residues in the Binding Site
int *BSReNu ------------ # rec residue in the Binding Site
int NGridx ----------- Tot # of grid points along x
int NGridy ----------- Tot # of grid points along y
int NGridz ----------- Tot # of grid points along z
int *PnxminBS ---------- grid point (along x) where the BS starts
int *PnyminBS ---------- grid point (along y) where the BS starts
int *PnzminBS ---------- grid point (along z) where the BS starts
int *PnxmaxBS ---------- grid point (along x) where the BS ends
int *PnymaxBS ---------- grid point (along y) where the BS ends
int *PnzmaxBS ---------- grid point (along z) where the BS ends
int *ResFirstAtom ---------- ResFirstAtom[n] = atom number of the first atom
                             belonging to residue n
int *NAtom_per_Res --------- NAtom_per_Res[n] = # of atoms of residue n
###########################################*/
{
  int iat,ires,Previous,Nres,Natom,nxlow,nylow,nzlow,nxhi,nyhi,nzhi;

/* Make two arrays that given the number of a residue, give the number
   of the first atom of that residue and the number of atoms belonging
   to this residue */

  Natom = 0;
  Previous = 0;
  Nres = 0;
  for (iat=1;iat<=ReAtNu;iat++) {
    if (ReResN[iat] == Previous+1 ) {
      if (Nres != 0 )
        NAtom_per_Res[Nres] = Natom;
      ResFirstAtom[++Nres] = iat;
      Previous++;
      Natom = 0;
    }
    Natom++;
  }
  NAtom_per_Res[Nres] = Natom;


/*  Calculate the extremes in the grid space (*PnxminBS,*PnxmaxBS,...) of
    the cube containing the residues of the binding site */

  iat = ResFirstAtom[BSReNu[1]];
  *PnxminBS = (ReCoor[iat][1] - ReRadOut[iat] - Min.x) / GrSiSo + 1;
  *PnyminBS = (ReCoor[iat][2] - ReRadOut[iat] - Min.y) / GrSiSo + 1;
  *PnzminBS = (ReCoor[iat][3] - ReRadOut[iat] - Min.z) / GrSiSo + 1;
  *PnxmaxBS = (ReCoor[iat][1] + ReRadOut[iat] - Min.x) / GrSiSo + 1;
  *PnymaxBS = (ReCoor[iat][2] + ReRadOut[iat] - Min.y) / GrSiSo + 1;
  *PnzmaxBS = (ReCoor[iat][3] + ReRadOut[iat] - Min.z) / GrSiSo + 1;

  for (ires=1;ires<=BSResN;ires++) {
    for (iat=ResFirstAtom[BSReNu[ires]];
         iat<=ResFirstAtom[BSReNu[ires]]+NAtom_per_Res[BSReNu[ires]]-1;
         iat++) {

	/* exception handling dey
	   this conditions are meant when e.g. the input file was not
	   setup properly */
	if(ResFirstAtom[BSReNu[ires]]+NAtom_per_Res[BSReNu[ires]]-1 <0)
	{
	    fprintf(stderr,"WARNING, caught exception while setting up binding site data, exiting ! ");
	    exit(12);
	}

      nxlow = (ReCoor[iat][1] - ReRadOut[iat] - Min.x - GrInSo) / GrSiSo + 1;
      nylow = (ReCoor[iat][2] - ReRadOut[iat] - Min.y - GrInSo) / GrSiSo + 1;
      nzlow = (ReCoor[iat][3] - ReRadOut[iat] - Min.z - GrInSo) / GrSiSo + 1;
      nxhi = (ReCoor[iat][1] + ReRadOut[iat] - Min.x + GrInSo) / GrSiSo + 1;
      nyhi = (ReCoor[iat][2] + ReRadOut[iat] - Min.y + GrInSo) / GrSiSo + 1;
      nzhi = (ReCoor[iat][3] + ReRadOut[iat] - Min.z + GrInSo) / GrSiSo + 1;
      *PnxminBS = (nxlow < *PnxminBS ) ? nxlow : *PnxminBS;
      *PnyminBS = (nylow < *PnyminBS ) ? nylow : *PnyminBS;
      *PnzminBS = (nzlow < *PnzminBS ) ? nzlow : *PnzminBS;
      *PnxmaxBS = (nxhi > *PnxmaxBS ) ? nxhi : *PnxmaxBS;
      *PnymaxBS = (nyhi > *PnymaxBS ) ? nyhi : *PnymaxBS;
      *PnzmaxBS = (nzhi > *PnzmaxBS ) ? nzhi : *PnzmaxBS;
    }
  }

  *PnxminBS = (*PnxminBS > 1) ? *PnxminBS : 1;
  *PnyminBS = (*PnyminBS > 1) ? *PnyminBS : 1;
  *PnzminBS = (*PnzminBS > 1) ? *PnzminBS : 1;
  *PnxmaxBS = (*PnxmaxBS < NGridx) ? *PnxmaxBS : NGridx;
  *PnymaxBS = (*PnymaxBS < NGridy) ? *PnymaxBS : NGridy;
  *PnzmaxBS = (*PnzmaxBS < NGridz) ? *PnzmaxBS : NGridz;

  if (ires == BSResN+1)
    return 1;
  else
    return 0;
}

int vdW_Surf(int ReAtNu,double **ReCoor,double *ReVdWE_sr,double *ReVdWR,
             double Sphere_apol,int BSResN,int *BSReNu,int *NAtom_per_Res,
             int *ResFirstAtom,struct point *surfpt_re,int *nsurf_re,
             int *pointsrf_re,double *vdWSurfPt)
/*##########################################
Calculate the vdW interaction of a sphere rolled over the rec surface
with the rec itself
###########################################*/

/*##########################################
int ReAtNu -------------- Tot # rec atoms
double **ReCoor ---------- Rec coordinates
double *ReVdWE_sr -------- Rec vdW energies
double *ReVdWR ----------- Rec vdW radii
double Sphere_apol ------ Dimension of the probe sphere used to generate SAS3
int BSResN -------------- Tot number of rec residues in the Binding Site
int *BSReNu ------------- # rec residue in the Binding Site
int *NAtom_per_Res ------ NAtom_per_Res[n] = # of atoms of residue n
int *ResFirstAtom ------- ResFirstAtom[n] = atom number of the first atom
                          belonging to residue n
struct point *surfpt_re -  Coor of points over SAS3
double *vdWSurfPt ------- vdWSurfPt[n] = vdW of a sphere placed on the
                          nth point of the SAS3
###########################################*/
{
  int i,iat,jat,ires;
  double r2,r6,r12,vdWMin,vdWMax,ProbeVdWE,SumRad,SumRad6;

/* Loop over S&R surface points, place a probe sphere over each of them
   and calculate the vdW interaction of the probe with the receptor */

  vdWMin=99999999.;
  vdWMax=-99999999.;
  ProbeVdWE = 0.3;  /* The square root of 0.0903 (carbon) */

  for (ires=1;ires<=BSResN;ires++) {
    for (iat=ResFirstAtom[BSReNu[ires]];iat<=ResFirstAtom[BSReNu[ires]]+
         NAtom_per_Res[BSReNu[ires]]-1;iat++) {
      for (i=pointsrf_re[iat];i<=pointsrf_re[iat]+nsurf_re[iat] - 1 ; i++) {
        for (jat=1;jat<=ReAtNu;jat++) {
          r2 = (surfpt_re[i].x - ReCoor[jat][1])*
               (surfpt_re[i].x - ReCoor[jat][1]) +
               (surfpt_re[i].y - ReCoor[jat][2])*
               (surfpt_re[i].y - ReCoor[jat][2]) +
               (surfpt_re[i].z - ReCoor[jat][3])*
               (surfpt_re[i].z - ReCoor[jat][3]);
          r6 = r2*r2*r2;
          r12 = r6*r6;
          SumRad=ReVdWR[jat]+Sphere_apol;
          SumRad6=SumRad*SumRad*SumRad*SumRad*SumRad*SumRad;
          vdWSurfPt[i] += ReVdWE_sr[jat]*ProbeVdWE *
                          SumRad6 * ( (SumRad6/r12)-(2/r6) );
        }
        vdWMin = (vdWSurfPt[i] < vdWMin) ? vdWSurfPt[i] : vdWMin;
        vdWMax = (vdWSurfPt[i] > vdWMax) ? vdWSurfPt[i] : vdWMax;
      }
    }
  }

  if (ires == BSResN+1 )
    return 1;
  else
    return 0;
}

int Desol_Surf(int ReAtNu,struct point Min,double Sphere_apol,
               double GrSiSo,int NGridx,int NGridy,int NGridz,int BSResN,
               int *BSReNu,int *NAtom_per_Res,int *ResFirstAtom,
               struct point *surfpt_re,int *nsurf_re,int *pointsrf_re,
               double ***DeltaPrDeso,double *ReSurf_apol,double NCutapolRatio,/*int NCutapol,*/
               double ScaleDeso,int *isurf_More_apol,
               int *PNapol_Vect_re,int *Nsurfpt_re_apol_BS,
	             int *iatsrf_re_apol,int distrPointBSNumb,double **distrPointBS,
	             double angle_rmin,double angle_rmax,double mult_fact_rmin,
	             double mult_fact_rmax,double **ReCoor,FILE *FPaOut)
/*##########################################
Calculate the desolvation operated by a sphere rolled over the rec surface
###########################################*/

/*##########################################
int ReAtNu -------------- Tot # rec atoms
struct point Min -------- Min coor of the grid box
double Sphere_apol ------ Dimension of the probe sphere used to generate SAS3
double GrSiSo ----------- Size of the 3D grid used for cont. electrostatics
int  NGridx ------------- Tot # of grid points along x
int  NGridy ------------- Tot # of grid points along y
int  NGridz ------------- Tot # of grid points along z
int BSResN -------------- Tot number of rec residues in the Binding Site
int *BSReNu ------------- # rec residue in the Binding Site
int *NAtom_per_Res ------ NAtom_per_Res[n] = # of atoms of residue n
int *ResFirstAtom ------- ResFirstAtom[n] = atom number of the first atom
                          belonging to residue n
struct point *surfpt_re -  Coor of points over SAS3
int *nsurf_re ----------- nsurf_re[n] = amount of SAS3 surface points
                          generated by atom n
int *pointsrf_re -------- pointsrf_re[n] = first rec SAS3 point (in the
                         list surfpt_re) that is generated by atom n
double ***DeltaPrDeso --- Elec rec desolvation due to the occupation of a grid point
double *ReSurf_apol ----- ReSurf_apol[n] = vdW + Desolvation of a sphere
                          placed on the nth point of the SAS3
int NCutapol ------------ Amount of most hydrophobic rec points used to seed
                         apolar fragments
int NCutapolRatio ------ Ratio of most hydrophobic rec points used to seed
                         apolar fragments
double ScaleDeso -------- Scaling factor for vdW interactions in the calculation
                         of ReSurf_apol[n]
int *isurf_More_apol ---- isurf_More_apol[n] = rec SAS3 point with the nth
                          highest hydrophobicity
int *PNapol_Vect_re ----- # of most hydrophobic points used to seed the apolar
                          vectors. It is an input value, but if the value chosen
                          is smaller than the tot # of hydrophobic points, it is
                          fixed to the tot #
int *Nsurfpt_re_apol_BS - Tot # of points over the rec SAS3
###########################################
   distrPointBSNumb  number of points in the binding site used for the
                     reduction of rec polar and apolar vectors using an
		     angle criterion
   distrPointBS  distribution of points in the binding site used for
                 the reduction of rec polar and apolar vectors using an
		 angle criterion
   angle_rmin  angle cutoff for the reduction of vectors
   angle_rmin_rad  angle cutoff for the reduction of vectors in radian
   angle_rmax  angle cutoff for the reduction of vectors
   angle_rmax_rad  angle cutoff for the reduction of vectors in radian
   mult_fact_rmin  multiplicative factor for the reduction of vectors
   mult_fact_rmax  multiplicative factor for the reduction of vectors
   keptVect_angle  apolar vectors on receptor (0 if rejected, 1 if kept)
                   in case of angle cutoff criterion
   keptVect_angle_numb  number of rec apolar vectors kept after angle
                        cutoff criterion applied
   vect_angle_stored  array for storing angle for each vector */
{
  int i,j,iat,ires,ix,iy,iz,ixmin,iymin,izmin,ixmax,iymax,izmax,*RankArr,
      *BSnumbering,Ntot,numbApolVectBS,closestPoint,
      *keptVect_angle,keptVect_angle_numb,numbApolVectBSCounter,
      helpCount,maxCount,NCutapol;/*new NCutapol == internal variable*/
  double Sphere_apol2,xtemp,x2temp,ytemp,xy2temp,ztemp,r2,DesoMin,
         DesoMax,DesoMid,scale,*ReSurf_apolBS,deso;
  double minSqDist,maxSqDist,vecPointSqDist,
	distSqClosest,vecPointAngle,minSqDistCutoff,maxSqDistCutoff,
	angleCutoff,angle_rmin_rad,angle_rmax_rad,*vect_angle_stored;
  FILE *FilePa;

  Sphere_apol2 = (Sphere_apol+.3) * (Sphere_apol+.3);

/* Loop over S&R surface points, place a probe sphere over each of them
   and calculate the desolvation of the protein if such a sphere would be
   removed */

  *Nsurfpt_re_apol_BS=0;
  *PNapol_Vect_re=0;
  DesoMin=999999999999.;
  DesoMax=-99999999999.;
  Ntot=0;
  for (ires=1;ires<=BSResN;ires++) {
    for (iat=ResFirstAtom[BSReNu[ires]];iat<=ResFirstAtom[BSReNu[ires]]+
         NAtom_per_Res[BSReNu[ires]]-1;iat++) {
      for (i=pointsrf_re[iat];i<=pointsrf_re[iat]+nsurf_re[iat] - 1 ;
           i++,Ntot++) {
        (*Nsurfpt_re_apol_BS)++;

/* Calculate the extremes, in the grid frame, of the cube containing the
   probe sphere. They are integer! */

        ixmin = (surfpt_re[i].x - Sphere_apol - Min.x) / GrSiSo + 1;
        iymin = (surfpt_re[i].y - Sphere_apol - Min.y) / GrSiSo + 1;
        izmin = (surfpt_re[i].z - Sphere_apol - Min.z) / GrSiSo + 1;
        ixmax = (surfpt_re[i].x + Sphere_apol - Min.x) / GrSiSo + 1;
        iymax = (surfpt_re[i].y + Sphere_apol - Min.y) / GrSiSo + 1;
        izmax = (surfpt_re[i].z + Sphere_apol - Min.z) / GrSiSo + 1;
        ixmin = (ixmin > 1) ? ixmin : 1;
        iymin = (iymin > 1) ? iymin : 1;
        izmin = (izmin > 1) ? izmin : 1;
        ixmax = (ixmax < NGridx) ? ixmax : NGridx;
        iymax = (iymax < NGridy) ? iymax : NGridy;
        izmax = (izmax < NGridz) ? izmax : NGridz;

/* xtemp,ytemp,ztemp are cartesian coordinates of the grid pt ix,iy,iz,
   in a frame with the origin in the suface point surfpt[i] */


        xtemp = GrSiSo * (ixmin - 1.5) + Min.x - surfpt_re[i].x;
        for (ix=ixmin;ix<=ixmax;ix++) {
          xtemp += GrSiSo;
          x2temp = xtemp*xtemp;
          ytemp = GrSiSo * (iymin - 1.5) + Min.y - surfpt_re[i].y;
          for (iy=iymin;iy<=iymax;iy++) {
            ytemp += GrSiSo;
            xy2temp = ytemp*ytemp + x2temp;

/*  Check if we are already out of the sphere */

            if (Sphere_apol2 > xy2temp) {
              ztemp = GrSiSo * (izmin - 1.5) + Min.z - surfpt_re[i].z;
              for (iz=izmin;iz<=izmax;iz++) {
                ztemp += GrSiSo;
                r2 = ztemp * ztemp + xy2temp;

/* Check if our grid point is inside the sphere. If yes change GridMat */

                if (Sphere_apol2 > r2)
                  ReSurf_apol[i] += ScaleDeso * DeltaPrDeso[ix][iy][iz];
              }
            }
          }
        }
        DesoMin = (ReSurf_apol[i] < DesoMin) ? ReSurf_apol[i] : DesoMin;
        DesoMax = (ReSurf_apol[i] > DesoMax) ? ReSurf_apol[i] : DesoMax;
      }
    }
  }

/* Write the SAS points of the binding site in FilePa */
  #ifdef ENABLE_MPI
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  if(myrank == MASTERRANK){
  #endif
  FilePa=fopen("./outputs/apolar_rec.mol2","w");
  fprintf(FilePa,"# TRIPOS MOL2 file generated by WITNOTP\n");
  fprintf(FilePa,"\n");
  fprintf(FilePa,"@<TRIPOS>MOLECULE\n");
  fprintf(FilePa,"apolar_rec\n");
  fprintf(FilePa,"%d 0 0 0 0\n",Ntot);
  fprintf(FilePa,"****\n");
  fprintf(FilePa,"USER_CHARGES\n");
  fprintf(FilePa,"INVALID_CHARGES\n");
  fprintf(FilePa,"\n");
  fprintf(FilePa,"@<TRIPOS>ATOM\n");
  j=1;
  for (ires=1;ires<=BSResN;ires++)
    for (iat=ResFirstAtom[BSReNu[ires]];iat<=ResFirstAtom[BSReNu[ires]]+
         NAtom_per_Res[BSReNu[ires]]-1;iat++)
      for (i=pointsrf_re[iat];i<=pointsrf_re[iat]+nsurf_re[iat] - 1 ; i++,j++)

        fprintf(FilePa,"%-5d %-3d%11.4f%11.4f%11.4f ***** 1    ****  0.000\n",
                j,j,surfpt_re[i].x,surfpt_re[i].y,surfpt_re[i].z);

  fclose(FilePa);
  #ifdef ENABLE_MPI
  }
  #endif




/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */



/* Compute the number of apolar vectors on the binding site surface */
  numbApolVectBS=0;
  for (ires=1;ires<=BSResN;ires++)
    for (iat=ResFirstAtom[BSReNu[ires]];iat<=ResFirstAtom[BSReNu[ires]]+
         NAtom_per_Res[BSReNu[ires]]-1;iat++)
      for (i=pointsrf_re[iat];i<=pointsrf_re[iat]+nsurf_re[iat]-1;i++)
        numbApolVectBS++;

/* Initialisation */
  vect_angle_stored=dvector(1,numbApolVectBS);
  keptVect_angle=ivector(1,numbApolVectBS);
  if (distrPointBSNumb>0)
  {
    for (i=1;i<=numbApolVectBS;i++)
      keptVect_angle[i]=0;
    keptVect_angle_numb=0;
  }
  else
  {
    for (i=1;i<=numbApolVectBS;i++)
      keptVect_angle[i]=1;
    keptVect_angle_numb=numbApolVectBS;
  }
  angle_rmin_rad=angle_rmin*3.1415927/180.0;
  angle_rmax_rad=angle_rmax*3.1415927/180.0;



  if (distrPointBSNumb>0)
  {
/* -------------------------------------- */
/*  reduction of rec apolar vectors using  */
/*  angle criterion                       */
/* -------------------------------------- */
/*
   The minimal (minDist) and maximal (maxDist) distances
   between the vectors and the points in the binding site
   (as defined in the SEED input file) are evaluated.
   A vector is discarded if the angle between the vector
   and its closest point in the binding site is larger than
   a cutoff angle value.
   The cutoff angle value follows the following distribution:
   -  angle_rmin  if distance <= (mult_fact_rmin*minDist)
   -  angle_rmax  if distance >= (mult_fact_rmax*maxDist)
   -  linear dependence (range between angle_rmin and angle_rmax)
      for other distances
*/

/* find minimal and maximal distances between vectors and
   points in the binding site */
    minSqDist=1000000.0;
    maxSqDist=0.0;

    for (ires=1;ires<=BSResN;ires++)
      for (iat=ResFirstAtom[BSReNu[ires]];iat<=ResFirstAtom[BSReNu[ires]]+
           NAtom_per_Res[BSReNu[ires]]-1;iat++)
        for (i=pointsrf_re[iat];i<=pointsrf_re[iat]+nsurf_re[iat]-1;i++)
        {

          for (j=1;j<=distrPointBSNumb;j++)
          {

	    vecPointSqDist=DistSq(surfpt_re[i].x,
	                          surfpt_re[i].y,
				  surfpt_re[i].z,
                                  distrPointBS[j][1],distrPointBS[j][2],
		       	          distrPointBS[j][3]);

            if (vecPointSqDist<minSqDist)
	      minSqDist=vecPointSqDist;
  	    if (vecPointSqDist>maxSqDist)
              maxSqDist=vecPointSqDist;

          }

        }

    minSqDistCutoff=(mult_fact_rmin*mult_fact_rmin)*minSqDist;
    maxSqDistCutoff=(mult_fact_rmax*mult_fact_rmax)*maxSqDist;

    fprintf(FPaOut,"Minimal and maximal distances between receptor ");
    fprintf(FPaOut,"apolar vectors \nand points in the binding site: ");
    fprintf(FPaOut,"%11.4f%11.4f\n\n",sqrtf(minSqDist),sqrtf(maxSqDist));

    if (minSqDistCutoff>=maxSqDistCutoff)
    {
      fprintf(FPaOut,"WARNING Parameters for reducing vectors ");
      fprintf(FPaOut,"(angle criterion) not appropriate\n");
      fclose(FPaOut);
      printf("Program exits\n");
      exit(0);
    }


/* compute angle for each vector with its closest point and
   discard vector if angle cutoff criterion not satisfied */
    numbApolVectBSCounter=0;
    for (ires=1;ires<=BSResN;ires++)
      for (iat=ResFirstAtom[BSReNu[ires]];iat<=ResFirstAtom[BSReNu[ires]]+
           NAtom_per_Res[BSReNu[ires]]-1;iat++)
        for (i=pointsrf_re[iat];i<=pointsrf_re[iat]+nsurf_re[iat]-1;i++)
        {

          numbApolVectBSCounter++;

/* find the closest point */
          closestPoint=-1;
          distSqClosest=1000000.0;

          for (j=1;j<=distrPointBSNumb;j++)
          {

	    vecPointSqDist=DistSq(surfpt_re[i].x,
	                          surfpt_re[i].y,
				  surfpt_re[i].z,
                                  distrPointBS[j][1],distrPointBS[j][2],
		       	          distrPointBS[j][3]);

            if (vecPointSqDist<distSqClosest)
	    {
	      distSqClosest=vecPointSqDist;
	      closestPoint=j;
	    }

          }

/* compute the angle */
	  vecPointAngle=PlaAng(ReCoor[iatsrf_re_apol[i]][1],
	                       ReCoor[iatsrf_re_apol[i]][2],
			       ReCoor[iatsrf_re_apol[i]][3],
			       surfpt_re[i].x,surfpt_re[i].y,surfpt_re[i].z,
                               distrPointBS[closestPoint][1],
			       distrPointBS[closestPoint][2],
		       	       distrPointBS[closestPoint][3]);

          vect_angle_stored[numbApolVectBSCounter]=vecPointAngle/3.1415927*180.0;

/* find angle cutoff and apply criterion to keep or discard the vector */
          if (distSqClosest<=minSqDistCutoff)
            angleCutoff=angle_rmin_rad;
          else if (distSqClosest>=maxSqDistCutoff)
            angleCutoff=angle_rmax_rad;
          else
            angleCutoff=angle_rmin_rad+
	                ((angle_rmax_rad-angle_rmin_rad)/
	                (sqrt(maxSqDistCutoff)-sqrt(minSqDistCutoff)))*
                        (sqrt(distSqClosest)-sqrt(minSqDistCutoff));

          if (vecPointAngle<=angleCutoff)
          {
            keptVect_angle[numbApolVectBSCounter]=1;
 	    keptVect_angle_numb++;
          }

        }

/* Print the kept vectors after angle cutoff criterion in a file */
    #ifdef ENABLE_MPI
    if(myrank == MASTERRANK){
    #endif
    FilePa=fopen("./outputs/apolar_rec_reduc_angle.mol2","w");
    fprintf(FilePa,"# TRIPOS MOL2 file generated by WITNOTP\n");
    fprintf(FilePa,"\n");
    fprintf(FilePa,"@<TRIPOS>MOLECULE\n");
    fprintf(FilePa,"apolar_rec_reduc_angle\n");
    fprintf(FilePa,"%d 0 0 0 0\n",keptVect_angle_numb);
    fprintf(FilePa,"****\n");
    fprintf(FilePa,"USER_CHARGES\n");
    fprintf(FilePa,"INVALID_CHARGES\n");
    fprintf(FilePa,"\n");
    fprintf(FilePa,"@<TRIPOS>ATOM\n");

    j=0;
    numbApolVectBSCounter=0;
    for (ires=1;ires<=BSResN;ires++)
      for (iat=ResFirstAtom[BSReNu[ires]];iat<=ResFirstAtom[BSReNu[ires]]+
           NAtom_per_Res[BSReNu[ires]]-1;iat++)
        for (i=pointsrf_re[iat];i<=pointsrf_re[iat]+nsurf_re[iat]-1;i++)
        {
          numbApolVectBSCounter++;
          if (keptVect_angle[numbApolVectBSCounter])
          {
            j=j+1;
            fprintf(FilePa,"%-5d %-3d%11.4f%11.4f%11.4f ***** 1    ****  %11.4f\n",
                    j,j,surfpt_re[i].x,surfpt_re[i].y,surfpt_re[i].z,
	            vect_angle_stored[numbApolVectBSCounter]);
          }
        }
    fclose(FilePa);
    #ifdef ENABLE_MPI
    }
    #endif

  }



/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */





/* A slightly complicated scheme for selecting the NCutapol points with
   the lowest desolvation */

/* Fill the array ReSurf_apolBS that contains only the surface points in
   the binding site. (ReSurf_apol contains also the ones out of the BS)
   BSnumbering gives the correspondance between the numbering of the surf
   points in the BS and the one over all the protein */
  ReSurf_apolBS = dvector(1,*Nsurfpt_re_apol_BS);
  BSnumbering = ivector(1,*Nsurfpt_re_apol_BS);
  *Nsurfpt_re_apol_BS=0;
  for (ires=1;ires<=BSResN;ires++)
    for (iat=ResFirstAtom[BSReNu[ires]];iat<=ResFirstAtom[BSReNu[ires]]+
         NAtom_per_Res[BSReNu[ires]]-1;iat++)
      for (i=pointsrf_re[iat];i<=pointsrf_re[iat]+nsurf_re[iat] - 1 ; i++) {
        ReSurf_apolBS[++(*Nsurfpt_re_apol_BS)] = ReSurf_apol[i];
        BSnumbering[(*Nsurfpt_re_apol_BS)] = i;
      }

/* Sort the array ReSurf_apolBS */
  RankArr=ivector(1,*Nsurfpt_re_apol_BS);
  for (i=1;i<=*Nsurfpt_re_apol_BS;i++)
    RankArr[i] = i;
  DSort(*Nsurfpt_re_apol_BS,RankArr,ReSurf_apolBS);



/* This part has been replaced by the following part to
   take into account the angle criterion selection of the
   apolar vectors of the binding site */
/*
   Get the NCutapol more hydrophobic surf points
  *PNapol_Vect_re = 0;
  for (i=1;i<=NCutapol && i<=(*Nsurfpt_re_apol_BS);i++)
    isurf_More_apol[++(*PNapol_Vect_re)] = BSnumbering[RankArr[i]];

   If NCutapol is bigger than the # surf points in
   the BS (*Nsurfpt_re_apol_BS) it is lowered to this value
  if (NCutapol<*Nsurfpt_re_apol_BS)
    *PNapol_Vect_re = NCutapol;
  else
    *PNapol_Vect_re = *Nsurfpt_re_apol_BS;
*/


  /* calculate the maximal number of apolar points
     as a sum from NcutapolRatio(user-def.) and Nsurfpt_re_apol_BS
     (=Tot # of points over the rec SAS3)
  */
  if (distrPointBSNumb>0)
      NCutapol=(NCutapolRatio * (keptVect_angle_numb));/*dey new */
  else
      NCutapol=(NCutapolRatio * (*Nsurfpt_re_apol_BS));/*dey new*/


/* Keep only the apolar vectors with the lowest energies (vdW+desolvation)
   taking into account the angle cutoff criterion selection */
  if ((numbApolVectBS>=NCutapol)&&
      (keptVect_angle_numb>=NCutapol))
  {
    maxCount=NCutapol;
  }
  else
  {
    maxCount=keptVect_angle_numb;
  }

  helpCount=0;
  for (i=1;i<=numbApolVectBS;i++)
  {
    if ((keptVect_angle[RankArr[i]]==1)&&
        (helpCount<maxCount))
    {
      helpCount++;
      isurf_More_apol[helpCount]=BSnumbering[RankArr[i]];
    }
  }



  *PNapol_Vect_re=helpCount;





  scale = 50. / *Nsurfpt_re_apol_BS ;
  deso = 0.;
  DesoMid = (DesoMax - DesoMin)/2. + DesoMin;
  #ifdef ENABLE_MPI
  if(myrank == MASTERRANK){
  #endif
  FilePa=fopen("outputs/sas_apolar.pdb","w");
  fprintf(FilePa,"GRASP PDB FILE\n");
  fprintf(FilePa,"FORMAT NUMBER=1\n");
  for (i=1;i<=(*Nsurfpt_re_apol_BS);i++,deso+=scale)
    fprintf(FilePa,"ATOM %6d  H   SURF    1    %8.3f%8.3f%8.3f %6.5f %6.5f\n",
            i,surfpt_re[BSnumbering[RankArr[i]]].x,
            surfpt_re[BSnumbering[RankArr[i]]].y,
            surfpt_re[BSnumbering[RankArr[i]]].z,deso,
            ReSurf_apol[BSnumbering[RankArr[i]]]);

  j=i;
  for (iat=1;iat<=ReAtNu;iat++)
    for (i=pointsrf_re[iat];i<=pointsrf_re[iat]+nsurf_re[iat] - 1;i++,j++)
      if (ReSurf_apol[i] == 0.)
        fprintf(FilePa,"ATOM %6d  H   SURF    1    %8.3f%8.3f%8.3f %6.5f 25.0\n",
                j,surfpt_re[i].x,surfpt_re[i].y,surfpt_re[i].z,DesoMid);

  fprintf(FilePa,"END\n");
  fclose(FilePa);
  #ifdef ENABLE_MPI
  }
  #endif

  free_ivector(BSnumbering,1,*Nsurfpt_re_apol_BS);
  free_dvector(ReSurf_apolBS,1,*Nsurfpt_re_apol_BS);
  free_ivector(RankArr,1,*Nsurfpt_re_apol_BS);

  free_dvector(vect_angle_stored,1,numbApolVectBS);
  free_ivector(keptVect_angle,1,numbApolVectBS);

  if (ires == BSResN+1 )
    return 1;
  else
    return 0;
}

int Calc_D_uhbd(int ReAtNu,double *ReRad,struct point Min,
                struct point Max,double *RePaCh,double WaMoRa,int NPtSphere,
                double DielRe,double DielWa,double GrSiSo,double GrInSo,
                double pi4,int NGridx,int NGridy,int NGridz,int nxminBS,
                int nyminBS,int nzminBS,int nxmaxBS,int nymaxBS,
                int nzmaxBS,double *XGrid,double *YGrid,double *ZGrid,
                char ***GridMat,double UnitVol,double corr_re_desofd,
                char *RecFilPDB,char *FDexe,char *FDdir,char *DesoMapAcc,
                char *DesoMapFile,double ***DeltaPrDeso)
/*##########################################
Calculate the electric displacement from a finite difference approach (UHBD)
###########################################*/

/*##########################################
int ReAtNu -------------- Tot # rec atoms
double *ReRad ----------- Rec charge radii
struct point Min -------- Min coor of the grid box
struct point Max -------- Max coor of the grid box
double *RePaCh ----------- Rec partial charges
double WaMoRa ----------- Radius of the water molecule
int NPtSphere ----------- Amount of points placed over each atom to generate
                          the SAS (built in order to obtain the volume
                          enclosed by the MS)
double DielRe ----------- Dielectric constant of the rec (and of the fragments)
double DielWa ----------- Dielectric constant of the water
double GrSiSo ----------- Size of the 3D grid used for cont. electrostatics
double GrInSo ----------- Margin left along each dimension (positive and
                          negative) of the rec for building the 3D grid box
double pi4 -------------- 4 * greekpi
int  NGridx ------------- Tot # of grid points along x
int  NGridy ------------- Tot # of grid points along y
int  NGridz ------------- Tot # of grid points along z
int nxminBS ------------- grid point (along x) where the BS starts
int nyminBS ------------- grid point (along y) where the BS starts
int nzminBS ------------- grid point (along z) where the BS starts
int nxmaxBS ------------- grid point (along x) where the BS ends
int nymaxBS ------------- grid point (along y) where the BS ends
int nzmaxBS ------------- grid point (along z) where the BS ends
double *XGrid ----------- X coor of the grid points
double *YGrid ----------- Y coor of the grid points
double *ZGrid ----------- Z coor of the grid points
char ***GridMat --------- Matrix telling if a grid point is occupied by the
                          rec (o), empty (e), or if it belongs to the interface
                          between SAS and MS (s)
double UnitVol ---------- Volume of the grid element for cont. elec.
double corr_re_desofd --- Correction factor for slow rec elec desolvation
                          (fd approx.)
char *RecFilPDB --------- name of the PDB file of the rec
                          (by default is ./outputs/receptor_uhbd.pdb)
char *FDexe ------------- Name of the UHBD executable file
char *FDdir ------------- Directory where to put the outputs of UHBD
char *DesoMapAcc -------- Access to deolvation map file
                          (w = calculate and write / n = calculate / r = read)
char *DesoMapFile ------- Name of the desolvation map file
double ***DeltaPrDeso --- Elec rec desolvation due to the occupation of a grid point
###########################################*/
{
  int ix,iy,iz,iat,nxcoarse,nycoarse,nzcoarse,nxfine,nyfine,nzfine;
  double xcen,ycen,zcen,Kdesol,Dx,Dy,Dz;
  double ***phi,***eps;
  char UxStr[300];
  FILE *FilChk,*FileOut;
  int dummy = 0; //clangini

    Kdesol = - corr_re_desofd * UnitVol *
             ( 1./DielWa - 1./DielRe ) / (2. * pi4 * 332.0716);

/* Calculate center and dimensions of the grid */

    xcen = (XGrid[nxmaxBS]-XGrid[nxminBS])/2. + XGrid[nxminBS];
    ycen = (YGrid[nymaxBS]-YGrid[nyminBS])/2. + YGrid[nyminBS];
    zcen = (ZGrid[nzmaxBS]-ZGrid[nzminBS])/2. + ZGrid[nzminBS];
    nxcoarse = (Max.x - Min.x + 20.*2 - 2. * GrInSo) / 2. + 1;
    nycoarse = (Max.y - Min.y + 20.*2 - 2. * GrInSo) / 2. + 1;
    nzcoarse = (Max.z - Min.z + 20.*2 - 2. * GrInSo) / 2. + 1;
    nxfine = nxmaxBS-nxminBS+3;
    nyfine = nymaxBS-nyminBS+3;
    nzfine = nzmaxBS-nzminBS+3;

/*-------------Create the script for UHBD that will be used to calculate the
   x(y,z)-shifted potentials and dielectric grid--------------------*/

    sprintf(UxStr,"%s",FDdir); /* prev "%s\0" */

    if ((FilChk=fopen(UxStr,"r"))==NULL) {
      fclose(FilChk);
      sprintf(UxStr,"mkdir %s",FDdir); /* prev "mkdir %s\0" */
      dummy = system(UxStr);
    }
    else
      fclose(FilChk);

    sprintf(UxStr,"%s%s%s%s%s","echo read mol1 file \\\'",RecFilPDB,
          "\\\' pdb end \\\\n >",FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s","echo assign  so  =  1  end >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s","echo assign  gx  =  XCENTER   end >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s","echo assign  gy  =  YCENTER   end >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s","echo assign  gz  =  ZCENTER   end >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%f%s%s%s","echo assign  sdi = ",DielWa," end >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%f%s%s%s","echo assign  eps = ",DielRe," end >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s","echo assign  ion =   0.0  end >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s","echo assign strn =   0.0  end >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s","echo assign   te = 298.0  end >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s","echo assign mits = 600    end >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%f%s%s%s","echo assign wrad = ",WaMoRa," end >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%d%s%s%s","echo assign ndots = ",NPtSphere," end >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s","echo assign  bc  =  2    end >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s","echo assign  dh  =  0.0  end >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s%s%s","echo \\\\n stream \\\'",FDdir,
          "/par.str\\\' end \\\\n >>",FDdir,"/uhbd_script\0");
    dummy = system(UxStr);

    sprintf(UxStr,"%s%s%s","echo ! ------- Energy total ------- >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s","echo elec calc mol1 ions \\$ion rion \\$strn >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s","echo solver \\$so bcflag \\$bc dhrad \\$dh >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s","echo pdie \\$eps sdie \\$sdi grid 2.0 >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%d %d %d%s%s%s","echo dime ",nxcoarse,nycoarse,nzcoarse,
          " >>",FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s","echo gcen \\$gx \\$gy \\$gz maxits \\$mits nmap >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s","echo \\$wrad nsph \\$ndots temp \\$te nosmooth end >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);

    sprintf(UxStr,"%s%s%s","echo print elec phisave mol1 end \\\\n >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);

    sprintf(UxStr,"%s%s%s","echo ! now focussing first step and only one >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s","echo elec calc mol1 ions \\$ion rion \\$strn >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s","echo solver \\$so bcflag 4 samsrf pdie \\$eps >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%f%s%s%s","echo sdie \\$sdi grid ",GrSiSo," >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%d %d %d%s%s%s","echo dime ",nxfine,nyfine,nzfine," >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s","echo gcen \\$gx \\$gy \\$gz maxits \\$mits >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s","echo nmap \\$wrad nsph \\$ndots temp \\$te >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s","echo nosmooth end >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s","echo print elec phisave mol1 end \\\\n >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s%s%s","echo write grid EPSILON ascii file \\\'",
          FDdir,"/eps.dat\\\' noexclude end >>",FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s%s%s","echo write grid phi ascii file \\\'",
          FDdir,"/phi.dat\\\' noexclude end >>",FDdir,"/uhbd_script\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s%s%s","echo \\\\n stop >>",
          FDdir,"/uhbd_script\0");
    dummy = system(UxStr);

/*------------Create the parameter file for UHBD--------------------*/

    sprintf(UxStr,"%s%s",FDdir,"/par.str\0");
    FileOut=fopen(UxStr,"w");
    for (iat=1;iat<=ReAtNu;iat++) {
      fprintf(FileOut,"%s %d %f %s",
              "edit mol1 charge atnum",iat,RePaCh[iat],"end\n");
      fprintf(FileOut,"%s %d %f %s",
              "edit mol1 radius atnum",iat,ReRad[iat],"end\n");
    }
    fclose(FileOut);

/*------------Create the UHBD script for the x-shifted potential--------*/

    sprintf(UxStr,"%s%f%s%s%s%f%s%f%s%s%s",
            "sed s/\"XCENTER\"/\"",xcen + GrSiSo/2.,
            "\"/g ",FDdir,"/uhbd_script | sed s/\"YCENTER\"/\"",ycen,
            "\"/g | sed s/\"ZCENTER\"/\"",zcen,
            "\"/g | sed s/\"EPSILON\"/\"epsi\"/g > ",FDdir,"/uhbd_x.inp \0");
    dummy = system(UxStr);

/*------------Run UHBD for the x-shifted potential--------*/

    sprintf(UxStr,"%s%s",FDdir,"/phi.dat");

    if ((FilChk=fopen(UxStr,"r"))!=NULL) {
      fclose(FilChk);
      sprintf(UxStr,"%s %s%s","rm",FDdir,"/phi.dat");
      dummy = system(UxStr);
    }
    else
      fclose(FilChk);

    sprintf(UxStr,"%s%s",FDdir,"/eps.dat");

    if ((FilChk=fopen(UxStr,"r"))!=NULL) {
      fclose(FilChk);
      sprintf(UxStr,"%s %s%s","rm",FDdir,"/eps.dat");
      dummy = system(UxStr);
    }
    else
      fclose(FilChk);

    sprintf(UxStr,"%s<%s%s%s%s",FDexe,FDdir,"/uhbd_x.inp > ",
            FDdir,"/uhbd_x.out\0");
    dummy = system(UxStr);

/*------------Read the x-shifted potential and dielectric maps--------*/
    phi = d3tensor(nxminBS-1,nxmaxBS+1,nyminBS-1,nymaxBS+1,nzminBS-1,nzmaxBS+1);
    eps = d3tensor(nxminBS-1,nxmaxBS+1,nyminBS-1,nymaxBS+1,nzminBS-1,nzmaxBS+1);

/*-----------Read the x-shifted potential-----------*/
    printf("phii\n");
    sprintf(UxStr,"%s%s",FDdir,"/phi.dat\0");
    Read_uhbd_map(UxStr,nxminBS,nyminBS,nzminBS,nxmaxBS,nymaxBS,nzmaxBS,phi);

/*-----------Read the x-shifted dielectric constant-----------*/
    printf("epsii\n");
    sprintf(UxStr,"%s%s",FDdir,"/eps.dat\0");
    Read_uhbd_map(UxStr,nxminBS,nyminBS,nzminBS,nxmaxBS,nymaxBS,nzmaxBS,eps);

    for (ix=nxminBS;ix<=nxmaxBS;ix++)
      for (iy=nyminBS;iy<=nymaxBS;iy++)
        for (iz=nzminBS;iz<=nzmaxBS;iz++) {
          Dx = - eps[ix-1][iy][iz] *
                (phi[ix][iy][iz]-phi[ix-1][iy][iz]) / GrSiSo;
          if (GridMat[ix][iy][iz] != 'o')
            DeltaPrDeso[ix][iy][iz] += Kdesol * (Dx * Dx);
        }
/*------------Create the UHBD script for the y-shifted potential--------*/

    sprintf(UxStr,"%s%f%s%s%s%f%s%f%s%s%s",
            "sed s/\"XCENTER\"/\"",xcen,"\"/g ",FDdir,
            "/uhbd_script | sed s/\"YCENTER\"/\"",ycen + GrSiSo/2.,
            "\"/g | sed s/\"ZCENTER\"/\"",zcen,
            "\"/g | sed s/\"EPSILON\"/\"epsj\"/g > ",FDdir,"/uhbd_y.inp \0");
    dummy = system(UxStr);

/*------------Run UHBD for the y-shifted potential--------*/

    sprintf(UxStr,"%s %s%s %s%s","rm",FDdir,"/phi.dat",FDdir,"/eps.dat\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s<%s%s%s%s",FDexe,FDdir,"/uhbd_y.inp > ",
            FDdir,"/uhbd_y.out\0");
    dummy = system(UxStr);

/*------------Read the y-shifted potential and dielectric maps--------*/

/*-----------Read the y-shifted potential-----------*/
    printf("phij\n");
    sprintf(UxStr,"%s%s",FDdir,"/phi.dat\0");
    Read_uhbd_map(UxStr,nxminBS,nyminBS,nzminBS,nxmaxBS,nymaxBS,nzmaxBS,phi);

/*-----------Read the y-shifted dielectric constant-----------*/
    printf("epsjj\n");
    sprintf(UxStr,"%s%s",FDdir,"/eps.dat\0");
    Read_uhbd_map(UxStr,nxminBS,nyminBS,nzminBS,nxmaxBS,nymaxBS,nzmaxBS,eps);

    for (ix=nxminBS;ix<=nxmaxBS;ix++)
      for (iy=nyminBS;iy<=nymaxBS;iy++)
        for (iz=nzminBS;iz<=nzmaxBS;iz++) {
          Dy = - eps[ix][iy-1][iz] *
                (phi[ix][iy][iz]-phi[ix][iy-1][iz]) / GrSiSo;
          if (GridMat[ix][iy][iz] != 'o')
            DeltaPrDeso[ix][iy][iz] += Kdesol * ( Dy * Dy );
        }

/*------------Create the UHBD script for the z-shifted potential--------*/

    sprintf(UxStr,"%s%f%s%s%s%f%s%f%s%s%s",
            "sed s/\"XCENTER\"/\"",xcen,"\"/g ",FDdir,
            "/uhbd_script | sed s/\"YCENTER\"/\"",ycen,
            "\"/g | sed s/\"ZCENTER\"/\"",zcen + GrSiSo/2.,
            "\"/g | sed s/\"EPSILON\"/\"epsk\"/g > ",FDdir,"/uhbd_z.inp \0");
    dummy = system(UxStr);

/*------------Run UHBD for the z-shifted potential--------*/

    sprintf(UxStr,"%s %s%s %s%s","rm",FDdir,"/phi.dat",FDdir,"/eps.dat\0");
    dummy = system(UxStr);
    sprintf(UxStr,"%s<%s%s%s%s",FDexe,FDdir,"/uhbd_z.inp > ",
            FDdir,"/uhbd_z.out\0");
    dummy = system(UxStr);

/*------------Read the z-shifted potential and dielectric maps--------*/

/*-----------Read the z-shifted potential-----------*/
    printf("phik\n");
    sprintf(UxStr,"%s%s",FDdir,"/phi.dat\0");
    Read_uhbd_map(UxStr,nxminBS,nyminBS,nzminBS,nxmaxBS,nymaxBS,nzmaxBS,phi);

/*-----------Read the z-shifted dielectric constant-----------*/
    printf("epskk\n");
    sprintf(UxStr,"%s%s",FDdir,"/eps.dat\0");
    Read_uhbd_map(UxStr,nxminBS,nyminBS,nzminBS,nxmaxBS,nymaxBS,nzmaxBS,eps);

    for (ix=nxminBS;ix<=nxmaxBS;ix++)
      for (iy=nyminBS;iy<=nymaxBS;iy++)
        for (iz=nzminBS;iz<=nzmaxBS;iz++) {
          Dz = - eps[ix][iy][iz-1] *
                (phi[ix][iy][iz]-phi[ix][iy][iz-1]) / GrSiSo;
          if (GridMat[ix][iy][iz] != 'o')
            DeltaPrDeso[ix][iy][iz] += Kdesol * ( Dz * Dz );
        }

    sprintf(UxStr,"%s %s%s %s%s","rm",FDdir,"/phi.dat",FDdir,"/eps.dat\0");
    dummy = system(UxStr);


    free_d3tensor(phi,nxminBS-1,nxmaxBS+1,nyminBS-1,nymaxBS+1,
                  nzminBS-1,nzmaxBS+1);
    free_d3tensor(eps,nxminBS-1,nxmaxBS+1,nyminBS-1,nymaxBS+1,
                  nzminBS-1,nzmaxBS+1);


  if (ix == nxmaxBS+1)
    return 1;
  else
    return 0;
}

int Read_uhbd_map(char *MapFile,int ibeginx,int ibeginy,int ibeginz,
                  int iendx,int iendy,int iendz,double ***arr)
/*#####################################################################
Read a map (eps or phi) created by UHBD (peculiar format...)
#######################################################################*/

/*#####################################################################
char *MapFile ------- name of the file to be read
int ibeginxS -------- grid point (along x) where the BS starts
int ibeginyS -------- grid point (along y) where the BS starts
int ibeginzS -------- grid point (along z) where the BS starts
int iendx  S -------- grid point (along x) where the BS ends
int iendy  S -------- grid point (along y) where the BS ends
int iendz  S -------- grid point (along z) where the BS ends
double ***arr -------- quantity that is read (eps or phi)
######################################################################*/
{
  int ix,iy,iz,nleft;
  char StrLin[_STRLENGTH];
  FILE *ReaFile;

  ReaFile=fopen(MapFile,"r");
  fgets_wrapper(StrLin,_STRLENGTH,ReaFile);
  fgets_wrapper(StrLin,_STRLENGTH,ReaFile);
  fgets_wrapper(StrLin,_STRLENGTH,ReaFile);
  fgets_wrapper(StrLin,_STRLENGTH,ReaFile);
  fgets_wrapper(StrLin,_STRLENGTH,ReaFile);
  for (iz=ibeginz-1;iz<=iendz+1;iz++) {
    fgets_wrapper(StrLin,_STRLENGTH,ReaFile);
    nleft = ibeginx-1;
    for (iy=ibeginy-1;iy<=iendy+1;iy++) {
      for (ix=nleft;ix<=iendx+1;ix+=6) {
        fgets_wrapper(StrLin,_STRLENGTH,ReaFile);
        if (ix+5 <= iendx+1) {
          nleft = ibeginx-1;
          sscanf(StrLin,"%lf%lf%lf%lf%lf%lf",&arr[ix][iy][iz],&arr[ix+1][iy][iz],
                                       &arr[ix+2][iy][iz],&arr[ix+3][iy][iz],
                                       &arr[ix+4][iy][iz],&arr[ix+5][iy][iz]);
        }
        else if ( iendx+1 == ix+4 ) {
          nleft = ibeginx;
          if ( iy < iendy+1)
            sscanf(StrLin,"%lf%lf%lf%lf%lf%lf",&arr[ix][iy][iz],&arr[ix+1][iy][iz],
                                       &arr[ix+2][iy][iz],&arr[ix+3][iy][iz],
                                       &arr[ix+4][iy][iz],
                                       &arr[ibeginx-1][iy+1][iz]);
          else
            sscanf(StrLin,"%lf%lf%lf%lf%lf",&arr[ix][iy][iz],&arr[ix+1][iy][iz],
                                       &arr[ix+2][iy][iz],&arr[ix+3][iy][iz],
                                       &arr[ix+4][iy][iz]);
        }
        else if ( iendx+1 == ix+3 ) {
          nleft = ibeginx+1;
          if ( iy < iendy+1)
            sscanf(StrLin,"%lf%lf%lf%lf%lf%lf",&arr[ix][iy][iz],&arr[ix+1][iy][iz],
                                       &arr[ix+2][iy][iz],&arr[ix+3][iy][iz],
                                       &arr[ibeginx-1][iy+1][iz],
                                       &arr[ibeginx][iy+1][iz]);
          else
            sscanf(StrLin,"%lf%lf%lf%lf",&arr[ix][iy][iz],&arr[ix+1][iy][iz],
                                       &arr[ix+2][iy][iz],&arr[ix+3][iy][iz]);
        }
        else if ( iendx+1 == ix+2 )  {
          nleft = ibeginx+2;
          if ( iy < iendy+1)
            sscanf(StrLin,"%lf%lf%lf%lf%lf%lf",&arr[ix][iy][iz],&arr[ix+1][iy][iz],
                                       &arr[ix+2][iy][iz],
                                       &arr[ibeginx-1][iy+1][iz],
                                       &arr[ibeginx][iy+1][iz],
                                       &arr[ibeginx+1][iy+1][iz]);
          else
            sscanf(StrLin,"%lf%lf%lf",&arr[ix][iy][iz],&arr[ix+1][iy][iz],
                                       &arr[ix+2][iy][iz]);
        }
        else if ( iendx+1 == ix+1 ) {
          nleft = ibeginx+3;
          if ( iy < iendy+1)
            sscanf(StrLin,"%lf%lf%lf%lf%lf%lf",&arr[ix][iy][iz],&arr[ix+1][iy][iz],
                                       &arr[ibeginx-1][iy+1][iz],
                                       &arr[ibeginx][iy+1][iz],
                                       &arr[ibeginx+1][iy+1][iz],
                                       &arr[ibeginx+2][iy+1][iz]);
          else
            sscanf(StrLin,"%lf%lf",&arr[ix][iy][iz],&arr[ix+1][iy][iz]);
        }
        else if ( iendx+1 == ix ) {
          nleft = ibeginx+4;
          if ( iy < iendy+1)
            sscanf(StrLin,"%lf%lf%lf%lf%lf%lf",&arr[ix][iy][iz],
                                       &arr[ibeginx-1][iy+1][iz],
                                       &arr[ibeginx][iy+1][iz],
                                       &arr[ibeginx+1][iy+1][iz],
                                       &arr[ibeginx+2][iy+1][iz],
                                       &arr[ibeginx+3][iy+1][iz]);
          else
            sscanf(StrLin,"%lf",&arr[ix][iy][iz]);
        }
      }
    }
  }
  fclose(ReaFile);

  return 1;
}

int Calc_D_Coul(int ReAtNu,double **ReCoor,double *ReRad2,struct point Min,
                struct point Max,double *RePaCh,double DielRe,double DielWa,
                double GrSiSo,double pi4,int nxminBS,int nyminBS,
                int nzminBS,int nxmaxBS,int nymaxBS,int nzmaxBS,
                double *XGrid,double *YGrid,double *ZGrid,char ***GridMat,
                double UnitVol,double corr_re_desoco,double ***DeltaPrDeso)
/*##########################################
Calculate the electric displacement with the Coulomb approximation
###########################################*/

/*##########################################
int ReAtNu -------------- Tot # rec (frag) atoms
double **ReCoor ---------- Rec (frag) coordinates
double *ReRad2 ---------- (Rec (frag) charge radii)^2
struct point Min -------- Min coor of the grid box
struct point Max -------- Max coor of the grid box
double *RePaCh ----------- Rec (frag) partial charges
double DielRe ----------- Dielectric constant of the rec (and of the fragments)
double DielWa ----------- Dielectric constant of the water
double GrSiSo ----------- Size of the 3D grid used for cont. electrostatics
double pi4 -------------- 4 * greekpi
int nxminBS ------------- starting grid point (along x) to calculate desolvation
int nyminBS ------------- starting grid point (along y) to calculate desolvation
int nzminBS ------------- starting grid point (along z) to calculate desolvation
int nxmaxBS ------------- ending grid point (along x) to calculate deoslvation
int nymaxBS ------------- ending grid point (along y) to calculate deoslvation
int nzmaxBS ------------- ending grid point (along z) to calculate deoslvation
double *XGrid ----------- X coor of the grid points
double *YGrid ----------- Y coor of the grid points
double *ZGrid ----------- Z coor of the grid points
char ***GridMat --------- Matrix telling if a grid point is occupied by the (o),
                          empty (e), or if it belongs to the interface
                          between SAS and MS (s)
double UnitVol ---------- Volume of the grid element for cont. elec.
double corr_re_desoco --- Correction factor for slow rec (frag) elec desolvation
                          (Coulomb approx.)
double ***DeltaPrDeso --- Elec rec (frag) desolvation due to the occupation
                          of a grid point
###########################################*/
{
  int ix,iy,iz,ixc,iyc,izc,ixmin,iymin,izmin,ixmax,iymax,izmax,iat,
     /*unused variables : nx,ny,nz,*/ boxL;
  double r,r2,vx,vy,vz,cutoff,cutoff2,***Dx,***Dy,***Dz,Kdesol;
/*unused variables :  char StrLin[60]; */
/*unused variables : FILE *WriFile,*ReaFile;*/

  Dx = d3tensor(nxminBS,nxmaxBS,nyminBS,nymaxBS,nzminBS,nzmaxBS);
  Dy = d3tensor(nxminBS,nxmaxBS,nyminBS,nymaxBS,nzminBS,nzmaxBS);
  Dz = d3tensor(nxminBS,nxmaxBS,nyminBS,nymaxBS,nzminBS,nzmaxBS);
  for (ix=nxminBS;ix<=nxmaxBS;ix++)
    for (iy=nyminBS;iy<=nymaxBS;iy++)
      for (iz=nzminBS;iz<=nzmaxBS;iz++) {
        Dx[ix][iy][iz] = 0.0000000000000000;
        Dy[ix][iy][iz] = 0.0000000000000000;
        Dz[ix][iy][iz] = 0.0000000000000000;
      }

/* Calculate the electric displacement generated by all the charges
   on every BS grid point. The electric displacement of a charge is added
   only if the distance between the charge and the grid point is lower than
   a "cutoff" distance (much faster and still very precise) */
  cutoff = 11.;
  cutoff2 = cutoff * cutoff;
  boxL = (cutoff + 0.00001) / GrSiSo + 1;
  for (iat=1;iat<=ReAtNu;iat++) {
    if (RePaCh[iat] != 0. ) {
      ixc = (ReCoor[iat][1] - Min.x) / GrSiSo + 1;
      iyc = (ReCoor[iat][2] - Min.y) / GrSiSo + 1;
      izc = (ReCoor[iat][3] - Min.z) / GrSiSo + 1;
      ixmin = ixc-boxL;
      iymin = iyc-boxL;
      izmin = izc-boxL;
      ixmax = ixc+boxL;
      iymax = iyc+boxL;
      izmax = izc+boxL;
      ixmin = (ixmin > nxminBS) ? ixmin : nxminBS;
      iymin = (iymin > nyminBS) ? iymin : nyminBS;
      izmin = (izmin > nzminBS) ? izmin : nzminBS;
      ixmax = (ixmax < nxmaxBS) ? ixmax : nxmaxBS;
      iymax = (iymax < nymaxBS) ? iymax : nymaxBS;
      izmax = (izmax < nzmaxBS) ? izmax : nzmaxBS;
      for (ix=ixmin;ix<=ixmax;ix++) {
        for (iy=iymin;iy<=iymax;iy++) {
          for (iz=izmin;iz<=izmax;iz++) {
            if (GridMat[ix][iy][iz] != 'o') {
              r2 = (XGrid[ix]-ReCoor[iat][1]) *
                   (XGrid[ix]-ReCoor[iat][1]) +
                   (YGrid[iy]-ReCoor[iat][2]) *
                   (YGrid[iy]-ReCoor[iat][2]) +
                   (ZGrid[iz]-ReCoor[iat][3]) *
                   (ZGrid[iz]-ReCoor[iat][3]);
              if ( r2 > ReRad2[iat] && r2 < cutoff2) {
                r = sqrt(r2);
                vx = (XGrid[ix]-ReCoor[iat][1]) / r;
                vy = (YGrid[iy]-ReCoor[iat][2]) / r;
                vz = (ZGrid[iz]-ReCoor[iat][3]) / r;
                Dx[ix][iy][iz] += RePaCh[iat] * vx / r2;
                Dy[ix][iy][iz] += RePaCh[iat] * vy / r2;
                Dz[ix][iy][iz] += RePaCh[iat] * vz / r2;
              }
            }
          }
        }
      }
    }
  }

/* From the electric displacement calculate the desolvation
   of every grid point */
  Kdesol = - corr_re_desoco * UnitVol * 332.0716 *
           ( 1./DielWa - 1./DielRe ) / (2. * pi4);
  for (ix=nxminBS;ix<=nxmaxBS;ix++)
    for (iy=nyminBS;iy<=nymaxBS;iy++)
      for (iz=nzminBS;iz<=nzmaxBS;iz++)
        if (GridMat[ix][iy][iz] != 'o')
          DeltaPrDeso[ix][iy][iz] = Kdesol *
                                  ( Dx[ix][iy][iz] * Dx[ix][iy][iz] +
                                    Dy[ix][iy][iz] * Dy[ix][iy][iz] +
                                    Dz[ix][iy][iz] * Dz[ix][iy][iz] ); // eq. (5) of SEED3.3.6 manual. clangini
/* Free memory */
  free_d3tensor(Dx,nxminBS,nxmaxBS,nyminBS,nymaxBS,nzminBS,nzmaxBS);
  free_d3tensor(Dy,nxminBS,nxmaxBS,nyminBS,nymaxBS,nzminBS,nzmaxBS);
  free_d3tensor(Dz,nxminBS,nxmaxBS,nyminBS,nymaxBS,nzminBS,nzmaxBS);

  if (ix == nxmaxBS+1)
    return 1;
  else
    return 0;
}

int Read_Write_Desol_Map(double WaMoRa,int NPtSphere,double DielRe,
                         double DielWa,double GrSiSo,double GrInSo,
                         int NGridx,int NGridy,int NGridz,int nxminBS,
                         int nyminBS,int nzminBS,int nxmaxBS,int nymaxBS,
                         int nzmaxBS,double *XGrid,double *YGrid,
                         double *ZGrid,char ***GridMat,double UnitVol,
                         char *ReDesoAlg,char *DesoMapAcc,
                         char *DesoMapFile,double ***DeltaPrDeso)
/*##########################################
Read or write a desolvation map
###########################################*/

/*##########################################
double WaMoRa ----------- Radius of the water molecule
int NPtSphere ----------- Amount of points placed over each atom to generate
                          the SAS (built in order to obtain the volume
                          enclosed by the MS)
double DielRe ----------- Dielectric constant of the rec (and of the fragments)
double DielWa ----------- Dielectric constant of the water
double GrSiSo ----------- Size of the 3D grid used for cont. electrostatics
double GrInSo ----------- Margin left along each dimension (positive and
                          negative) of the rec for building the 3D grid box
int  NGridx ------------- Tot # of grid points along x
int  NGridy ------------- Tot # of grid points along y
int  NGridz ------------- Tot # of grid points along z
int nxminBS ------------- grid point (along x) where the BS starts
int nyminBS ------------- grid point (along y) where the BS starts
int nzminBS ------------- grid point (along z) where the BS starts
int nxmaxBS ------------- grid point (along x) where the BS ends
int nymaxBS ------------- grid point (along y) where the BS ends
int nzmaxBS ------------- grid point (along z) where the BS ends
double *XGrid ----------- X coor of the grid points
double *YGrid ----------- Y coor of the grid points
double *ZGrid ----------- Z coor of the grid points
char ***GridMat --------- Matrix telling if a grid point is occupied by the
                          rec (o), empty (e), or if it belongs to the interface
                          between SAS and MS (s)
double UnitVol ---------- Volume of the grid element for cont. elec.
char *ReDesoAlg --------- Type of algorithm to calculate rec desolvation
                         (fd = finite diff. / co = coulomb approx.)
char *DesoMapAcc -------- Access to deolvation map file
                          (w = calculate and write / n = calculate / r = read)
char *DesoMapFile ------- Name of the desolvation map file
double ***DeltaPrDeso --- Elec desolvation due to the occupation of a grid point
###########################################*/
{
  int ix,iy,iz,NPt,Nx,Ny,Nz,nminx,nminy,nminz,nmaxx,nmaxy,nmaxz;
  double Wat,Size,Incr,epsre,epswa,Xfirst,Yfirst,Zfirst;
  char StrLin[_STRLENGTH],Alg[3];
  FILE *WriFile,*ReaFile;

  if (DesoMapAcc[0] == 'w') {
    #ifdef ENABLE_MPI
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    if (myrank == MASTERRANK){
    #endif
/* Write the D vector on the output file */
#ifndef NOWRITE
    WriFile=fopen(DesoMapFile,"w");
    fprintf(WriFile,"%f %d %f %f %f %f %s\n",WaMoRa,NPtSphere,GrSiSo,GrInSo,
                                          DielRe,DielWa,ReDesoAlg);
    fprintf(WriFile,"%d %d %d %f %f %f\n",NGridx,NGridy,NGridz,
                                          XGrid[1],YGrid[1],ZGrid[1]);
    fprintf(WriFile,"%d %d %d %d %d %d\n",nxminBS,nyminBS,nzminBS,
                                          nxmaxBS,nymaxBS,nzmaxBS);
    for (ix=nxminBS;ix<=nxmaxBS;ix++)
      for (iy=nyminBS;iy<=nymaxBS;iy++)
        for (iz=nzminBS;iz<=nzmaxBS;iz++)
          if (GridMat[ix][iy][iz] != 'o' )
            fprintf(WriFile,"%.12f\n",DeltaPrDeso[ix][iy][iz]);
    fclose(WriFile);
#endif
    #ifdef ENABLE_MPI
    }
    #endif
  }
  else if (DesoMapAcc[0] == 'r') {
/* Read the D vector from the input file */
    ReaFile=fopen(DesoMapFile,"r");
    fgets_wrapper(StrLin,_STRLENGTH,ReaFile);
    sscanf(StrLin,"%lf%d%lf%lf%lf%lf%s",&Wat,&NPt,&Size,&Incr,
                                        &epsre,&epswa,Alg);
    fgets_wrapper(StrLin,_STRLENGTH,ReaFile);
    sscanf(StrLin,"%d%d%d%lf%lf%lf",&Nx,&Ny,&Nz,&Xfirst,&Yfirst,&Zfirst);
    fgets_wrapper(StrLin,_STRLENGTH,ReaFile);
    sscanf(StrLin,"%d%d%d%d%d%d",&nminx,&nminy,&nminz,&nmaxx,&nmaxy,&nmaxz);
    if (fabs(Wat - WaMoRa) > 0.000001 || NPt != NPtSphere ||
        fabs(Size - GrSiSo) > 0.000001 || fabs(Incr - GrInSo) > 0.000001 ||
        fabs(epsre - DielRe) > 0.000001 || fabs(epswa - DielWa) > 0.000001 ||
        strcmp(ReDesoAlg,Alg) != 0 || Nx != NGridx || Ny != NGridy ||
        Nz != NGridz || fabs(Xfirst - XGrid[1]) > 0.000001 ||
        fabs(Yfirst - YGrid[1]) > 0.000001 ||
        fabs(Zfirst - ZGrid[1]) > 0.000001 ) {
      printf("\n");
      printf("Error: The Desolvation map you're trying to read has not\n");
      printf("been created with the same parameter you have specified\n");
      printf("in the input file:\n");
      printf("Actual water molecule radius, # points per sphere, grid size, grid increase:\n");
      printf("%f %d %f %f\n",WaMoRa,NPtSphere,GrSiSo,GrInSo);
      printf("Molecule radius, # points per sphere, grid size, grid increase you're trying to read:\n");
      printf("%f %d %f %f\n",Wat,NPt,Size,Incr);
      printf("Actual algoritm for the evaluation of the receptor desolvation:\n");
      printf("%s\n",ReDesoAlg);
      printf("Algoritm for the evaluation of the receptor desolvation in the file map you're trying to read:\n");
      printf("%s\n",Alg);
      printf("Actual solute and solvent dielectric constants:\n");
      printf("%f %f\n",DielRe,DielWa);
      printf("Solute and solvent dielectric constants in the file map you're trying to read:\n");
      printf("%f %f\n",epsre,epswa);
      printf("Actual # grid points in the box containing the receptor:\n");
      printf("%d %d %d\n",NGridx,NGridy,NGridz);
      printf("# grid points in the box containing the receptor in the file map you're trying to read:\n");
      printf("%d %d %d\n",Nx,Ny,Nz);
      printf("Actual coordinates of the grid point 1,1,1:\n");
      printf("%f %f %f\n",XGrid[1],YGrid[1],ZGrid[1]);
      printf("Coordinates of the grid point 1,1,1 in the file map you're trying to read:\n");
      printf("%f %f %f\n",Xfirst,Yfirst,Zfirst);
      printf("\n");
      printf("Program exits\n");
      exit(0);
    }
    else if ( nminx > nxminBS || nminy > nyminBS || nminz > nzminBS ||
              nmaxx < nxmaxBS || nmaxy < nymaxBS || nmaxz < nzmaxBS ) {
      printf("\n");
      printf("Error: The Desolvation map you're trying does not cover all \n");
      printf("the binding site.\n");
      printf("\n");
      printf("Program exits\n");
      exit(0);
    }
    for (ix=nminx;ix<=nmaxx;ix++)
      for (iy=nminy;iy<=nmaxy;iy++)
        for (iz=nminz;iz<=nmaxz;iz++)
          if (GridMat[ix][iy][iz] != 'o' ) {
            fgets_wrapper(StrLin,_STRLENGTH,ReaFile);
            sscanf(StrLin,"%lf",&(DeltaPrDeso)[ix][iy][iz]);
          }
    fclose(ReaFile);
  }
  else
    {
      if(DesoMapAcc[0] != 'n')
	printf("Undefined desolvation map access mode %s\n",DesoMapAcc);
    ix = 0;
    }

  if (ix == nxmaxBS+1 || ix == nmaxx+1)
    return 1;
  else
    return 0;
}

int Fast_Desol_Surf(struct point Min,double WaMoRa,double GrSiSo,
                   int NStartGridx,int NStartGridy,int NStartGridz,
                   int NGridx,int NGridy,int NGridz,int Nsurfpt,
                   struct point *surfpt,double ***DeltaPrDeso,
                   double corr_fast_deso,double *ReSurf_deso)
/*##########################################
Calculate the desolvation operated by a sphere rolled over the solute surface
This quantity will be used for a fast and approximate evaluation of the
elec desolvation of the rec and of the frag
###########################################*/

/*##########################################
struct point Min -------- Min coor of the grid box
double WaMoRa ----------- Radius of the water molecule
double GrSiSo ----------- Size of the 3D grid used for cont. electrostatics
int NStartGridx --------- First grid point along x (1)
int NStartGridy --------- First grid point along y (1)
int NStartGridz --------- First grid point along z (1)
int  NGridx ------------- Tot # of grid points along x
int  NGridy ------------- Tot # of grid points along y
int  NGridz ------------- Tot # of grid points along z
int Nsurfpt ------------- Tot # of points over SAS2
struct point *surfpt ---- Coor of points over SAS2
double ***DeltaPrDeso --- Elec desolvation due to the occupation of a grid point
double corr_fast_deso --- Correction factor for fast elec desolvation
double **ReSurf_deso ---- ReSurf_deso[n] = Desolvation of a sphere placed
                          on the nth point of the SAS2
###########################################*/
{
  int i,ix,iy,iz,ixmin,iymin,izmin,ixmax,iymax,izmax;
  double WaMoRa2,xtemp,x2temp,ytemp,xy2temp,ztemp,r2;

  WaMoRa2 = WaMoRa * WaMoRa;

  for (i=1;i<=Nsurfpt;i++) {
    ReSurf_deso[i] = 0.;

/* Calculate the extremes, in the grid frame, of the cube containing the
   probe sphere. They are integer! */

    ixmin = (surfpt[i].x - WaMoRa - Min.x) / GrSiSo + 1;
    iymin = (surfpt[i].y - WaMoRa - Min.y) / GrSiSo + 1;
    izmin = (surfpt[i].z - WaMoRa - Min.z) / GrSiSo + 1;
    ixmax = (surfpt[i].x + WaMoRa - Min.x) / GrSiSo + 1;
    iymax = (surfpt[i].y + WaMoRa - Min.y) / GrSiSo + 1;
    izmax = (surfpt[i].z + WaMoRa - Min.z) / GrSiSo + 1;
    ixmin = (ixmin > NStartGridx) ? ixmin : NStartGridx;
    iymin = (iymin > NStartGridy) ? iymin : NStartGridy;
    izmin = (izmin > NStartGridz) ? izmin : NStartGridz;
    ixmax = (ixmax < NGridx) ? ixmax : NGridx;
    iymax = (iymax < NGridy) ? iymax : NGridy;
    izmax = (izmax < NGridz) ? izmax : NGridz;

/* xtemp,ytemp,ztemp are cartesian coordinates of the grid pt ix,iy,iz,
   in a frame with the origin in the suface point surfpt[i] */

    xtemp = GrSiSo * (ixmin - 1.5) + Min.x - surfpt[i].x;
    for (ix=ixmin;ix<=ixmax;ix++) {
      xtemp += GrSiSo;
      x2temp = xtemp*xtemp;
      ytemp = GrSiSo * (iymin - 1.5) + Min.y - surfpt[i].y;
      for (iy=iymin;iy<=iymax;iy++) {
        ytemp += GrSiSo;
        xy2temp = ytemp*ytemp + x2temp;

/*  Check if we are already out of the sphere */

        if (WaMoRa2 > xy2temp) {
          ztemp = GrSiSo * (izmin - 1.5) + Min.z - surfpt[i].z;
          for (iz=izmin;iz<=izmax;iz++) {
            ztemp += GrSiSo;
            r2 = ztemp * ztemp + xy2temp;

/* Check if our grid point is inside the sphere. If yes change GridMat */

            if (WaMoRa2 > r2){
              ReSurf_deso[i] += DeltaPrDeso[ix][iy][iz];
            }
          }
        }
      }
    }
    ReSurf_deso[i] *= corr_fast_deso;
  }

  if (i == Nsurfpt+1 )
    return 1;
  else
    return 0;
}
