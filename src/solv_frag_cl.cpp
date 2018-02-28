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
#include <math.h>
#include <time.h>
#include "nrutil.h"
#include "funct.h"
#include <iomanip> // added by clangini
#include <limits>  // added by clangini
#include "solv_frag_cl.h"

//using namespace std;

int FragSolvEn_cl(int FrAtNu,double **FrCoor,double *FrPaCh,
               double *FrVdWR, double *FrRad, double *FrRad2, double *FrRadOut,
               double *FrRadOut2, double **Frdist2,int Nsurfpt_fr,
               struct point *surfpt_fr_orig,double WaMoRa,double GrSiSo,
               double Ksolv,double pi4,double *PFrSolvEn,char *EmpCorrB,FILE*FPaOut)
/*##########################################
Calculate the solvation energy of the isolated frag
according to the GB formula
###########################################*/

/*##########################################
int FrAtNu -------------- Tot # frag atoms
double **FrCoor ---------- Frag coordinates in the original location
double *FrPaCh ----------- Frag partial charges
double *FrVdWR ----------- Frag vdW radii
double *FrRadOut -------- Frag charge radii + WaMoRa
double *FrRadOut2 ------- (Frag charge radii + WaMoRa)^2
double **Frdist2 --------- Squared interatomic frag distances
int Nsurfpt_fr ---------- Tot # of points over the frag SAS1
struct point *surfpt_fr_orig - Coor of points over frag SAS1 (obtained from
                               FrCoor coordinates)
double WaMoRa ----------- Radius of the water molecule
double GrSiSo ----------- Size of the 3D grid used for cont. electrostatics
double Ksolv ------------ Constant
double pi4 -------------- 4 * greekpi
double *PFrSolvEn ------- Frag solvation energy
###########################################*/
{
  int iat,ix,iy,iz,nn,NGridx,NGridy,NGridz,NGrid;
  double *FrSelfVol,*FrEffRad,FrSelfEn,FrIntEn,UnitVol,*XGrid,*YGrid,*ZGrid,
         *FrSelfVol_corrB;
  struct point Min,Max;
  char ***FrGridMat;
/*##########################################
int iat,ix,iy,iz,nn ----- Multipurpose variables
int NGridx -------------- Tot # of frag grid points along x
int NGridy -------------- Tot # of frag grid points along y
int NGridz -------------- Tot # of frag grid points along z
double *FrSelfVol ------- FrSelfVol[iat] = integral of 1/r^4 over the frag
                          volume for atom iat
double *FrEffRad -------- FrEffRad[n] = effective radius of frag atom n
double FrSelfEn --------- Total self energy of the isolated frag
double FrIntEn ---------- Total intramolecular interaction energy
                          of the isolated frag
double UnitVol ---------- Volume of the grid element for cont. elec.
double *XGrid ----------- X coor of the grid points
double *YGrid ----------- Y coor of the grid points
double *ZGrid ----------- Z coor of the grid points
struct point Min -------- Min coor of the frag grid box
struct point Max -------- Max coor of the frag grid box
char ***FrGridMat ------- Matrix telling if a grid point is occupied by the
                          frag (o), empty (e), or if it belongs to the interface
                          between SAS and MS (s)
###########################################*/

  /* Get the # of grid pts needed to cover the frag molecule along the
     three axis (no grid increase) */
  nn = get_Grid_Dim(FrAtNu,FrCoor,FrVdWR,WaMoRa,GrSiSo,0.,&NGridx,&NGridy,
                    &NGridz,&NGrid,&Min,&Max,&UnitVol); // instead of 0.

  XGrid=dvector(1,NGridx);
  YGrid=dvector(1,NGridy);
  ZGrid=dvector(1,NGridz);
  /* Get the cartesian coor for the grid points */
  nn = Coor_Grid_Pts(GrSiSo,NGridx,NGridy,NGridz,Min,XGrid,YGrid,ZGrid);

  FrSelfVol=dvector(1,FrAtNu);
  FrSelfVol_corrB=dvector(1,FrAtNu);
  FrEffRad=dvector(1,FrAtNu);
  FrGridMat = c3tensor(1,NGridx, 1, NGridy, 1, NGridz);
  for (ix = 1; ix <= NGridx; ix++)
    for (iy = 1; iy <= NGridy; iy++)
      for (iz = 1; iz <= NGridz; iz++)
        FrGridMat[ix][iy][iz] = 'e';
/* Make the map (FrGridMat) of the 3D grid points occupied by the volume
   enclosed by the frag SAS  */
  nn = SAS_Volume_cl(FrAtNu,FrCoor,FrRadOut,FrRadOut2,Min,GrSiSo,
                  1,1,1,NGridx,NGridy,NGridz,FrGridMat);
  for (iat=1;iat<=FrAtNu;iat++)
  {
    FrSelfVol[iat] = 0.;
    FrSelfVol_corrB[iat] = 0.;
  }
/* Place a sphere of radius WaMoRa on every SAS surface
   grid point (*surfpt_fr_orig) and mark as empty all the volume grid
   points falling inside the sphere.
   This corresponds to building the MS */
  nn = Excl_Grid_cl(FrAtNu,FrCoor,Min,FrRadOut,FrRadOut2,WaMoRa,GrSiSo,
                 1,1,1,NGridx,NGridy,NGridz,UnitVol,
                 FrGridMat,Nsurfpt_fr,surfpt_fr_orig,FrSelfVol,FrSelfVol_corrB,
		             EmpCorrB);
/* Calculate the frag self energy */
  nn = Get_Self_En_Fr_cl(FrAtNu,FrCoor,FrPaCh,FrRad,FrRad2,
                      XGrid,YGrid,ZGrid,UnitVol,Ksolv,pi4,FrGridMat,
                      1,1,1,NGridx,NGridy,NGridz,FrSelfVol,FrEffRad,&FrSelfEn,
		                  FrSelfVol_corrB,EmpCorrB,FPaOut);
  // for (int jj=1;jj<=FrAtNu;jj++){
  //   std::cout << "FrSelfVol_corrB[" << jj << "]= " << FrSelfVol_corrB[jj] << '\n';
  //   std::cout << "FrSelfVol[" << jj << "]= " << FrSelfVol[jj] << '\n';
  // }

/* Calculate the frag interaction energy */
  nn = GB_int_fr(FrAtNu,Frdist2,FrPaCh,FrEffRad,Ksolv,&FrIntEn);

  *PFrSolvEn = FrSelfEn + FrIntEn;
/*  printf("%lf\n",*PFrSolvEn); */

/* Free memory */
  free_dvector(XGrid,1,NGridx);
  free_dvector(YGrid,1,NGridy);
  free_dvector(ZGrid,1,NGridz);
  free_dvector(FrEffRad,1,FrAtNu);
  free_dvector(FrSelfVol,1,FrAtNu);
  free_c3tensor(FrGridMat,1,NGridx,1,NGridy,1,NGridz);
  free_dvector(FrSelfVol_corrB,1,FrAtNu);

  if ( iat == FrAtNu+1 )
    return 1;
  else
    return 0;
}

int SAS_Volume_cl(int ReAtNu,double **ReCoor,double *ReRadOut,
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

    if (ixmax > NGridx || iymax > NGridy || izmax > NGridz){
      std::cout << "ERROR in determining solute occupancy grid" << std::endl;
    }
    // ixmin = (ixmin > NStartGridx) ? ixmin : NStartGridx;
    // iymin = (iymin > NStartGridy) ? iymin : NStartGridy;
    // izmin = (izmin > NStartGridz) ? izmin : NStartGridz;
    // ixmax = (ixmax < NGridx) ? ixmax : NGridx;
    // iymax = (iymax < NGridy) ? iymax : NGridy;
    // izmax = (izmax < NGridz) ? izmax : NGridz;

    for (ix=ixmin;ix <= ixmax;ix++) {
      xtemp = GrSiSo * (ix - 0.5) + Min.x - ReCoor[iat][1];
      x2temp = xtemp*xtemp;
      for (iy=iymin;iy <= iymax;iy++) {
        ytemp = GrSiSo * (iy - 0.5) + Min.y - ReCoor[iat][2];
        xy2temp = ytemp*ytemp + x2temp;
        if (ReRadOut2[iat] > xy2temp) {
          for (iz=izmin; iz <= izmax;iz++) {
            ztemp = GrSiSo * (iz - 0.5) + Min.z - ReCoor[iat][3];
            if (GridMat[ix][iy][iz] == 'e') {
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


int Excl_Grid_cl(int ReAtNu,double **ReCoor,struct point Min,double *ReRadOut,
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

    // if (ixmin < 1 || iymin < 1 || izmin < 1){
    //   std::cout << "ERROR in determining interface between MS and SAS (min)" << std::endl;
    //   std::cout << "ix, iy, iz = " << ixmin << ", " << iymin << ", " << izmin << std::endl;
    //   std::cout << "(surfpt[i].x - WaMoRa - Min.x) = " << (surfpt[i].x - WaMoRa - Min.x) << std::endl;
    // }
    // if (ixmax > NGridx || iymax > NGridy || izmax > NGridz){
    //   std::cout << "ERROR in determining interface between MS and SAS (max)" << std::endl;
    // }

    ixmin = (ixmin > NStartGridx) ? ixmin : NStartGridx;
    iymin = (iymin > NStartGridy) ? iymin : NStartGridy;
    izmin = (izmin > NStartGridz) ? izmin : NStartGridz;
    ixmax = (ixmax < NGridx) ? ixmax : NGridx;
    iymax = (iymax < NGridy) ? iymax : NGridy;
    izmax = (izmax < NGridz) ? izmax : NGridz;

/* xtemp,ytemp,ztemp are cartesian coordinates of the grid pt ix,iy,iz,
   in a frame with the origin in the suface point surfpt[i] */

    for (ix = ixmin; ix <= ixmax; ix++) {
      xtemp = GrSiSo * (ix - 0.5) + Min.x - surfpt[i].x;
      x2temp = xtemp * xtemp;
      for (iy = iymin; iy <= iymax; iy++) {
        ytemp = GrSiSo * (iy - 0.5) + Min.y - surfpt[i].y;
        xy2temp = ytemp*ytemp + x2temp;

/*  Check if we are already out of the sphere */

        if (WaMoRa2 > xy2temp) {
          for (iz = izmin; iz <= izmax; iz++) {
            ztemp = GrSiSo * (iz - 0.5) + Min.z - surfpt[i].z;
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

  // //int jumppoint; //debug
  // for (iat=1;iat<=ReAtNu;iat++) {
  //   //jumppoint = 0;
  //   ixmin = ( (ReCoor[iat][1] - ReRadOut[iat] - Min.x) / GrSiSo + 1 );
  //   iymin = ( (ReCoor[iat][2] - ReRadOut[iat] - Min.y) / GrSiSo + 1 );
  //   izmin = ( (ReCoor[iat][3] - ReRadOut[iat] - Min.z) / GrSiSo + 1 );
  //   ixmax = ( (ReCoor[iat][1] + ReRadOut[iat] - Min.x) / GrSiSo + 1 );
  //   iymax = ( (ReCoor[iat][2] + ReRadOut[iat] - Min.y) / GrSiSo + 1 );
  //   izmax = ( (ReCoor[iat][3] + ReRadOut[iat] - Min.z) / GrSiSo + 1 );
  //   ixmin = (ixmin > NStartGridx) ? ixmin : NStartGridx;
  //   iymin = (iymin > NStartGridy) ? iymin : NStartGridy;
  //   izmin = (izmin > NStartGridz) ? izmin : NStartGridz;
  //   ixmax = (ixmax < NGridx) ? ixmax : NGridx;
  //   iymax = (iymax < NGridy) ? iymax : NGridy;
  //   izmax = (izmax < NGridz) ? izmax : NGridz;
  //
  //   xtemp = GrSiSo * (ixmin - 1.5) + Min.x - ReCoor[iat][1];
  //   for (ix=ixmin;ix<=ixmax;ix++) {
  //     xtemp += GrSiSo;
  //     x2temp = xtemp*xtemp;
  //     ytemp = GrSiSo * (iymin - 1.5) + Min.y - ReCoor[iat][2];
  //     for (iy=iymin;iy<=iymax;iy++) {
  //       ytemp += GrSiSo;
  //       xy2temp = ytemp*ytemp + x2temp;
  //       if (ReRadOut2[iat] > xy2temp) {
  //         ztemp = GrSiSo * (izmin - 1.5) + Min.z - ReCoor[iat][3];
  //         for (iz=izmin;iz<=izmax;iz++) {
  //           ztemp += GrSiSo;
  //           r2 = ztemp * ztemp + xy2temp;
  //           if (ReRadOut2[iat] > r2) {
  //             if (GridMat[ix][iy][iz] == 's') { // 's'
  //               r4 = r2 * r2;
  //               SelfVol[iat] -= UnitVol/r4;
	// 	            if (EmpCorrB[0]=='y')
	// 	              SelfVol_corrB[iat] -= UnitVol/(r4*sqrt(r2));
  //               /*std::cout << "here subtracting volume for atom " << iat << "\n";*/
  //               //jumppoint++;
  //             }
  //           }
  //         }
  //       }
  //     }
  //   }
  //   //std::cout << "Subtracted surf point for atom " << iat << " = " << jumppoint <<"\n";
  // }
  if (iat == ReAtNu+1 )
    return 1;
  else
    return 0;
}

int Get_Self_En_Fr_cl(int FrAtNu,double **RoSFCo,double *FrPaCh,double *FrRad,
                   double *FrRad2,double *XGrid,double *YGrid,double *ZGrid,
                   double UnitVol,double Ksolv,double pi4,char ***FrGridMat,
                   int nxminFr,int nyminFr,int nzminFr,
                   int nxmaxFr,int nymaxFr,int nzmaxFr,
                   double *FrSelfVol,double *FrEffRad,double *PFrSelfEn,
		               double *FrSelfVol_corrB,char *EmpCorrB,FILE * FPaOut)
/*##########################################
Calculate the frag self energy
###########################################*/

/*##########################################
int FrAtNu -------------- Tot # frag atoms
double **RoSFCo ---------- Frag coordinates in the seeded conformation
double *FrPaCh ----------- Frag partial charges
double *FrRadOut -------- Frag charge radii + WaMoRa
double *FrRadOut2 ------- (Frag charge radii + WaMoRa)^2
double *XGrid ----------- X coor of the grid points
double *YGrid ----------- Y coor of the grid points
double *ZGrid ----------- Z coor of the grid points
double UnitVol ---------- Volume of the grid element for cont. elec.
double Ksolv ------------ Constant
double pi4 -------------- 4 * greekpi
char ***FrGridMat ------- Small submatrix of GridMat around the frag telling the
                          volume occupied by the frag in the bound conformation
int nxminFr ------------- grid point (along x) where the frag starts
int nyminFr ------------- grid point (along y) where the frag starts
int nzminFr ------------- grid point (along z) where the frag starts
int nxmaxFr ------------- grid point (along x) where the frag ends
int nymaxFr ------------- grid point (along y) where the frag ends
int nzmaxFr ------------- grid point (along z) where the frag ends
double *FrSelfVol ------- FrSelfVol[iat] = integral of 1/r^4 over the solute
                          volume for frag atom iat
double *FrEffRad -------- FrEffRad[n] = effective radius of frag atom n
double *PFrSelfEn ------- Tot frag self-energy
###########################################*/
{
  int ix,iy,iz,iat;
  double r2,r4;

  ix=nxminFr;

  *PFrSelfEn = 0.;
  for (iat=1;iat<=FrAtNu;iat++) {
    if (FrPaCh[iat] != 0. ) {
      for (ix = nxminFr; ix <= nxmaxFr; ix++) {
        for (iy = nyminFr; iy <= nymaxFr; iy++) {
          for (iz = nzminFr; iz <= nzmaxFr; iz++) {
            if ( FrGridMat[ix][iy][iz] == 'o' ) { // only integrate on 'o' (not 's' or 'e')
              r2 = (XGrid[ix] - RoSFCo[iat][1]) *
                   (XGrid[ix] - RoSFCo[iat][1]) +
                   (YGrid[iy] - RoSFCo[iat][2]) *
                   (YGrid[iy] - RoSFCo[iat][2]) +
                   (ZGrid[iz] - RoSFCo[iat][3]) *
                   (ZGrid[iz] - RoSFCo[iat][3]);
              if ( r2 > FrRad2[iat]) {
                r4 = r2 * r2;
                FrSelfVol[iat] += UnitVol / r4;
		            if (EmpCorrB[0]=='y')
		              FrSelfVol_corrB[iat] += UnitVol / (r4*sqrt(r2));
              }
            }
          }
        }
      }

      //std::cout << "FrSelfVol[" << iat << "]/pi4 = " << FrSelfVol[iat]/pi4 << "\n";
      //std::cout << "Calculated Born radius[" << iat << "] = " << 1. / ( 1./(FrRadOut[iat]) - FrSelfVol[iat]/pi4 ) << "\n";
      //std::cout << "My hypothesis Born radius[" << iat << "] = " << 1. / ( 1./(FrRadOut[iat]-3.5) + FrSelfVol[iat]/pi4 ) << "\n";
      //std::cout << "FrRadOut[" << iat << "] = " << FrRadOut[iat] << "\n" ;

      if (EmpCorrB[0]!='y')
        FrEffRad[iat] = 1. / ( 1./FrRad[iat] - FrSelfVol[iat]/pi4 );
      else
      {

        FrEffRad[iat] = 1./( (-1.*(1./FrRad[iat] - FrSelfVol[iat]/pi4))
      		       + 3.0*sqrt( (1./(2.*FrRad2[iat])) - (FrSelfVol_corrB[iat]/pi4) ) )
                 + 0.215;
        std::cout << "=====" << "\n";
        std::cout << "1/2R^2 = " << 1./(2.*FrRad2[iat]) << " FrSelfVol_corrB[iat]/pi4 " << (FrSelfVol_corrB[iat]/pi4) << std::endl;
        std::cout << "FrEffRad_corrB[" << iat << "] = " << FrEffRad[iat] << std::endl;
        std::cout << "FrEffRad[" << iat << "] = " << 1. / ( 1./FrRad[iat] - FrSelfVol[iat]/pi4 ) << std::endl;
        std::cout << "Lower bound (Born radius) = " << FrRad[iat] << std::endl;

        /*
          Dey exception handling :
          in rare cases the expression :
          (1./(2.*ReRadOut[iat]*ReRadOut[iat])) - ((*SelfVol_corrB)[iat]/pi4)
          can become < 0 -> the sqrt() function cannot be evaluated, which leads
          to "nan" values or the effective born radius is smaller than 0

        */
        if(FrEffRad[iat]<=0 || isnan(FrEffRad[iat]))
        {
      #ifndef NOWARN
          fprintf(FPaOut,"WARNING could not calculate empirically-corrected effective born radius of fragment atom %d, using standard approach\n",iat);
      #endif
          //std::cout << "Fragment effective radius is 0 or nan" << std::endl;
          FrEffRad[iat] = 1. / ( 1./FrRad[iat] - FrSelfVol[iat]/pi4 );
        }

      }

      *PFrSelfEn += Ksolv * FrPaCh[iat] * FrPaCh[iat] / (2. * FrEffRad[iat]);

    }
    //std::cout << "Jumped surf point for atom " << iat << " = " << jumppoint <<"\n";
  }

  if (ix == FrAtNu+1 )
    return 1;
  else
    return 0;
}
