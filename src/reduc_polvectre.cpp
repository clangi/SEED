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
#include <math.h>
#include "nrutil.h"
#include "funct.h"
#include <stdlib.h> /*dey for gcc 4.0 compatibility*/

void Reduc_polvectre(int ReAtNu,double **ReCoor,double *ReVdWR,double *ReVdWE_sr,
                     int ReAcNu,int ReDoNu,int *ReDATy,int *ReDAAt,int *ReHydN,
                     double **ReVeCo,double RedRPV_rp,double RedRPV_nkvRatio, /*new double RedRPV_nkvRatio,*/
                     int *PolVect_rec_01,FILE *FPaOut,
		     int distrPointBSNumb,double **distrPointBS,double angle_rmin,
		     double angle_rmax,double mult_fact_rmin,double mult_fact_rmax,
		     double Sphere_apol)
/* This function reduces the number of polar vectors of the receptor by
   evaluating van der Waals energies of a probe located on these vectors :
   RedRPV_rp  van der Waals radius of the probe
   RedRPV_nkv  number of kept vectors (maximal number given by the user)
   RedRPV_nkvRatio ratio of kept vectors (maximal number given by the user)
                   after the angle criterion is applied
   PolVect_rec_01  polar vectors of the receptor (0 if it is rejected, 1 if it
                   is kept)
   ProbCoor  coordinates of the probe
   VdWEn_prob  van der Waals energy of the probe
   VdWEn_prob_arr_in  van der Waals energies of the probe
                      (ordering before sorting)
   VdWEn_prob_arr_in2  van der Waals energies of the probe
                       (ordering before sorting,then modified with cutoff)
   RPSqDi  squared distance between a receptor atom and the probe atom
   SumRad  sum of the van der Waals radii
   RedRPV_ep_sr  square root of the van der Waals energy parameter of the probe
   SortArray  array used in the sorting step
   NumKeVe  actual number of kept vectors after the algorithm
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
   keptVect_angle  polar vectors on receptor (0 if rejected, 1 if kept)
                   in case of angle cutoff criterion
   keptVect_angle_numb  number of rec polar vectors kept after angle
                        cutoff criterion applied
   vect_angle_stored  array for storing angle for each vector
   ReVeCo_elongated  in such a way the rec polar vectors are on a kind
                     of SAS */
{
  FILE *FilePa;
  int i,j,*SortArray,NumKeVe,closestPoint,*keptVect_angle,
      keptVect_angle_numb,helpCount,maxCount,RedRPV_nkv; /*new RedRPV_nkv == internal variable*/
  double ProbCoor[4],VdWEn_prob,RPSqDi,RPSqDi_p3,SumRad,SumRad_p6,
        *VdWEn_prob_arr,RedRPV_ep_sr,/*unused variable : ScaleFac,*/*VdWEn_prob_arr_in,
        *VdWEn_prob_arr_in2,minSqDist,maxSqDist,vecPointSqDist,
	distSqClosest,vecPointAngle,minSqDistCutoff,maxSqDistCutoff,
	angleCutoff,angle_rmin_rad,angle_rmax_rad,*vect_angle_stored,
	ReVeCo_elongated[4];



/* Initialisation */
  vect_angle_stored=dvector(1,ReAcNu+ReDoNu);
  keptVect_angle=ivector(1,ReAcNu+ReDoNu);
  if (distrPointBSNumb>0)
  {
    for (i=1;i<=(ReAcNu+ReDoNu);i++)
      keptVect_angle[i]=0;
    keptVect_angle_numb=0;
  }
  else
  {
    for (i=1;i<=(ReAcNu+ReDoNu);i++)
      keptVect_angle[i]=1;
    keptVect_angle_numb=ReAcNu+ReDoNu;
  }
  angle_rmin_rad=angle_rmin*3.1415927/180.0;
  angle_rmax_rad=angle_rmax*3.1415927/180.0;



  if (distrPointBSNumb>0)
  {
/* -------------------------------------- */
/*  reduction of rec polar vectors using  */
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

    for (i=1;i<=(ReAcNu+ReDoNu);i++)
    {

      for (j=1;j<=distrPointBSNumb;j++)
      {

        if (ReDATy[i]==0)
	  PoCoVe(ReVeCo[i][1],ReVeCo[i][2],ReVeCo[i][3],
	         ReVeCo[i][4],ReVeCo[i][5],ReVeCo[i][6],
                 ReVdWR[ReDAAt[i]]+Sphere_apol,
		 &(ReVeCo_elongated[1]),
		 &(ReVeCo_elongated[2]),
		 &(ReVeCo_elongated[3]));
	else
	  PoCoVe(ReVeCo[i][4],ReVeCo[i][5],ReVeCo[i][6],
	         ReVeCo[i][1],ReVeCo[i][2],ReVeCo[i][3],
                 ReVdWR[ReDAAt[i]]+Sphere_apol,
		 &(ReVeCo_elongated[1]),
		 &(ReVeCo_elongated[2]),
		 &(ReVeCo_elongated[3]));

	vecPointSqDist=DistSq(ReVeCo_elongated[1],
	                      ReVeCo_elongated[2],
			      ReVeCo_elongated[3],
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
    fprintf(FPaOut,"polar vectors \nand points in the binding site: ");
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
    for (i=1;i<=(ReAcNu+ReDoNu);i++)
    {

/* find the closest point */
      closestPoint=-1;
      distSqClosest=1000000.0;

      for (j=1;j<=distrPointBSNumb;j++)
      {

        if (ReDATy[i]==0)
	  PoCoVe(ReVeCo[i][1],ReVeCo[i][2],ReVeCo[i][3],
	         ReVeCo[i][4],ReVeCo[i][5],ReVeCo[i][6],
                 ReVdWR[ReDAAt[i]]+Sphere_apol,
		 &(ReVeCo_elongated[1]),
		 &(ReVeCo_elongated[2]),
		 &(ReVeCo_elongated[3]));
	else
	  PoCoVe(ReVeCo[i][4],ReVeCo[i][5],ReVeCo[i][6],
	         ReVeCo[i][1],ReVeCo[i][2],ReVeCo[i][3],
                 ReVdWR[ReDAAt[i]]+Sphere_apol,
		 &(ReVeCo_elongated[1]),
		 &(ReVeCo_elongated[2]),
		 &(ReVeCo_elongated[3]));

	vecPointSqDist=DistSq(ReVeCo_elongated[1],
	                      ReVeCo_elongated[2],
			      ReVeCo_elongated[3],
                              distrPointBS[j][1],distrPointBS[j][2],
		       	      distrPointBS[j][3]);

        if (vecPointSqDist<distSqClosest)
	{
	  distSqClosest=vecPointSqDist;
	  closestPoint=j;
	}

      }

      if(closestPoint<0) /* this was not checked previously */
	continue;
      /* compute the angle */
      if (ReDATy[i]==0) /* -> error closestPoint == -1 is not checked !!!! */
	{

	vecPointAngle=PlaAng(ReVeCo[i][1],ReVeCo[i][2],ReVeCo[i][3],
	                     ReVeCo[i][4],ReVeCo[i][5],ReVeCo[i][6],
                             distrPointBS[closestPoint][1],
			     distrPointBS[closestPoint][2],
		       	     distrPointBS[closestPoint][3]);
	}
      else
	vecPointAngle=PlaAng(ReVeCo[i][4],ReVeCo[i][5],ReVeCo[i][6],
	                     ReVeCo[i][1],ReVeCo[i][2],ReVeCo[i][3],
                             distrPointBS[closestPoint][1],
			     distrPointBS[closestPoint][2],
		   	     distrPointBS[closestPoint][3]);

      vect_angle_stored[i]=vecPointAngle/3.1415927*180.0;

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
        keptVect_angle[i]=1;
	keptVect_angle_numb++;
      }

    }


/* Print the kept vectors after angle cutoff criterion in a file */
    FilePa=fopen("./outputs/polar_rec_reduc_angle.mol2","w");
    fprintf(FilePa,"# TRIPOS MOL2 file generated by WITNOTP\n");
    fprintf(FilePa,"\n");
    fprintf(FilePa,"@<TRIPOS>MOLECULE\n");
    fprintf(FilePa,"polar_rec_reduc_angle\n");
    fprintf(FilePa,"%d 0 0 0 0\n",keptVect_angle_numb);
    fprintf(FilePa,"****\n");
    fprintf(FilePa,"USER_CHARGES\n");
    fprintf(FilePa,"INVALID_CHARGES\n");
    fprintf(FilePa,"\n");
    fprintf(FilePa,"@<TRIPOS>ATOM\n");

    j=0;
    for (i=1;i<=(ReAcNu+ReDoNu);i++)
    {
      if (keptVect_angle[i])
      {
        j=j+1;
        if (ReDATy[i]==0)
          fprintf(FilePa,"%-5d %-3d%11.4f%11.4f%11.4f ***** 1    ****  %11.4f\n",
                  j,j,ReVeCo[i][4],ReVeCo[i][5],ReVeCo[i][6],
		  vect_angle_stored[i]);
        else
          fprintf(FilePa,"%-5d %-3d%11.4f%11.4f%11.4f ***** 1    ****  %11.4f\n",
                  j,j,ReVeCo[i][1],ReVeCo[i][2],ReVeCo[i][3],
		  vect_angle_stored[i]);
      }
    }
    fclose(FilePa);


  }



/*keptVect_angle_numb*/

/* -------------------------------------- */
/*  reduction of rec polar vectors using  */
/*  van der Waals energy                  */
/* -------------------------------------- */

/* Allocate the memory and initialization */
  SortArray=ivector(1,ReAcNu+ReDoNu);
  VdWEn_prob_arr=dvector(1,ReAcNu+ReDoNu);
  VdWEn_prob_arr_in=dvector(1,ReAcNu+ReDoNu);
  VdWEn_prob_arr_in2=dvector(1,ReAcNu+ReDoNu);
  RedRPV_ep_sr=sqrtf(0.09);
  for (i=1;i<=(ReAcNu+ReDoNu);i++)
    PolVect_rec_01[i]=0;

/* Compute the van der Waals energy of the probe for each polar vector
   position */
  for (i=1;i<=(ReAcNu+ReDoNu);i++) {

/* Find the coordinates of the probe */
    if (ReDATy[i]==0)
      PoCoVe(ReVeCo[i][1],ReVeCo[i][2],ReVeCo[i][3],ReVeCo[i][4],ReVeCo[i][5],
             ReVeCo[i][6],ReVdWR[ReDAAt[i]]+RedRPV_rp,&ProbCoor[1],
             &ProbCoor[2],&ProbCoor[3]);
    else
      PoCoVe(ReVeCo[i][4],ReVeCo[i][5],ReVeCo[i][6],ReVeCo[i][1],ReVeCo[i][2],
             ReVeCo[i][3],ReVdWR[ReDAAt[i]]+RedRPV_rp,&ProbCoor[1],
             &ProbCoor[2],&ProbCoor[3]);

/* Compute the van der Waals energy of this probe. One doesn't take into
   account the receptor H involved in the current polar vector */
    VdWEn_prob=0.0;
    for (j=1;j<=ReAtNu;j++) {

      if (j!=ReHydN[i]) {

        RPSqDi=DistSq(ReCoor[j][1],ReCoor[j][2],ReCoor[j][3],
                      ProbCoor[1],ProbCoor[2],ProbCoor[3]);

        if (RPSqDi>(1.e-4)) {

          RPSqDi_p3=RPSqDi*RPSqDi*RPSqDi;
          SumRad=ReVdWR[j]+RedRPV_rp;
          SumRad_p6=SumRad*SumRad*SumRad*SumRad*SumRad*SumRad;

          VdWEn_prob=VdWEn_prob+ReVdWE_sr[j]*RedRPV_ep_sr*SumRad_p6*
                     ((SumRad_p6/(RPSqDi_p3*RPSqDi_p3))-(2/RPSqDi_p3));

        }
        else
          VdWEn_prob=VdWEn_prob+(1.e+12);

      }

    }

    VdWEn_prob_arr[i]=VdWEn_prob;
    VdWEn_prob_arr_in[i]=VdWEn_prob;
    VdWEn_prob_arr_in2[i]=VdWEn_prob;

  }

/* Sort the van der Waals energies of the probe */
  for (i=1;i<=(ReAcNu+ReDoNu);i++)
    SortArray[i]=i;

  Sort(ReAcNu+ReDoNu,SortArray,VdWEn_prob_arr);

/* calculate total number of requested polar vectors from user defined ratio : */
  if (distrPointBSNumb>0)
      RedRPV_nkv=(RedRPV_nkvRatio*(keptVect_angle_numb)); /* new  DEBUG CHECK REMOVE ???? */
  else
      RedRPV_nkv=(RedRPV_nkvRatio*(ReAcNu+ReDoNu)); /* new */

/* Keep only the polar vectors with the best van der Waals probe energies
   taking into account the angle cutoff criterion selection */
  if (((ReAcNu+ReDoNu)>=RedRPV_nkv)&&
      (keptVect_angle_numb>=RedRPV_nkv))
  {
    maxCount=RedRPV_nkv;
  }
  else
  {
    maxCount=keptVect_angle_numb;
  }

  helpCount=0;
  for (i=1;i<=(ReAcNu+ReDoNu);i++)
  {
    if ((keptVect_angle[SortArray[i]]==1)&&
        (helpCount<maxCount))
    {
      helpCount++;
      PolVect_rec_01[SortArray[i]]=1;
    }
  }

/* Count the actual number of kept vectors after the algorithm */
  NumKeVe=0;
  for (i=1;i<=(ReAcNu+ReDoNu);i++)
    NumKeVe=NumKeVe+PolVect_rec_01[i];

  if(NumKeVe==0)
      fprintf(stderr,"WARNING, number of kept polar vectors == 0 !\n");

  fprintf(FPaOut,"Total number of polar vectors generated on the ");
  fprintf(FPaOut,"receptor binding site : %d\n",ReAcNu+ReDoNu);
  fprintf(FPaOut,"Number of kept polar vectors for seeding of ");
  fprintf(FPaOut,"polar fragments : %d\n\n",NumKeVe);

/* Print the kept vectors in a file */
  FilePa=fopen("./outputs/polar_rec_reduc.mol2","w");
  fprintf(FilePa,"# TRIPOS MOL2 file generated by WITNOTP\n");
  fprintf(FilePa,"\n");
  fprintf(FilePa,"@<TRIPOS>MOLECULE\n");
  fprintf(FilePa,"polar_rec_reduc\n");
  fprintf(FilePa,"%d 0 0 0 0\n",NumKeVe);
  fprintf(FilePa,"****\n");
  fprintf(FilePa,"USER_CHARGES\n");
  fprintf(FilePa,"INVALID_CHARGES\n");
  fprintf(FilePa,"\n");
  fprintf(FilePa,"@<TRIPOS>ATOM\n");

  j=0;
  for (i=1;i<=(ReAcNu+ReDoNu);i++) {
    if (PolVect_rec_01[i]) {
      j=j+1;
/*
        fprintf(FilePa,"%-5d %-3d%11.4f%11.4f%11.4f ***** 1    ****  0.000\n",
                j*2-1,j*2-1,ReVeCo[i][1],ReVeCo[i][2],ReVeCo[i][3]);
        fprintf(FilePa,"%-5d %-3d%11.4f%11.4f%11.4f ***** 1    ****  0.000\n",
                j*2,j*2,ReVeCo[i][4],ReVeCo[i][5],ReVeCo[i][6]);
*/
      if (ReDATy[i]==0)
        fprintf(FilePa,"%-5d %-3d%11.4f%11.4f%11.4f ***** 1    ****  0.000\n",
                j,j,ReVeCo[i][4],ReVeCo[i][5],ReVeCo[i][6]);
      else
        fprintf(FilePa,"%-5d %-3d%11.4f%11.4f%11.4f ***** 1    ****  0.000\n",
                j,j,ReVeCo[i][1],ReVeCo[i][2],ReVeCo[i][3]);
    }
  }
  fclose(FilePa);



/* Print a file containing all the polar vectors (only the vector point which
   is not on the atom) with the corresponding van der Waals probe energy.
   The van der Waals energies which are larger than 100.0 are set to 100.0 */

/*
  FilePa=fopen("./outputs/polar_rec_map.mol2","w");
  fprintf(FilePa,"# TRIPOS MOL2 file generated by WITNOTP\n");
  fprintf(FilePa,"\n");
  fprintf(FilePa,"@<TRIPOS>MOLECULE\n");
  fprintf(FilePa,"polar_rec_map\n");
  fprintf(FilePa,"%d 0 %d 0 0\n",ReAcNu+ReDoNu,ReAcNu+ReDoNu);
  fprintf(FilePa,"****\n");
  fprintf(FilePa,"USER_CHARGES\n");
  fprintf(FilePa,"INVALID_CHARGES\n");
  fprintf(FilePa,"\n");
  fprintf(FilePa,"@<TRIPOS>ATOM\n");

  ScaleFac=50.0/(100.0-VdWEn_prob_arr[1]);
  for (i=1;i<=(ReAcNu+ReDoNu);i++) {
    if (VdWEn_prob_arr_in2[i]>100.0)
      VdWEn_prob_arr_in2[i]=100.0;
  }

  for (i=1;i<=(ReAcNu+ReDoNu);i++) {
    if (ReDATy[i]==0)
      fprintf(FilePa,"%-5d H  %11.4f%11.4f%11.4f H %4d   M%-4d  %7.4f\n",
              i,ReVeCo[i][4],ReVeCo[i][5],ReVeCo[i][6],i,i,
              ScaleFac*(VdWEn_prob_arr_in2[i]-VdWEn_prob_arr[1]));
    else
      fprintf(FilePa,"%-5d H  %11.4f%11.4f%11.4f H %4d   M%-4d  %7.4f\n",
              i,ReVeCo[i][1],ReVeCo[i][2],ReVeCo[i][3],i,i,
              ScaleFac*(VdWEn_prob_arr_in2[i]-VdWEn_prob_arr[1]));
  }

  for (i=1;i<=(ReAcNu+ReDoNu);i++)
    fprintf(FilePa,"# %d  %7.4f\n",i,VdWEn_prob_arr_in[i]);
  fclose(FilePa);
*/



/* Deallocate the memory */
  free_ivector(SortArray,1,ReAcNu+ReDoNu);
  free_dvector(VdWEn_prob_arr,1,ReAcNu+ReDoNu);
  free_dvector(VdWEn_prob_arr_in,1,ReAcNu+ReDoNu);
  free_dvector(VdWEn_prob_arr_in2,1,ReAcNu+ReDoNu);
  free_ivector(keptVect_angle,1,ReAcNu+ReDoNu);
  free_dvector(vect_angle_stored,1,ReAcNu+ReDoNu);

}
