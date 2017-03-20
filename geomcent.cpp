#include <stdio.h>
#include "nrutil.h"
#include "funct.h"

#ifndef _STRLENGTH
#define _STRLENGTH 500
#endif

void GeomCenter_FFLD(int CurFra,char **FrFiNa_out,int FrAtNu,int *FrAtEl_nu,
                     int NuSdClKe,int *ReprSdClAr,float *To_s_ro,
                     float ***FrCoPo,char **ResN_fr,char *FragNa,
                     int gc_reprke,float gc_cutclus,float gc_endifclus,
                     float gc_weighneg,float gc_weighpos,int gc_maxwrite)
/* This function clusters the geometrical centers of the postprocessed 
   representative positions (of the conserved second clusters). They 
   are supposed to be ranked according to increasing energy. The 
   resulting geometrical centers are written out in a file that can be 
   used by the software FFLD :
   GeomCentAr       array of geometrical centers (of the original positions)
   GeomCentList     list of geometrical centers (for FFLD)
   GeomCentEner     energy for each geom. center (for FFLD)
   NumNonHy         number of non-hydrogen atoms in the fragment
   distSqu          squared distance
   ClusCutSqDi      cutoff (squared distance) for the clustering
   ClustListLabel   list of cluster labels
   clusCount        cluster counter
   NumClusReprKe    number of cluster representatives to keep
   totNumbGeomCent  total number of geometrical centers (for FFLD)
   AllowEnerDif     allowed energy difference in a cluster
   maxNuClusMemb    maximal number of cluster members
   minEnClus        minimal energy in cluster
   maxEnClus        maximal energy in cluster
   numbClusMemb     number of cluster members
   listClustMemb    list of cluster members
   countGeomCent    count geometrical centers
   countNegEner     count cluster members with a negative energy
   totClusEner      total cluster energy
   multiplFac       multiplicative factor
   weightPosEn      weight for positive energies
   weightNegEn      weight for negative energies
   weightFactPos    weight factor for positive energies
   weightFactNeg    weight factor for negative energies
   sortFloatAr      array of float numbers to be sorted
   sortIntAr        array of integer numbers for sorting 
   gc_reprke        number of cluster representatives to keep
   gc_cutclus       cutoff for the clustering
   gc_endifclus     allowed energy difference in a cluster
   gc_weighneg      weight for negative energies
   gc_weighpos      weight for positive energies
   gc_maxwrite      maximal number of geometrical centers to write */
{

  FILE *FilePath;
  char WriPat[_STRLENGTH],AtomName[50];
  float **GeomCentAr,sum_x,sum_y,sum_z,distSqu,ClusCutSqDi,**GeomCentList,
        AllowEnerDif,minEn,minEnClus,maxEnClus,totClusEner,multiplFac,
        weightPosEn,weightNegEn,weightFactPos,weightFactNeg,*GeomCentEner,
        *sortFloatAr;
  int i1,i2,j,NumNonHy,*ClustListLabel,clusCount,NumClusReprKe,exLoop,
      totNumbGeomCent,counter,boolVar,maxNuClusMemb,numbClusMemb,
      countGeomCent,*listClustMemb,countNegEner,*sortIntAr,list1,list2;

  /* uninitialized value fix */
  minEn     =  10e10;
  minEnClus =  10e10;
  maxEnClus = -10e10;

/* initialization of variables */
  NumClusReprKe = gc_reprke;
  ClusCutSqDi = gc_cutclus*gc_cutclus;
  AllowEnerDif = gc_endifclus;
  weightNegEn = gc_weighneg;
  weightPosEn = gc_weighpos;

/* compute the geometrical center of every position. Hydrogen 
   atoms are not taken into account */
  GeomCentAr=matrix(1,NuSdClKe,1,3);

  NumNonHy = 0;
  for (i1=1;i1<=FrAtNu;i1++)
  {
    if (FrAtEl_nu[i1] != 1) {NumNonHy++;}
  }

  for (i1=1;i1<=NuSdClKe;i1++)
  {
    sum_x = 0.0;
    sum_y = 0.0;
    sum_z = 0.0;
    for (i2=1;i2<=FrAtNu;i2++)
    {
      if (FrAtEl_nu[i2] != 1)
      {
        sum_x += FrCoPo[ReprSdClAr[i1]][i2][1];
        sum_y += FrCoPo[ReprSdClAr[i1]][i2][2];
        sum_z += FrCoPo[ReprSdClAr[i1]][i2][3];
      }
    }
    GeomCentAr[i1][1] = sum_x/((float) NumNonHy);
    GeomCentAr[i1][2] = sum_y/((float) NumNonHy);
    GeomCentAr[i1][3] = sum_z/((float) NumNonHy);
  }

/* cluster the geometrical centers */
  ClustListLabel=ivector(1,NuSdClKe);
  for (i1=1;i1<=NuSdClKe;i1++) {ClustListLabel[i1] = 0;}

  clusCount = 0;

  for (i1=1;i1<=NuSdClKe;i1++)
  {

    if (ClustListLabel[i1] == 0)
    {

      clusCount++;
      ClustListLabel[i1] = clusCount;

      for (i2=1;i2<=NuSdClKe;i2++)
      {

        if (ClustListLabel[i2] == 0)
        {
          distSqu = DistSq(GeomCentAr[i1][1],GeomCentAr[i1][2],
                           GeomCentAr[i1][3],
                           GeomCentAr[i2][1],GeomCentAr[i2][2],
                           GeomCentAr[i2][3]);
          if (distSqu<ClusCutSqDi) {ClustListLabel[i2] = clusCount;}
        }

      }

    }

  }

/* reassignment of NumClusReprKe */
  if (NumClusReprKe>clusCount)
  {
    NumClusReprKe = clusCount;
  }

/* label (with -1) the first NumClusReprKe cluster representatives */
  for (i1=1;i1<=NumClusReprKe;i1++)
  {
    exLoop = 0;
    for (i2=1;((i2<=NuSdClKe)&&(!exLoop));i2++)    
    {
      if (ClustListLabel[i2] == i1)
      {
        ClustListLabel[i2] = -1;
        exLoop = 1;
      }
    }  
  }

/* cluster again the geometrical centers without the first 
   NumClusReprKe cluster representatives */
  for (i1=1;i1<=NuSdClKe;i1++)
  {
    if (ClustListLabel[i1]!=(-1)) {ClustListLabel[i1] = 0;}
  }

  clusCount = 0;

  for (i1=1;i1<=NuSdClKe;i1++)
  {

    if (ClustListLabel[i1] == 0)
    {

      clusCount++;
      ClustListLabel[i1] = clusCount;

      for (i2=1;i2<=NuSdClKe;i2++)
      {

        if (ClustListLabel[i2] == 0)
        {
          distSqu = DistSq(GeomCentAr[i1][1],GeomCentAr[i1][2],
                           GeomCentAr[i1][3],
                           GeomCentAr[i2][1],GeomCentAr[i2][2],
                           GeomCentAr[i2][3]);
          if (distSqu<ClusCutSqDi) {ClustListLabel[i2] = clusCount;}
        }

      }

    }

  }

/* label (with -2) cluster members with an energy larger than 
   AllowEnerDif with respect to the best cluster energy. 
   Compute the maximal number of cluster members */
  maxNuClusMemb = 0;

  for (i1=1;i1<=clusCount;i1++)
  {

    boolVar = 1;
    counter = 0;

    for (i2=1;i2<=NuSdClKe;i2++)
    {

      if (ClustListLabel[i2] == i1)
      {

        if (boolVar)
        {
          minEn = To_s_ro[ReprSdClAr[i2]];
          boolVar = 0;
          counter++;
        }
        else
        {
          if (To_s_ro[ReprSdClAr[i2]] > (minEn+AllowEnerDif))
          {
            ClustListLabel[i2] = -2;
          }
          else
          {
            counter++;
          }
        }

      }

    }

    if (counter>maxNuClusMemb) {maxNuClusMemb = counter;}

  }

/* define GeomCentList and GeomCentEner */
/* store the geometrical centers and energies of the first 
   NumClusReprKe cluster representatives */
  totNumbGeomCent = NumClusReprKe+clusCount;
  GeomCentList=matrix(1,totNumbGeomCent,1,3);
  GeomCentEner=vector(1,totNumbGeomCent);

  countGeomCent = 0;
  for (i1=1;i1<=NuSdClKe;i1++)
  {
    if (ClustListLabel[i1]==(-1))
    {
      countGeomCent++;
      GeomCentList[countGeomCent][1] = GeomCentAr[i1][1];
      GeomCentList[countGeomCent][2] = GeomCentAr[i1][2];
      GeomCentList[countGeomCent][3] = GeomCentAr[i1][3];
      GeomCentEner[countGeomCent] = To_s_ro[ReprSdClAr[i1]];
    }
  }

/* store the following geometrical centers */
  listClustMemb=ivector(1,maxNuClusMemb);

  for (i1=1;i1<=clusCount;i1++)
  {

    numbClusMemb = 0;
    countNegEner = 0;
    for (j=1;j<=maxNuClusMemb;j++) {listClustMemb[j] = 0;}

    for (i2=1;i2<=NuSdClKe;i2++)
    {
      if (ClustListLabel[i2] == i1)
      {

        numbClusMemb++;
        listClustMemb[numbClusMemb] = i2;

        if (To_s_ro[ReprSdClAr[i2]] < 0) {countNegEner++;}

      }
    }

    countGeomCent++;
    GeomCentList[countGeomCent][1] = 0.0;
    GeomCentList[countGeomCent][2] = 0.0;
    GeomCentList[countGeomCent][3] = 0.0;
    GeomCentEner[countGeomCent] = 0.0;

    if ((countNegEner == 0)||(countNegEner == numbClusMemb))
    {
      weightFactNeg = 1.0;
      weightFactPos = 1.0;
    }
    else
    {
      weightFactNeg = weightNegEn;
      weightFactPos = weightPosEn;
    }

/* subset of positions with a negative energy */
    if (countNegEner)
    {
      minEnClus = To_s_ro[ReprSdClAr[listClustMemb[1]]];
      maxEnClus = To_s_ro[ReprSdClAr[listClustMemb[countNegEner]]];
    }

    totClusEner = 0.0;
    for (j=1;j<=countNegEner;j++)
    {totClusEner += To_s_ro[ReprSdClAr[listClustMemb[j]]];}

    for (j=1;j<=countNegEner;j++)
    {
      multiplFac = To_s_ro[ReprSdClAr[listClustMemb[j]]]/totClusEner;
      multiplFac *= weightFactNeg;

      GeomCentList[countGeomCent][1] += 
                       (multiplFac * GeomCentAr[listClustMemb[j]][1]);
      GeomCentList[countGeomCent][2] += 
                       (multiplFac * GeomCentAr[listClustMemb[j]][2]);
      GeomCentList[countGeomCent][3] += 
                       (multiplFac * GeomCentAr[listClustMemb[j]][3]);
      GeomCentEner[countGeomCent] += 
                       (multiplFac * To_s_ro[ReprSdClAr[listClustMemb[j]]]);
    }

/* subset of positions with a positive energy */
    if ((countNegEner+1)<=numbClusMemb)
    {
      minEnClus = To_s_ro[ReprSdClAr[listClustMemb[countNegEner+1]]];
      maxEnClus = To_s_ro[ReprSdClAr[listClustMemb[numbClusMemb]]];
    }

    totClusEner = 0.0;
    for (j=countNegEner+1;j<=numbClusMemb;j++)
    {
      totClusEner += (To_s_ro[ReprSdClAr[listClustMemb[j]]]
                      - (minEnClus + maxEnClus));
    }

    for (j=countNegEner+1;j<=numbClusMemb;j++)
    {
      multiplFac = (To_s_ro[ReprSdClAr[listClustMemb[j]]]
                    - (minEnClus + maxEnClus)) / totClusEner;
      multiplFac *= weightFactPos;

      GeomCentList[countGeomCent][1] += 
                       (multiplFac * GeomCentAr[listClustMemb[j]][1]);
      GeomCentList[countGeomCent][2] += 
                       (multiplFac * GeomCentAr[listClustMemb[j]][2]);
      GeomCentList[countGeomCent][3] += 
                       (multiplFac * GeomCentAr[listClustMemb[j]][3]);
      GeomCentEner[countGeomCent] += 
                       (multiplFac * To_s_ro[ReprSdClAr[listClustMemb[j]]]);
    }

  }

/* sort "geometrical center energies" without the first 
   NumClusReprKe cluster representatives */
  sortFloatAr=vector(1,clusCount);
  sortIntAr=ivector(1,clusCount);
  for (i1=1;i1<=clusCount;i1++)
  {
    sortFloatAr[i1] = GeomCentEner[i1+NumClusReprKe];
    sortIntAr[i1] = i1+NumClusReprKe;
  }
  Sort(clusCount,sortIntAr,sortFloatAr);

/* write out geometrical centers in a mol2 file */
  sprintf(WriPat,"%s%s%s","./outputs/",&FrFiNa_out[CurFra][1],
          "_geomcent.mol2\0");
  FilePath=fopen(WriPat,"w");

/* adjust the number of geometrical centers to write according 
   to gc_maxwrite; 
   list1 refers to the list of cluster representatives to keep, 
         the maximal value would be NumClusReprKe; 
   list2 refers to the list of following geometrical centers, 
         the maximal value would be (countGeomCent-NumClusReprKe) */
/* Modification 06.08.2004 NM */
  list1=0;
  list2=0;
  if (gc_maxwrite<=NumClusReprKe)
  {
    list1=gc_maxwrite;
    list2=0;
  }
  else if ( (gc_maxwrite>NumClusReprKe) && 
            (gc_maxwrite<=countGeomCent) )
  {
    list1=NumClusReprKe;
    list2=gc_maxwrite-NumClusReprKe;
  }
  else if (gc_maxwrite>countGeomCent)
  {
    list1=NumClusReprKe;
    list2=countGeomCent-NumClusReprKe;
  }

  fprintf(FilePath,"# TRIPOS MOL2 file generated by SEED\n\n");

  fprintf(FilePath,"@<TRIPOS>MOLECULE\n");
  fprintf(FilePath,"%s\n",FragNa);
  fprintf(FilePath," %d 0 0 0 0\n",(list1+list2));
  fprintf(FilePath,"****\n");
  fprintf(FilePath,"USER_CHARGES\n");
  fprintf(FilePath,"INVALID_CHARGES\n\n");

  fprintf(FilePath,"@<TRIPOS>ATOM\n");

  counter=0;

/* first NumClusReprKe geometrical centers */
  for (j=1;j<=list1;j++)
  {
    counter++;
    sprintf(AtomName,"ZN%d",j); /* prev "ZN%d\0" */
    fprintf(FilePath,"%7d %-8s%10.4f%10.4f%10.4f     MZN 1 %s%14.5f\n",
            counter,AtomName,GeomCentList[j][1],GeomCentList[j][2],
            GeomCentList[j][3],&ResN_fr[CurFra][1],GeomCentEner[j]);
  }

/* following geometrical centers */
  for (j=1;j<=list2;j++)
  {
    counter++;
    i1 = sortIntAr[j];
    sprintf(AtomName,"CL%d",j); /* prev "CL%d\0" */
    fprintf(FilePath,"%7d %-8s%10.4f%10.4f%10.4f     XCL 1 %s%14.5f\n",
            counter,AtomName,GeomCentList[i1][1],GeomCentList[i1][2],
            GeomCentList[i1][3],&ResN_fr[CurFra][1],GeomCentEner[i1]);
  }

  fclose(FilePath);

  free_ivector(sortIntAr,1,clusCount);
  free_vector(sortFloatAr,1,clusCount);
  free_ivector(listClustMemb,1,maxNuClusMemb);
  free_vector(GeomCentEner,1,totNumbGeomCent);
  free_matrix(GeomCentList,1,totNumbGeomCent,1,3);
  free_ivector(ClustListLabel,1,NuSdClKe);
  free_matrix(GeomCentAr,1,NuSdClKe,1,3);

}
