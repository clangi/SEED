#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"
#include "funct.h"

float MaxVector(float *Vector, int VectorStart, int VectorEnd)
/* Find the maximum of Array[*] between ArrayStart and ArrayEnd (float) */
{
  int i;
  float max;

  max=Vector[VectorStart];

  for (i=VectorStart;i<=VectorEnd;i++) 
    max = (Vector[i]>max) ? Vector[i] : max;

  return max;
}

double MaxDVector(double *Vector, int VectorStart, int VectorEnd)
/* Find the maximum of Array[*] between ArrayStart and ArrayEnd (double) */
{
  int i;
  double max;

  max=Vector[VectorStart];

  for (i=VectorStart;i<=VectorEnd;i++) 
    max = (Vector[i]>max) ? Vector[i] : max;

  return max;
}

double MinVector(double *Vector, int VectorStart, int VectorEnd)
/* Find the maximum of Array[*] between ArrayStart and ArrayEnd */
{
  int i;
  double min;

  min=Vector[VectorStart];

  for (i=VectorStart;i<=VectorEnd;i++)
    min = (Vector[i]<min) ? Vector[i] : min;

  return min;
}

int getColumnFrom2DArray(float **Array2D, int ColumnNum, 
                         int start, int end, double *Column)
/* Get column "ColumnNum" from the matrix Array2D. Start from 
   line "start" and go up to line "end" */
{
  int i;

  for (i=start;i<=end;i++)
    Column[i]=Array2D[i][ColumnNum];

  if (i == end) 
    return 1;
  else
    return 0;
}

int SumVectors(double *vect1, float *vect2, int start, int end, double *sum)
/* Sum vect1 and vect2 between "start" and "end" and put the sum 
   in the vector sum */
{
  int i;

  for (i=start;i<=end;i++)
    sum[i]=vect1[i]+vect2[i];

  if (i == end) 
    return 1;
  else
    return 0;
}

int SubstractVectors(double *vect1, float *vect2, int start, int end, 
    double *sum)
/* Sum vect2 from vect1 between "start" and "end" and put the result 
   in the vector sum */
{
  int i;

  for (i=start;i<=end;i++)
    sum[i]=vect1[i]-vect2[i];

  if (i == end) 
    return 1;
  else
    return 0;
}

struct point c_per_pt(double c, struct point pt)
/* Multiplicate the struct point pt for a constant c */
{
  struct point result;

  result.x = c * pt.x;
  result.y = c * pt.y;
  result.z = c * pt.z;

  return result;
}

struct point pt_plus_pt(struct point pt1, struct point pt2)
/* Sum the 2 struct point pt1 and pt2 and return the result */
{
  struct point result;

  result.x = pt1.x + pt2.x;
  result.y = pt1.y + pt2.y;
  result.z = pt1.z + pt2.z;

  return result;
}

struct point pt_minus_pt(struct point pt1, struct point pt2)
/* Sum the struct point pt2 from struct point pt1 and return the result */
{
  struct point result;

  result.x = pt1.x - pt2.x;
  result.y = pt1.y - pt2.y;
  result.z = pt1.z - pt2.z;

  return result;
}

double pt_scal_pt(struct point pt1, struct point pt2)
/* Get the scalar product of the 2 struct point pt1 and pt2 
   and return the result */
{
  return (pt1.x * pt2.x + pt1.y * pt2.y +  pt1.z * pt2.z);
}

void pt_eq_pt(struct point *Ppt1, struct point pt2)
/* Make struct point pt1 equal to struct point pt2 */
{
  Ppt1->x = pt2.x;
  Ppt1->y = pt2.y;
  Ppt1->z = pt2.z;
}

void DSort(int N,int *IndArr,double *SorArr)
/* This function sorts an array of double SorArr and gives also the new order 
   in the array IndArr. The result goes from the smallest to the largest. */
{
  int m,k,i,j,l,nn,tempin,loopva;
  float aln2i,lognb2,temp;

  aln2i=1.0/0.69314718;
  lognb2= (float) ((int) (logf(N)*aln2i+5.e-7));

  m=N;
  for (nn=1;nn<=lognb2;nn++) {
    m=m/2;
    k=N-m;
    for (j=1;j<=k;j++) {
      i=j;
      loopva=1;
      while (loopva) {
        loopva=0;
        l=i+m;
        if (SorArr[l]<SorArr[i]) {
          temp=SorArr[i];
          SorArr[i]=SorArr[l];
          SorArr[l]=temp;
          tempin=IndArr[i];
          IndArr[i]=IndArr[l];
          IndArr[l]=tempin;
          i=i-m;
          if (i>=1) 
            loopva=1;
        }
      }
    }
  }

}
