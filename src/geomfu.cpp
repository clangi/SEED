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
#include "funct.h"

double VeNorm(double a1,double a2,double a3)
/* This function returns the norm of the vector A */
{
  double norm;

  norm=sqrtf(a1*a1+a2*a2+a3*a3);

  return norm;
}



/*void NormVe(float *a1,float *a2,float *a3)*/
/* This function normalizes the vector A */
/*{
  float norm;

  norm=VeNorm(*a1,*a2,*a3);

  *a1=*a1/norm;
  *a2=*a2/norm;
  *a3=*a3/norm;
}*/

void NormVe(double *a1, double *a2, double *a3)
/* Overloaded NormVe for doubles. clangini */
{
  double norm;
  norm = sqrtf((*a1)*(*a1)+(*a2)*(*a2)+(*a3)*(*a3));

  *a1=*a1/norm;
  *a2=*a2/norm;
  *a3=*a3/norm;
}


void PoCoVe(double a1,double a2,double a3,double b1,double b2,double b3,double len,
            double *c1,double *c2,double *c3)
/* This function computes the coordinates of C which is along the vector AB
   and at a distance of len from A */
{
  double vec1,vec2,vec3;

  vec1=b1-a1;
  vec2=b2-a2;
  vec3=b3-a3;

  NormVe(&vec1,&vec2,&vec3);

  *c1=a1+len*vec1;
  *c2=a2+len*vec2;
  *c3=a3+len*vec3;
}



void VectPr(double a1,double a2,double a3,double b1,double b2,double b3,double *c1,
            double *c2,double *c3)
/* This function computes the vector product C from A and B (C=A^B) */
{
  *c1=a2*b3-a3*b2;
  *c2=a3*b1-a1*b3;
  *c3=a1*b2-a2*b1;
}



double PlaAng(double o1,double o2,double o3,double a1,double a2,double a3,double b1,
             double b2,double b3)
/* This function returns the angle between the vectors OA and OB. The angle
   value ranges from 0 to PI */
{
  double vecA1,vecA2,vecA3,vecB1,vecB2,vecB3,normA,normB,ScaPro,StorVa;

  vecA1=a1-o1;
  vecA2=a2-o2;
  vecA3=a3-o3;

  vecB1=b1-o1;
  vecB2=b2-o2;
  vecB3=b3-o3;

  normA=VeNorm(vecA1,vecA2,vecA3);
  normB=VeNorm(vecB1,vecB2,vecB3);

  ScaPro=vecA1*vecB1+vecA2*vecB2+vecA3*vecB3;

  StorVa=ScaPro/(normA*normB);

  if (StorVa>1.0) StorVa=1.0;
  if (StorVa<(-1.0)) StorVa=-1.0;

  return acosf(StorVa);
}



void RoArVe(double a1,double a2,double a3,double vec1,double vec2,double vec3,
            double ang,double *b1,double *b2,double *b3)
/* This function rotates the point A around the vector Vec through origin for
   an angle ang in radian. The result is given in B */
{
  double v1,v2,v3,rtl1,rtl2,rtl3,sinl,cosl,rl,help1,help2,help3;

  NormVe(&vec1,&vec2,&vec3);
  sinl=sinf(ang);
  cosl=1.0-cosf(ang);
  help1=a1;
  help2=a2;
  help3=a3;
  VectPr(vec1,vec2,vec3,help1,help2,help3,&v1,&v2,&v3);
  rl=a1*vec1+a2*vec2+a3*vec3;
  rtl1=rl*vec1;
  rtl2=rl*vec2;
  rtl3=rl*vec3;

  *b1=a1+v1*sinl+(rtl1-a1)*cosl;
  *b2=a2+v2*sinl+(rtl2-a2)*cosl;
  *b3=a3+v3*sinl+(rtl3-a3)*cosl;
}



void AxisVe(double a1,double a2,double a3,double b1,double b2,double b3,
            double c1,double c2,double c3,double *d1,double *d2,double *d3)
/* This function computes the normed vector D normal to the plane ABC */
{
  double ba1,ba2,ba3,bc1,bc2,bc3; /* ,norm; */

  ba1=a1-b1;
  ba2=a2-b2;
  ba3=a3-b3;

  bc1=c1-b1;
  bc2=c2-b2;
  bc3=c3-b3;

  VectPr(ba1,ba2,ba3,bc1,bc2,bc3,d1,d2,d3);

  NormVe(d1,d2,d3);
}



void RotPla(double a1,double a2,double a3,double b1,double b2,double b3,double c1,
            double c2,double c3,double angl,double *d1,double *d2,double *d3)
/* This function computes the coordinates of the point D which results from
   the rotation of C with angle angl in radian around B in the ABC plane */
{
  double bc1,bc2,bc3,arot1,arot2,arot3,bc1ret,bc2ret,bc3ret;

  AxisVe(a1,a2,a3,b1,b2,b3,c1,c2,c3,&arot1,&arot2,&arot3);

  bc1=c1-b1;
  bc2=c2-b2;
  bc3=c3-b3;

  RoArVe(bc1,bc2,bc3,arot1,arot2,arot3,angl,&bc1ret,&bc2ret,&bc3ret);

  *d1=bc1ret+b1;
  *d2=bc2ret+b2;
  *d3=bc3ret+b3;
}



double DihAng(double a1,double a2,double a3,double b1,double b2,double b3,
             double c1,double c2,double c3,double d1,double d2,double d3)
/* This function returns the dihedral angle between the ABC and BCD planes.
   The dihedral angle ranges from 0 to PI */
{
  double ba1,ba2,ba3,bc1,bc2,bc3,VePrA1,VePrA2,VePrA3,cb1,cb2,cb3,cd1,cd2,cd3,
        VePrB1,VePrB2,VePrB3,ScaPro,normA,normB,StorVa;

  ba1=a1-b1;
  ba2=a2-b2;
  ba3=a3-b3;

  bc1=c1-b1;
  bc2=c2-b2;
  bc3=c3-b3;

  VectPr(ba1,ba2,ba3,bc1,bc2,bc3,&VePrA1,&VePrA2,&VePrA3);
  normA=VeNorm(VePrA1,VePrA2,VePrA3);

  cb1=b1-c1;
  cb2=b2-c2;
  cb3=b3-c3;

  cd1=d1-c1;
  cd2=d2-c2;
  cd3=d3-c3;

  VectPr(cb1,cb2,cb3,cd1,cd2,cd3,&VePrB1,&VePrB2,&VePrB3);
  normB=VeNorm(VePrB1,VePrB2,VePrB3);

  ScaPro=VePrA1*VePrB1+VePrA2*VePrB2+VePrA3*VePrB3;

  StorVa=ScaPro/(normA*normB);

  if (StorVa>1.0) StorVa=1.0;
  if (StorVa<(-1.0)) StorVa=-1.0;

  return acosf(StorVa);
}



//inline float DistSq(float a1,float a2,float a3,float b1,float b2,float b3)
/* This function returns the squared distance between A and B */
//{
  /* float DisSqu; */

    /* DisSqu= */
//  return (b1-a1)*(b1-a1)+(b2-a2)*(b2-a2)+(b3-a3)*(b3-a3);

  /* return DisSqu; */

//}
/*inline function has to be defined in header file-> moved to funct.h clangini*/


double TrProd(double ox,double oy,double oz,double ax,double ay,double az,
             double bx,double by,double bz,double cx,double cy,double cz)
/* This function returns the triple scalar product of the vectors oa, ob
   and oc */
{
  double vec1x,vec1y,vec1z,vec2x,vec2y,vec2z,vec3x,vec3y,vec3z,norm1,norm2,
        norm3,TripPr;

  vec1x=ax-ox;
  vec1y=ay-oy;
  vec1z=az-oz;

  vec2x=bx-ox;
  vec2y=by-oy;
  vec2z=bz-oz;

  vec3x=cx-ox;
  vec3y=cy-oy;
  vec3z=cz-oz;

  norm1=VeNorm(vec1x,vec1y,vec1z);
  norm2=VeNorm(vec2x,vec2y,vec2z);
  norm3=VeNorm(vec3x,vec3y,vec3z);

  TripPr=(1.0/(norm1*norm2*norm3))*(vec1x*(vec2y*vec3z-vec2z*vec3y)
                -vec2x*(vec1y*vec3z-vec1z*vec3y)
                +vec3x*(vec1y*vec2z-vec1z*vec2y));

  return fabsf(TripPr);

}
