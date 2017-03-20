#include <stdio.h>
#include <math.h>
#include "funct.h"
//#include <iomanip> //clangini

float VWFrEn(int FrAtNu,float **RoSFCo,float *BSMinC,float VWGrIn,float VWGrSi,
             float *FrVdWR_p3,float *FrVdWE_sr,int *VWGPoN,float ***VWGrRP_at,
             float ***VWGrRP_re,FILE *FPaOut,
	     int *print)
/* This function evaluates the van der Waals interaction energy between the
   current fragment and the receptor whose effect is precalculated in a grid :
   GrCo_i  "integer coordinate" along the x axis of the grid point just "below"
           the x coordinate of the current fragment atom
   GrCo_j  "integer coordinate" along the y axis of the grid point just "below"
           the y coordinate of the current fragment atom
   GrCo_k  "integer coordinate" along the z axis of the grid point just "below"
           the z coordinate of the current fragment atom
   Weig_t  weight factor for the computation of the van der Waals energy
   Weig_u  weight factor for the computation of the van der Waals energy
   Weig_v  weight factor for the computation of the van der Waals energy
   VWEnEv  evaluation of the van der Waals interaction energy */
{
  int i,GrCo_i,GrCo_j,GrCo_k;
  float Weig_t,Weig_u,Weig_v,VWEnEv,hlp1,hlp2,hlp3;

  //float at_term, rep_term;//clangini
  //at_term = 0.0; //clangini
  //rep_term = 0.0; //clangini
  //double Weig_t_double, Weig_u_double, Weig_v_double; // clangini
  //double VWEnEv_double = 0.0; // clangini

  VWEnEv=0.0;

  for (i=1;i<=FrAtNu;i++) {

    hlp1 = (RoSFCo[i][1]-(BSMinC[1]-VWGrIn))/VWGrSi;//necessary to fix loss of
    hlp2 = (RoSFCo[i][2]-(BSMinC[2]-VWGrIn))/VWGrSi;//precision bug. clangini
    hlp3 = (RoSFCo[i][3]-(BSMinC[3]-VWGrIn))/VWGrSi;//(thx to A.Vitalis)
    GrCo_i=ffloor(hlp1)+1;
    GrCo_j=ffloor(hlp2)+1;
    GrCo_k=ffloor(hlp3)+1;

    if ((GrCo_i>=1)&&(GrCo_i<VWGPoN[1])&&(GrCo_j>=1)&&(GrCo_j<VWGPoN[2])&&
        (GrCo_k>=1)&&(GrCo_k<VWGPoN[3])) {
      //Energy on the grid is computed with a trilinear interpolation. clangini
      //Weig_t=(RoSFCo[i][1]-(BSMinC[1]-VWGrIn+(GrCo_i-1)*VWGrSi))/VWGrSi;
      //Weig_u=(RoSFCo[i][2]-(BSMinC[2]-VWGrIn+(GrCo_j-1)*VWGrSi))/VWGrSi;
      //Weig_v=(RoSFCo[i][3]-(BSMinC[3]-VWGrIn+(GrCo_k-1)*VWGrSi))/VWGrSi;
      Weig_t=hlp1 - float(GrCo_i-1);
      Weig_u=hlp2 - float(GrCo_j-1);
      Weig_v=hlp3 - float(GrCo_k-1);

      //clangini START
      if((Weig_t>1)||(Weig_u>1)||(Weig_v>1)){
        std::cout<<"Here value larger than one"<<std::endl;
        std::cout<<"Weig_tuv: "<<Weig_t<<", "<<Weig_u<<", "<<Weig_v<<std::endl;
      }
      if((Weig_t<0)||(Weig_u<0)||(Weig_v<0)){
       std::cout<<"Here value smaller than zero"<<std::endl;
       std::cout<<"Weig_tuv: "<<Weig_t<<", "<<Weig_u<<", "<<Weig_v<<std::endl;
      }
      //clangini END

      VWEnEv=VWEnEv + FrVdWE_sr[i]*FrVdWR_p3[i] * (//attractive term. clangini
             (1-Weig_t)*(1-Weig_u)*(1-Weig_v)*
              VWGrRP_at[GrCo_i][GrCo_j][GrCo_k] +
             Weig_t*(1-Weig_u)*(1-Weig_v)*VWGrRP_at[GrCo_i+1][GrCo_j][GrCo_k] +
             Weig_t*Weig_u*(1-Weig_v)*VWGrRP_at[GrCo_i+1][GrCo_j+1][GrCo_k] +
             (1-Weig_t)*Weig_u*(1-Weig_v)*VWGrRP_at[GrCo_i][GrCo_j+1][GrCo_k] +
             (1-Weig_t)*(1-Weig_u)*Weig_v*VWGrRP_at[GrCo_i][GrCo_j][GrCo_k+1] +
             Weig_t*(1-Weig_u)*Weig_v*VWGrRP_at[GrCo_i+1][GrCo_j][GrCo_k+1] +
             Weig_t*Weig_u*Weig_v*VWGrRP_at[GrCo_i+1][GrCo_j+1][GrCo_k+1] +
             (1-Weig_t)*Weig_u*Weig_v*VWGrRP_at[GrCo_i][GrCo_j+1][GrCo_k+1] )
        + FrVdWE_sr[i]*FrVdWR_p3[i]*FrVdWR_p3[i] * ( //repulsive term. clangini
             (1-Weig_t)*(1-Weig_u)*(1-Weig_v)*
              VWGrRP_re[GrCo_i][GrCo_j][GrCo_k] +
             Weig_t*(1-Weig_u)*(1-Weig_v)*VWGrRP_re[GrCo_i+1][GrCo_j][GrCo_k] +
             Weig_t*Weig_u*(1-Weig_v)*VWGrRP_re[GrCo_i+1][GrCo_j+1][GrCo_k] +
             (1-Weig_t)*Weig_u*(1-Weig_v)*VWGrRP_re[GrCo_i][GrCo_j+1][GrCo_k] +
             (1-Weig_t)*(1-Weig_u)*Weig_v*VWGrRP_re[GrCo_i][GrCo_j][GrCo_k+1] +
             Weig_t*(1-Weig_u)*Weig_v*VWGrRP_re[GrCo_i+1][GrCo_j][GrCo_k+1] +
             Weig_t*Weig_u*Weig_v*VWGrRP_re[GrCo_i+1][GrCo_j+1][GrCo_k+1] +
             (1-Weig_t)*Weig_u*Weig_v*VWGrRP_re[GrCo_i][GrCo_j+1][GrCo_k+1] ) ;
    }
#ifndef NOWARN
    else if(print==0) {
      ++(*print);
      /* can produce a lot of output in some cases */
      /* fprintf(FPaOut,"WARNING One fragment is not completely in the "); */
      fprintf(FPaOut,"WARNING at least one fragment is not completely in the ");
      fprintf(FPaOut,"van der Waals energy grid\n\n");
      return 1e6; /* New -> give penalty */
    }
#else
    {
      return 1e6; /* New -> give penalty */
    }
#endif
  }

  // //clangini START
  // if(VWEnEv == -10390967296.00){
  //   //std:cout<<"GrCo_i: "<<GrCo_i<<std::endl;
  //   //std::cout <<"Weig_u_double"
  //   VWEnEv=0.0;
  //   for (i=1;i<=FrAtNu;i++) {
  //     std::cout<<"Atom "<<i<<std::endl;
  //     GrCo_i=ffloor((RoSFCo[i][1]-(BSMinC[1]-VWGrIn))/VWGrSi)+1;
  //     GrCo_j=ffloor((RoSFCo[i][2]-(BSMinC[2]-VWGrIn))/VWGrSi)+1;
  //     GrCo_k=ffloor((RoSFCo[i][3]-(BSMinC[3]-VWGrIn))/VWGrSi)+1;
  //     //std::cout<<"GrCo_ijk: "<<GrCo_i<<", "<<GrCo_j<<", "<<GrCo_k<<std::endl;
  //     if ((GrCo_i>=1)&&(GrCo_i<VWGPoN[1])&&(GrCo_j>=1)&&(GrCo_j<VWGPoN[2])&&
  //         (GrCo_k>=1)&&(GrCo_k<VWGPoN[3])) {
  //       //Energy on the grid is computed with a trilinear interpolation. clangini
  //       Weig_t=(RoSFCo[i][1]-(BSMinC[1]-VWGrIn+(GrCo_i-1)*VWGrSi))/VWGrSi;
  //       Weig_u=(RoSFCo[i][2]-(BSMinC[2]-VWGrIn+(GrCo_j-1)*VWGrSi))/VWGrSi;
  //       Weig_v=(RoSFCo[i][3]-(BSMinC[3]-VWGrIn+(GrCo_k-1)*VWGrSi))/VWGrSi;
  //       //Weig_v=(static_cast<double>(RoSFCo[i][3])-(static_cast<double>(BSMinC[3])-static_cast<double>(VWGrIn)
  //       //        +static_cast<double>((GrCo_k-1))*static_cast<double>(VWGrSi)))/static_cast<double>(VWGrSi);
  //       std::cout<<"GrCo_ijk: "<<GrCo_i<<", "<<GrCo_j<<", "<<GrCo_k<<std::endl;
  //       std::cout<<"Weig_tuv: "<<Weig_t<<", "<<Weig_u<<", "<<Weig_v<<std::endl;
  //       if(i==15){
  //         std::cout<<"VWGrIn: "<<VWGrIn<<std::endl;
  //         std::cout<<"VWGrSi: "<<VWGrSi<<std::endl;
  //         std::cout<<"RoSFCo[15][3]: "<<RoSFCo[i][3]<<std::endl;
  //         std::cout<<"BSMinC[3]: "<<BSMinC[3]<<std::endl;
  //         std::cout<<"Intermed: "<<(RoSFCo[i][3]-(BSMinC[3]-VWGrIn+(GrCo_k-1)*VWGrSi))<<std::endl;
  //         std::cout<<"Intermed 2: "<<(RoSFCo[i][3]-BSMinC[3])<<std::endl;
  //         std::cout<<"Intermed 3: "<<(-VWGrIn+(GrCo_k-1)*VWGrSi)<<std::endl;
  //         std::cout<<"Intermed 4: "<<(BSMinC[3]-VWGrIn+(GrCo_k-1)*VWGrSi)<<std::endl;
  //         std::cout<<"Intermed 5: "<< (((RoSFCo[i][3])-(BSMinC[3]-VWGrIn+(GrCo_k-1)*VWGrSi))-
  //                                 (-(RoSFCo[i][3])+(BSMinC[3]-VWGrIn+(GrCo_k-1)*VWGrSi)))/2<<std::endl;
  //         //Weig_v = 0.0;
  //       }
  //       at_term = (//attractive term. clangini
  //              (1-Weig_t)*(1-Weig_u)*(1-Weig_v)*
  //               VWGrRP_at[GrCo_i][GrCo_j][GrCo_k] +
  //              Weig_t*(1-Weig_u)*(1-Weig_v)*VWGrRP_at[GrCo_i+1][GrCo_j][GrCo_k] +
  //              Weig_t*Weig_u*(1-Weig_v)*VWGrRP_at[GrCo_i+1][GrCo_j+1][GrCo_k] +
  //              (1-Weig_t)*Weig_u*(1-Weig_v)*VWGrRP_at[GrCo_i][GrCo_j+1][GrCo_k] +
  //              (1-Weig_t)*(1-Weig_u)*Weig_v*VWGrRP_at[GrCo_i][GrCo_j][GrCo_k+1] +
  //              Weig_t*(1-Weig_u)*Weig_v*VWGrRP_at[GrCo_i+1][GrCo_j][GrCo_k+1] +
  //              Weig_t*Weig_u*Weig_v*VWGrRP_at[GrCo_i+1][GrCo_j+1][GrCo_k+1] +
  //              (1-Weig_t)*Weig_u*Weig_v*VWGrRP_at[GrCo_i][GrCo_j+1][GrCo_k+1] );
  //       rep_term = (//repulsive term. clangini
  //                   (1-Weig_t)*(1-Weig_u)*(1-Weig_v)*
  //                   VWGrRP_re[GrCo_i][GrCo_j][GrCo_k] +
  //                   Weig_t*(1-Weig_u)*(1-Weig_v)*VWGrRP_re[GrCo_i+1][GrCo_j][GrCo_k] +
  //                   Weig_t*Weig_u*(1-Weig_v)*VWGrRP_re[GrCo_i+1][GrCo_j+1][GrCo_k] +
  //                   (1-Weig_t)*Weig_u*(1-Weig_v)*VWGrRP_re[GrCo_i][GrCo_j+1][GrCo_k] +
  //                   (1-Weig_t)*(1-Weig_u)*Weig_v*VWGrRP_re[GrCo_i][GrCo_j][GrCo_k+1] +
  //                   Weig_t*(1-Weig_u)*Weig_v*VWGrRP_re[GrCo_i+1][GrCo_j][GrCo_k+1] +
  //                   Weig_t*Weig_u*Weig_v*VWGrRP_re[GrCo_i+1][GrCo_j+1][GrCo_k+1] +
  //                   (1-Weig_t)*Weig_u*Weig_v*VWGrRP_re[GrCo_i][GrCo_j+1][GrCo_k+1] );
  //       VWEnEv=VWEnEv + FrVdWE_sr[i]*FrVdWR_p3[i] * at_term
  //                     + FrVdWE_sr[i]*FrVdWR_p3[i]*FrVdWR_p3[i] * rep_term ;
  //
  //     std::cout<<"at_term: "<<at_term<<std::endl;
  //     std::cout<<"rep_term: "<<rep_term<<std::endl;
  //     std::cout<<VWGrRP_re[GrCo_i][GrCo_j][GrCo_k]<<std::endl;
  //     std::cout<<VWGrRP_re[GrCo_i][GrCo_j][GrCo_k+1]<<std::endl;
  //     std::cout<<VWGrRP_re[GrCo_i][GrCo_j+1][GrCo_k]<<std::endl;
  //     std::cout<<VWGrRP_re[GrCo_i][GrCo_j+1][GrCo_k+1]<<std::endl;
  //     std::cout<<VWGrRP_re[GrCo_i+1][GrCo_j][GrCo_k]<<std::endl;
  //     std::cout<<VWGrRP_re[GrCo_i+1][GrCo_j][GrCo_k+1]<<std::endl;
  //     std::cout<<VWGrRP_re[GrCo_i+1][GrCo_j+1][GrCo_k]<<std::endl;
  //     std::cout<<VWGrRP_re[GrCo_i+1][GrCo_j+1][GrCo_k+1]<<std::endl;
  //     std::cout<<"Running energy: "<<VWEnEv<<std::endl;
  //     std::cout<<"----------------"<<std::endl;
  //     //std::cout<<"en double: "<<VWEnEv_double<<std::endl;
  //    }
  //  }
  //}
  //clangini END

  return VWEnEv;

}
