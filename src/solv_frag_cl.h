#ifndef _SOLV_FRAG_CL_H
#define _SOLV_FRAG_CL_H

int FragSolvEn_cl(int FrAtNu,double **FrCoor,double *FrPaCh,
               double *FrVdWR, double *FrRad, double *FrRad2, double *FrRadOut,
               double *FrRadOut2, double **Frdist2,int Nsurfpt_fr,
               struct point *surfpt_fr_orig,double WaMoRa,double GrSiSo,
               double Ksolv,double pi4,double *PFrSolvEn,char *EmpCorrB,FILE*FPaOut);
int SAS_Volume_cl(int ReAtNu,double **ReCoor,double *ReRadOut,
              double *ReRadOut2,struct point Min,
              double GrSiSo,int NStartGridx,int NStartGridy,
              int NStartGridz,int NGridx,int NGridy,int NGridz,
              char ***GridMat);
int Excl_Grid_cl(int ReAtNu,double **ReCoor,struct point Min,double *ReRadOut,
              double *ReRadOut2,double WaMoRa,double GrSiSo,
              int NStartGridx,int NStartGridy,int NStartGridz,
              int NGridx,int NGridy,int NGridz,
              double UnitVol,char ***GridMat,int Nsurfpt,
              struct point *surfpt,double *SelfVol,double *SelfVol_corrB,
	            char *EmpCorrB);
int Get_Self_En_Fr_cl(int FrAtNu,double **RoSFCo,double *FrPaCh,double *FrRad,
                   double *FrRad2,double *XGrid,double *YGrid,double *ZGrid,
                   double UnitVol,double Ksolv,double pi4,char ***FrGridMat,
                   int nxminFr,int nyminFr,int nzminFr,
                   int nxmaxFr,int nymaxFr,int nzmaxFr,
                   double *FrSelfVol,double *FrEffRad,double *PFrSelfEn,
		               double *FrSelfVol_corrB,char *EmpCorrB,FILE * FPaOut);

#endif
