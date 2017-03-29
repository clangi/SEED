#include <iostream>
#include <cstdio>
#include "quaternion.h"
#define sqrtf sqrt
#define sinf sin
#define cosf cos
#define fabsf fabs
#define ffloor floor
#define acosf acos
#define expf exp
#define modff modf
#define logf log

/*void ReInFi(char *InpFil,char *RecFil,int *BSResN,int **BSReNu,
            int *FragNu,char ***FrFiNa,char *TREFiP,double *SphAng,
            int *SphPoN,int *NuRoAx,double *VdWFaB,double *CoDieV,int *CoDieP,
            double *CoGrIn,double *CoGrSi,char *OutFil,double *BuEvFa,
            double **FrMaEn,double *PsSpRa,double *GrSiCu_en,
            int *FiNuMa,double *GrInSo,double *GrSiSo,double *WaMoRa,
            int *NPtSphere,double *DielWa,double *DielRe,char *ReDesoAlg,
            char *DesoMapAcc,char *DesoMapFile,char *FDexe,char *FDdir,
            double *ReSurfDens_apol,double *PtDensFr,double *Sphere_apol,
            double *NCutapolRatio*//*int *NCutapol*//*,double *ScaleDeso,
            double *ScaleVDW,double **SimWei,double *SimExp,
            double *SimCut,double **FrMaEn_sd,double *SimExp_sd,double *SimCut_sd,
            int *BSMeNu,int **BSMeAN,double ***BSMeVE,char *CoGrAcc,
            char *CoGrFile,char *EvalEn,char *Solv_typ,
            char *SpPoCh_opt,double *SpPoCh_cent,double *SpPoCh_rad,
            double *SFDeso_fr,double *SFDeso_re,double *SFVWEn,double *SFIntElec,
            int *NuClusMem,double *RedRPV_rp,double *RedRPV_nkvDens,double *ScMaBump, *//*int *RedRPV_nkv, dey new*/
            /*double *MuFaVdWCoff_ap,int *NuLiEnClus,char ***ApPoChoi,
            double *VWGrIn,double *VWGrSi,double *BumpFaCut,char *VWGrAcc,
            char *VWGrFile,int *MaxPosClus,int *PrintLev,
            int *NumbAT,char ***AtTyAr,int **AtENAr,double **VdWRad,
            double **VdWEne,double ***BLAtTy,int *distrPointBSNumb,
	    double ***distrPointBS,double *angle_rmin,double *angle_rmax,
	    double *mult_fact_rmin,double *mult_fact_rmax,char *EmpCorrB,
            char *gc_opt,int *gc_reprke,double *gc_cutclus,double *gc_endifclus,
            double *gc_weighneg,double *gc_weighpos,int *gc_maxwrite,
            char *write_pproc_opt,char *write_pproc_chm_opt,int *CorrFiNumb);*/
void ReInFi(char *InpFil,char *RecFil,int *BSResN,int **BSReNu,
            char *FrFiNa,char *TREFiP,double *SphAng,
            int *SphPoN,int *NuRoAx,double *VdWFaB,double *CoDieV,int *CoDieP,
            double *CoGrIn,double *CoGrSi,char *OutFil,double *BuEvFa,
            double *FrMaEn,double *PsSpRa,double *GrSiCu_en,
            int *FiNuMa,double *GrInSo,double *GrSiSo,double *WaMoRa,
            int *NPtSphere,double *DielWa,double *DielRe,char *ReDesoAlg,
            char *DesoMapAcc,char *DesoMapFile,char *FDexe,char *FDdir,
            double *ReSurfDens_apol,double *PtDensFr,double *Sphere_apol,
            double *NCutapolRatio/*int *NCutapol*/,double *ScaleDeso,
            double *ScaleVDW,double **SimWei,double *SimExp,
            double *SimCut,double *FrMaEn_sd,double *SimExp_sd,double *SimCut_sd,
            int *BSMeNu,int **BSMeAN,double ***BSMeVE,char *CoGrAcc,
            char *CoGrFile,char *EvalEn,char *Solv_typ,
            char *SpPoCh_opt,double *SpPoCh_cent,double *SpPoCh_rad,
            double *SFDeso_fr,double *SFDeso_re,double *SFVWEn,double *SFIntElec,
            int *NuClusMem,int *NuPosMem,double *RedRPV_rp,double *RedRPV_nkvDens,double *ScMaBump, /*int *RedRPV_nkv, dey new*/
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
            char *write_sumtab_opt,double **AtWei);

void ReReFi_mol2(char *RecFil,int *ReAtNu,int *ReBdNu,int *ReReNu,
                 char ***ReAtEl,double ***ReCoor,char ***ReAtTy,int **ReResN,
                 double **RePaCh,int ***ReBdAr);
/*void ReFrFi_mol2(int CurFra,char **FrFiNa,char *FragNa,int *FrAtNu,int *FrBdNu,
                 char ***FrAtEl,double ***FrCoor,char ***FrAtTy,double **FrPaCh,
                 int ***FrBdAr,char ***FrBdTy,char *FrSubN,char *FrSubC,
                 int *FrCoNu); Redefined. clangini*/
int count_heavy_atom(int *FrAtEl_nu,int FrAtNu);
double molecular_weight(int *FrAtEl_nu,int FrAtNu,double *AtWei);
bool set_align_ref(double **FrCoor,int FrAtNu,double **FrAlRef,double **FrAlSet);
void struct_align(double **FrCoor,int FrAtNu, double **FrAlRef,
                double **FrAlSet,int AlAtNu);
template<class T>
void align_3D(int n,double **FrAlSet,double **FrAlRef,Quaternion<T> **qrot,
                double *tvec,double *cto,double *ctn);
void CentroidCalc(double **Coord,int nAt,double *centroid);

double **convert2pp_return(int irl, int irh,int icl, int ich,double *first,
                            double *mat_pp_rows[]);
void print_pp(int irl, int irh,int icl, int ich,double **mat);
/* CLANGINI 2016 */
int ReFrFi_mol2(std::istream *inStream, std::streampos *strPos,
                int *SkiFra,int *CurFraTot,char *FragNa,std::string &FragNa_str,
                int *FrAtNu,int *FrBdNu,
                char ***FrAtEl,double ***FrCoor,char ***FrAtTy,  char ***FrSyAtTy,double **FrPaCh,
                int ***FrBdAr,char ***FrBdTy,char *FrSubN,char *FrSubC,
                int *FrCoNu,char ***SubNa, std::string &AlTySp);
void ReFrFi_mol2_syb(int CurFra,char **FrFiNa,char *FragNa,int *FrAtNu,int *FrBdNu,
                 char ***FrAtEl,double ***FrCoor,char ***FrAtTy, char ***FrSyAtTy,double **FrPaCh,
                 int ***FrBdAr,char ***FrBdTy,char *FrSubN,char *FrSubC,
                 int *FrCoNu);
void write_mol2_clus_pproc_syb(int CurFra,int NuPosSdCl,int FrAtNu,char **FrAtTy, char **FrSyAtTy,
                           double ***FrCoPo,int *FrPosAr_sort,char **ResN_fr,
                           char **FrFiNa_out,int FrBdNu,char **FrAtEl,
                           double *FrPaCh,int **FrBdAr,char **FrBdTy,
                           char *FragNa,int *SdClusAr_sort,double *To_s_ro);
void append_pose_to_mol2_old(FILE *FilePa,char *FragNa,int FrAtNu,int FrBdNu,int imol,
												 char **FrAtEl,double ***FrCoPo,int *FrPosAr_sort,char **FrSyAtTy,
												 int CurFra,char **ResN_fr,int **FrBdAr,char **FrBdTy,
												 int *SdClusAr_sort,double *To_s_ro,double *FrPaCh);
void write_mol2_clus_pproc_separate(int CurFra,int NuPosSdCl,int FrAtNu,char **FrAtTy,
													 char **FrSyAtTy,
                           double ***FrCoPo,int *FrPosAr_sort,char **ResN_fr,
                           char **FrFiNa_out,int FrBdNu,char **FrAtEl,
                           double *FrPaCh,int **FrBdAr,char **FrBdTy,
                           char *FragNa,int *SdClusAr_sort,double *To_s_ro);
void append_to_mol2(int CurFra,int NuPosSdCl,int FrAtNu,char **FrAtTy,
                    char **FrSyAtTy,double ***FrCoPo,int *FrPosAr_sort,
                    char **SubNa,char *FrFiNa_out,int FrBdNu,char **FrAtEl,
                    double *FrPaCh,int **FrBdAr,char **FrBdTy,
                    char *FragNa,int *SdClusAr_sort,double *To_s_ro, char *WriPat);
void append_pose_to_mol2(FILE *FilePa,char *FragNa,/*int FragNa_count,*/
                         int FrAtNu,int FrBdNu,int imol,
                         char **FrAtEl,double ***FrCoPo,int Fr_nu,char **FrSyAtTy,
                         char **FrAtTy,int CurFra,int **FrBdAr,char **FrBdTy,
                         int SdClu,double To_s_ro,double *FrPaCh,char **SubNa,
                         std::string const& AlTySp);
#ifdef DEBUG_CL
void append_pose_to_mol2_double(FILE *FilePa,char *FragNa,int FrAtNu,int FrBdNu,int imol,
                char **FrAtEl,double **FrCoPo,int Fr_nu,char **FrSyAtTy,
                char **FrAtTy,int CurFra,int **FrBdAr,char **FrBdTy,
                int SdClu,double To_s_ro,double *FrPaCh,char **SubNa,
                std::string const&AlTySp);
#endif
void append_pose_to_mol2(FILE *FilePa,char *FragNa,int FrAtNu,int FrBdNu,int imol,
 char **FrAtEl,double **FrCoPo,int Fr_nu,char **FrSyAtTy,
 char **FrAtTy,int CurFra,int **FrBdAr,char **FrBdTy,
 int SdClu,double To_s_ro,double *FrPaCh,char **SubNa,
 std::string const&AlTySp); // Overloaded: FrCoPo passed
 //as **double instead of ***double. clangini

 inline void fgets_wrapper(char * str, int num, FILE * stream ){
   /* This is a wrapper function for standard fgets,in order to
      correclty ignore the return value */
   if(fgets(str,num,stream )!=NULL){
     return;
   }
   else {
     std::cout << "Could not read input string." << std::endl;
   }
 }

/* CLANGINI 2016 END */
void FiAlHy(int FrAtNu,int FrBdNu,int *FrAtEl_nu,int **FrBdAr,int *AliHyd,
            int *NonAliHy);
double VeNorm(double a1,double a2,double a3);
void NormVe(double *a1,double *a2,double *a3);
/*clangini START*/
/* overloaded NormVe for double */
void NormVe(double *a1, double *a2, double *a3);
/*clangini END*/
void PoCoVe(double a1,double a2,double a3,double b1,double b2,double b3,double len,
            double *c1,double *c2,double *c3);
void VectPr(double a1,double a2,double a3,double b1,double b2,double b3,double *c1,
            double *c2,double *c3);
double PlaAng(double o1,double o2,double o3,double a1,double a2,double a3,double b1,
             double b2,double b3);
void RoArVe(double a1,double a2,double a3,double vec1,double vec2,double vec3,
            double ang,double *b1,double *b2,double *b3);
void AxisVe(double a1,double a2,double a3,double b1,double b2,double b3,
            double c1,double c2,double c3,double *d1,double *d2,double *d3);
void RotPla(double a1,double a2,double a3,double b1,double b2,double b3,double c1,
            double c2,double c3,double angl,double *d1,double *d2,double *d3);
double DihAng(double a1,double a2,double a3,double b1,double b2,double b3,
             double c1,double c2,double c3,double d1,double d2,double d3);
inline double DistSq(double a1,double a2,double a3,double b1,double b2,double b3){
  return (b1-a1)*(b1-a1)+(b2-a2)*(b2-a2)+(b3-a3)*(b3-a3);
}
void SeedFr(int ReCuVe,double **ReVeCo,int FrCuVe,double **FrVeCo,int FrAtNu,
            double **FrCoor,int *ReDATy,double **SeFrCo,int *FrAtoTyp_nu,
            int *FrDAAt,int *ReAtoTyp_nu,int *ReDAAt,double **BLAtTy,
            FILE *FPaOut);
void write_mol2(char *WriPat,char *WriNam,int xxAtNu,int xxBdNu,char **xxAtEl,
                double **xxCoor,char **xxAtTy,double *xxPaCh,int **xxBdAr,
                char **xxBdTy,char *FrSubN,char *FrSubC,char **ResN_fr,
                int CurFra);
void RoSeFr(int ReCuVe,int *RexxAt,double **ReCoor,int FrCuVe,int *FrxxAt,
            int FrAtNu,double **SeFrCo,double AnglRo,double **RoSFCo);
int AssiRE(int NumbAT,char **AtTyAr,double *VdWRad,double *VdWEne,int xxAtNu,
	   char **xxAtTy,double *xxVdWR,double *xxVdWE,FILE *FPaOut,
	   int *AtENAr,int *xxAtEl_nu,int *xxAtoTyp_nu);
void VdWRaM(int NumbAT,double *VdWRad,double *LaVdWR);
void ReMMCo(int ReAtNu,double **ReCoor,double *ReMaxC,double *ReMinC);
void ReLAIC(int ReAtNu,double **ReCoor,double LaVdWR,double *ReMaxC,double *ReMinC,
            int *CubNum,int ****CubFAI,int ****CubLAI,int **CubLiA);
int CoNuBu(int *CubNum,int ***CubFAI,int ***CubLAI,int *CubLiA,double LaVdWR,
           double *ReMaxC,double *ReMinC,double **ReCoor,double *ReVdWR, int FrAtNu,
           double **RoSFCo,double *FrVdWR,double VdWFaB,int BumpMa,int RecHyN,
           int FraHyN,double BuEvFa);
void BSMMCo(double **ReCoor,int BSAtNu,int *BSAtLi,double *BSMaxC,double *BSMinC);
void CoGReP(int ReAtNu,double **ReCoor,double *RePaCh,double CoDieV,
            int CoDieP,double CoGrIn,double CoGrSi,double *BSMinC,int *CoGPoN,
            double ***CoGrRP);
double CoFrEn(int FrAtNu,double *FrPaCh,double **RoSFCo,double *BSMinC,
             double CoGrIn,double CoGrSi,int *CoGPoN,double ***CoGrRP,
             FILE *FPaOut,int * print);
void SqRoEn(int xxAtNu,double *xxVdWE,double *xxVdWE_sr);
void ReLAIC_en(int ReReNu,double **ReCoor,double GrSiCu_en,double *ReMinC,
               int *CubNum_en,int ***CubFAI_en,int ***CubLAI_en,int *CubLiA_en,
               int *AtReprRes,FILE *FPaOut);
void PsSpMa(double PsSpRa,double GrSiCu_en,int PsSpNC,int ***PsSphe);
void PsSpEE(int FrAtNu,int ReAtNu,double *ReVdWE_sr,double *FrVdWE_sr,
            double *ReVdWR,double *FrVdWR,double *VWEnEv_ps,double **SDFrRe_ps);
/* deprecated */
/* void ReaOut(int CurFra,int SFWrNu,int *Index_ro,double *VW_f_ro,double *VW_s_ro, */
/*             double *In_f_ro,double *In_s_ro,double *Dr_f_ro,double *Dr_s_ro, */
/*             double *Df_f_ro,double *Df_s_ro,double *To_f_ro,double *To_s_ro); */
void Sort(int N,int *IndArr,double *SorArr);
void write_charmm(int CurFra,int SFWrNu,int FrAtNu,int *Index_ro,
                  double *Coulo_ro,double *Vande_ro,double *TotEn_ro,
                  char **FrAtTy,double ***FrCoPo,char **ResN_fr,
                  char **FrFiNa_out);
void MakSpV(int DAType,int DoAcAt,int HyOrZe,double SphAng,int SphPoN,
            int *DANumb,int *ReDATy_L,int *ReDAAt_L,int *ReHydN_L,
            double **ReVeCo_L);
void Simila(int FrNu_1,int FrNu_2,int FrAtNu,int *FrAtEl_nu,double ***FrCoPo,
            double **SimWei,double SimExp,double *FraSim);
void write_chm_clus(int CurFra,int SFWrNu,int FrAtNu,int *AliHyd,
                    double *Coulo_ro,double *Vande_ro,double *TotEn_ro,
                    char **FrAtTy,double ***FrCoPo,int *ClusLi_sd,
                    int *ClusLi_sd_01,int NonAliHy,char **ResN_fr,
                    char **FrFiNa_out);
void MakBSAtList(int ReAtNu,int *ReResN,int BSResN,int *BSReNu,int *BSAtNu,
                 int **BSAtLi);
double TrProd(double ox,double oy,double oz,double ax,double ay,double az,
             double bx,double by,double bz,double cx,double cy,double cz);
void HybStAt(int xxAtNu,int *xxAtEl_nu,double **xxCoor,int xxBdNu,int **xxBdAr,
             int *HybxxAt,FILE *FPaOut);
void FrAcDo(int FrAtNu,int *FrAtEl_nu,double **FrCoor,int *HybFrAt,int FrBdNu,
            int **FrBdAr,int *FrAcNu,int *FrDoNu,int **FrDATy,int **FrDAAt,
            int **FrHydN,double ***FrVeCo,FILE *FPaOut,char *OutFil);
void ReAcDo(int BSAtNu,int *BSAtLi,int *ReAtEl_nu,double **ReCoor,
            int *HybReAt,int ReBdNu,int **ReBdAr,double SphAng,int SphPoN,
            int *ReAcNu,int *ReDoNu,int **ReDATy,int **ReDAAt,int **ReHydN,
            double ***ReVeCo,FILE *FPaOut,char *OutFil,int BSMeNu,int *BSMeAN,
            double **BSMeVE);
void SeedFr_ap(int ReCuVe,double **apol_Vect_re,int *ReApAt,int FrCuVe,
               double **apol_Vect_fr,int *FrApAt,double **FrCoor,double *ReVdWR,
               double *FrVdWR,int FrAtNu,FILE *FPaOut,double **SeFrCo);
void Sort_all(int FragNu,int *SFWrNu_ar,char **ResN_fr,char *Solv_typ);
void SqDisFrRe_ps(int FrAtNu,double **RoSFCo,double **ReCoor,double *ReMinC,
            double GrSiCu_en,int *CubNum_en,int ***CubFAI_en,int ***CubLAI_en,
            int *CubLiA_en,int PsSpNC,int ***PsSphe,double **SDFrRe_ps,
            int ReAtNu,double PsSpRa,double *RePaCh,int ReReNu,int *AtReprRes,
            int *FiAtRes,int *LaAtRes,double *TotChaRes,int NuChResEn,
            int *LiChResEn,double **SDFrRe_ps_elec,double **ChFrRe_ps_elec);
/* void ExtOutNam(int FragNu,char **FrFiNa,char **FrFiNa_out); clangini */
void ExtOutNam(char *FrFiNa,char *FrFiNa_out); /* clangini */
void FindFrSym(int FrAtNu,double **FrCoor,char **FrAtTy,int *UndisAt_fr);
void write_chm_clus_reduc(int CurFra,int SFWrNu,int FrAtNu,int *AliHyd,
                          double *Coulo_ro,double *Vande_ro,double *TotEn_ro,
                          char **FrAtTy,double ***FrCoPo,int *ClusLi_sd,
                          int *ClusLi_sd_01_reduc,int NonAliHy,char **ResN_fr,
                          char **FrFiNa_out);
void Wr_pdb_uhbd(int ReAtNu,double **ReCoor,char *RecFilPDB);
void Reduc_polvectre(int ReAtNu,double **ReCoor,double *ReVdWR,double *ReVdWE_sr,
                     int ReAcNu,int ReDoNu,int *ReDATy,int *ReDAAt,int *ReHydN,
                     double **ReVeCo,double RedRPV_rp,double RedRPV_nkvRatio, /* RedRPV_nkv */
                     int *PolVect_rec_01,FILE *FPaOut,
		     int distrPointBSNumb,double **distrPointBS,double angle_rmin,
		     double angle_rmax,double mult_fact_rmin,double mult_fact_rmax,
		     double Sphere_apol);
void SkipComLin(FILE *FilePa,char *StrLin);
void FeatRes(int ReAtNu,double **ReCoor,int ReReNu,int *ReResN,double *RePaCh,
             int *FiAtRes,int *LaAtRes,int *AtReprRes,double *TotChaRes,
             FILE *FPaOut);
void ListChRes(int BSResN,int *BSReNu,int ReReNu,double **ReCoor,int *AtReprRes,
               double *TotChaRes,double PsSpRa,FILE *FPaOut,int *NuChResEn,
               int *LiChResEn);
void CheckRESN(char **FrFiNa,int FragNu,char **ResN_fr);
void ReFrFi_coo_mol2(int CurFra,int Ind_num_cn,char **FrFiNa,int FrAtNu,
                     double **FrCoor);
void VWGReP(int ReAtNu,double **ReCoor,double *ReVdWR_p3,double *ReVdWE_sr,
            double VWGrIn,double VWGrSi,double *BSMinC,int *VWGPoN,
            double *ReMaxC,double *ReMinC,
            double ***VWGrRP_at,double ***VWGrRP_re);
double VWFrEn(int FrAtNu,double **RoSFCo,double *BSMinC,double VWGrIn,double VWGrSi,
             double *FrVdWR_p3,double *FrVdWE_sr,int *VWGPoN,double ***VWGrRP_at,
             double ***VWGrRP_re,FILE *FPaOut,int * print);
void write_chm_clus_pproc(int CurFra,int NuPosSdCl,int FrAtNu,int *AliHyd,
                          double *In_s_ro,double *VW_s_ro,double *To_s_ro,
                          char **FrAtTy,double ***FrCoPo,int *FrPosAr_sort,
                          int NonAliHy,char **ResN_fr,char **FrFiNa_out);
void write_chm_clus_pprocr(int CurFra,int NuSdClKe,int FrAtNu,char **FrAtTy,
                           double ***FrCoPo,char **ResN_fr,char **FrFiNa_out,
                           int *ReprSdClAr);
void write_mol2_clus_pproc(int CurFra,int NuPosSdCl,int FrAtNu,char **FrAtTy,
                           double ***FrCoPo,int *FrPosAr_sort,char **ResN_fr,
                           char **FrFiNa_out,int FrBdNu,char **FrAtEl,
                           double *FrPaCh,int **FrBdAr,char **FrBdTy,
                           char *FragNa,int *SdClusAr_sort,double *To_s_ro);
void GeomCenter_FFLD(int CurFra,char **FrFiNa_out,int FrAtNu,int *FrAtEl_nu,
                     int NuSdClKe,int *ReprSdClAr,double *To_s_ro,
                     double ***FrCoPo,char **ResN_fr,char *FragNa,
                     int gc_reprke,double gc_cutclus,double gc_endifclus,
                     double gc_weighneg,double gc_weighpos,int gc_maxwrite);



struct point {
  double x;
  double y;
  double z;
};
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
               double ReSurfDens_apol,double Sphere_apol,double NCutapolRatio, /*nCutapol old*/
               double ScaleDeso,double ScaleVDW,
               double ***apol_Vect_re,int **ReApAt,int *Nsurfpt_re_apol_BS,
               int *PNapol_Vect_re,double **SelfVol,
               double *PKelec,double *PKsolv,double *PUnitVol,double pi4,
               double corr_re_desoco,double corr_re_desofd,double corr_fast_deso,
	       int distrPointBSNumb,double **distrPointBS,double angle_rmin,
	       double angle_rmax,double mult_fact_rmin,double mult_fact_rmax,
	       FILE *FPaOut,double **SelfVol_corrB,char *EmpCorrB);
int getColumnFrom2DArray(double **Array2D, int ColumnNum,
                         int start, int end, double *Column);
int SumVectors(double *vect1, double *vect2, int start, int end, double *sum);
int SubstractVectors(double *vect1, double *vect2, int start, int end, double *sum);
double MaxVector(double *Vector, int VectorStart, int VectorEnd);
double MaxDVector(double *Vector, int VectorStart, int VectorEnd);
double MinVector(double *Vector, int VectorStart, int VectorEnd);
struct point c_per_pt(double c, struct point pt);
struct point pt_plus_pt(struct point pt1, struct point pt2);
struct point pt_minus_pt(struct point pt1, struct point pt2);
double pt_scal_pt(struct point pt1, struct point pt2);
void pt_eq_pt(struct point *Ppt1, struct point pt2);
void DSort(int N,int *IndArr,double *SorArr);
int mk_Grid_Re(int ReAtNu,double **ReCoor,double *ReVdWR,double WaMoRa,
                  double GrSiSo,int NGridx,int NGridy,int NGridz,
                  double **XGrid,double **YGrid,double **ZGrid,double UnitVol);
int get_Grid_Dim(int ReAtNu,double **ReCoor,double *ReVdWR,double WaMoRa,
                 double GrSiSo,double GrInSo,int *PNGridx,int *PNGridy,
                 int *PNGridz,int *PNGrid,struct point *PMin,
                 struct point *PMax,double *PUnitVol);
int Coor_Grid_Pts(double GrSiSo,int NGridx,int NGridy,int NGridz,
                  struct point Min,double *XGrid,double *YGrid,double *ZGrid);
int get_Ch_Rad(int ReAtNu,double **ReCoor,double *ReVdWR,double WaMoRa,
               double *ReRad,double *ReRad2,double *ReRadOut,
               double *ReRadOut2,double *PRmin,double *PRmax);
int get_Ch_Rad_Fr(int ReAtNu,double **ReCoor,double *ReVdWR,double WaMoRa,
                  double *ReRad,double *ReRad2,double *ReRadOut,
                  double *ReRadOut2,double **dist2,double *PRmin,double *PRmax);
struct point *SpherePoints(int *PNPtSphere);
int Excl_Grid(int ReAtNu,double **ReCoor,struct point Min,double *ReRadOut,
              double *ReRadOut2,double WaMoRa,double GrSiSo,
              int NStartGridx,int NStartGridy,int NStartGridz,
              int NGridx,int NGridy,int NGridz,
              double UnitVol,char ***GridMat,int Nsurfpt,
              struct point *surfpt,double *SelfVol,double *SelfVol_corrB,
	      char *EmpCorrB);
int Get_Self_Vol(int ReAtNu,double **ReCoor,double *ReRadOut2,
                 double GrSiSo,double *XGrid,double *YGrid,double *ZGrid,
                 int NStartGridx,int NStartGridy,int NStartGridz,
                 int NGridx,int NGridy,int NGridz,
                 double UnitVol,char ***GridMat,double *SelfVol,
		 double *SelfVol_corrB,char *EmpCorrB);
int SAS_Volume(int ReAtNu,double **ReCoor,double *ReRadOut,
               double *ReRadOut2,struct point Min,
               double GrSiSo,int NStartGridx,int NStartGridy,
               int NStartGridz,int NGridx,int NGridy,int NGridz,
               char ***GridMat);
int GB_int(int ReAtNu,double **ReCoor,double *RePaCh,double *EffRad,
           double Ksolv,double *PIntEnTot);
int Read_uhbd_map(char *MapFile,int ibeginx,int ibeginy,int ibeginz,
                  int iendx,int iendy,int iendz,double ***arr);
int Calc_D_uhbd(int ReAtNu,double *ReRad,struct point Min,
                struct point Max,double *RePaCh,double WaMoRa,int NPtSphere,
                double DielRe,double DielWa,double GrSiSo,double GrInSo,
                double pi4,int NGridx,int NGridy,int NGridz,int nxminBS,
                int nyminBS,int nzminBS,int nxmaxBS,int nymaxBS,
                int nzmaxBS,double *XGrid,double *YGrid,double *ZGrid,
                char ***GridMat,double UnitVol,double corr_re_desofd,
                char *RecFilPDB,char *FDexe,char *FDdir,char *DesoMapAcc,
                char *DesoMapFile,double ***DeltaPrDeso);
int Calc_D_Coul(int ReAtNu,double **ReCoor,double *ReRad2,struct point Min,
                struct point Max,double *RePaCh,double DielRe,double DielWa,
                double GrSiSo,double pi4,int nxminBS,int nyminBS,
                int nzminBS,int nxmaxBS,int nymaxBS,int nzmaxBS,
                double *XGrid,double *YGrid,double *ZGrid,char ***GridMat,
                double UnitVol,double corr_re_desoco,double ***DeltaPrDeso);
int Read_Write_Desol_Map(double WaMoRa,int NPtSphere,double DielRe,
                         double DielWa,double GrSiSo,double GrInSo,
                         int NGridx,int NGridy,int NGridz,int nxminBS,
                         int nyminBS,int nzminBS,int nxmaxBS,int nymaxBS,
                         int nzmaxBS,double *XGrid,double *YGrid,
                         double *ZGrid,char ***GridMat,double UnitVol,
                         char *ReDesoAlg,char *DesoMapAcc,
                         char *DesoMapFile,double ***DeltaPrDeso);
int Get_BS_Info(int ReAtNu,double **ReCoor,double *ReRadOut,
                struct point Min,double GrSiSo,double GrInSo,int *ReResN,
                int BSResN,int *BSReNu,int NGridx,
                int NGridy,int NGridz,int *PnxminBS,int *PnyminBS,
                int *PnzminBS,int *PnxmaxBS,int *PnymaxBS,
                int *PnzmaxBS,int *ResFirstAtom,int *NAtom_per_Res);
int Fast_Desol_Surf(struct point Min,double WaMoRa,double GrSiSo,
                   int NStartGridx,int NStartGridy,int NStartGridz,
                   int NGridx,int NGridy,int NGridz,int Nsurfpt,
                   struct point *surfpt,double ***DeltaPrDeso,
                   double corr_fast_deso,double *ReSurf_deso);
int Desol_Surf(int ReAtNu,struct point Min,double Sphere_apol,
               double GrSiSo,int NGridx,int NGridy,int NGridz,int BSResN,
               int *BSReNu,int *NAtom_per_Res,int *ResFirstAtom,
               struct point *surfpt_re,int *nsurf_re,int *pointsrf_re,
               double ***DeltaPrDeso,double *ReSurf_apol,double NCutapol,/*int NCutapol,*/
               double ScaleDeso,int *isurf_More_apol,
               int *PNapol_Vect_re,int *Nsurfpt_re_apol_BS,
	       int *iatsrf_re_apol,int distrPointBSNumb,double **distrPointBS,
	       double angle_rmin,double angle_rmax,double mult_fact_rmin,
	       double mult_fact_rmax,double **ReCoor,FILE *FPaOut);
int vdW_Surf(int ReAtNu,double **ReCoor,double *ReVdWE_sr,double *ReVdWR,
             double Sphere_apol,int BSResN,int *BSReNu,int *NAtom_per_Res,
             int *ResFirstAtom,struct point *surfpt_re,int *nsurf_re,
             int *pointsrf_re,double *vdWSurfPt);
int Int_Rec_q(struct point Min,struct point Max,double GrSiSo,double GrInSo,
              double *XGrid,double *YGrid,double *ZGrid,
              int NGridx,int NGridy,int NGridz,char ***GridMat,
              double ***Dx,double ***Dy,double ***Dz,double Kcoul,int nxminBS,
              int nyminBS,int nzminBS,int nxmaxBS,int nymaxBS,int nzmaxBS,
              double ***UnitIntBS);
int FragSolvEn(int FrAtNu,double **FrCoor,double *FrPaCh,
               double *FrVdWR,double *FrRadOut,
               double *FrRadOut2,double **Frdist2,int Nsurfpt_fr,
               struct point *surfpt_fr_orig,double WaMoRa,double GrSiSo,
               double Ksolv,double pi4,double *PFrSolvEn,char *EmpCorrB,FILE*FPaOut);
int FragDesoSurf(int FrAtNu,double **FrCoor,double *FrPaCh,
                 double *FrVdWR,double *FrRad2,double *FrRadOut,
                 double *FrRadOut2,int Nsurfpt_fr,
                 struct point *surfpt_fr_orig,double WaMoRa,double GrSiSo,
                 double DielRe,double DielWa,double pi4,
                 double corr_fast_deso,double corr_re_desoco,
                 double **FrSurf_deso);
void ElecFrag(int ReAtNu,double **ReCoor,double *RePaCh,
              double **RePaCh_Fr,double *ReRad,
              double *ReRad2,double *ReRadOut,double *ReRadOut2,
              struct point *surfpt_re,int *nsurf_re,
              int *pointsrf_re,double *ReSelfVol,int FrAtNu,double **RoSFCo,
              double **FrCoor,double *FrPaCh,double *FrRad,double *FrRad2,
              double *FrRadOut,double *FrRadOut2,double **Frdist2,double **dist2,
              double *FrMinC_orig,double *FrMaxC_orig,
              double *PFrSolvEn,int Nsurfpt_fr,struct point *surfpt_fr_orig,
              int *nsurf_fr,int *pointsrf_fr,struct point *surfpt_ex,
              double Tr[4],double U1[4][4],double U2[4][4],double WaMoRa,
              double GrSiSo,int NPtSphere,struct point Min,
              struct point Max,double *XGrid,double *YGrid,double *ZGrid,
              int NGridx,int NGridy,int NGridz,char ***GridMat,
              double ***DeltaPrDeso,
              double Kelec,double Ksolv,double UnitVol,double pi4,
              int nxminBS,int nyminBS,int nzminBS,int nxmaxBS,
              int nymaxBS,int nzmaxBS,double corr_scrint,
              double corr_fr_deso,double *PReDesoElec,
              double *PReFrIntElec,double *PFrDesoElec,double *ReSelfVol_corrB,
	      char *EmpCorrB,FILE *FPaOut);
void ElecFragFast(int ReAtNu,double **ReCoor,double *ReRadOut,double *ReRadOut2,
                 struct point *surfpt_re,int *nsurf_re,double *ReSurf_deso,
                 int *pointsrf_re,int ***nlist_re,int **list_re,int nstep_re,
                 struct point cbmid_re,int midg_re,double rscale_re,
                 int FrAtNu,double **RoSFCo,double **FrCoor,double *FrRad,
                 double *FrRadOut,double *FrRadOut2,double *FrMinC_orig,
                 double *FrMaxC_orig,double FrRmax,int Nsurfpt_fr,
                 struct point *surfpt_fr_orig,struct point *surfpt_fr,
                 int *nsurf_fr,double *FrSurf_deso,int *pointsrf_fr,
                 double Tr[4],double U1[4][4],double U2[4][4],
                 double WaMoRa,double rslop,int MaxNeigh,double *PReDesoElec,
		 double *PFrDesoElec);
int get_frag_dim(int FrAtNu,double **RoSFCo,double *FrRadOut,double GrSiSo,
                 double WaMoRa,int NGridx,int NGridy,int NGridz,
                 struct point Min,int *PnxminFr,
                 int *PnyminFr,int *PnzminFr,int *PnxmaxFr,
                 int *PnymaxFr,int *PnzmaxFr,int *PFrOut);
int get_Rec_neigh(int ReAtNu,int FrAtNu,double **dist2,double **dist,
                  double *ReRad,double *FrRad,double *RePaCh,double WaMoRa,
                  int *nsurf_re,int *PNNeigh1,int *NeighList1,int *PNNeigh2,
                  int *NeighList2,int *PNNeigh3,int *NeighList3);
int Mov_surf(int Nsurfpt_fr,struct point *surfpt_fr_orig,double **FrCoor,
             double *FrMinC_orig,double *FrMaxC_orig,double Tr[4],double U1[4][4],
             double U2[4][4],struct point *surfpt_fr,double *FrMinC,
             double *FrMaxC);
int join_surf(double **ReCoor,double *ReRadOut2,struct point *surfpt_re,
              int *nsurf_re,int *pointsrf_re,int FrAtNu,double **RoSFCo,
              double *FrRadOut2,struct point *surfpt_fr,int *nsurf_fr,
              int *pointsrf_fr,int NNeigh,int *NeighList,int *PNsurfpt_ex,
              struct point *surfpt_ex);
int SAS_Volume_Fr(int ExFrAtNu,double **ExRoSFCo,double *ExFrRadOut,
                  double *ExFrRadOut2,struct point Min,
                  double GrSiSo,int NStartGridx,int NStartGridy,
                  int NStartGridz,int NGridx,int NGridy,int NGridz,
                  char ***FrGridMat,char ***GridMat,int FrOut);
int Excl_Grid_Fr(struct point Min,double WaMoRa,double GrSiSo,
                 int nxminFr,int nyminFr,int nzminFr,
                 int nxmaxFr,int nymaxFr,int nzmaxFr,
                 char ***FrGridMat,int Nsurfpt,struct point *surfpt);
int Get_Self_Vol_Re(int ReAtNu,double **ReCoor,int NFrNeigh,int *FrNeighList,
                    double *ReRad2,
                    double GrSiSo,double *XGrid,double *YGrid,double *ZGrid,
                    int NStartGridx,int NStartGridy,int NStartGridz,
                    int NGridx,int NGridy,int NGridz,
                    double UnitVol,char ***GridMat,double *SelfVol,
		    double *SelfVol_corrB,char *EmpCorrB);
int Get_Self_Vol_Re2(int FrAtNu,double **RoSFCo,double *FrPaCh,int NNeigh3,
                     int *NeighList3,
                    double *FrRad2,double *FrRadOut,double *FrRadOut2,
                    double GrSiSo,double *XGrid,double *YGrid,double *ZGrid,
                    int NGridx,int NGridy,int NGridz,struct point Min,
                    double UnitVol,char ***GridMat,char ***FrGridMat,
                    int nxminFr,int nyminFr,int nzminFr,
                    int nxmaxFr,int nymaxFr,int nzmaxFr,double *SelfVol);
int Get_Self_Vol_Fr(int FrAtNu,double **RoSFCo,double *FrPaCh,
                    double *FrRad2,double *FrRadOut,double *FrRadOut2,
                    double GrSiSo,double *XGrid,double *YGrid,double *ZGrid,
                    struct point Min,double UnitVol,int NGridx,int NGridy,
                    int NGridz,char ***GridMat,int nxminFr,int nyminFr,
                    int nzminFr,int nxmaxFr,int nymaxFr,int nzmaxFr,
                    char ***FrGridMat,int FrOut,int nxmin_big,int nymin_big,
                    int nzmin_big,int nxmax_big,int nymax_big,int nzmax_big,
                    double *SelfVol,double *SelfVol_corrB,char *EmpCorrB);
int Get_Self_En_Fr(int FrAtNu,double **RoSFCo,double *FrPaCh,double *FrRadOut,
                   double *FrRadOut2,double *XGrid,double *YGrid,double *ZGrid,
                   double UnitVol,double Ksolv,double pi4,char ***FrGridMat,
                   int nxminFr,int nyminFr,int nzminFr,
                   int nxmaxFr,int nymaxFr,int nzmaxFr,
                   double *FrSelfVol,double *FrEffRad,double *PFrSelfEn,
		   double *FrSelfVol_corrB,char *EmpCorrB,FILE* FPaOut);
int screened_int(double **RePaCh_Fr,double *ReEffRad,
                 int NNeigh3,int *NeighList3,int FrAtNu,
                 double *FrPaCh,double *FrEffRad,double **dist2,double **dist,
                 double Kelec,double Ksolv,double *PReFrIntElec);
int screened_int_uhbd(int FrAtNu,double **RoSFCo,double *FrPaCh,double ***phi,
                      double *XGrid,double *YGrid,double *ZGrid,
                      struct point Min,double GrSiSo,double *PReFrIntElec);
int GB_int_fr(int FrAtNu,double **Frdist2,double *FrPaCh,
              double *EffRad,double Ksolv,double *PIntEnTot);
int GB_int_re(int FrAtNu,double **ReCoor,double *FrPaCh,
              double *EffRad,double Ksolv,double *PIntEnTot);
int Surf_Grid_unif_dens(int ReAtNu,double **ReCoor,double *ReMaxC,double *ReMinC,
                        double *ReRad,double *ReRadOut2,double WaMoRa,
                        double ReSurfDens_apol,double pi4,int *PNsurfpt,
                        struct point *surfpt,int *iatsrf,int *nsurf,
                        int *pointsrf);
int Surf_Grid_unif_dens_neighlist(int ReAtNu,double **ReCoor,double *ReMinC,
                                  double *ReRad,double *ReRadOut2,
                                  double WaMoRa,double ReSurfDens_deso,
                                  double pi4,int *PNsurfpt,
                                  struct point *surfpt,int *iatsrf,int *nsurf,
                                  int *pointsrf,struct point len,double slen,
                                  int ***nlist,int **list,
                                  int nstep,struct point *Pcbmid,int *Pmidg,
                                  double rscale);
int Surf_Grid(int ReAtNu,double **ReCoor,double *ReMaxC,double *ReMinC,
               double *ReRad,double WaMoRa,int NPtSphere,int *PNsurfpt,
               struct point *surfpt,int *iatsrf,int *nsurf,int *pointsrf);
int Rot_Tran(int FrAtNu,double **FrCoor,double **RoSFCo,double Tr[4],
                  double U1[4][4],double U2[4][4]);
int Fr_neighlist(int FrAtNu,double **RoSFCo,double *FrMinC,double *FrRad,
                 double WaMoRa,struct point len,double slen,
                 int ***nlist,int **list,int nstep,struct point *Pcbmid,
                 int *Pmidg,double rscale);
int Fast_Desol(int AtNu,double **Coor,double *RadOut,double *RadOut2,
               struct point *surfpt,int *nsurf,int *pointsrf,int ***nlist,
               int **list,int nstep,struct point cbmid,int midg,double rscale,
               double *Surf_deso,double *PDesoElec);
int Fast_Desol_Receptor(int AtNu,double **Coor,double *RadOut,double *RadOut2,
			struct point *surfpt,int *nsurf,int *pointsrf,int ***nlist,
			int **list,int nstep,struct point cbmid,int midg,double rscale,
			double *Surf_deso,double *PDesoElec,
			double * minC,double *maxC, double ** xcoor,double maxrad,double WaMoRa);
int CheckFile(char*name,char type);
int CheckMol2File(char*name);
