#define sqrtf sqrt
#define sinf sin
#define cosf cos
#define fabsf fabs
#define ffloor floor
#define acosf acos
#define expf exp
#define modff modf
#define logf log

void ReInFi(char *InpFil,char *RecFil,int *BSResN,int **BSReNu,
            int *FragNu,char ***FrFiNa,char *TREFiP,float *SphAng,
            int *SphPoN,int *NuRoAx,float *VdWFaB,float *CoDieV,int *CoDieP,
            float *CoGrIn,float *CoGrSi,char *OutFil,float *BuEvFa,
            float **FrMaEn,float *PsSpRa,float *GrSiCu_en,
            int *FiNuMa,double *GrInSo,double *GrSiSo,double *WaMoRa,
            int *NPtSphere,double *DielWa,double *DielRe,char *ReDesoAlg,
            char *DesoMapAcc,char *DesoMapFile,char *FDexe,char *FDdir,
            double *ReSurfDens_apol,double *PtDensFr,double *Sphere_apol,
            float *NCutapolRatio/*int *NCutapol*/,double *ScaleDeso,
            double *ScaleVDW,float **SimWei,float *SimExp,
            float *SimCut,float **FrMaEn_sd,float *SimExp_sd,float *SimCut_sd,
            int *BSMeNu,int **BSMeAN,float ***BSMeVE,char *CoGrAcc,
            char *CoGrFile,char *EvalEn,char *Solv_typ,
            char *SpPoCh_opt,float *SpPoCh_cent,float *SpPoCh_rad,
            float *SFDeso_fr,float *SFDeso_re,float *SFVWEn,float *SFIntElec,
            int *NuClusMem,float *RedRPV_rp,float *RedRPV_nkvDens,float *ScMaBump, /*int *RedRPV_nkv, dey new*/
            float *MuFaVdWCoff_ap,int *NuLiEnClus,char ***ApPoChoi,
            float *VWGrIn,float *VWGrSi,float *BumpFaCut,char *VWGrAcc,
            char *VWGrFile,int *MaxPosClus,int *PrintLev,
            int *NumbAT,char ***AtTyAr,int **AtENAr,float **VdWRad,
            float **VdWEne,float ***BLAtTy,int *distrPointBSNumb,
	    float ***distrPointBS,float *angle_rmin,float *angle_rmax,
	    float *mult_fact_rmin,float *mult_fact_rmax,char *EmpCorrB,
            char *gc_opt,int *gc_reprke,float *gc_cutclus,float *gc_endifclus,
            float *gc_weighneg,float *gc_weighpos,int *gc_maxwrite,
            char *write_pproc_opt,char *write_pproc_chm_opt,int *CorrFiNumb);
void ReReFi_mol2(char *RecFil,int *ReAtNu,int *ReBdNu,int *ReReNu,
                 char ***ReAtEl,float ***ReCoor,char ***ReAtTy,int **ReResN,
                 float **RePaCh,int ***ReBdAr);
void ReFrFi_mol2(int CurFra,char **FrFiNa,char *FragNa,int *FrAtNu,int *FrBdNu,
                 char ***FrAtEl,float ***FrCoor,char ***FrAtTy,float **FrPaCh,
                 int ***FrBdAr,char ***FrBdTy,char *FrSubN,char *FrSubC,
                 int *FrCoNu);
/* C.LANGINI 2016 */
void ReFrFi_mol2_syb(int CurFra,char **FrFiNa,char *FragNa,int *FrAtNu,int *FrBdNu,
                 char ***FrAtEl,float ***FrCoor,char ***FrAtTy, char ***FrSyAtTy,float **FrPaCh,
                 int ***FrBdAr,char ***FrBdTy,char *FrSubN,char *FrSubC,
                 int *FrCoNu);
void write_mol2_clus_pproc_syb(int CurFra,int NuPosSdCl,int FrAtNu,char **FrAtTy, char **FrSyAtTy,
                           float ***FrCoPo,int *FrPosAr_sort,char **ResN_fr,
                           char **FrFiNa_out,int FrBdNu,char **FrAtEl,
                           float *FrPaCh,int **FrBdAr,char **FrBdTy,
                           char *FragNa,int *SdClusAr_sort,float *To_s_ro);
void append_pose_to_mol2(FILE *FilePa,char *FragNa,int FrAtNu,int FrBdNu,int imol,
                         char **FrAtEl,float ***FrCoPo,int *FrPosAr_sort,char **FrSyAtTy,
                         int CurFra,char **ResN_fr,int **FrBdAr,char **FrBdTy,
                         int *SdClusAr_sort,float *To_s_ro,float *FrPaCh);
void write_mol2_clus_pproc_separate(int CurFra,int NuPosSdCl,int FrAtNu,char **FrAtTy,
                                    char **FrSyAtTy,float ***FrCoPo,int *FrPosAr_sort,
                                    char **ResN_fr,char **FrFiNa_out,int FrBdNu,char **FrAtEl,
                                    float *FrPaCh,int **FrBdAr,char **FrBdTy,
                                    char *FragNa,int *SdClusAr_sort,float *To_s_ro);
void ReFrFi_mol2_cpp(std::istream *inStream, std::streampos *strPos,
		 int CurFra,char **FrFiNa,char *FragNa,int *FrAtNu,int *FrBdNu,
                 char ***FrAtEl,float ***FrCoor,char ***FrAtTy,  char ***FrSyAtTy,float **FrPaCh,
                 int ***FrBdAr,char ***FrBdTy,char *FrSubN,char *FrSubC,
                 int *FrCoNu)
/* C.LANGINI 2016 END */
void FiAlHy(int FrAtNu,int FrBdNu,int *FrAtEl_nu,int **FrBdAr,int *AliHyd,
            int *NonAliHy);
float VeNorm(float a1,float a2,float a3);
void NormVe(float *a1,float *a2,float *a3);
void PoCoVe(float a1,float a2,float a3,float b1,float b2,float b3,float len,
            float *c1,float *c2,float *c3);
void VectPr(float a1,float a2,float a3,float b1,float b2,float b3,float *c1,
            float *c2,float *c3);
float PlaAng(float o1,float o2,float o3,float a1,float a2,float a3,float b1,
             float b2,float b3);
void RoArVe(float a1,float a2,float a3,float vec1,float vec2,float vec3,
            float ang,float *b1,float *b2,float *b3);
void AxisVe(float a1,float a2,float a3,float b1,float b2,float b3,
            float c1,float c2,float c3,float *d1,float *d2,float *d3);
void RotPla(float a1,float a2,float a3,float b1,float b2,float b3,float c1,
            float c2,float c3,float angl,float *d1,float *d2,float *d3);
float DihAng(float a1,float a2,float a3,float b1,float b2,float b3,
             float c1,float c2,float c3,float d1,float d2,float d3);
float DistSq(float a1,float a2,float a3,float b1,float b2,float b3);
void SeedFr(int ReCuVe,float **ReVeCo,int FrCuVe,float **FrVeCo,int FrAtNu,
            float **FrCoor,int *ReDATy,float **SeFrCo,int *FrAtoTyp_nu,
            int *FrDAAt,int *ReAtoTyp_nu,int *ReDAAt,float **BLAtTy,
            FILE *FPaOut);
void write_mol2(char *WriPat,char *WriNam,int xxAtNu,int xxBdNu,char **xxAtEl,
                float **xxCoor,char **xxAtTy,float *xxPaCh,int **xxBdAr,
                char **xxBdTy,char *FrSubN,char *FrSubC,char **ResN_fr,
                int CurFra);
void RoSeFr(int ReCuVe,int *RexxAt,float **ReCoor,int FrCuVe,int *FrxxAt,
            int FrAtNu,float **SeFrCo,float AnglRo,float **RoSFCo);
int AssiRE(int NumbAT,char **AtTyAr,float *VdWRad,float *VdWEne,int xxAtNu,
	   char **xxAtTy,float *xxVdWR,float *xxVdWE,FILE *FPaOut,
	   int *AtENAr,int *xxAtEl_nu,int *xxAtoTyp_nu);
void VdWRaM(int NumbAT,float *VdWRad,float *LaVdWR);
void ReMMCo(int ReAtNu,float **ReCoor,float *ReMaxC,float *ReMinC);
void ReLAIC(int ReAtNu,float **ReCoor,float LaVdWR,float *ReMaxC,float *ReMinC,
            int *CubNum,int ****CubFAI,int ****CubLAI,int **CubLiA);
int CoNuBu(int *CubNum,int ***CubFAI,int ***CubLAI,int *CubLiA,float LaVdWR,
           float *ReMaxC,float *ReMinC,float **ReCoor,float *ReVdWR, int FrAtNu,
           float **RoSFCo,float *FrVdWR,float VdWFaB,int BumpMa,int RecHyN,
           int FraHyN,float BuEvFa);
void BSMMCo(float **ReCoor,int BSAtNu,int *BSAtLi,float *BSMaxC,float *BSMinC);
void CoGReP(int ReAtNu,float **ReCoor,float *RePaCh,float CoDieV,
            int CoDieP,float CoGrIn,float CoGrSi,float *BSMinC,int *CoGPoN,
            float ***CoGrRP);
float CoFrEn(int FrAtNu,float *FrPaCh,float **RoSFCo,float *BSMinC,
             float CoGrIn,float CoGrSi,int *CoGPoN,float ***CoGrRP,
             FILE *FPaOut,int * print);
void SqRoEn(int xxAtNu,float *xxVdWE,float *xxVdWE_sr);
void ReLAIC_en(int ReReNu,float **ReCoor,float GrSiCu_en,float *ReMinC,
               int *CubNum_en,int ***CubFAI_en,int ***CubLAI_en,int *CubLiA_en,
               int *AtReprRes,FILE *FPaOut);
void PsSpMa(float PsSpRa,float GrSiCu_en,int PsSpNC,int ***PsSphe);
void PsSpEE(int FrAtNu,int ReAtNu,float *ReVdWE_sr,float *FrVdWE_sr,
            float *ReVdWR,float *FrVdWR,float *VWEnEv_ps,float **SDFrRe_ps);
/* deprecated */
/* void ReaOut(int CurFra,int SFWrNu,int *Index_ro,float *VW_f_ro,float *VW_s_ro, */
/*             float *In_f_ro,float *In_s_ro,float *Dr_f_ro,float *Dr_s_ro, */
/*             float *Df_f_ro,float *Df_s_ro,float *To_f_ro,float *To_s_ro); */
void Sort(int N,int *IndArr,float *SorArr);
void write_charmm(int CurFra,int SFWrNu,int FrAtNu,int *Index_ro,
                  float *Coulo_ro,float *Vande_ro,float *TotEn_ro,
                  char **FrAtTy,float ***FrCoPo,char **ResN_fr,
                  char **FrFiNa_out);
void MakSpV(int DAType,int DoAcAt,int HyOrZe,float SphAng,int SphPoN,
            int *DANumb,int *ReDATy_L,int *ReDAAt_L,int *ReHydN_L,
            float **ReVeCo_L);
void Simila(int FrNu_1,int FrNu_2,int FrAtNu,int *FrAtEl_nu,float ***FrCoPo,
            float **SimWei,float SimExp,float *FraSim);
void write_chm_clus(int CurFra,int SFWrNu,int FrAtNu,int *AliHyd,
                    float *Coulo_ro,float *Vande_ro,float *TotEn_ro,
                    char **FrAtTy,float ***FrCoPo,int *ClusLi_sd,
                    int *ClusLi_sd_01,int NonAliHy,char **ResN_fr,
                    char **FrFiNa_out);
void MakBSAtList(int ReAtNu,int *ReResN,int BSResN,int *BSReNu,int *BSAtNu,
                 int **BSAtLi);
float TrProd(float ox,float oy,float oz,float ax,float ay,float az,
             float bx,float by,float bz,float cx,float cy,float cz);
void HybStAt(int xxAtNu,int *xxAtEl_nu,float **xxCoor,int xxBdNu,int **xxBdAr,
             int *HybxxAt,FILE *FPaOut);
void FrAcDo(int FrAtNu,int *FrAtEl_nu,float **FrCoor,int *HybFrAt,int FrBdNu,
            int **FrBdAr,int *FrAcNu,int *FrDoNu,int **FrDATy,int **FrDAAt,
            int **FrHydN,float ***FrVeCo,FILE *FPaOut,char *OutFil);
void ReAcDo(int BSAtNu,int *BSAtLi,int *ReAtEl_nu,float **ReCoor,
            int *HybReAt,int ReBdNu,int **ReBdAr,float SphAng,int SphPoN,
            int *ReAcNu,int *ReDoNu,int **ReDATy,int **ReDAAt,int **ReHydN,
            float ***ReVeCo,FILE *FPaOut,char *OutFil,int BSMeNu,int *BSMeAN,
            float **BSMeVE);
void SeedFr_ap(int ReCuVe,float **apol_Vect_re,int *ReApAt,int FrCuVe,
               float **apol_Vect_fr,int *FrApAt,float **FrCoor,float *ReVdWR,
               float *FrVdWR,int FrAtNu,FILE *FPaOut,float **SeFrCo);
void Sort_all(int FragNu,int *SFWrNu_ar,char **ResN_fr,char *Solv_typ);
void SqDisFrRe_ps(int FrAtNu,float **RoSFCo,float **ReCoor,float *ReMinC,
            float GrSiCu_en,int *CubNum_en,int ***CubFAI_en,int ***CubLAI_en,
            int *CubLiA_en,int PsSpNC,int ***PsSphe,float **SDFrRe_ps,
            int ReAtNu,float PsSpRa,float *RePaCh,int ReReNu,int *AtReprRes,
            int *FiAtRes,int *LaAtRes,float *TotChaRes,int NuChResEn,
            int *LiChResEn,float **SDFrRe_ps_elec,float **ChFrRe_ps_elec);
void ExtOutNam(int FragNu,char **FrFiNa,char **FrFiNa_out);
void FindFrSym(int FrAtNu,float **FrCoor,char **FrAtTy,int *UndisAt_fr);
void write_chm_clus_reduc(int CurFra,int SFWrNu,int FrAtNu,int *AliHyd,
                          float *Coulo_ro,float *Vande_ro,float *TotEn_ro,
                          char **FrAtTy,float ***FrCoPo,int *ClusLi_sd,
                          int *ClusLi_sd_01_reduc,int NonAliHy,char **ResN_fr,
                          char **FrFiNa_out);
void Wr_pdb_uhbd(int ReAtNu,float **ReCoor,char *RecFilPDB);
void Reduc_polvectre(int ReAtNu,float **ReCoor,float *ReVdWR,float *ReVdWE_sr,
                     int ReAcNu,int ReDoNu,int *ReDATy,int *ReDAAt,int *ReHydN,
                     float **ReVeCo,float RedRPV_rp,float RedRPV_nkvRatio, /* RedRPV_nkv */
                     int *PolVect_rec_01,FILE *FPaOut,
		     int distrPointBSNumb,float **distrPointBS,float angle_rmin,
		     float angle_rmax,float mult_fact_rmin,float mult_fact_rmax,
		     double Sphere_apol);
void SkipComLin(FILE *FilePa,char *StrLin);
void FeatRes(int ReAtNu,float **ReCoor,int ReReNu,int *ReResN,float *RePaCh,
             int *FiAtRes,int *LaAtRes,int *AtReprRes,float *TotChaRes,
             FILE *FPaOut);
void ListChRes(int BSResN,int *BSReNu,int ReReNu,float **ReCoor,int *AtReprRes,
               float *TotChaRes,float PsSpRa,FILE *FPaOut,int *NuChResEn,
               int *LiChResEn);
void CheckRESN(char **FrFiNa,int FragNu,char **ResN_fr);
void ReFrFi_coo_mol2(int CurFra,int Ind_num_cn,char **FrFiNa,int FrAtNu,
                     float **FrCoor);
void VWGReP(int ReAtNu,float **ReCoor,float *ReVdWR_p3,float *ReVdWE_sr,
            float VWGrIn,float VWGrSi,float *BSMinC,int *VWGPoN,
            float *ReMaxC,float *ReMinC,
            float ***VWGrRP_at,float ***VWGrRP_re);
float VWFrEn(int FrAtNu,float **RoSFCo,float *BSMinC,float VWGrIn,float VWGrSi,
             float *FrVdWR_p3,float *FrVdWE_sr,int *VWGPoN,float ***VWGrRP_at,
             float ***VWGrRP_re,FILE *FPaOut,int * print);
void write_chm_clus_pproc(int CurFra,int NuPosSdCl,int FrAtNu,int *AliHyd,
                          float *In_s_ro,float *VW_s_ro,float *To_s_ro,
                          char **FrAtTy,float ***FrCoPo,int *FrPosAr_sort,
                          int NonAliHy,char **ResN_fr,char **FrFiNa_out);
void write_chm_clus_pprocr(int CurFra,int NuSdClKe,int FrAtNu,char **FrAtTy,
                           float ***FrCoPo,char **ResN_fr,char **FrFiNa_out,
                           int *ReprSdClAr);
void write_mol2_clus_pproc(int CurFra,int NuPosSdCl,int FrAtNu,char **FrAtTy,
                           float ***FrCoPo,int *FrPosAr_sort,char **ResN_fr,
                           char **FrFiNa_out,int FrBdNu,char **FrAtEl,
                           float *FrPaCh,int **FrBdAr,char **FrBdTy,
                           char *FragNa,int *SdClusAr_sort,float *To_s_ro);
void GeomCenter_FFLD(int CurFra,char **FrFiNa_out,int FrAtNu,int *FrAtEl_nu,
                     int NuSdClKe,int *ReprSdClAr,float *To_s_ro,
                     float ***FrCoPo,char **ResN_fr,char *FragNa,
                     int gc_reprke,float gc_cutclus,float gc_endifclus,
                     float gc_weighneg,float gc_weighpos,int gc_maxwrite);



struct point {
  double x;
  double y;
  double z;
};
void Solvation(int ReAtNu,float **ReCoor,float *ReVdWE_sr,float *ReVdWR,
               double *ReRad,double *ReRad2,double *ReRadOut,double *ReRadOut2,
               float *ReMaxC,float *ReMinC,float *RePaCh,double DielRe,
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
               double ReSurfDens_apol,double Sphere_apol,float NCutapolRatio, /*nCutapol old*/
               double ScaleDeso,double ScaleVDW,
               float ***apol_Vect_re,int **ReApAt,int *Nsurfpt_re_apol_BS,
               int *PNapol_Vect_re,double **SelfVol,
               double *PKelec,double *PKsolv,double *PUnitVol,double pi4,
               double corr_re_desoco,double corr_re_desofd,double corr_fast_deso,
	       int distrPointBSNumb,float **distrPointBS,float angle_rmin,
	       float angle_rmax,float mult_fact_rmin,float mult_fact_rmax,
	       FILE *FPaOut,double **SelfVol_corrB,char *EmpCorrB);
int getColumnFrom2DArray(float **Array2D, int ColumnNum,
                         int start, int end, double *Column);
int SumVectors(double *vect1, float *vect2, int start, int end, double *sum);
int SubstractVectors(double *vect1, float *vect2, int start, int end, double *sum);
float MaxVector(float *Vector, int VectorStart, int VectorEnd);
double MaxDVector(double *Vector, int VectorStart, int VectorEnd);
double MinVector(double *Vector, int VectorStart, int VectorEnd);
struct point c_per_pt(double c, struct point pt);
struct point pt_plus_pt(struct point pt1, struct point pt2);
struct point pt_minus_pt(struct point pt1, struct point pt2);
double pt_scal_pt(struct point pt1, struct point pt2);
void pt_eq_pt(struct point *Ppt1, struct point pt2);
void DSort(int N,int *IndArr,double *SorArr);
int mk_Grid_Re(int ReAtNu,float **ReCoor,float *ReVdWR,double WaMoRa,
                  double GrSiSo,int NGridx,int NGridy,int NGridz,
                  double **XGrid,double **YGrid,double **ZGrid,double UnitVol);
int get_Grid_Dim(int ReAtNu,float **ReCoor,float *ReVdWR,double WaMoRa,
                 double GrSiSo,double GrInSo,int *PNGridx,int *PNGridy,
                 int *PNGridz,int *PNGrid,struct point *PMin,
                 struct point *PMax,double *PUnitVol);
int Coor_Grid_Pts(double GrSiSo,int NGridx,int NGridy,int NGridz,
                  struct point Min,double *XGrid,double *YGrid,double *ZGrid);
int get_Ch_Rad(int ReAtNu,float **ReCoor,float *ReVdWR,double WaMoRa,
               double *ReRad,double *ReRad2,double *ReRadOut,
               double *ReRadOut2,double *PRmin,double *PRmax);
int get_Ch_Rad_Fr(int ReAtNu,float **ReCoor,float *ReVdWR,double WaMoRa,
                  double *ReRad,double *ReRad2,double *ReRadOut,
                  double *ReRadOut2,float **dist2,double *PRmin,double *PRmax);
struct point *SpherePoints(int *PNPtSphere);
int Excl_Grid(int ReAtNu,float **ReCoor,struct point Min,double *ReRadOut,
              double *ReRadOut2,double WaMoRa,double GrSiSo,
              int NStartGridx,int NStartGridy,int NStartGridz,
              int NGridx,int NGridy,int NGridz,
              double UnitVol,char ***GridMat,int Nsurfpt,
              struct point *surfpt,double *SelfVol,double *SelfVol_corrB,
	      char *EmpCorrB);
int Get_Self_Vol(int ReAtNu,float **ReCoor,double *ReRadOut2,
                 double GrSiSo,double *XGrid,double *YGrid,double *ZGrid,
                 int NStartGridx,int NStartGridy,int NStartGridz,
                 int NGridx,int NGridy,int NGridz,
                 double UnitVol,char ***GridMat,double *SelfVol,
		 double *SelfVol_corrB,char *EmpCorrB);
int SAS_Volume(int ReAtNu,float **ReCoor,double *ReRadOut,
               double *ReRadOut2,struct point Min,
               double GrSiSo,int NStartGridx,int NStartGridy,
               int NStartGridz,int NGridx,int NGridy,int NGridz,
               char ***GridMat);
int GB_int(int ReAtNu,float **ReCoor,float *RePaCh,double *EffRad,
           double Ksolv,double *PIntEnTot);
int Read_uhbd_map(char *MapFile,int ibeginx,int ibeginy,int ibeginz,
                  int iendx,int iendy,int iendz,float ***arr);
int Calc_D_uhbd(int ReAtNu,double *ReRad,struct point Min,
                struct point Max,float *RePaCh,double WaMoRa,int NPtSphere,
                double DielRe,double DielWa,double GrSiSo,double GrInSo,
                double pi4,int NGridx,int NGridy,int NGridz,int nxminBS,
                int nyminBS,int nzminBS,int nxmaxBS,int nymaxBS,
                int nzmaxBS,double *XGrid,double *YGrid,double *ZGrid,
                char ***GridMat,double UnitVol,double corr_re_desofd,
                char *RecFilPDB,char *FDexe,char *FDdir,char *DesoMapAcc,
                char *DesoMapFile,double ***DeltaPrDeso);
int Calc_D_Coul(int ReAtNu,float **ReCoor,double *ReRad2,struct point Min,
                struct point Max,float *RePaCh,double DielRe,double DielWa,
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
int Get_BS_Info(int ReAtNu,float **ReCoor,double *ReRadOut,
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
               double ***DeltaPrDeso,double *ReSurf_apol,float NCutapol,/*int NCutapol,*/
               double ScaleDeso,int *isurf_More_apol,
               int *PNapol_Vect_re,int *Nsurfpt_re_apol_BS,
	       int *iatsrf_re_apol,int distrPointBSNumb,float **distrPointBS,
	       float angle_rmin,float angle_rmax,float mult_fact_rmin,
	       float mult_fact_rmax,float **ReCoor,FILE *FPaOut);
int vdW_Surf(int ReAtNu,float **ReCoor,float *ReVdWE_sr,float *ReVdWR,
             double Sphere_apol,int BSResN,int *BSReNu,int *NAtom_per_Res,
             int *ResFirstAtom,struct point *surfpt_re,int *nsurf_re,
             int *pointsrf_re,double *vdWSurfPt);
int Int_Rec_q(struct point Min,struct point Max,double GrSiSo,double GrInSo,
              double *XGrid,double *YGrid,double *ZGrid,
              int NGridx,int NGridy,int NGridz,char ***GridMat,
              double ***Dx,double ***Dy,double ***Dz,double Kcoul,int nxminBS,
              int nyminBS,int nzminBS,int nxmaxBS,int nymaxBS,int nzmaxBS,
              double ***UnitIntBS);
int FragSolvEn(int FrAtNu,float **FrCoor,float *FrPaCh,
               float *FrVdWR,double *FrRadOut,
               double *FrRadOut2,float **Frdist2,int Nsurfpt_fr,
               struct point *surfpt_fr_orig,double WaMoRa,double GrSiSo,
               double Ksolv,double pi4,double *PFrSolvEn,char *EmpCorrB,FILE*FPaOut);
int FragDesoSurf(int FrAtNu,float **FrCoor,float *FrPaCh,
                 float *FrVdWR,double *FrRad2,double *FrRadOut,
                 double *FrRadOut2,int Nsurfpt_fr,
                 struct point *surfpt_fr_orig,double WaMoRa,double GrSiSo,
                 double DielRe,double DielWa,double pi4,
                 double corr_fast_deso,double corr_re_desoco,
                 double **FrSurf_deso);
void ElecFrag(int ReAtNu,float **ReCoor,float *RePaCh,
              float **RePaCh_Fr,double *ReRad,
              double *ReRad2,double *ReRadOut,double *ReRadOut2,
              struct point *surfpt_re,int *nsurf_re,
              int *pointsrf_re,double *ReSelfVol,int FrAtNu,float **RoSFCo,
              float **FrCoor,float *FrPaCh,double *FrRad,double *FrRad2,
              double *FrRadOut,double *FrRadOut2,float **Frdist2,float **dist2,
              float *FrMinC_orig,float *FrMaxC_orig,
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
void ElecFragFast(int ReAtNu,float **ReCoor,double *ReRadOut,double *ReRadOut2,
                 struct point *surfpt_re,int *nsurf_re,double *ReSurf_deso,
                 int *pointsrf_re,int ***nlist_re,int **list_re,int nstep_re,
                 struct point cbmid_re,int midg_re,double rscale_re,
                 int FrAtNu,float **RoSFCo,float **FrCoor,double *FrRad,
                 double *FrRadOut,double *FrRadOut2,float *FrMinC_orig,
                 float *FrMaxC_orig,double FrRmax,int Nsurfpt_fr,
                 struct point *surfpt_fr_orig,struct point *surfpt_fr,
                 int *nsurf_fr,double *FrSurf_deso,int *pointsrf_fr,
                 double Tr[4],double U1[4][4],double U2[4][4],
                 double WaMoRa,double rslop,int MaxNeigh,double *PReDesoElec,
		 double *PFrDesoElec);
int get_frag_dim(int FrAtNu,float **RoSFCo,double *FrRadOut,double GrSiSo,
                 double WaMoRa,int NGridx,int NGridy,int NGridz,
                 struct point Min,int *PnxminFr,
                 int *PnyminFr,int *PnzminFr,int *PnxmaxFr,
                 int *PnymaxFr,int *PnzmaxFr,int *PFrOut);
int get_Rec_neigh(int ReAtNu,int FrAtNu,float **dist2,float **dist,
                  double *ReRad,double *FrRad,float *RePaCh,double WaMoRa,
                  int *nsurf_re,int *PNNeigh1,int *NeighList1,int *PNNeigh2,
                  int *NeighList2,int *PNNeigh3,int *NeighList3);
int Mov_surf(int Nsurfpt_fr,struct point *surfpt_fr_orig,float **FrCoor,
             float *FrMinC_orig,float *FrMaxC_orig,double Tr[4],double U1[4][4],
             double U2[4][4],struct point *surfpt_fr,float *FrMinC,
             float *FrMaxC);
int join_surf(float **ReCoor,double *ReRadOut2,struct point *surfpt_re,
              int *nsurf_re,int *pointsrf_re,int FrAtNu,float **RoSFCo,
              double *FrRadOut2,struct point *surfpt_fr,int *nsurf_fr,
              int *pointsrf_fr,int NNeigh,int *NeighList,int *PNsurfpt_ex,
              struct point *surfpt_ex);
int SAS_Volume_Fr(int ExFrAtNu,float **ExRoSFCo,double *ExFrRadOut,
                  double *ExFrRadOut2,struct point Min,
                  double GrSiSo,int NStartGridx,int NStartGridy,
                  int NStartGridz,int NGridx,int NGridy,int NGridz,
                  char ***FrGridMat,char ***GridMat,int FrOut);
int Excl_Grid_Fr(struct point Min,double WaMoRa,double GrSiSo,
                 int nxminFr,int nyminFr,int nzminFr,
                 int nxmaxFr,int nymaxFr,int nzmaxFr,
                 char ***FrGridMat,int Nsurfpt,struct point *surfpt);
int Get_Self_Vol_Re(int ReAtNu,float **ReCoor,int NFrNeigh,int *FrNeighList,
                    double *ReRad2,
                    double GrSiSo,double *XGrid,double *YGrid,double *ZGrid,
                    int NStartGridx,int NStartGridy,int NStartGridz,
                    int NGridx,int NGridy,int NGridz,
                    double UnitVol,char ***GridMat,double *SelfVol,
		    double *SelfVol_corrB,char *EmpCorrB);
int Get_Self_Vol_Re2(int FrAtNu,float **RoSFCo,float *FrPaCh,int NNeigh3,
                     int *NeighList3,
                    double *FrRad2,double *FrRadOut,double *FrRadOut2,
                    double GrSiSo,double *XGrid,double *YGrid,double *ZGrid,
                    int NGridx,int NGridy,int NGridz,struct point Min,
                    double UnitVol,char ***GridMat,char ***FrGridMat,
                    int nxminFr,int nyminFr,int nzminFr,
                    int nxmaxFr,int nymaxFr,int nzmaxFr,double *SelfVol);
int Get_Self_Vol_Fr(int FrAtNu,float **RoSFCo,float *FrPaCh,
                    double *FrRad2,double *FrRadOut,double *FrRadOut2,
                    double GrSiSo,double *XGrid,double *YGrid,double *ZGrid,
                    struct point Min,double UnitVol,int NGridx,int NGridy,
                    int NGridz,char ***GridMat,int nxminFr,int nyminFr,
                    int nzminFr,int nxmaxFr,int nymaxFr,int nzmaxFr,
                    char ***FrGridMat,int FrOut,int nxmin_big,int nymin_big,
                    int nzmin_big,int nxmax_big,int nymax_big,int nzmax_big,
                    double *SelfVol,double *SelfVol_corrB,char *EmpCorrB);
int Get_Self_En_Fr(int FrAtNu,float **RoSFCo,float *FrPaCh,double *FrRadOut,
                   double *FrRadOut2,double *XGrid,double *YGrid,double *ZGrid,
                   double UnitVol,double Ksolv,double pi4,char ***FrGridMat,
                   int nxminFr,int nyminFr,int nzminFr,
                   int nxmaxFr,int nymaxFr,int nzmaxFr,
                   double *FrSelfVol,double *FrEffRad,double *PFrSelfEn,
		   double *FrSelfVol_corrB,char *EmpCorrB,FILE* FPaOut);
int screened_int(float **RePaCh_Fr,double *ReEffRad,
                 int NNeigh3,int *NeighList3,int FrAtNu,
                 float *FrPaCh,double *FrEffRad,float **dist2,float **dist,
                 double Kelec,double Ksolv,double *PReFrIntElec);
int screened_int_uhbd(int FrAtNu,float **RoSFCo,float *FrPaCh,float ***phi,
                      double *XGrid,double *YGrid,double *ZGrid,
                      struct point Min,double GrSiSo,double *PReFrIntElec);
int GB_int_fr(int FrAtNu,float **Frdist2,float *FrPaCh,
              double *EffRad,double Ksolv,double *PIntEnTot);
int GB_int_re(int FrAtNu,float **ReCoor,float *FrPaCh,
              double *EffRad,double Ksolv,double *PIntEnTot);
int Surf_Grid_unif_dens(int ReAtNu,float **ReCoor,float *ReMaxC,float *ReMinC,
                        double *ReRad,double *ReRadOut2,double WaMoRa,
                        double ReSurfDens_apol,double pi4,int *PNsurfpt,
                        struct point *surfpt,int *iatsrf,int *nsurf,
                        int *pointsrf);
int Surf_Grid_unif_dens_neighlist(int ReAtNu,float **ReCoor,float *ReMinC,
                                  double *ReRad,double *ReRadOut2,
                                  double WaMoRa,double ReSurfDens_deso,
                                  double pi4,int *PNsurfpt,
                                  struct point *surfpt,int *iatsrf,int *nsurf,
                                  int *pointsrf,struct point len,double slen,
                                  int ***nlist,int **list,
                                  int nstep,struct point *Pcbmid,int *Pmidg,
                                  double rscale);
int Surf_Grid(int ReAtNu,float **ReCoor,float *ReMaxC,float *ReMinC,
               double *ReRad,double WaMoRa,int NPtSphere,int *PNsurfpt,
               struct point *surfpt,int *iatsrf,int *nsurf,int *pointsrf);
int Rot_Tran(int FrAtNu,float **FrCoor,float **RoSFCo,double Tr[4],
                  double U1[4][4],double U2[4][4]);
int Fr_neighlist(int FrAtNu,float **RoSFCo,float *FrMinC,double *FrRad,
                 double WaMoRa,struct point len,double slen,
                 int ***nlist,int **list,int nstep,struct point *Pcbmid,
                 int *Pmidg,double rscale);
int Fast_Desol(int AtNu,float **Coor,double *RadOut,double *RadOut2,
               struct point *surfpt,int *nsurf,int *pointsrf,int ***nlist,
               int **list,int nstep,struct point cbmid,int midg,double rscale,
               double *Surf_deso,double *PDesoElec);
int Fast_Desol_Receptor(int AtNu,float **Coor,double *RadOut,double *RadOut2,
			struct point *surfpt,int *nsurf,int *pointsrf,int ***nlist,
			int **list,int nstep,struct point cbmid,int midg,double rscale,
			double *Surf_deso,double *PDesoElec,
			float * minC,float *maxC, float ** xcoor,float maxrad,float WaMoRa);
int CheckFile(char*name,char type);
int CheckMol2File(char*name);
