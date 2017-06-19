#include <stdio.h>

void Wr_pdb_uhbd(int ReAtNu,double **ReCoor,char *RecFilPDB)
/* This function writes a kind of PDB format to be used by UHBD :
   RecFilPDB  name of the receptor file (PDB format for UHBD) */
{
  FILE *FilePa_1;
  int i;

/* Create RecFilPDB */
  sprintf(RecFilPDB,"%s%s","./outputs/","receptor_uhbd.pdb\0");
  FilePa_1=fopen(RecFilPDB,"w");

  for (i=1;i<=ReAtNu;i++) {
    fprintf(FilePa_1,"ATOM %6d X XXX 1 %12.3f%12.3f%12.3f\n",
            i,ReCoor[i][1],ReCoor[i][2],ReCoor[i][3]);
  }

  fprintf(FilePa_1,"END\n");

  fclose(FilePa_1);

}
