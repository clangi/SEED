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
