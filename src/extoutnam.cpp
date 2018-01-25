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
#include <string.h>
#include <iostream>
/* void ExtOutNam(int FragNu,char **FrFiNa,char **FrFiNa_out) clangini*/
/* This function extracts the output name for each fragment type from the
   fragments paths :
   FrFiNa_out  output name for each fragment type
   LengStr  length of the fragment file path (FrFiNa)
   FirCha  pointer on the first character of the name
   LasCha  pointer on the last character of the name */
/*{
  char *FirCha,*LasCha;
  int i,LengStr,j;

  for (i=1;i<=FragNu;i++) {

    LengStr=strlen(&FrFiNa[i][1]);

    FirCha=strrchr(&FrFiNa[i][1],'/');
    if (FirCha)
      FirCha=FirCha+1;
    else
      FirCha=&FrFiNa[i][1];

    LasCha=strrchr(&FrFiNa[i][1],'.');
    if (LasCha) {
      if (LasCha>FirCha)
        LasCha=LasCha-1;
      else
        LasCha=&FrFiNa[i][1]+LengStr-1;
    }
    else
      LasCha=&FrFiNa[i][1]+LengStr-1;

    for (j=1;j<=(LasCha-FirCha+1);j++)
      FrFiNa_out[i][j]=FrFiNa[i][FirCha-&FrFiNa[i][1]+j];
    FrFiNa_out[i][LasCha-FirCha+2]='\0';

  }

} clangini*/

/* clangini 2016 New version of ExtOutNam: now we have only one input file */
void ExtOutNam(char *FrFiNa,char *FrFiNa_out)
/* This function extracts the output name for each fragment type from the
   fragments paths :
   FrFiNa_out  output name for each fragment type
   LengStr  length of the fragment file path (FrFiNa)
   FirCha  pointer on the first character of the name
   LasCha  pointer on the last character of the name */
{
  char *FirCha,*LasCha;
  int /*i,*/LengStr,j;

  LengStr=strlen(&FrFiNa[0]);
  //std::cout << "LengStr: " << LengStr << std::endl;

  FirCha=strrchr(&FrFiNa[0],'/');
  if (FirCha)
    FirCha=FirCha+1;
  else
    FirCha=&FrFiNa[0];

  //std::cout << "FirCha: " << FirCha << std::endl;

  LasCha=strrchr(&FrFiNa[0],'.');
  if (LasCha) {
    if (LasCha>FirCha)
      LasCha=LasCha-1;
    else
      LasCha=&FrFiNa[0]+LengStr-1;
  }
  else
    LasCha=&FrFiNa[0]+LengStr-1;

  //std::cout << "LasCha-FirCha+1: " << (LasCha-FirCha+1) << std::endl;

   for (j=0;j<=(LasCha-FirCha);j++)
    FrFiNa_out[j]=FrFiNa[FirCha-&FrFiNa[0]+j];
  FrFiNa_out[LasCha-FirCha+2]='\0';

  //std::cout << "FrFiNa_out: " << FrFiNa_out <<"aaa"<< std::endl;

}
/* clangini 2016 end */
