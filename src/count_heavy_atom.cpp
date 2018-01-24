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
#include "nrutil.h"

/* This function calculates the HAC (Heavy Atom Count) clangini*/
int count_heavy_atom(int *FrAtEl_nu,int FrAtNu)
{
  int HeAtCo, i;

  HeAtCo = 0;
  for (i = 1; i <= FrAtNu; i++){
    if (FrAtEl_nu[i] > 1 )
      HeAtCo++;
  }
  return HeAtCo;
}

/* This function compute the Molecular Weight (MW) clangini*/
double molecular_weight(int *FrAtEl_nu,int FrAtNu,double *AtWei)
{
  int i;
  double MolWei = 0.0;
  for(i = 1; i <= FrAtNu; i++){
    MolWei += AtWei[FrAtEl_nu[i]];
  }
  return MolWei;
}
