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
