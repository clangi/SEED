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
#include <stdlib.h>
#include "funct.h"
#include "nrutil.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#ifndef _STRLENGTH
#define _STRLENGTH 500
#endif
/* The code for the mol2 reader was partially inspired by RDKIT http://www.rdkit.org/ */
// /* Older version, rewritten (see below) by clangini 2017 to read ALT_TYPE_SET */
// void ReReFi_mol2(char *RecFil,int *ReAtNu,int *ReBdNu,int *ReReNu,
//                  char ***ReAtEl,double ***ReCoor,char ***ReAtTy,int **ReResN,
//                  double **RePaCh,int ***ReBdAr)
// /* This function reads the data of the receptor file (RecFil) in the mol2
//    format :
//    ReAtNu  number of receptor atoms
//    ReBdNu  number of receptor bonds
//    ReReNu  number of receptor residues
//    ReAtEl  receptor atoms elements
//    ReCoor  receptor atoms coordinates
//    ReAtTy  receptor atoms types
//    ReResN  receptor residues numbers
//    RePaCh  receptor partial charges
//    ReBdAr  receptor bonds array */
// {
//   FILE *FilePa;
//   char StrLin[_STRLENGTH],**ReAtEl_L,**ReAtTy_L,dummy2[10];
//   double **ReCoor_L,*RePaCh_L;
//   int i,dummy1,*ReResN_L;
//
//   FilePa=fopen(RecFil,"r");
//
// /* Read ReAtNu ReBdNu ReReNu */
//   fgets_wrapper(StrLin,_STRLENGTH,FilePa);
//   for (i=0;strncmp(StrLin,"@<TRIPOS>MOLECULE",17);i++)
//     fgets_wrapper(StrLin,_STRLENGTH,FilePa);
//   fgets_wrapper(StrLin,_STRLENGTH,FilePa);
//   fgets_wrapper(StrLin,_STRLENGTH,FilePa);
//   sscanf(StrLin,"%d%d%d",ReAtNu,ReBdNu,ReReNu);
//
// /* Read ReAtEl ReCoor ReAtTy ReResN RePaCh */
//   ReAtEl_L=cmatrix(1,*ReAtNu,1,5);
//   ReCoor_L=dmatrix(1,*ReAtNu,1,3);
//   ReAtTy_L=cmatrix(1,*ReAtNu,1,7);
//   ReResN_L=ivector(1,*ReAtNu);
//   RePaCh_L=dvector(1,*ReAtNu);
//   *ReAtEl=ReAtEl_L;
//   *ReCoor=ReCoor_L;
//   *ReAtTy=ReAtTy_L;
//   *ReResN=ReResN_L;
//   *RePaCh=RePaCh_L;
//   fgets_wrapper(StrLin,_STRLENGTH,FilePa);
//   for (i=0;strncmp(StrLin,"@<TRIPOS>ATOM",13);i++)
//     fgets_wrapper(StrLin,_STRLENGTH,FilePa);
//   for (i=1;i<=*ReAtNu;i++) {
//     fgets_wrapper(StrLin,_STRLENGTH,FilePa);
//     sscanf(StrLin,"%d%s%lf%lf%lf%s%d%s%lf",&dummy1,&ReAtEl_L[i][1],
//            &ReCoor_L[i][1],&ReCoor_L[i][2],&ReCoor_L[i][3],&ReAtTy_L[i][1],
//            &ReResN_L[i],dummy2,&RePaCh_L[i]);
//   }
//
// /* Read ReBdAr */
//   *ReBdAr=imatrix(1,*ReBdNu,1,2);
//   fgets_wrapper(StrLin,_STRLENGTH,FilePa);
//   for (i=0;strncmp(StrLin,"@<TRIPOS>BOND",13);i++)
//     fgets_wrapper(StrLin,_STRLENGTH,FilePa);
//   for (i=1;i<=*ReBdNu;i++) {
//     fgets_wrapper(StrLin,_STRLENGTH,FilePa);
//     sscanf(StrLin,"%d%d%d%s",&dummy1,&((*ReBdAr)[i][1]),&((*ReBdAr)[i][2]),
//            dummy2);
//   }
//
//   /* clangini 2017
//   add check for the correct number of residues */
//   int ReReNu_count = 0;
//   int ReRe_idx;
//   fgets_wrapper(StrLin,_STRLENGTH,FilePa);
//   for (i=0;strncmp(StrLin,"@<TRIPOS>SUBSTRUCTURE",13);i++)
//     fgets_wrapper(StrLin,_STRLENGTH,FilePa);
//   while ((fgets(StrLin,_STRLENGTH,FilePa )!= NULL) &&
//                               !feof(FilePa) && !ferror(FilePa) ){
//       sscanf(StrLin,"%d",&ReRe_idx);
//       if (ReRe_idx == (ReReNu_count+1))
//         ++ReReNu_count;
//       else {
//         break;
//       }
//   }
//   std::cout << ReReNu_count << "  " << *ReReNu << std::endl;
//   if (ReReNu_count != *ReReNu){
//     std::cerr << "WARNING, number of residues in the receptor mol2 file (" <<
//     ReReNu_count << "), does not match the number declared in the " <<
//     "@<TRIPOS>MOLECULE record (" << *ReReNu << ")" << std::endl;
//     std::cerr << "Program exits" << std::endl;
//     exit(12);
//   }
//   /* clangini 2017 end of check for correct number of residues */
//
//   fclose(FilePa);
// }
/* End of older version,
  rewritten (see below) by clangini 2017 to read ALT_TYPE_SET */

using namespace std;
/* New version: rewritten by clangini 2017 to read ALT_TYPE_SET */
void ReReFi_mol2(char *RecFil,int *ReAtNu,int *ReBdNu,int *ReReNu,
                 char ***ReAtEl,double ***ReCoor,char ***ReAtTy,int **ReResN,
                 double **RePaCh,int ***ReBdAr)
/* This function reads the data of the receptor file (RecFil) in the mol2
   format :
   ReAtNu  number of receptor atoms
   ReBdNu  number of receptor bonds
   ReReNu  number of receptor residues
   ReAtEl  receptor atoms elements
   ReCoor  receptor atoms coordinates
   ReAtTy  receptor atoms types
   ReResN  receptor residues numbers
   RePaCh  receptor partial charges
   ReBdAr  receptor bonds array */
{
  typedef boost::tokenizer< boost::char_separator<char> > tokenizer;
  boost::char_separator<char> sep(" \t\n");

  std::ifstream inStream;
  std::string StrLin, AlTySp, firstToken, ResNa_tmp, ResNa;
  std::size_t found;
  bool AtNu_flag;
  int CuAtNu, AtCount, ReReNu_prev;
  char **ReAtEl_L,**ReAtTy_L;//dummy2[10];
  int ReReNu_count;
  //char dummy_AtTy[10];
  double **ReCoor_L,*RePaCh_L;
  int i,dummy1,*ReResN_L;

  // Open file stream
  inStream.open(RecFil, std::ios_base::binary); /* open input stream for the fragment */
  if (inStream.fail()) {
    fprintf(stderr,"Cannot open receptor input file %s\nProgram exits!\n",RecFil);
    exit(13);
  }
  // Initialize tokenizer:
  std::getline(inStream, StrLin);
  boost::trim(StrLin);
	tokenizer tokens(StrLin,sep); //Initialize tokenizer
	tokenizer::const_iterator itItem = tokens.begin();

  // Move to molecule section
  while (!inStream.eof() && StrLin != "@<TRIPOS>MOLECULE"){
    std::getline(inStream, StrLin);
    boost::trim(StrLin);
  }
  if (inStream.eof()){
    std::cerr << "No Receptor found.\nProgram exits!\n" << std::endl;
    exit(13);
  }
  std::getline(inStream,StrLin);
  std::getline(inStream,StrLin);
  std::stringstream(StrLin) >> (*ReAtNu) >> (*ReBdNu) >> dummy1;

  /* Read ReAtEl ReCoor ReAtTy ReResN RePaCh */
  ReAtEl_L=cmatrix(1,*ReAtNu,1,5);
  ReCoor_L=dmatrix(1,*ReAtNu,1,3);
  ReAtTy_L=cmatrix(1,*ReAtNu,1,7); //Now read from ALT_TYPE
  ReResN_L=ivector(1,*ReAtNu);
  RePaCh_L=dvector(1,*ReAtNu);
  *ReAtEl=ReAtEl_L;
  *ReCoor=ReCoor_L;
  *ReAtTy=ReAtTy_L;
  *ReResN=ReResN_L;
  *RePaCh=RePaCh_L;
  // Move to atom section:
  while (!inStream.eof() && StrLin != "@<TRIPOS>ATOM"){
    std::getline(inStream,StrLin);
    boost::trim(StrLin);
  }
  if (inStream.eof()){
    std::cerr << "No Atoms in receptor found.\nProgram exits!\n" << std::endl;
    exit(13);
  }
  /* Read atom section */
  ReReNu_count = 0; //Initialize residue number counter
  ReReNu_prev = -1;
  for (i=1; i<=(*ReAtNu); i++ ){
    std::getline(inStream,StrLin);
    //std::cout << "atom line "<<i<<": "<< StrLin << std::endl;
    tokens.assign(StrLin, sep);
    itItem = tokens.begin();
    ++itItem;
    //std::cout << "atom element "<<i<<": "<< *itItem << std::endl;
    //FrAtEl_L[i][1] = (*itItem).c_str();
    strcpy(&ReAtEl_L[i][1],(*itItem).c_str());
    ++itItem;
    //std::cout << "X-coordinate "<<i<<": "<< *itItem << std::endl;
    ReCoor_L[i][1] = boost::lexical_cast<double> (*itItem);
    //std::cout << "X-coordinate "<<i<<": "<<FrCoor_L[i][1]<< std::endl;
    ++itItem;
    //std::cout << "Y-coordinate "<<i<<": "<< *itItem << std::endl;
    ReCoor_L[i][2] = boost::lexical_cast<double> (*itItem);
    //std::cout << "Y-coordinate "<<i<<": "<<FrCoor_L[i][2]<< std::endl;
    ++itItem;
    //std::cout << "Z-coordinate "<<i<<": "<< *itItem << std::endl;
    ReCoor_L[i][3] = boost::lexical_cast<double> (*itItem);
    //std::cout << "Z-coordinate "<<i<<": "<<FrCoor_L[i][3]<< std::endl;
    ++itItem; ++itItem;
    ReResN_L[i] = boost::lexical_cast<int> (*itItem);
    ++itItem;
    ResNa_tmp = *itItem;
    if ((ResNa_tmp != ResNa) || (ReResN_L[i] > ReReNu_prev)){ // change res if res name or res num change
      ++ReReNu_count;
      ReReNu_prev = ReResN_L[i];
      ResNa = ResNa_tmp;
    }
    ++itItem;
    //std::cout << "Partial charge "<<i<<": "<< (*itItem)<<"ebbasta" << std::endl;
    RePaCh_L[i] = boost::lexical_cast<double> (*itItem);
    //std::cout << "Charge "<<i<<": "<<FrPaCh_L[i]<< std::endl;
  }
  *ReReNu = ReReNu_count;
  cout << "\n\tNumber of Residues in receptor: " << (*ReReNu) << endl;

  while (!inStream.eof() && StrLin != "@<TRIPOS>BOND"){
    std::getline(inStream,StrLin);
    boost::trim(StrLin);
  }
  if (inStream.eof()){
    std::cerr << "No Bonds in receptor found.\nProgram exits!\n" << std::endl;
    exit(13);
  }
  // Read bond section
  *ReBdAr=imatrix(1,*ReBdNu,1,2);
  for (i = 1; i <= *ReBdNu; i++) {
    std::getline(inStream,StrLin);
    tokens.assign(StrLin,sep);
    itItem = tokens.begin();
    ++itItem; // skip bond id
    (*ReBdAr)[i][1] = boost::lexical_cast<int>(*itItem);
    ++itItem;
    (*ReBdAr)[i][2] = boost::lexical_cast<int>(*itItem);
  }
  //cout << "Receptor bonds read" << endl;

  // Move to ALT_TYPE section
  while (!inStream.eof() && StrLin != "@<TRIPOS>ALT_TYPE"){
    std::getline(inStream,StrLin);
    boost::trim(StrLin);
  }
  if (inStream.eof()){
    std::cerr << "No Atom Types in receptor found.\nProgram exits!\n"
              << std::endl;
    exit(13);
  }

  std::getline(inStream,StrLin);
  boost::trim(StrLin);
  found = StrLin.find("ALT_TYPE_SET");
  if (found == std::string::npos){
    std::cerr << "No standard ALT_TYPE_SET signature find for receptor"
                  << "\nProgram exits"<< endl;
    exit(13);
  }
  AlTySp = StrLin.substr(0,(found-1)); // Save the ALT_TYPE_SET name (for example CHARMM)
  //std::cout <<"Alternative atom type specification is: "<<AlTySp<<std::endl;
  std::getline(inStream,StrLin);
  tokens.assign(StrLin, sep);
  itItem = tokens.begin();
  firstToken = *(itItem);
  if (firstToken != AlTySp){
    std::cerr << "Names of alternative atom type set do not coincide"
              << "for receptor\nProgram exits!" << endl;
    exit(13);
  }
  ++itItem; // skip the alternative atom type set name
  AtCount = 0;
  AtNu_flag = false;
  while (AtCount < *ReAtNu){
    if (*itItem != "\\"){
      if (!AtNu_flag){
        CuAtNu =  boost::lexical_cast<int>(*itItem); // Current atom number
        AtNu_flag = true;
        ++itItem;
      } else {
        //FrAtTy_L[CuAtNu][1] = (*itItem).c_str();
        //std::cout<<"Atom type for atom "<<AtCount+1<<": "<<(*itItem)<<std::endl;
        strcpy(&ReAtTy_L[CuAtNu][1],(*itItem).c_str());
        //std::cout<<"Atom type for atom "<<CuAtNu<<": "<<(*FrAtTy)[CuAtNu][1]<<std::endl;
        //std::cout<<"Atom type for atom "<<AtCount+1<<": "<<FrAtTy_L[CuAtNu][1]<<std::endl;
        ++AtCount;
        ++itItem;
        AtNu_flag = false;
      }
    } else {
      //StrLin = getline(inStream);
      std::getline(inStream,StrLin);
      //tokenizer tokens(StrLin, sep);
      tokens.assign(StrLin, sep);
      //tokenizer::const_iterator itItem = tokens.begin();
      itItem = tokens.begin();
    }
  }
  //cout << "Receptor Atom Types read" << endl;
  inStream.close();

// /* Read ReAtNu ReBdNu ReReNu */
//   fgets_wrapper(StrLin,_STRLENGTH,FilePa);
//   for (i=0;strncmp(StrLin,"@<TRIPOS>MOLECULE",17);i++)
//     fgets_wrapper(StrLin,_STRLENGTH,FilePa);
//   fgets_wrapper(StrLin,_STRLENGTH,FilePa);
//   fgets_wrapper(StrLin,_STRLENGTH,FilePa);
//   //sscanf(StrLin,"%d%d%d",ReAtNu,ReBdNu,ReReNu);
//   sscanf(StrLin,"%d%d%d",ReAtNu,ReBdNu,&dummy1); //No longer read ReReNu
//                                                  //from here
//
// /* Read ReAtEl ReCoor ReAtTy ReResN RePaCh */
//   ReAtEl_L=cmatrix(1,*ReAtNu,1,5);
//   ReCoor_L=dmatrix(1,*ReAtNu,1,3);
//   //ReAtTy_L=cmatrix(1,*ReAtNu,1,7); Now read from ALT_TYPE
//   ReResN_L=ivector(1,*ReAtNu);
//   RePaCh_L=dvector(1,*ReAtNu);
//   *ReAtEl=ReAtEl_L;
//   *ReCoor=ReCoor_L;
//   //*ReAtTy=ReAtTy_L;
//   *ReResN=ReResN_L;
//   *RePaCh=RePaCh_L;
//   fgets_wrapper(StrLin,_STRLENGTH,FilePa);
//   for (i=0;strncmp(StrLin,"@<TRIPOS>ATOM",13);i++)
//     fgets_wrapper(StrLin,_STRLENGTH,FilePa);
//
//   ReReNu_count = 0; //Initialize residue number counter
//   for (i=1;i<=*ReAtNu;i++) {
//     fgets_wrapper(StrLin,_STRLENGTH,FilePa);
//     sscanf(StrLin,"%d%s%lf%lf%lf%s%d%s%lf",&dummy1,&ReAtEl_L[i][1],
//            &ReCoor_L[i][1],&ReCoor_L[i][2],&ReCoor_L[i][3],dummy_AtTy,
//            &ReResN_L[i],ResNa_tmp,&RePaCh_L[i]);
//     if (strncmp(ResNa_tmp,ResNa) != 0){
//       ReReNu_count++;
//       strcpy(ResNa,ResNa_tmp);
//     }
//   }
//   *ReReNu = ReReNu_count;
//   cout << "Number of Residues: " << (*ReReNu) << endl;
//
// /* Read ReBdAr */
//   *ReBdAr=imatrix(1,*ReBdNu,1,2);
//   fgets_wrapper(StrLin,_STRLENGTH,FilePa);
//   for (i=0;strncmp(StrLin,"@<TRIPOS>BOND",13);i++)
//     fgets_wrapper(StrLin,_STRLENGTH,FilePa);
//   for (i=1;i<=*ReBdNu;i++) {
//     fgets_wrapper(StrLin,_STRLENGTH,FilePa);
//     sscanf(StrLin,"%d%d%d%s",&dummy1,&((*ReBdAr)[i][1]),&((*ReBdAr)[i][2]),
//            dummy2);
//   }

  /* clangini 2017
  add check for the correct number of residues */
  /*int ReReNu_count = 0;
  int ReRe_idx;
  fgets_wrapper(StrLin,_STRLENGTH,FilePa);
  for (i=0;strncmp(StrLin,"@<TRIPOS>SUBSTRUCTURE",13);i++)
    fgets_wrapper(StrLin,_STRLENGTH,FilePa);
  while ((fgets(StrLin,_STRLENGTH,FilePa )!= NULL) &&
                              !feof(FilePa) && !ferror(FilePa) ){
      sscanf(StrLin,"%d",&ReRe_idx);
      if (ReRe_idx == (ReReNu_count+1))
        ++ReReNu_count;
      else {
        break;
      }
  }
  std::cout << ReReNu_count << "  " << *ReReNu << std::endl;
  if (ReReNu_count != *ReReNu){
    std::cerr << "WARNING, number of residues in the receptor mol2 file (" <<
    ReReNu_count << "), does not match the number declared in the " <<
    "@<TRIPOS>MOLECULE record (" << *ReReNu << ")" << std::endl;
    std::cerr << "Program exits" << std::endl;
    exit(12);
  }*/
  /* clangini 2017 end of check for correct number of residues */

  //fclose(FilePa);
}
