#!/bin/bash
# clangini 26/07/2017
# script for splitting a concatenated mol2 file library.mol2 into N files
# Usage: bash split_library.sh library.mol2 N

library_fn_full=$1
npart=$2

# library_path=${library%/*}
# library_fn_full=${library##*/}
library_fn_noext=${library_fn_full%.*}
extension=${library_fn_full##*.}

# echo $library_path
echo $library_fn_full
echo $library_fn_noext
echo $extension


nmol=$(grep "@<TRIPOS>MOLECULE" $library_fn_full | wc -l )
nmol_part=$((($nmol + $npart - 1)/$npart))
echo $nmol
echo $nmol_part

awk -v fname=${library_fn_noext} -v pat="@<TRIPOS>MOLECULE" -v ext=$extension -v maxmol=$nmol_part '
BEGIN {
  i = 0
  filenum = 1
}
{
  if ($0 ~ pat){
    i++
  }
  if (i <= maxmol){
    print > fname"_part"filenum"."ext
  } else {
    filenum++
    print > fname"_part"filenum"."ext
    i = 1
  }

}
' ${library_fn_full}
