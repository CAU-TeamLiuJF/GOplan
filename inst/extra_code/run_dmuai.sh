#!/bin/bash
# Date: 2021-04-08
# Author: liwn
# Description: Run dmu_ai in parallel

if [  $# -eq 0 ]
  then
    echo "Usage: /full/path/of/run_dmuai your_parameter_card(dmuai).DIR"
	exit 1
  else
    name=$1
fi

if [ -f "$1.DIR" ]; then
  echo Starting DMU with $1.DIR as directive file
  ./dmu1 < $1.DIR > $1.lst
  
  if [ -f MODINF ]; then
    ./dmuai >> $1.lst
        if [ $? -eq 0 ]
        then
          printf "program %-27s OK \n" $1   >> ../run_ex.log
        else
          echo
          echo "${1} DMUai failed - Check output files"
          echo "${1} DMUai failed - Check output files" >> ../run_ex.log
        fi
  else
    echo
    echo "${1} DMUai failed - Check output files" >> ../run_ex.log
    echo DMUAI Not started due to errors in DMU1
  fi
  
  [ -f SOL ] && mv SOL $1.SOL
  [ -f PAROUT ] && mv PAROUT $1.PAROUT
  [ -f PAROUT_STD ] && mv PAROUT_STD $1.PAROUT_STD
  [ -f LLIK ] && mv LLIK $1.LLIK
  [ -f RESIDUAL ] && mv RESIDUAL $1.RESIDUAL
  
  [ -f CODE_TABLE ] && rm CODE_TABLE
  [ -f DMU1.dir ] && rm DMU1.dir
  [ -f DMUAI.dir ] && rm DMUAI.dir
  [ -f DMU_LOG ] && rm DMU_LOG
  [ -f DUMMY ] && rm DUMMY
  [ -f FSPAKWK ] && rm FSPAKWK
  [ -f Latest_parm ] && rm Latest_parm
  [ -f LEVAL ] && rm LEVAL
  [ -f MODINF ] && rm MODINF
  [ -f PARIN ] && rm PARIN
  [ -f RCDATA_I ] && rm RCDATA_I
  [ -f RCDATA_R ] && rm RCDATA_R
  [ -f PARIN ] && rm PARIN
  
  [ -f AINV* ] && rm AINV*
  [ -f COR* ] && rm COR*
  [ -f DOM* ] && rm DOM*
  [ -f IBD* ] && rm IBD*
  [ -f PEDFILE* ] && rm PEDFILE*
  rm -f fort.* 1
    
  if [ -f INBREED ]; then
    nchar=`cat INBREED | wc -c `
    if [ $nchar == 0 ]; then
      rm INBREED
    else
      mv INBREED $1.INBREED
    fi
  fi
  
else
  echo
  echo File $1.DIR not in current directory
fi
