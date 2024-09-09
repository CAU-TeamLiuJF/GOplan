#!/bin/bash

# macOS script based on original bat script

dmu1="./dmu1"   # Assuming these are executables on macOS
dmuai="./dmuai"

export MKL_NUM_THREADS=2

if [ "$1" == "-np" ]; then
  shift
  MKL_NUM_THREADS=$1
  shift
fi

FN="$1.DIR"
if [ -f "$FN" ]; then
  echo "Starting DMU with $1.DIR as directive file"
  if [ "$MKL_NUM_THREADS" -gt 1 ]; then
    echo "Running on $MKL_NUM_THREADS CPU/Cores"
  fi

  $dmu1 < "$1.DIR" > "$1.lst"

  if [ -e MODINF ]; then
    $dmuai >> "$1.lst"
  else
    echo "DMUAI Not started due to errors in DMU1"
  fi

  [ -e "$1.SOL" ] && rm "$1.SOL"
  [ -e SOL ] && mv SOL "$1.SOL"

  [ -e "$1.PAROUT" ] && rm "$1.PAROUT"
  [ -e PAROUT ] && mv PAROUT "$1.PAROUT"

  [ -e "$1.PAROUT_STD" ] && rm "$1.PAROUT_STD"
  [ -e PAROUT_STD ] && mv PAROUT_STD "$1.PAROUT_STD"

  [ -e "$1.LLIK" ] && rm "$1.LLIK"
  [ -e LLIK ] && mv LLIK "$1.LLIK"

  if [ -e RESIDUAL ]; then
    [ -e "$1.RESIDUAL" ] && rm "$1.RESIDUAL"
    mv RESIDUAL "$1.RESIDUAL"
  fi

  [ -e "$1.INBREED" ] && rm "$1.INBREED"

  if [ -e INBREED ]; then
    if [ ! -s INBREED ]; then
      rm INBREED
    else
      mv INBREED "$1.INBREED"
    fi
  fi

  rm -f CODE_TABLE DMU1.dir DMUAI.dir DMU_LOG DUMMY FSPAKWK Latest_parm LEVAL MODINF PARIN RCDATA_I RCDATA_R

  # Deleting matching files
  rm -f AINV* COR* DOM* IBD* PEDFILE* fort.*
else
  echo "File $1.DIR not in current directory"
fi
