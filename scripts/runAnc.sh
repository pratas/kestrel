#!/bin/bash
GET_GOOSE=1;
GET_KESTREL=1;
RUN_KESTREL=1;
#==============================================================================
# GET GOOSE
if [[ "$GET_GOOSE" -eq "1" ]]; then
  rm -fr goose/ goose-*
  git clone https://github.com/pratas/goose.git
  cd goose/src/
  make
  cp goose-* ../../
  cd ../../
fi
#==============================================================================
# GET KESTREL
if [[ "$GET_KESTREL" -eq "1" ]]; then
  rm -rf kestrel
  git clone https://github.com/pratas/kestrel.git
  cd kestrel/src/
  cmake .
  make
  cp KESTREL ../../
  cd ../../
fi
#==============================================================================
# RUN KESTREL
if [[ "$RUN_KESTREL" -eq "1" ]]; then
  ./goose-extractreadbypattern olonicum < ~/DATABASE/DB/DB.fa > PP.fa
  (time ./KESTREL -v -n 10 -F -m 13:200:1:4/10 -m 11:100:1:0/0 -g 0.9 -o filtered_pp.fq -t 0.5 PP.fa ~/DATABASE/DENISOVA/DENI ) &> REPORT;
fi
#==============================================================================
