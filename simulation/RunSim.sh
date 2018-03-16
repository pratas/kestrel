#!/bin/bash
GET_GOOSE=1;
GET_GARGAMMEL=1;
GET_KESTREL=1;
GET_FALCON=1;
GET_BWA=1;
GET_BOWTIE=1;
#
SIMULATE=1;
#
RUN=1;
#
PLOT=1;
#==============================================================================
# GET GOOSE
if [[ "$GET_GOOSE" -eq "1" ]]; then
  rm -fr goose/ goose-*
  git clone https://github.com/pratas/goose.git
  cd goose/src/
  make
  cd ../../
  cp goose/src/goose-* .
  cp goose/scripts/*.sh .
fi
#==============================================================================
# GET FALCON
if [[ "$GET_GECO" -eq "1" ]]; then
  git clone https://github.com/pratas/geco.git
  cd geco/src/
  cmake .
  make
  cp GeCo ../../
  cd ../../
fi
#==============================================================================
# GET GARGAMMEL
if [[ "$GET_GARGAMMEL" -eq "1" ]]; then
  rm -fr gargammel/
  git clone --depth 1 https://github.com/grenaud/gargammel.git
  cd gargammel/
  make
  cd ..
fi
#==============================================================================
# SIMULATE
if [[ "$SIMULATE" -eq "1" ]]; then
  cd gargammel/
  mkdir data
  mkdir data/bact
  make bacterialex
  cp -v bactDBexample/clovis/fasta/* data/bact/
  mkdir data/endo
  cd data/endo/
  wget https://raw.githubusercontent.com/pratas/goose/master/scripts/GetHumanFastaReference.sh
  chmod +x GetHumanFastaReference.sh
  ./GetHumanFastaReference.sh
  rm -f GetHumanFastaReference.sh
  samtools faidx HS.fa
  cd ../../
    
fi
