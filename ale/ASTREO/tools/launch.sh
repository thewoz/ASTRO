#!/bin/bash

TLE='sentinel2A_02072019.TLE'
SATNUM=40697

NFRAME=9
OBS1=urbe
OBS2=colle

# erase files
rm *.$SATNUM *.$SATNUM.* *.$OBS1 *.$OBS2

# compiling
./compile.sh

# SGP4 propagation
./viewsat.x $TLE

# synthetic photos
./photosat.x
./photostar.x

# sort stars
./sort.x $NFRAME 'photostar.'${SATNUM}'.'${OBS1}'.dat'
./sort.x $NFRAME 'photostar.rot.'${SATNUM}'.'${OBS1}'.dat'
./sort.x $NFRAME 'photostar.'${SATNUM}'.'${OBS2}'.dat'
./sort.x $NFRAME 'photostar.rot.'${SATNUM}'.'${OBS2}'.dat'

# rename sorted files
mv 'photostar.'${SATNUM}'.'$OBS1'.dat_sorted' 'photostar.'${SATNUM}'.'$OBS1'.dat'
mv 'photostar.rot.'${SATNUM}'.'$OBS1'.dat_sorted' 'photostar.rot.'${SATNUM}'.'$OBS1'.dat'
mv 'photostar.'${SATNUM}'.'$OBS2'.dat_sorted' 'photostar.'${SATNUM}'.'$OBS2'.dat'
mv 'photostar.rot.'${SATNUM}'.'$OBS2'.dat_sorted' 'photostar.rot.'${SATNUM}'.'$OBS2'.dat'
