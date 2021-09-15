#!/bin/bash

# local, Pb-Pb, 2011, real data, EMCal

list="localaod-11-pbpb-real.txt" # path to the file with paths to input AOD files
isPbPb=1 # 1 if input are Pb-Pb collisions, 0 otherwise
year=2011 # year when the real data was taken
isMC=0 # 1 if simulated (Monte Carlo) input, 0 if real data
nFiles="$1" # number of input AOD files to process

bash clean.sh
aliroot -b -q "runLocal.C(\"$list\", $isPbPb, \"$year\", $isMC, $nFiles)" 1>stdout.log 2>stderr.log && echo "Done" || echo "Error"
