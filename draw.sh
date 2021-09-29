#!/bin/bash

# $1 train output file
# $2 directory with efficieny files
# $3 directory with feed-down files
# $4 modifier of the branch name
# $5 mode: 0 uncorrected, 1 simulated, 2 corrected

AnalysisOutput="$1"
PathEff="$2"
PathFD="$3"
Flag="$4"
Mode="$5"

EffK="$PathEff/OutputK0s.root"
EffL="$PathEff/OutputLambda.root"
EffAL="$PathEff/OutputALambda.root"
FDL="$PathFD/OutputFDLambda.root"
FDAL="$PathFD/OutputFDALambda.root"

# Create a directory named after the flag
#if [ ! -z $Flag ] # string is not empty
#then
##rm -rf $Flag
#mkdir -p $Flag
#fi

LogK0S="stdouterr_K0s.log"
LogLambda="stdouterr_Lambda.log"
LogALambda="stdouterr_ALambda.log"
LogRatios="stdouterr_Ratio.log"

date +"%Y-%m-%d %H-%M-%S"
# K0S
echo "Processing K0S"
root -b -q "compileAndRunDrawResults.C(\"$AnalysisOutput\", 0, \"$EffK\", \"\", 0, \"$Flag\", $Mode)" > "$LogK0S" 2>&1
date +"%Y-%m-%d %H-%M-%S"
# Lambda
echo "Processing Lambda"
root -b -q "compileAndRunDrawResults.C(\"$AnalysisOutput\", 1, \"$EffL\", \"$FDL\", 0, \"$Flag\", $Mode)" > "$LogLambda" 2>&1
date +"%Y-%m-%d %H-%M-%S"
# anti-Lambda
echo "Processing anti-Lambda"
root -b -q "compileAndRunDrawResults.C(\"$AnalysisOutput\", 2, \"$EffAL\", \"$FDAL\", 0, \"$Flag\", $Mode)" > "$LogALambda" 2>&1
date +"%Y-%m-%d %H-%M-%S"
# Lambda/K0S
echo "Processing ratios"
root -b -q "compileAndRunDrawResults.C(\"\", 0, \"\", \"\", 1, \"$Flag\", $Mode)" > "$LogRatios" 2>&1
date +"%Y-%m-%d %H-%M-%S"

bash sort-output.sh

grep "MassPtSigma:" "$LogK0S" > FitK0s.log
grep "MassPtSigma:" "$LogLambda" > FitLambda.log
grep "Statistics:" "$LogK0S" > StatK0s.log
grep "Statistics:" "$LogLambda" > StatLambda.log

echo "Errors in $LogK0S:"
grep [Ee]rror "$LogK0S"
echo "Errors in $LogLambda:"
grep [Ee]rror "$LogLambda"
echo "Errors in $LogALambda:"
grep [Ee]rror "$LogALambda"
echo "Errors in "$LogRatios":"
grep [Ee]rror "$LogRatios"

echo "Deleting"
find . -empty -delete

exit 0

