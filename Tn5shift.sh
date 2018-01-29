#!/bin/bash

#this scriptt takes in one command line argument which is the base of the bed filename 
# (includes path, but without .bed extension)

baseName=$1
echo $baseName
inputBed=`echo ${baseName}.bed`
outputBed=`echo ${baseName}.Tn5shifted.bed`

awk -F'\t' 'BEGIN {OFS = FS} {if($6 == "+") { $2 = $2 + 4; $3 = $2 } \
else if($6 == "-") {$3 = $3 - 5; $2= $3 } print $0 }' $inputBed > $outputBed
