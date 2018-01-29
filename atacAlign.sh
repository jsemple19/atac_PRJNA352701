#!/bin/bash

#BSUB -L /bin/bash
#BSUB -o atacAlign-output.txt
#BSUB -e atacAlign-error.txt
#BSUB -J atacAlign
#BSUB -u jennifer.semple@izb.unibe.ch
#BSUB -N
#BSUB -n 3
#BSUB -R "span[ptile=


module add UHTS/Quality_control/fastqc/0.11.5      #fastqc
module add UHTS/Quality_control/cutadapt/1.13     #cutadapt
module add UHTS/Aligner/bwa/0.7.15                 #bwa
module add UHTS/Analysis/samtools/1.4             #samtools
module add UHTS/Analysis/picard-tools/2.9.0        #picard.jar
module add UHTS/Quality_control/qualimap/2.2.1    #qualimap.jar
module add UHTS/Analysis/BEDTools/2.26.0        #bedToBam
module add UHTS/Analysis/macs/2.1.1.20160309;    #macs2

#before running script, make directory on scratch
dataDir=/scratch/cluster/weekly/jsemple/public/atac
scriptDir=/home/jsemple/publicData/atac_PRJNA352701
mkdir -p $dataDir
#go to the scratch directory and copy script from home directory and run from there
cd $dataDir
cp ${scriptDir}/Makefile
cp ${scriptDir}/Tn5shift.sh

make -j 3 all

