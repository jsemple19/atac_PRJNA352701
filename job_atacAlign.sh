#!/bin/bash

#BSUB -L /bin/bash
#BSUB -o atacAlign-output__%I.txt
#BSUB -e atacAlign-error__%I.txt
#BSUB -J atacAlign
#BSUB -u jennifer.semple@izb.unibe.ch
#BSUB -N
##BSUB -n 3
##BSUB -R "span[ptile=4]"
##BSUB –R "rusage[mem=4096]"
##BSUB -M 4194304
#BSUB -J array[1-3]


module add UHTS/Analysis/sratoolkit/2.8.2.1;	#fastq-dump
module add UHTS/Quality_control/fastqc/0.11.5      #fastqc
module add UHTS/Quality_control/cutadapt/1.13     #cutadapt
module add UHTS/Aligner/bwa/0.7.15                 #bwa
module add UHTS/Analysis/samtools/1.4             #samtools
module add UHTS/Analysis/picard-tools/2.9.0        #picard.jar
module add UHTS/Quality_control/qualimap/2.2.1    #qualimap.jar
module add UHTS/Analysis/BEDTools/2.26.0        #bedToBam
module add UHTS/Analysis/macs/2.1.1.20160309;    #macs2

#before running script, make directory on scratch
dataDir=/scratch/cluster/weekly/jsemple/publicData/atac_PRJNA352701
scriptDir=/home/jsemple/publicData/atac_PRJNA352701
mkdir -p $HOME/.ncbi
echo '/repository/user/main/public/root = "/scratch/cluster/weekly/jsemple/publicData/"' > ${HOME}/.ncbi/user-settings.mkfg
#go to the scratch directory and copy script from home directory and run from there
cd $dataDir
cp -u ${scriptDir}/Makefile $dataDir/
cp -u ${scriptDir}/Tn5shift.sh $dataDir/

SRRlist=(SRR5000676 SRR5000679 SRR5000682)
let i=${LSB_JOBINDEX}-1
SRRfile=${SRRlist[$i]}

make SRR=${SRRfile} all
