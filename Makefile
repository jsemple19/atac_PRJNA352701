#!/usr/bin/make -f

########### VARIABLES #########
SRR := SRR5000676 SRR5000679 SRR5000682

#make list of fastq files based on SRR accession
fastqFiles := $(addsuffix .fastq.gz,$(addprefix fastq/,$(SRR)))
OBJECTS := \
	$(addsuffix _1_fastqc.html, $(addprefix fastQC/,$(SRR))) \
	$(addsuffix _2_fastqc.html, $(addprefix fastQC/,$(SRR))) \
	$(addsuffix _cutadapt.txt, $(addprefix cutadapt/report_,$(SRR))) \
 	$(addsuffix _flagstats.txt, $(addprefix aln/report_,$(SRR))) \
 	$(addsuffix _stats.txt, $(addprefix aln/report_,$(SRR))) \
 	$(addsuffix _picard.txt, $(addprefix aln/report_,$(SRR))) \
	$(addsuffix _model.r, $(addprefix macs2peaks/,$(SRR))) \
 	$(addsuffix _peaks.narrowPeak, $(addprefix macs2peaks/,$(SRR))) \
 	$(addsuffix _peaks.xls, $(addprefix macs2peaks/,$(SRR))) \
 	$(addsuffix _summits.bed, $(addprefix macs2peaks/,$(SRR))) 


# list of intermediate files to be deleted at end of run
intermediateFiles := $(addprefix aln/, $(addsuffix .sam, $(SRR))) \
	$(addprefix aln/, $(addsuffix .noDup.bam, $(SRR))) \
	$(addprefix aln/, $(addsuffix .sorted.bam, $(SRR))) \
	$(addprefix cutadapt/, $(addsuffix _1.fastq.gz, $(SRR))) \
	$(addprefix cutadapt/, $(addsuffix _2.fastq.gz, $(SRR))) 

# list secondary files to keep
secondaryFiles := $(addprefix aln/, $(addsuffix .noMito.bam, $(SRR)))

#OBJECTS := $(addsuffix _fastqc.html,$(addprefix fastqc/,$(SRR))) 
genomefile:=/Users/semple/Documents/MeisterLab/GenomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.fa



########### RECIPES ##########
all:  $(fastqFiles) $(OBJECTS)

.PHONY: all clean cleanall
.INTERMEDIATE: $(intermediateFiles)
.SECONDARY: $(secondaryFiles)

cleanall: 
	rm -f $(OBJECTS)
	rm -f $(intermediateFiles)
	rm -f $(secondaryFiles)
	rm -f fastq/*
	
clean:
	rm -f $(intermediateFiles)

#download reads from SRA
$(fastqFiles):
	mkdir -p fastq
	fastq-dump -O fastq -X 1000 --split-files --gzip $(SRR)
	echo "Data downloaded:" > fastq/logfile
	date >>fastq/logfile
	echo $(SRR) >>fastq/logfile
	
# run fastqc on downloaded sequences
fastQC/%_fastqc.html: fastq/%.fastq.gz
	mkdir -p fastQC
	fastqc $^ -o fastQC 

# use cutadapt to trim nextera sequence in paired end mode
cutadapt/%_1.fastq.gz cutadapt/%_2.fastq.gz cutadapt/report_%_cutadapt.txt: fastq/%_1.fastq.gz fastq/%_2.fastq.gz
	mkdir -p cutadapt
	cutadapt -a CTGTCTCTTATA -A CTGTCTCTTATA -q 10 -m 25 -o cutadapt/$*_1.fastq.gz -p cutadapt/$*_2.fastq.gz \
	fastq/$*_1.fastq.gz fastq/$*_2.fastq.gz > cutadapt/report_$*_cutadapt.txt

# use bwa-mem to align sequences
aln/%.sam: $(genomefile) cutadapt/%_1.fastq.gz cutadapt/%_2.fastq.gz
	mkdir -p aln
	bwa mem -t 2 $^ > $@
	
# use samtools to convert to bam, sort and clean
aln/%.sorted.bam: aln/%.sam
	samtools view -q 10 -F2304 -u $^ | samtools sort -o $@

# 	get alignment stats
aln/report_%_flagstats.txt: aln/%.sorted.bam
	samtools flagstat $^ > $@

aln/report_%_stats.txt: aln/%.sorted.bam
	samtools stats $^ > $@

# find duplicates with picard  (path to picard should be set in $PICARD variable in .bash_profile or in session)
aln/%.noDup.bam aln/report_%_picard.txt: aln/%.sorted.bam $(PICARD)
	java -Xmx5g -jar $(PICARD) MarkDuplicates I=aln/$*.sorted.bam O=aln/$*.noDup.bam M=aln/report_$*_picard.txt

# 	remove mitochondrial reads
aln/%.noMito.bam: aln/%.noDup.bam
	samtools view -h -F1024 $^ | grep -v -e '\tMtDNA\t' | samtools view -b -> $@

# call peaks with macs2peaks (need to activate conda python 2.7 emvironment. you need to invoke bash shell and )
macs2peaks/%_model.r macs2peaks/%_peaks.narrowPeak macs2peaks/%_peaks.xls macs2peaks/%_summits.bed: aln/%.noMito.bam
	mkdir -p ./macs2peaks
	( bash -c "source ${HOME}/anaconda/bin/activate py27; \
		macs2 callpeak --keep-dup all --format BAMPE -t aln/$*.noMito.bam -n macs2peaks/$* ")

