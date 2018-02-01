#!/usr/bin/make -f
## required software: fastqc, cutadapt, bwa, samtools, picard, qualimap, macs2peaks

########### VARIABLES #########
#SRR := SRR5000676 SRR5000679 SRR5000682

PICARD:=/software/UHTS/Analysis/picard-tools/2.9.0/bin/picard.jar

# make list of input fastq files based on SRR accession
fastqFiles := $(addsuffix .fastq.gz,$(addprefix fastq/,$(SRR)))

# make a list of objects based on SRR accession
OBJECTS := \
	$(addsuffix _1_fastqc.html, $(addprefix fastQC/,$(SRR))) \
	$(addsuffix _2_fastqc.html, $(addprefix fastQC/,$(SRR))) \
	$(addsuffix _cutadapt.txt, $(addprefix cutadapt/report_,$(SRR))) \
 	$(addsuffix _flagstats.txt, $(addprefix aln/report_,$(SRR))) \
 	$(addsuffix _stats.txt, $(addprefix aln/report_,$(SRR))) \
 	$(addsuffix _picard.txt, $(addprefix aln/report_,$(SRR))) \
 	$(addsuffix _peaks.narrowPeak, $(addprefix macs2peaks/,$(SRR))) \
 	$(addsuffix _peaks.xls, $(addprefix macs2peaks/,$(SRR))) \
 	$(addsuffix _summits.bed, $(addprefix macs2peaks/,$(SRR))) \
 	$(addsuffix _insert_size_metrics.txt, $(addprefix aln/report_picard_,$(SRR))) \
 	$(addsuffix _insert_size_histogram.pdf, $(addprefix aln/report_picard_,$(SRR))) \
 	$(addsuffix .pdf, $(addprefix aln/report_qualimap_,$(SRR))) 

# list of intermediate files to be deleted at end of run
intermediateFiles := \
	$(addsuffix _1.fastq.gz, $(addprefix fastq/,$(SRR))) \
        $(addsuffix _2.fastq.gz, $(addprefix fastq/,$(SRR))) \
	$(addprefix aln/, $(addsuffix .sam, $(SRR))) \
	$(addprefix aln/, $(addsuffix .noDup.bam, $(SRR))) \
	$(addprefix aln/, $(addsuffix .sorted.bam, $(SRR))) \
	$(addprefix cutadapt/, $(addsuffix _1.fastq.gz, $(SRR))) \
	$(addprefix cutadapt/, $(addsuffix _2.fastq.gz, $(SRR))) \
	$(addprefix aln/, $(addsuffix .bed, $(SRR))) \
	$(addprefix aln/, $(addsuffix .Tn5shifted.bed, $(SRR)))  

# list secondary files to keep
secondaryFiles := $(addprefix aln/, $(addsuffix .noMito.bam, $(SRR)))

# reference geneme file location
genomefile:=/home/jsemple/publicData/genomeVer/sequence/c_elegans.PRJNA13758.WS250.genomic.fa

########### RECIPES ##########
all:  $(fastqFiles) $(OBJECTS)

.PHONY: all clean cleanall
.INTERMEDIATE: $(intermediateFiles)
.SECONDARY: $(secondaryFiles)

cleanall: 
	rm -f $(OBJECTS)
	rm -f $(intermediateFiles)
	rm -f $(secondaryFiles)
	
clean:
	rm -f $(intermediateFiles)

#download reads from SRA
$(fastqFiles):
	mkdir -p fastq
	fastq-dump -O fastq --split-files --gzip $(SRR)
	echo "Data downloaded:" >> fastq/logfile
	date >>fastq/logfile
	echo $(SRR) >>fastq/logfile_${SRR}_.txt
	touch $@
	
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
	samtools view -q 30 -F 1804 -u $^ | samtools sort -o $@

# 	get alignment stats
aln/report_%_flagstats.txt: aln/%.sorted.bam
	samtools flagstat $^ > $@

aln/report_%_stats.txt: aln/%.sorted.bam
	samtools stats $^ > $@

# find duplicates with picard  (path to picard should be set in $PICARD variable in .bash_profile or in session)
aln/%.noDup.bam aln/report_%_picard.txt: aln/%.sorted.bam $(PICARD)
	java -Xmx5g -jar ${PICARD}  MarkDuplicates I=aln/$*.sorted.bam O=aln/$*.noDup.bam M=aln/report_$*_picard.txt

# Get insert size statistics and plots with picard and qualimap (set path to $QUALIMAP in .bash_profile)
aln/report_picard_%_insert_size_metrics.txt aln/report_picard_%_insert_size_histogram.pdf \
aln/report_qualimap_%.pdf: aln/%.noDup.bam 
	java -jar ${PICARD}  CollectInsertSizeMetrics I=aln/$*.noDup.bam \
       O=aln/report_picard_$*_insert_size_metrics.txt \
       H=aln/report_picard_$*_insert_size_histogram.pdf
	qualimap bamqc -bam aln/$*.noDup.bam -c -outdir aln -outfile report_qualimap_$*.pdf -outformat PDF


# 	remove duplicate and mitochondrial reads
aln/%.noMito.bam: aln/%.noDup.bam
	samtools view -h -F 1024 $^ | grep -v -e '\tMtDNA\t' | samtools view -b -> $@
	
# convert bam to bed and shift by Tn5 footprint using Tn5shift.sh script (awk within file was
# not running possibly for some reason)
aln/%.Tn5shifted.bed: aln/%.noMito.bam
	bedtools bamtobed -i aln/$*.noMito.bam > aln/$*.bed
	./Tn5shift.sh aln/$*

# call peaks with macs2peaks (need to activate conda python 2.7 emvironment. you need to invoke bash shell and )
macs2peaks/%_peaks.narrowPeak macs2peaks/%_peaks.xls macs2peaks/%_summits.bed: aln/%.Tn5shifted.bed
	mkdir -p ./macs2peaks
	macs2 callpeak -t aln/$*.Tn5shifted.bed -f BED -n macs2peaks/$* \
	-g 9e7 --nomodel --extsize 150 --shift -75 -B --keep-dup all --call-summits
