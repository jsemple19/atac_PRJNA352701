# atac_PRJNA352701
The makefile downloads data from given SRS accessions to the SRA for L3 worms and aligns them to the genome.

These scripts are adapted to run as an array job on the vital-it cluster.
The job is launched with job_atacAlign.sh. Each SRR file is downloaded and processed in a separate job in the array and processed by the makefile.

