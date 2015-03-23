misc_scripts
============

A collection of scripts I've made to do various processing of sequencing data.

**demultiplex_reads.sh** - Demultiplexes samples from non-demultiplexed MiSeq runs.
Relies on Qiime backend and a set of fixed mapping files to operate.

**Per_Sample_Raw_Reads.sh** - Creates raw read files for each sample of an amplicon sequencing experiment.
Relies on fastq_filter.py

**fastq_filter.py** - Custom script to create a filtered set of fastq sequences based on an input seqID list. Output files are automatically gzip compressed. Based upon and uses some code derived from the filter_fasta.py command of QIIME. 
