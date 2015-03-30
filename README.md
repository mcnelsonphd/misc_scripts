misc_scripts
============

A collection of scripts I've made to do various processing of sequencing data.

**demultiplex_reads.sh** - Demultiplexes samples from non-demultiplexed MiSeq runs.
Relies on Qiime backend and a set of fixed mapping files to operate.

**demultiplex_reads.py** - Demultiplexes single samples from a non-demultiplexed MiSeq run.
No outside dependencies outside of standard python2.7.
Contains a built-in dictionary of NexteraXT, TruSeq6 and TruSeq8 index sequences. 

**Per_Sample_Raw_Reads.sh** - Creates raw read files for each sample of an amplicon sequencing experiment.
Relies on dual_fastq_filter.py

**dual_fastq_filter.py** - Custom script to create a paired set of filtered fastq sequences based on an input seqID list. Output files are automatically gzip compressed. 

**fastq_filter.py** - Custom script to create a filtered set of fastq sequences based on an input seqID list. Output files are automatically gzip compressed. Based upon and uses some code derived from the filter_fasta.py command of QIIME. 
