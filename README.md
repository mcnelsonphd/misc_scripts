misc_scripts
============

A collection of scripts I've made to do various processing of sequencing data.

**demultiplex_reads.sh** - Demultiplexes samples from non-demultiplexed MiSeq runs. Requires the script to be called from the directory where the raw Undetermined files are located and outputs to the results to individual directories within the current directory.
NOTE: Relies on Qiime backend and a set of fixed mapping files to operate.

Example usage:
	`demultiplex_reads.sh MyGenomeA N701 MyGenomeB N706`

**demultiplex_reads.py** - Demultiplexes single samples from a non-demultiplexed MiSeq run according to a provided Illumina index name, writing the output to a newly created directory named according to the user.  
No outside dependencies outside of standard python2.7.  
Contains a built-in dictionary of NexteraXT, TruSeq6 and TruSeq8 index sequences. 

Example usage:
	`demultiplex_reads.py -f Undetermined_R1.fastq.gz -r Undetermined_R2.fastq.gz -b Undetermined_I1.fastq.gz -o MyGenome_reads -i N701 -n MyGenome`

**Per_Sample_Raw_Reads.sh** - Creates raw read files for each sample of an amplicon sequencing experiment.  
NOTE: Relies on dual\_fastq\_filter.py

Example usage:
	`Per_Sample_Raw_Reads.sh Map.txt seqs.fna Undetermined_R1.fastq.gz Undetermined_R2.fastq.gz`
	
**dual_fastq_filter.py** - Custom script to create a paired set of filtered fastq sequences based on an input seqID list. Output files are automatically gzip compressed. 

Example usage:
	`dual_fastq_filter.py -f Read1.fastq.gz -r Read2.fastq.gz -o Output1.fastq.gz -p Output2.fastq.gz -n Names.txt`

**fastq_filter.py** - Custom script to create a filtered set of fastq sequences based on an input seqID list. Output files are automatically gzip compressed. Based upon and uses some code derived from the filter_fasta.py command of QIIME. 

Example usage:
	`fastq_filter.py -i Input.fastq -o Output.fastq.gz -n Names.txt`
