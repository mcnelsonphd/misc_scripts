#Miscellaneous scripts

This is a collection of scripts that I've made to do various processing of sequencing data.

###demultiplex_reads.py
Demultiplexes single samples from a non-demultiplexed MiSeq run according to a provided Illumina index name, writing the output to a user specific directory which is created if needed, and naming the output files according to the ID provided by the user.  
- No outside dependencies outside of standard python2.7.  
- Contains a built-in dictionary of NexteraXT, TruSeq6 and TruSeq8 index sequences. 

Example usage:
```Python
demultiplex_reads.py -f Undetermined_R1.fastq.gz -r Undetermined_R2.fastq.gz -b Undetermined_I1.fastq.gz -o MyGenome_reads_folder -i N701 -n MyGenomeID
```


###Per_Sample_Raw_Reads.sh
Creates the raw read1/read2 files for each sample of an amplicon sequencing experiment that are needed for submission to one of the INSDC repositories. The raw files are gzip compressed and will be named according to the SampleIDs found in the user provided mapping file that was used for initial demultiplexing. Along with a log file, the script will also create a file containing the MD5 checksum values for each of the raw files that are necessary when submitting the reads to the repository.

Example usage:
```
Per_Sample_Raw_Reads.sh Map.txt seqs.fna Undetermined_R1.fastq.gz Undetermined_R2.fastq.gz
```
NOTE: Relies on dual\_fastq\_filter.py which must be locatable in the users $PATH.


###dual_fastq_filter.py
Custom script to create a paired set of filtered fastq sequences based on an input seqID list. Output files are automatically gzip compressed. Requires a standard python 2.7 install with no other outside dependencies.

Example usage:
```Python
dual_fastq_filter.py -f Read1.fastq.gz -r Read2.fastq.gz -o Output1.fastq.gz -p Output2.fastq.gz -n Names.txt
```

###fastq_filter.py
Custom script to create a filtered set of fastq sequences based on an input seqID list. Output files are automatically gzip compressed. Based upon and uses some code derived from the filter_fasta.py command of QIIME.  Requires a standard python 2.7 install with no other outside dependencies.

Example usage:  
`fastq_filter.py -i Input.fastq -o Output.fastq.gz -n Names.txt`
	
## Deprecated Scripts
The following scripts are no longer under development and are not recommended for use. 
###demultiplex_reads.sh
Demultiplexes samples from non-demultiplexed MiSeq runs. Requires the script to be called from the directory where the raw Undetermined files are located and outputs to the results to individual directories within the current directory.

Example usage:
```shell
demultiplex_reads.sh MyGenomeA N701 MyGenomeB N706
```
NOTE: Relies on QIIME for backend processing and a set of fixed mapping files to operate that are not provided.

