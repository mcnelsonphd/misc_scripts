#!/bin/bash

############################################################################################
#
# Per_sample_raw_reads.sh
# 
# Script for creating raw Read 1 and Read 2 files for each individual sample of a 16S 
# This script take a users QIIME mapping file, the demultiplexed sequences file and the two raw read files as arguments.
# A time-stamped log file of all steps that are conducted is created.
# This file also shows input file MD5 checksums as well as MD5 checksums for all output files.
#
# TODO: Figure out a way to run the filter_fasta.py command on R1 and R2 in parallel to increase speed.
#       This should be feasible with gnu parallel but should figure out how best to implement
#
# Created by Michael C. Nelson on 2014-09-09.
# Last revised: 2014-09-12
# Revision #: 2
# Copyright 2014 Michael C. Nelson and the University of Connecticut. All rights reserved.
# 
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 
#
############################################################################################

## pull in macqiime path
source /macqiime/configs/bash_profile.txt

# Define intitial variables
DATE=`date +%Y-%m-%d`
TIME=`date +%H:%M`
TM=`date +%Y%m%d-%H%M`
LOG=Log_$TM.txt
MAP2=MAP_$TM.txt

#Create log file
echo "Script executed on $DATE at $TIME using the command call: Per_sample_raw_reads.sh $*" | tee $LOG
echo ''
# Initiail sanity checking for script operation
# Are the correct number of arguments given?
if [ "$1" = "" ]; then
    echo '' | tee  -a $LOG
    echo 'ERROR: No arguments were passed to script.' | tee  -a $LOG
    echo 'USAGE: Per_sample_raw_reads.sh Map.txt seqs.fna Undetermined_R1.fastq.gz Undetermined_R2.fastq.gz' | tee -a $LOG
    echo '' | tee  -a $LOG
    exit 1
elif [ "$(($# % 4))" != 0 ]; then
    echo '' | tee  -a $LOG
    echo 'ERROR: Invalid number of arguments given.' | tee  -a $LOG
    echo 'USAGE: Per_sample_raw_reads.sh Map.txt seqs.fna Undetermined_R1.fastq.gz Undetermined_R2.fastq.gz' | tee  -a $LOG
    echo '' | tee  -a $LOG
    exit 1
fi

# Are the require input files in the current dir?
if [ ! -f $MAP ] && [ ! -f $seqsfna ] && [ ! -f $READ1 ] && [ ! -f $READ2 ]; then
    echo '' | tee  -a $LOG
    echo 'ERROR: Required input files could not be found.' | tee  -a $LOG
    echo 'Script must be executed in the directory containing the Undertermined Read 1, Read 2, and Index files.'  | tee  -a $LOG
    echo 'USAGE: Per_sample_raw_reads.sh Map.txt seqs.fna Undetermined_R1.fastq.gz Undetermined_R2.fastq.gz' | tee  -a $LOG
    echo '' | tee  -a $LOG
    exit 1
else
    MAP=$1
    MD5MAP=`md5sum $1`
    seqsfna=$2
    MD5SEQS=`md5sum $2`
    READ1=$3
    MD5READ1=`md5sum $3`
    READ2=$4
    MD5READ2=`md5sum $4`
    echo '' | tee -a $LOG
    echo "Using $MAP as the input mapping file." | tee -a $LOG
    echo $MD5MAP| tee -a $LOG
    # Changing any carriage returns to newlines in the mapping file... so above is a slight lie
    tr '\r' '\n' <$MAP>$MAP2
    echo "Using $seqsfna as the input seqs.fna file." | tee -a $LOG
    echo $MD5SEQS | tee -a $LOG
    echo "Using $READ1 as the raw Read 1 file." | tee -a $LOG
    echo $MD5READ1 | tee -a $LOG
    echo "Using $READ2 as the raw Read 2 file." | tee -a $LOG
    echo $MD5READ2 | tee -a $LOG
fi

# Check to see if GNU parallel is installed.
if hash parallel 2>/dev/null; then
    SMP=TRUE
fi

# Unzip the input raw read files, won't affect them if they're already unzipped but will throw a non-lethal gzip error
DATE=`date +%Y-%m-%d`                                         
TIME=`date +%H:%M`                                            
echo '' | tee -a $LOG
echo "$DATE $TIME: Decompressing the input raw read files."

if $SMP; then
    parallel gunzip ::: $READ1 $READ2
else
    gunzip $READ1
    gunzip $READ2
fi

# If the input files were gzipped, then we now need to capture the file name w/o the .gz extension
R1=`echo $READ1 | sed 's/.gz//'`
R2=`echo $READ2 | sed 's/.gz//'`

# Lets start actually doing something why don't we
line=1                                 # We start with line 1
total=`grep -c '^' $MAP2`              # Determine how many lines are actually in the file (safer than wc -l if map doesn't have a final newline character)
(( samples = $total - 1 ))             # Total number of samples should be num lines minus header line
echo '' | tee  -a $LOG
echo "There are $samples samples in your mapfile." | tee  -a $LOG
echo '' | tee  -a $LOG
DATE=`date +%Y-%m-%d`                                         
TIME=`date +%H:%M`                                            
echo "$DATE $TIME: Proceeding to demultiplex the raw reads into per-sample R1 and R2 files." | tee  -a $LOG

while [ $line -lt $total ] 		                                  # While the current line number is less than the total number of sample lines,
do                             	                                  # Do the following actions
    DATE=`date +%Y-%m-%d`                                         # Reset Date
    TIME=`date +%H:%M`                                            # Reset Time
    printf "$DATE $TIME   " | tee  -a $LOG                        # Print time stamp so user can track progress rate
    printf "Sample: $line   " | tee  -a $LOG 	                  # First we'll print the current sample number
    (( line++ )) 	                                              # Now we need to increase the line count to dissociate from the header line
    sampleID=`sed -n "$line{p;q;}" $MAP2 | cut -f1,1`             # Now we find out what the sample ID is
    names=$sampleID.txt                                           # Set an output file for the read names based on the sample ID
    names2=$sampleID2.txt                                         # Create name of second names folder in case we're running in parallel mode
    searchID=$sampleID\_
    printf "$sampleID	" | tee  -a $LOG                          # Print what the name of the names file is for each sample
    touch $names                                                  # Create the output file as empty
    count=`grep -c $searchID $seqsfna`                            # Check to see how many reads are in seqs.fna
    echo "$count seqs" | tee  -a $LOG                             # Print out how many sequences are present for the sample
    grep $searchID $seqsfna | tr -d '>' | cut -d\  -f2,2 > $names # Compile the list of SeqIDs for filter_fasta command
    RAW1=$sampleID"_R1.fastq"                                     # Define the Read1 output file
    RAW2=$sampleID"_R2.fastq"                                     # Define the Read2 output file
    if $SMP; then
        cp $names $names2
        parallel -N3 echo filter_fastq.py -f {1} -o {2} -s {3} ::: $R1 $RAW1 $names $R2 $RAW $names2
    else
        filter_fasta.py -f $R1 -o $RAW1 -s $names                     # Create Read1 raw read file using QIIME
        filter_fasta.py -f $R2 -o $RAW2 -s $names                     # Create Read2 raw read file using QIIME
    fi
    if $SMP; then                                                 # Now we need to compress the files b/c thats what the SRA wants.
        parallel gzip ::: $RAW1 $RAW2
    else
        gzip $RAW1
        gzip $RAW2
    fi
    rm $names                                                     # We no longer need the names file so let's get rid of it
done

# Cleanup phase: step 1, re-zip the input files to again save file space and delete the "cleaned" map file.
DATE=`date +%Y-%m-%d`                                         
TIME=`date +%H:%M`                                            
echo '' | tee -a $LOG
echo "$DATE $TIME: Recompressing the input raw read files."
rm $MAP2

if $SMP; then
    parallel gzip ::: $R1 $R2
else
    gzip $R1
    gzip $R2
fi

# Step 2, calculate md5 checksums for all of the raw read files. These are needed for SRA submissions and also just nice to have.
echo '' | tee -a $LOG
DATE=`date +%Y-%m-%d`                                         
TIME=`date +%H:%M`                                            
echo "$DATE $TIME: Calculating md5 checksum values for all sample files." | tee -a $LOG
md5sum *_R1.fastq.gz | tee -a $LOG | tee md5sums.txt     # This could be done in parallel, but, eh.
md5sum *_R2.fastq.gz | tee -a $LOG | tee -a md5sums.txt  # Note that if the files are decompressed and then recompressed using a different method
                                                         # (e.g. pigz) then the MD5 values will change.

DATE=`date +%Y-%m-%d`                                         
TIME=`date +%H:%M`                                            
echo '' | tee -a $LOG
echo '' | tee -a $LOG
echo "$DATE $TIME: Script is now finished." | tee -a $LOG
echo "Each sample should have a gzip compressed R1 and R2 read file that you will need to upload to the SRA." | tee -a $LOG
echo "md5 checksum values have been calculated for these files and can be found in the md5sums.txt file." | tee -a $LOG
