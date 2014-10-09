#!/bin/bash +x

## pull in macqiime path
source /macqiime/configs/bash_profile.txt

# Goal:
# User calls the script along with a series of command line arguments that correspond to the pattern of 'SampleID Index'
# The script will then execute a loop where it will read sets of arguments and process the files accordingly.
# The output will be in folders labeled according to the given SampleID and the demultiplexing will use the appropriate Index mapping file that has already been created
# A sanity check should be in place at the very beginning to ensure that the number of arguments passed to the script is even.

# Define script variables
DATE=`date +%Y-%m-%d`
TIME=`date +%H:%M`
TM=`date +%Y%m%d-%H%M`
LOG=demultiplex_reads_log_$TM.txt
MAP_LOC=/Volumes/DATA_RAID/Reference_Files/Index_Maps
READ1=$PWD/Undetermined_S0_L001_R1_001.fastq.gz
READ2=$PWD/Undetermined_S0_L001_R2_001.fastq.gz
INDEX=$PWD/Undetermined_S0_L001_I1_001.fastq.gz

echo "Script executed on $DATE at $TIME using the command call: demultiplex_reads.sh $*" > $LOG

# Initiail sanity checking for script operation
# Are the correct number of arguments given?
if [ "$1" = "" ]; then
    echo '' | tee -a $LOG
    echo 'ERROR: No arguments were passed to script.' | tee -a $LOG
    echo 'USAGE: demultiplex.sh SampleID1 Index1 SampleID2 Index2 ...' | tee -a $LOG
    echo '' | tee -a $LOG
    exit 1
elif [ "$(($# % 2))" != 0 ]; then
    echo '' | tee -a $LOG
    echo 'ERROR: Number of arguments given was not even, script must be passed pairs of sampleID/index.' | tee -a $LOG
    echo 'USAGE: demultiplex.sh SampleID1 Index1 SampleID2 Index2 ...' | tee -a $LOG
    echo '' | tee -a $LOG
    exit 1
fi
# Are the require input files in the current dir?
if [ ! -f Undetermined_S0_L001_R1_001.fastq.gz ] && [ ! -f Undetermined_S0_L001_R2_001.fastq.gz ] && [ ! -f Undetermined_S0_L001_I1_001.fastq.gz ]; then
    echo '' | tee -a $LOG
    echo 'ERROR: Required input files could not be found.' | tee -a $LOG
    echo 'Script must be executed in the directory containing the Undertermined Read 1, Read 2, and Index files.'  | tee -a $LOG
    echo 'USAGE: demultiplex.sh SampleID1 Index1 SampleID2 Index2 ...' | tee -a $LOG
    echo '' | tee -a $LOG
    exit 1
fi

## Demultiplexing code block, should never be necessary to change.
while [ $1 ]
do
    echo '' | tee -a $LOG
    echo "Processing sample $1 using index $2" | tee -a $LOG
    mkdir $1
    split_libraries_fastq.py -i $READ1 -o $1/R1/ -m $MAP_LOC/$2_Map.txt -b $INDEX --store_demultiplexed_fastq -r 250 -p 0 -n 250 -q 0 --barcode_type 12 --max_barcode_errors 0
    sed 's|^@X_[0-9]* |@|g' $1/R1/seqs.fastq | sed 's| orig_bc.*||g' >$1/$1_R1.fastq
    echo "Sample $1 Read 1 done!" | tee -a $LOG
    if [ -s $1/$1_R1.fastq ]; then
        split_libraries_fastq.py -i $READ2 -o $1/R2/ -m $MAP_LOC/$2_Map.txt -b $INDEX --store_demultiplexed_fastq -r 250 -p 0 -n 250 -q 0 --barcode_type 12 --max_barcode_errors 0
        sed 's|^@X_[0-9]* |@|g' $1/R2/seqs.fastq | sed 's| orig_bc.*||g' >$1/$1_R2.fastq
        echo "Sample $1 Read 2 done!" | tee -a $LOG
    else
        echo 'Sample $1 has no Read 1 reads, skipping Read 2 and proceeding to next sample.' | tee -a $LOG
    fi
    shift 2
done
