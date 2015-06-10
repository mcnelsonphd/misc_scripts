#!/bin/sh

# Demultiplex_reads_multi.sh
# 
# Script for demultiplexing multiple genomes as a batch
#
# Created by Michael C. Nelson on 2015-05-16.
# Last revised: 2015-06-10
# Revision #: 2
# Copyright 2015 Michael C. Nelson and the University of Connecticut. All rights reserved.
#
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 

#Set Script Name variable
SCRIPT=`basename ${BASH_SOURCE[0]}`

#Help function
function HELP {
    echo "$SCRIPT <Required Options> <GenomeID:IndexID pairs>"
    echo "Required arguments."
    echo "-f  The path to the Undetermined Read 1 file."
    echo "-r  The path to the Undetermined Read 2 file."
    echo "-b  The path to the Undetermined Index file."
    echo "Optional arguments."
    echo "-h  Displays this help message. No further functions are performed."\\n
    echo "Example usage:\\n$SCRIPT -f Undetermined_R1.fastq.gz -r Undetermined_R2.fastq.gz -b Undetermined_I1.fastq.gz Genome1ID N701 Genome2ID N702 ..."\\n
    exit 1
}

#Check the number of arguments. If none are passed, print help and exit.
if [ $# -eq 0 ]; then
    echo \\n"ERROR: No arguments given."\\n
    HELP
fi

### Start getopts code ###
while getopts :f:r:b:h FLAG; do
    case $FLAG in
        f) FREADS=$OPTARG;;
        r) RREADS=$OPTARG;;
        b) BREADS=$OPTARG;;
        h) HELP;;
        '?') # Unrecognized option
            echo \\n"ERROR: Option -$OPTARG not recognized. See help (-h)."\\n
            exit 2
        ;;
        :)  # No argument given for an option (f,r,b)
            echo "Option -$OPTARG requires an argument. See help (-h)."\\n
            exit 2
        ;;
    esac
done
shift $((OPTIND-1)) # This allows us to parse any remaining positional arguments like normal

### End getopts code ###
if [ $# -lt 2 ]; then
    echo \\n"ERROR: No GenomeID:IndexID pairs given."\\n
    HELP
fi

if [ $(( $# % 2 )) != 0 ]; then
    echo \\n"ERROR: Number of GenomeID:IndexID arguments is not even."\\n
    HELP
fi


# Check that the required options are given and that the files exist.
if [ ! $FREADS ]; then
    echo "ERROR: No Read 1 file was given."\\n
    HELP
elif [ ! -f $FREADS ]; then
    echo "ERROR: Read 1 file could not be found, please check that your filepath is correct."\\n
    exit 1
fi
if [ ! $RREADS ]; then
    echo "ERROR: No Read 2 file was given."\\n
    HELP
elif [ ! -f $RREADS ]; then
    echo "ERROR: Read 2 file could not be found, please check that your filepath is correct."\\n
    exit 1
fi
if [ ! $BREADS ]; then
    echo "ERROR: No Index sequence file was given."\\n
    HELP
elif [ ! -f $BREADS ]; then
    echo "ERROR: Index sequenc file could not be found, please check that your filepath is correct."\\n
    exit 1
fi

# Initialize timestamp, log and start time
DATE=`date +%Y-%m-%d`
TIME=`date +%H:%M:%S`
TM=`date +%Y%m%d-%H%M`
LOG=Demultiplex_reads_multi_log_$TM.txt
echo "Executing Demultiplex_reads_multi.sh on $DATE at $TIME "\\n | tee $LOG
START=`date +%s`
echo "Using $FREADS as the Read1 file.\\nUsing $RREADS as the Read2 file.\\nUsing $BREADS as the Index sequence file." | tee -a $LOG

## Demultiplexing code block, should never be necessary to change.
while [ $1 ]
    do
    echo '' | tee -a $LOG
    echo "Processing sample $1 using index $2" | tee -a $LOG
    demultiplex_reads.py -f $FREADS -r $RREADS -b $BREADS -o $1 -n $1 -i $2
    echo "Demultiplexing of sample $1 done!" | tee -a $LOG
    shift 2
done

END=`date +%s`
RUNTIME=$(( END - START ))
echo \\n'Demultiplexing has finished!' | tee -a $LOG
echo "This job took $RUNTIME seconds to complete" | tee -a $LOG



