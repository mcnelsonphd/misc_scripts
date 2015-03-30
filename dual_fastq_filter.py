#!/usr/bin/env python

__author__ = "Michael C. Nelson"
__created__ = "2015-03-29"
__copyright__ = "Copyright 2015, Michael C. Nelson/University of Connecticut"
__credits__ = ["Created using code lifted from the QIIME project (www.qiime.org)."]
__license__ = "GPL3"
__version__ = "1.0"
__revised__ = "2015-03-29"

import argparse
import gzip
from itertools import izip

class ParseError(Exception):
    pass

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--forward', required=True, help="The forward input fastq file to search through. [REQUIRED]", metavar='Input1.fastq')
parser.add_argument('-r', '--reverse', required=True, help="The reverse input fastq file to search through. [REQUIRED]", metavar='Input2.fastq')
parser.add_argument('-o', '--output1', required=True, help="The output fastq to create. [REQUIRED]", metavar='Output1.fastq.gz')
parser.add_argument('-p', '--output2', required=True, help="The output fastq to create. [REQUIRED]", metavar='Output2.fastq.gz')
parser.add_argument('-n', '--names', required=True, help="The file containing the read IDs to search for. [REQUIRED]", metavar='Names.txt')
args = parser.parse_args()

def check_seq_headers(headerr1, headerr2):
    """ Checks that the sequence headers are all equal to each other """
    header1 = headerr1.split(':')
    header2 = headerr2.split(':')
    for e1,e2 in zip(header1,header2):
        if e1.split(' ')[0] != e2.split(' ')[0]:
            return False
    return True


def fastqparser(data):
    if type(data) == str:
        data = open(data, 'rU')
    line_num = -1
    record = []
    for line in data:
        line_num += 1
        if line_num == 4:
            yield record[0][1:], record[1], record[3]
            line_num = 0
            record = []
        record.append(line.strip())
    if record:
        if record[0]:  # could be just an empty line at eof
            yield record[0][1:], record[1], record[3]
    if type(data) == file:
        data.close()


def filter_paired_fastq(fastq_read1_f, fastq_read2_f, seqs_to_keep):
    # seqs_to_keep_lookup = {}.fromkeys([seq_id.split()[0] for seq_id in seqs_to_keep])

    for read1_data, read2_data in izip(fastqparser(fastq_read1_f), fastqparser(fastq_read2_f)):
        if not check_seq_headers(read1_data[0], read2_data[0]):
            raise ParseError, "Headers of barcode and read do not match. Can't continue. Confirm that the barcode fastq and read fastq that you are passing match one another."
        else:
            header1 = read1_data[0]
            header2 = read2_data[0]
            seq1 = read1_data[1]
            seq2 = read2_data[1]
            qual1 = read1_data[2]
            qual2 = read2_data[2]
            if header1.split()[0] in seqs_to_keep:
                yield header1, seq1, qual1, header2, seq2,qual2


def main():
    fwd_reads_fp = args.forward
    rev_reads_fp = args.reverse
    output_r1 = args.output1
    output_r2 = args.output2
    output_r1_f = gzip.open(output_r1,'w')      # Going to automatically write to gzip compressed file.
    output_r2_f = gzip.open(output_r2,'w')      # Going to automatically write to gzip compressed file.

    if output_r1.endswith('.fastq'):
        print "Warning: Output files are automatically gzip compressed, adding correct extension to output file name."
        output_r1 += '.gz'
    if output_r2.endswith('.fastq'):
        print "Warning: Output files are automatically gzip compressed, adding correct extension to output file name."
        output_r2 += '.gz'
    NAMES = args.names
    keep_IDs = list([l.split()[0].strip() for l in open(NAMES, 'U')]) # Creates a list of seqIDs for seqs we want to keep

    # Check file extensions of read and index files to determine which method to use for opening file for processing
    if fwd_reads_fp.endswith('.gz'):            # If the input read file is gzipped ## sequence_read_fp == fwd_reads_fp
        fwd_read_f = gzip.open(fwd_reads_fp)    # Then gzopen it for reading
    else:                                       # Else
        fwd_read_f = open(fwd_reads_fp,'U')     # Open as a standard text file
    if rev_reads_fp.endswith('.gz'):
        rev_read_f = gzip.open(rev_reads_fp)
    else:
        rev_read_f = open(rev_reads_fp,'U')

    # Parse reads and return to a generator expression to reduce memory footprint:
    seq_generator = filter_paired_fastq(fwd_read_f, rev_read_f, keep_IDs)

    # Write the parsed data to both output files
    for header_f, sequence_f, quality_f, header_r, sequence_r, quality_r in seq_generator:
        output_r1_f.write('@%s\n%s\n+\n%s\n' % (header_f, sequence_f, quality_f)) # This will write out fastq output
        output_r2_f.write('@%s\n%s\n+\n%s\n' % (header_r, sequence_r, quality_r)) # This will write out fastq output

    output_r1_f.close()
    output_r2_f.close()

if __name__ == "__main__":
    main()
