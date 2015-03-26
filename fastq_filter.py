#!/usr/bin/env python

__author__ = "Michael C. Nelson"
__copyright__ = "Copyright 2015, Michael C. Nelson/University of Connecticut"
__credits__ = ["Created using code lifted from the QIIME project (www.qiime.org)."]
__license__ = "GPL3"
__version__ = "1.0"

import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("-i", '--input', required=True, help="The input fastq file to search through. [REQUIRED]", metavar='Input.fastq')
parser.add_argument('-o', '--output', required=True, help="The output fastq to create. [REQUIRED]", metavar='Output.fastq.gz')
parser.add_argument('-n', '--names', required=True, help="The file containing the read IDs to search for. [REQUIRED]", metavar='Names.txt')
args = parser.parse_args()

INFILE = args.input
OUTFILE = args.output
if OUTFILE.endswith('.fastq'):
   print "Warning: Output files are automatically gzip compressed, adding correct extension to output file name. "
   OUTFILE += '.gz'
NAMES = args.names

def get_ids(names_file):
    return set([l.split()[0].strip() for l in names_file if not l.startswith('#') and l])

def minimalfastqparser(data):
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


def filter_fastq(input, output, seqs_to_keep):
    input_seqs = minimalfastqparser(input)
    output_f = gzip.open(output, 'w')
    seqs_to_keep_lookup = {}.fromkeys([seq_id.split()[0] for seq_id in seqs_to_keep])

    for seq_id, seq, qual in input_seqs:
        if seq_id.split()[0] in seqs_to_keep_lookup:
            output_f.write('@%s\n%s\n+\n%s\n' % (seq_id, seq, qual))
    output_f.close()


def main():
    keep_IDs = get_ids(open(NAMES, 'U')) # Creates a dictionary of seqIDs for seqs we want to keep

    if INFILE.endswith('.fastq'):
        filter_fastq(INFILE, OUTFILE, keep_IDs)
    elif INFILE.endswith('.fastq.gz'):
        filter_fastq(INFILE, OUTFILE, keep_IDs)
    else:
        print "ERROR: Input file must be a non-compressed fastq file with the extension .fastq"


if __name__ == "__main__":
    main()
