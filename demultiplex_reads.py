#!/usr/bin/env python2.7
__author__ = 'Michael C. Nelson'
__copyright__ = "Copyright 2015, Michael C. Nelson/University of Connecticut"
__license__ = 'GPL3'
__version__ = '1'
__creationdate__ = '2015-03-23'
__lastmoddate__ = '2015-03-36'
__credits__ = ["Created using code lifted from the QIIME project (www.qiime.org)."]

"""
Demultiplex paired reads for sample MyGenome having index N701 for one MiSeq run and write results to directory GenomeA/,
demultiplex_reads.py -f Undetermined_read1.fastq.gz -r Undetermined_read2.fastq.gz -b Undetermined_I1.fastq.gz -o GenomeA/ -i N701 -n MyGenome
"""
import argparse
import gzip
import sys
from os import rename, makedirs
from os.path import isdir, exists
from itertools import izip

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = MyParser()
parser.add_argument('-f', '--fwd_reads', required=True, help="The Read 1 fastq file to search through. [REQUIRED]", metavar='Undetermined_R1.fastq.gz')
parser.add_argument('-r', '--rev_reads', required=True, help="The Read 2 fastq file to search through. [REQUIRED]", metavar='Undetermined_R2.fastq.gz')
parser.add_argument('-b', '--index_reads', required=True, help="The index fastq file to search against. [REQUIRED]", metavar='Undetermined_I1.fastq.gz')
parser.add_argument('-o', '--output_dir', required=True, help="The output directory to use, creating if not already present. [REQUIRED]", metavar='MyGenome_reads')
parser.add_argument('-i', '--index_name', required=True, help="The Illumina name of the index sequence. [REQUIRED]", metavar='N701')
parser.add_argument('-n', '--sample_id', required=True, help="Sample ID that should be used to name the output files. [REQUIRED]", metavar='MyGenome')
parser.add_argument('-c', '--custom_index', required=False, help="The sequence of a custom index primer not that is not part of the standard Illumina set.", metavar='ATCGATCG')
parser.add_argument('-l', '--list_indices', required=False, help="Print the built-in dictionary of Illumina indices and exit.")

# Defining a dict of the Illumina index sequences
Indices = {'N701': 'TAAGGCGA', 'N702': 'CGTACTAG', 'N703': 'AGGCAGAA', 'N704': 'TCCTGAGC', 'N705': 'GGACTCCT', 'N706': 'TAGGCATG', 'N707': 'CTCTCTAC', 'N708': 'CAGAGAGG', 'N709': 'GCTACGCT', 'N710': 'CGAGGCTG', 'N711': 'AAGAGGCA', 'N712': 'GTAGAGGA', 'AD001': 'ATCACG', 'AD002': 'CGATGT', 'AD003': 'TTAGGC', 'AD004': 'TGACCA', 'AD005': 'ACAGTG', 'AD006': 'GCCAAT', 'AD007': 'CAGATC', 'AD008': 'ACTTGA', 'AD009': 'GATCAG', 'AD010': 'TAGCTT', 'AD011': 'GGCTAC', 'AD012': 'CTTGTA', 'AD013': 'AGTCAA', 'AD014': 'AGTTCC'}

class ParseError(Exception):
    pass

def create_dir(dir_name, fail_on_exist=True):
    if exists(dir_name):
        if isdir(dir_name):
            if fail_on_exist:
                ror = "Directory already exists: %s" % dir_name
                return ror
        else:
            ror = "File with same name exists: %s" % dir_name
            return ror
    else:
        try:
            makedirs(dir_name)
        except OSError:
            ror = "Could not create output directory: %s. " % dir_name + "Check the permissions."
            return ror

    return 0

def check_map(infile):
    try:
        lines = open(infile,'U')
    except IOError:
        raise ParseError, ("A string was passed that doesn't refer to an accessible filepath.")

    # creating a throwaway function to strip dbl quotes and spaces from input lines
    strip_f = lambda x: x.replace('"','').strip()
    mapping_data = []
    headers = []

    for line in lines:                                             # Begin iterating over lines
        line = strip_f(line)                                       # Strip out witespace and dbl quotes
        if line.startswith('#'):                                   # If the input line begins with a hash,
            line = line[1:]                                        # The the line data is shifted by one character to strip the hash,
            headers = line.strip().split('\t')                     # Then strip whitespaces and split on tabs, putting line contents into the header list
        else:                                                      # Else, if the line doesn't begin with a hash, it's a sample line
            tmp_line = map(strip_f, line.split('\t'))              # Put the contents of the line into a temp list
            if len(tmp_line)<len(headers):                         # If the length of the temp list is less than the headder list
                tmp_line.extend(['']*(len(headers)-len(tmp_line))) # Pad out the temp list to equal the header list length
            mapping_data.append(tmp_line)                          # Now add the data from the tmp list to the mapping_data list

    id_map = {}                                                     # Create an empty dictionary instance called id_map
    for curr_data in mapping_data:                                  # For eash list in mapping_data
        id_map[curr_data[0]] = {}                                   # Assign the first value of the list to be a key with an empty dict associated
    for value in range(len(headers)):                               # For each list in the hds list
        for curr_data in mapping_data:                              # For each list in the mapping data
            id_map[curr_data[0]][headers[value]] = curr_data[value] #

    barcode_to_sample_id = {}                                       # Create an empty dictionary instance
    for sample_id, sample in id_map.items():                        # For each key:value pair in the id_map
        barcode_to_sample_id[sample['BarcodeSequence'].upper()] = sample_id

    return barcode_to_sample_id

def check_seq_headers(headerr1, headerr2, headeri1):
    """ Checks that the sequence headers are all equal to each other """
    header1 = headerr1.split(':')
    header2 = headerr2.split(':')
    headeri = headeri1.split(':')
    for e1,e2 in zip(header1,header2):
        if e1.split(' ')[0] != e2.split(' ')[0]:
            return False
    for l1,l2 in zip(header1,headeri):
        if l1.split(' ')[0] != l2.split(' ')[0]:
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
            yield record[0][1:], record[1], record[3] # yield the sequence header without the @ symbol, the DNA seqeunce, and quality sequence as strings
            line_num = 0
            record = []
        record.append(line.strip())
    if record:
        if record[0]:  # could be just an empty line at eof
            yield record[0][1:], record[1], record[3]
    if type(data) == file:
        data.close()

def format_log(input_sequence_count, sample_seq_count):
    """ Format the split libraries log """
    log_out = ["Demultiplexing results"]
    log_out.append("Total number of input sequences: %d" % input_sequence_count)
    log_out.append("Number of seqs for sample: %d" % sample_seq_count)
    return '\n'.join(log_out)

def parse_paired_reads(fastq_read1_f, fastq_read2_f, fastq_barcode_f, index, log_f):
    """parses out paired sequence reads according to given index sequence """
    # Define the index positions of the header, DNA sequence, and quality sequence of the fastqparer result.
    header_index = 0
    sequence_index = 1
    quality_index = 2
    sample_seq_count = 0

    barcode_length = len(index) # Setting length of barcode form the dict of index sequences

    # prep data for logging
    input_sequence_count = 0    # Keep track of how many sequences we're reading in
    for bc_data, read1_data, read2_data in izip(fastqparser(fastq_barcode_f), fastqparser(fastq_read1_f), fastqparser(fastq_read2_f)):
        # The fastqparser yields three things as a list: the header w/o the @, the DNA sequence, and quality sequence
        input_sequence_count += 1
        # Confirm match between barcode and read headers
        if not check_seq_headers(bc_data[header_index],read1_data[header_index], read2_data[header_index]):
            raise ParseError, ("Headers of barcode and read do not match. Can't continue. Confirm that the barcode fastq and read fastq that you are passing match one another.")
        else:
            header1 = read1_data[header_index]
            header2 = read2_data[header_index]                 # Set the read2 sequence header
            barcode = bc_data[sequence_index][:barcode_length] # Set the barcode to be the 1st n bases of the barcode string, with n = barcode length
            sequence1 = read1_data[sequence_index]             # Grab the read1 sequence string
            sequence2 = read2_data[sequence_index]             # Grad the read2 sequence string
            quality1 = read1_data[quality_index]               # Grab the read1 quality string
            quality2 = read2_data[quality_index]               # Grab the read2 quality string
            if barcode == index:
                # Returns a generator of the sequence data.
                yield header1, sequence1, quality1, header2, sequence2, quality2
                sample_seq_count +=1

    if log_f != None:
        log_str = format_log(input_sequence_count, sample_seq_count)
        log_f.write(log_str)

def main():
    opts = parser.parse_args()
    fwd_reads_fp = opts.fwd_reads                                           # Set the input R1 filepath
    rev_reads_fp = opts.rev_reads                                           # Set the input R2 filepath
    index_read_fp = opts.index_reads                                        # Set the index read filepath
    barcode = Indices[opts.index_name]                                      # Lookup what the barcode to be used for demultiplexing is
    sample = opts.sample_id                                                 # Set what the sample name is
    output_dir = opts.output_dir                                            # Set what the output directory will be called
    create_dir(output_dir)                                                  # Create the output directory safely
    output_r1_fp_temp = "%s/%s_R1.fastq.incomplete" % (output_dir,sample)   # Set the name and filepath of the temporary output R1 file
    output_r2_fp_temp = "%s/%s_R2.fastq.incomplete" % (output_dir,sample)   # Set the name and filepath of the temporary output R2 file
    output_r1_fp = '%s/%s_R1.fastq.gz' % (output_dir, sample)               # Set the name and filepath of the actual output R1 file
    output_r2_fp = '%s/%s_R2.fastq.gz' % (output_dir, sample)               # Set the name and filepath of the actual output R2 file
    output_r1_f = gzip.open(output_r1_fp_temp,'w')                          # Going to automatically write to gzip compressed file.
    output_r2_f = gzip.open(output_r2_fp_temp,'w')                          # Going to automatically write to gzip compressed file.
    log_fp = '%s/demultiplex_log.txt' % output_dir                          # Set the name and filepath of the log file
    log_f = open(log_fp,'w')                                                # Writing the log file as plain text (no compression).


    # Check file extensions of read and index files to determine which method to use for opening file for processing
    if fwd_reads_fp.endswith('.gz'):            # If the input read file is gzipped ## sequence_read_fp == fwd_reads_fp
        fwd_read_f = gzip.open(fwd_reads_fp)    # Then gzopen it for reading
    else:                                       # Else
        fwd_read_f = open(fwd_reads_fp,'U')     # Open as a standard text file
    if rev_reads_fp.endswith('.gz'):
        rev_read_f = gzip.open(rev_reads_fp)
    else:
        rev_read_f = open(rev_reads_fp,'U')
    if index_read_fp.endswith('.gz'):
        barcode_read_f = gzip.open(index_read_fp)
    else:
        barcode_read_f = open(index_read_fp,'U')

    # Parse reads and return to a generator expression to reduce memory footprint:
    seq_generator = parse_paired_reads(fwd_read_f, rev_read_f, barcode_read_f, barcode, log_f)

    # Write the parsed data to both output files
    for header_f, sequence_f, quality_f, header_r, sequence_r, quality_r in seq_generator:
        output_r1_f.write('@%s\n%s\n+\n%s\n' % (header_f, sequence_f, quality_f)) # This will write out fastq output
        output_r2_f.write('@%s\n%s\n+\n%s\n' % (header_r, sequence_r, quality_r)) # This will write out fastq output

    log_f.write('\n---\n\n')

    output_r1_f.close()
    rename(output_r1_fp_temp, output_r1_fp)
    output_r2_f.close()
    rename(output_r2_fp_temp, output_r2_fp)

if __name__ == "__main__":
    main()
