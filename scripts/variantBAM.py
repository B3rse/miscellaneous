#!/usr/bin/env python3

################################################
#
#   Insert SNV at position in bam file
#
#   Michele Berselli
#   berselli.michele@gmail.com
#
################################################

################################################
#   Libraries
################################################
import sys, argparse, os
import pysam as ps
import random
import re

################################################
#   Functions
################################################
def md_list(md_, pos_, read_):
    ''' '''
    md_list_, is_new = [], False
    for md_ in re.sub(r'([\\^]*[ACGT]+)[0]*', ' \\1 ', md_).split():
        try:
            int_md_ = int(md_)
            md_list_ += ['M'] * int_md_
        except: md_list_ += [md_]

    if md_list_[pos_] == 'M':
        md_list_[pos_] = read_.seq[pos_]
        is_new = True

    return md_list_, is_new

def md_from_list(md_list_):
    ''' '''
    md_str, c_M = '', 0
    for md_ in md_list_:
        if md_ == 'M':
            c_M += 1
        else:
            if c_M:
                md_str += str(c_M)
                c_M = 0
                md_str += md_
            else:
                if md_.startswith('^'):
                    md_str += md_
                else: md_str += '0' + md_
    if c_M:
        md_str += str(c_M)

    return md_str

def add_snv(read, pos, base):
    ''' '''
    a = ps.AlignedSegment()

    # Add snv to seq and set high base quality score
    qual_ = read.qual[:pos-read.reference_start] + 'A' + read.qual[pos-read.reference_start+1:]
    seq_ = read.seq[:pos-read.reference_start] + base + read.seq[pos-read.reference_start+1:]

    # Calculate tag MD
    MD = read.get_tag('MD')
    md_list_, is_new = md_list(MD, pos-read.reference_start, read)
    MD_ = md_from_list(md_list_)

    # Calculate tag NM
    NM_ = read.get_tag('NM')
    if is_new: NM_ += 1

    # Create modified read
    a.query_name = read.query_name
    a.query_sequence = seq_
    a.flag = read.flag
    a.reference_id = read.reference_id
    a.reference_start = read.reference_start
    a.mapping_quality = read.mapping_quality
    a.cigar = read.cigar
    a.next_reference_id = read.next_reference_id
    a.next_reference_start = read.next_reference_start
    a.template_length = read.template_length
    a.query_qualities = ps.qualitystring_to_array(qual_)
    a.tags = read.get_tags()
    a.set_tag('MD', MD_, 'Z')
    a.set_tag('NM', NM_, 'i')

    return a

def rnd_gen():
    '''return 0 or 1 with 50% chance for each'''
    return random.randrange(2)

def main(args):

    # Variables
    threads = int(args['threads']) if args['threads'] else 1
    is_snv = False

    # Check SNV to insert
    if args['snv']:
        is_snv = True
        chr_snv, pos_snv, base_snv, gt_snv = args['snv'].split(':')
        pos_snv = int(pos_snv)
    else:
        sys.exit('specify SNV to add')

    # Open bamfile using pysam
    samfile = ps.AlignmentFile(args['inputfile'], "rb", threads=threads)

    # Modify and write reads to stdout
    with ps.AlignmentFile('{0}_{1}.bam'.format(args['inputfile'].split('.')[0], args['snv'].replace(':', '_')), "wb", header=samfile.header, threads=threads) as outf:
        for i, read in enumerate(samfile):
            sys.stderr.write('\r writing read #: ' + str(i+1))
            a = read
            # Check if read to modify
            if chr_snv == read.reference_name:
                # indexes in read are 0-based
                if not a.is_unmapped:
                    if pos_snv-1 >= read.reference_start and pos_snv-1 <= read.reference_end:
                        # hom snv
                        if gt_snv == 'HOM':
                            a = add_snv(read, pos_snv-1, base_snv)
                        # het snv
                        elif gt_snv == 'HET':
                            if rnd_gen(): a = add_snv(read, pos_snv-1, base_snv)

            outf.write(a)
    sys.stderr.write('\n')

################################################
#   MAIN
################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Add read groups to a BAM file')

    parser.add_argument('-i','--inputfile', help='input aligned paired-ends BAM file', required=True)
    parser.add_argument('-t','--threads', help='number of threads to use for compression/decompression', required=False)
    parser.add_argument('--snv', help='SNV to insert in the format CHR:POS:BASE:<HET|HOM>', required=False)

    args = vars(parser.parse_args())

    main(args)
