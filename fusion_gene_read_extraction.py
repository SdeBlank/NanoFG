#!/usr/bin/python

import argparse
import vcf as pyvcf
import pysam
import sys

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('-b', '--bam', required=True, type=str, help='Input bam file')
parser.add_argument('-v', '--vcf', required=True, type=str, help='Input NanoSV vcf file')

args = parser.parse_args()

def create_fasta( chr, start, end, svid, exclude):
    print( chr, start, end )
    bamfile = pysam.AlignmentFile(args.bam, "rb" )
    fasta = open(svid+".fasta", 'a+')
    for read in bamfile.fetch(chr, start, end):
        if read.query_name in exclude or read.seq == None:
            continue
        fasta.write( ">"+read.query_name+"\n")
        fasta.write(read.seq+"\n")
        exclude.append(read.query_name)
    fasta.close()

vcf_reader = pyvcf.Reader(open(args.vcf, 'r'))
for record in vcf_reader:
    create_fasta(record.CHROM, record.POS-40000, record.POS, record.ID, record.INFO['REF_READ_IDS_1'])
    create_fasta(record.ALT[0].chr, record.ALT[0].pos, record.ALT[0].pos+1000, record.ID, record.INFO['REF_READ_IDS_2'])
    break
