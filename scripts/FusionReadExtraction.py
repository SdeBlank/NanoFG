#!/usr/bin/python

import argparse
import vcf as pyvcf
import pysam
import sys
import datetime
from EnsemblRestClient import EnsemblRestClient

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('-b', '--bam', required=True, type=str, help='Input bam file')
parser.add_argument('-v', '--vcf', required=True, type=str, help='Input NanoSV vcf file')
parser.add_argument('-o', '--output_dir', required=True, type=str, help='Output directory for fasta files')

args = parser.parse_args()

def get_gene_overlap( chr, pos, ori, bp ):
    SERVER='http://grch37.rest.ensembl.org'
    SPECIES="human"
    ENDPOINT="/overlap/region/"+SPECIES+"/"+str(chr)+":"+str(pos)+"-"+str(pos)
    HEADERS={"Content-Type" : "application/json"}
    PARAMS={"feature": "gene"}

    genes_data=EnsemblRestClient.perform_rest_action(SERVER, ENDPOINT, HEADERS, PARAMS)

    fusions = dict()

    for gene in genes_data:
        if gene["biotype"]=="protein_coding":
            fusion=[]
            if not ori and gene['strand'] == 1:
                if 'donor' not in fusions:
                    fusions['donor'] = dict()
                fusions['donor'][gene['id']+"\t"+bp] = gene['start']
            elif not ori and gene['strand'] == -1:
                if 'acceptor' not in fusions:
                    fusions['acceptor'] = dict()
                fusions['acceptor'][gene['id']+"\t"+bp] = gene['start']
            elif ori and gene['strand'] == 1:
                if 'acceptor' not in fusions:
                    fusions['acceptor'] = dict()
                fusions['acceptor'][gene['id']+"\t"+bp] = gene['end']
            elif ori and gene['strand'] == -1:
                if 'donor' not in fusions:
                    fusions['donor'] = dict()
                fusions['donor'][gene['id']+"\t"+bp] = gene['end']

    return( fusions )

def create_fasta( chr, start, end, svid, exclude ):
    if end < start:
        end, start = start, end
    bamfile = pysam.AlignmentFile(args.bam, "rb" )
    fasta = open(args.output_dir+"/"+svid+".fasta", 'a+')
    for read in bamfile.fetch(chr, start, end):
        #### Breakpoints that only have supplementary reads and not a primary read spanning the breakpoint will be excluded
        #### No effect when testing on the truthset in recall, but see if it ever happens in real sets
        if read.query_name in exclude or read.seq == None or read.is_supplementary:
            continue
        fasta.write( ">"+read.query_name+"\n")
        fasta.write(read.seq+"\n")
        #exclude.append(read.query_name)
    fasta.close()
    bamfile.close()

print("Start:", datetime.datetime.now())

EnsemblRestClient=EnsemblRestClient()
vcf_reader = pyvcf.Reader(open(args.vcf, 'r'))
for record in vcf_reader:
    if not isinstance(record.ALT[0], pyvcf.model._Breakend):
        #Possibility to add converter for different SV callers
        continue
    fusions={'donor':{}, 'acceptor':{}}
    bnd1_fusions = get_gene_overlap(record.CHROM, record.POS, record.ALT[0].orientation, '1')

    ### Skip next request if the first BND already falls outside of a gene
    if not fusions:
        continue
    bnd2_fusions=get_gene_overlap(record.ALT[0].chr, record.ALT[0].pos, record.ALT[0].remoteOrientation, '2')

    fusions.update(bnd1_fusions)
    if 'donor' in bnd2_fusions:
        fusions['donor'].update(bnd2_fusions['donor'])
    if 'acceptor' in bnd2_fusions:
        fusions['acceptor'].update(bnd2_fusions['acceptor'])

    good_fusion=False
    if 'donor' in fusions and 'acceptor' in fusions:
        largest_donor_size=0
        largest_acceptor_size=0
        REF_READ_IDS=record.INFO['REF_READ_IDS_1']+record.INFO['REF_READ_IDS_2']
        for donor in fusions['donor']:
            donor_gene, donor_bp = donor.split("\t")
            for acceptor in fusions['acceptor']:
                acceptor_gene, acceptor_bp = acceptor.split("\t")
                fusion = donor_gene+"_"+acceptor_gene
                if donor_gene != acceptor_gene:
                    if donor_bp == '1' and acceptor_bp == '2':
                        good_fusion=True
                        donor_size=abs(record.POS-fusions['donor'][donor])
                        if donor_size>largest_donor_size:
                            largest_donor_size=donor_size
                            donor_chr = record.CHROM
                            donor_start = fusions['donor'][donor]
                            donor_end = record.POS

                        acceptor_size=abs(record.ALT[0].pos-fusions['acceptor'][acceptor])
                        if acceptor_size>largest_acceptor_size:
                            largest_acceptor_size=acceptor_size
                            acceptor_chr = record.ALT[0].chr
                            acceptor_start = record.ALT[0].pos
                            acceptor_end = fusions['acceptor'][acceptor]

                    elif donor_bp == '2' and acceptor_bp == '1':
                        good_fusion=True
                        donor_size=abs(record.ALT[0].pos-fusions['donor'][donor])
                        if donor_size>largest_donor_size:
                            largest_donor_size=donor_size
                            donor_chr = record.ALT[0].chr
                            donor_start = fusions['donor'][donor]
                            donor_end = record.ALT[0].pos

                        acceptor_size=abs(record.POS-fusions['acceptor'][acceptor])

                        if acceptor_size>largest_acceptor_size:
                            largest_acceptor_size=acceptor_size
                            acceptor_chr = record.CHROM
                            acceptor_start = fusions['acceptor'][acceptor]
                            acceptor_end = record.POS
    if good_fusion:
        create_fasta(donor_chr, donor_start, donor_end, record.ID, REF_READ_IDS)
        create_fasta(acceptor_chr, acceptor_start, acceptor_end, record.ID, REF_READ_IDS)

print("End:", datetime.datetime.now())
    #create_fasta(record.ALT[0].chr, record.ALT[0].pos, record.ALT[0].pos+1000, record.ID, record.INFO['REF_READ_IDS_2'])
    #break
