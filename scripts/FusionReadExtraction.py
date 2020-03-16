#!/usr/bin/python

import argparse
import vcf as pyvcf
import pysam
import sys
import datetime
import copy
from pybiomart import Dataset
from EnsemblRestClient import EnsemblRestClient

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('-b', '--bam', required=True, type=str, help='Input bam file')
parser.add_argument('-v', '--vcf', required=True, type=str, help='Input NanoSV vcf file')
parser.add_argument('-o', '--output_dir', required=True, type=str, help='Output directory for fasta files')

args = parser.parse_args()

#############################################   Convert unknown ALT fields to bracket notations N[Chr:pos[   #############################################
def alt_convert( record ):
    orientation = None
    remoteOrientation = None
    if 'INS' in record.INFO['SVTYPE']:
        return( record )
    elif record.INFO['SVTYPE'] == 'DEL':
        CHR2=record.CHROM
        orientation = False
        remoteOrientation = True
    elif record.INFO['SVTYPE'] == 'DUP':
        CHR2=record.CHROM
        orientation = True
        remoteOrientation = False
    elif record.INFO['SVTYPE'] == 'INV':
        CHR2=record.CHROM
        strands=record.INFO['STRANDS'][0]
        if strands == "++":
            orientation = False
            remoteOrientation = False
        elif strands == "--":
            orientation = True
            remoteOrientation = True
    elif record.INFO['SVTYPE'] == 'TRA':
        CHR2=record.INFO['CHR2']
        strands=record.INFO['STRANDS'][0]
        if strands == "++":
            orientation = False
            remoteOrientation = False
        elif strands == "+-":
            orientation = False
            remoteOrientation = True
        elif strands == "-+":
            orientation = True
            remoteOrientation = False
        elif strands == "--":
            orientation = True
            remoteOrientation = True
    elif record.INFO['SVTYPE'] == 'INVDUP':
        CHR2=record.INFO['CHR2']
        strands=record.INFO['STRANDS'][0]
        if strands == "++":
            orientation = False
            remoteOrientation = False
        elif strands == "+-":
            orientation = False
            remoteOrientation = True
        elif strands == "-+":
            orientation = True
            remoteOrientation = False
        elif strands == "--":
            orientation = True
            remoteOrientation = True
    else:
        CHR2=record.INFO['CHR2']
        strands=record.INFO['STRANDS'][0]
        if strands == "++":
            orientation = False
            remoteOrientation = False
        elif strands == "+-":
            orientation = False
            remoteOrientation = True
        elif strands == "-+":
            orientation = True
            remoteOrientation = False
        elif strands == "--":
            orientation = True
            remoteOrientation = True
    if orientation is None or remoteOrientation is None:
        sys.exit("Error in alt_convert; Unknown ALT field")
    record.ALT = [ pyvcf.model._Breakend( CHR2, record.INFO['END'], orientation, remoteOrientation, record.REF, True ) ]
    return( record )

#############################################   Get the overlap of a breakpoint with genes from ensembl and extract information of the gene   #############################################
def get_gene_overlap( bnd1_chr, bnd1_pos, bnd1_ori, bnd2_chr, bnd2_pos, bnd2_ori, ensembl_gene_regions):
    overlap=[]
    for region in ensembl_gene_regions:
        if bnd1_chr==region['chromosome'] and bnd1_pos>region['start'] and bnd1_pos<region['end']:
            region["bnd"]="1"
            overlap.append(copy.deepcopy(region))
        if bnd2_chr==region['chromosome'] and bnd2_pos>region['start'] and bnd2_pos<region['end']:
            region["bnd"]="2"
            overlap.append(copy.deepcopy(region))

    fusions = dict()
    for gene in overlap:
        # if gene["biotype"]=="protein_coding":
        if gene['bnd']=="1":
            ori=bnd1_ori
            if not ori:
                if 'donor' not in fusions:
                    fusions['donor'] = dict()
                fusions['donor'][gene['id']+"\t"+gene['bnd']] = gene['start']
            elif ori:
                if 'donor' not in fusions:
                    fusions['donor'] = dict()
                fusions['donor'][gene['id']+"\t"+gene['bnd']] = gene['end']
        elif gene['bnd']=="2":
            ori=bnd2_ori
            if not ori:
                if 'acceptor' not in fusions:
                    fusions['acceptor'] = dict()
                fusions['acceptor'][gene['id']+"\t"+gene['bnd']] = gene['start']
            elif ori:
                if 'acceptor' not in fusions:
                    fusions['acceptor'] = dict()
                fusions['acceptor'][gene['id']+"\t"+gene['bnd']] = gene['end']
    return( fusions )

########################################   Obtain all reads from the bam file that support a breakpoint, excluding reference or non-overlapping reads   ########################################
def create_fasta( chr, start, end, svid, exclude, include ):
    if end < start:
        end, start = start, end
    bamfile = pysam.AlignmentFile(args.bam, "rb" )
    fasta = open(args.output_dir+"/"+svid+".fasta", 'a+')

    for read in bamfile.fetch(chr, start, end):
        '''
        Breakpoints that only have supplementary reads and not a primary read spanning the breakpoint will be excluded
        No effect when testing on the truthset in recall, but see if it ever happens in real sets
        '''
        if read.query_name in include and not read.seq == None and not read.is_supplementary and read.query_name not in exclude:
            fasta.write( ">"+svid+"."+read.query_name+"\n")
            fasta.write(read.seq+"\n")
            exclude.append(read.query_name)
        #### Uncomment to select all the reads supporting the breakpoint and adjacent regions of the gene to get full gene sequence back
        # if read.query_name in exclude or read.seq == None or read.is_supplementary:
        #     continue
        # fasta.write( ">"+svid+"."+read.query_name+"\n")
        # fasta.write(read.seq+"\n")

    fasta.close()
    bamfile.close()




########################################   Main code   ########################################

print("Start:", datetime.datetime.now())

EnsemblRestClient=EnsemblRestClient()

### DOWNLOAD BASIC GENE INFORMATION FROM ENSEMBL (ID, CHROMOSOME, POSITION, STRAND, BIOTYPE)
dataset = Dataset(name='hsapiens_gene_ensembl', host='http://grch37.ensembl.org')
ensembl_genes = dataset.query(attributes=['ensembl_gene_id', 'chromosome_name','start_position', 'end_position', 'strand', 'gene_biotype'],
                filters={'chromosome_name': ['1','2','3','4','5','6','7','8','9','10','11','12','13','14',
                                            '15','16','17','18','19','20','21','22','X','Y','MT']})
regions=[]

for line in ensembl_genes.iterrows():
    index, data = line
    columns=data.tolist()
    regions.append({'id':columns[0], 'chromosome':columns[1], 'start':int(columns[2]), 'end':int(columns[3]), 'strand':int(columns[4]), 'biotype':columns[5]})

vcf_reader = pyvcf.Reader(open(args.vcf, 'r'))

### DETERMINE IF VCF IS PRODUCED BY SNIFFLES OR NANOSV
if "source" in vcf_reader.metadata:
    if vcf_reader.metadata["source"][0].lower()=="sniffles":
        vcf_type="Sniffles"
elif "cmdline" in vcf_reader.metadata:
    if "nanosv" in vcf_reader.metadata["cmdline"][0].lower():
        vcf_type="NanoSV"

### GO THROUGH EVERY SV IN THE VCF AND DETERMINE IF THE SV PRODUCES A VALID FUSION BETWEEN 2 GENES
for record in vcf_reader:
    #CONVERT THE ALT FIELD IN THE VCF TO A BND 'N]]' structure'
    if not isinstance(record.ALT[0], pyvcf.model._Breakend):
        record = alt_convert(record)
    if not isinstance(record.ALT[0], pyvcf.model._Breakend):
        continue

    fusions={'donor':{}, 'acceptor':{}}

    #SNIFFLES DOES NOT SHOW A CORRECT BND STRUCTURE FOR ALL BREAKPOINTS. FOR THAT REASON, THE STRANDS VALUE IN THE INFO FIELD IS USED TO PRODUCE A CORRECT BND STRUCTURE
    if vcf_type=="Sniffles":
        if record.INFO["STRANDS"][0][0]=="+":
            strand1=False
        else:
            strand1=True
        if record.INFO["STRANDS"][0][1]=="+":
            strand2=False
        else:
            strand2=True

    if vcf_type=="NanoSV":
        strand1=record.ALT[0].orientation
        strand2=record.ALT[0].remoteOrientation

    fusions=get_gene_overlap(record.CHROM, record.POS, strand1, record.ALT[0].chr, record.ALT[0].pos, strand2, regions)

    #CHECK IF A BREAKPOINT CONTAIN THE CORRECT ORIENTATION TO PRODUCE A FUSION GENE
    good_fusion=False
    if 'donor' in fusions and 'acceptor' in fusions:
        largest_donor_size=0
        largest_acceptor_size=0
        if "ALT_READ_IDS" in record.INFO:
            SUPP_READ_IDS=record.INFO['ALT_READ_IDS']
            REF_READ_IDS=record.INFO['REF_READ_IDS_1']+record.INFO['REF_READ_IDS_2']
        elif "RNAMES" in record.INFO:
            SUPP_READ_IDS=record.INFO['RNAMES']
            REF_READ_IDS=[]
        # GO THROUGH ALL POSSIBLE FUSIONS FOR EACH BREAKPOINT TO CHECK THE ORIENTATION
        for donor in fusions['donor']:
            donor_gene, donor_bp = donor.split("\t")
            for acceptor in fusions['acceptor']:
                acceptor_gene, acceptor_bp = acceptor.split("\t")
                fusion = donor_gene+"_"+acceptor_gene
                if donor_gene != acceptor_gene:
                    good_fusion=True
                    if donor_bp == '1' and acceptor_bp == '2':
                    #     good_fusion=True
                    # donor_start = fusions['donor'][donor]
                    # if donor_strand=="1":
                    #     donor_end = record.POS

                    ### OPTIONAL CODE TO SELECT ALL THE READS FOR A WHOLE GENE (CURRENTLY NOT USED)
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
                        # good_fusion=True

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
        #EXTRACT ALL THE READS THAT SUPPORT THE BREAKPOINT FOR LATER REMAPPING WITH MORE ACCURATE PARAMETERS
        create_fasta(donor_chr, donor_start, donor_end, record.ID, REF_READ_IDS, SUPP_READ_IDS)
        create_fasta(acceptor_chr, acceptor_start, acceptor_end, record.ID, REF_READ_IDS, SUPP_READ_IDS)

print("End:", datetime.datetime.now())
