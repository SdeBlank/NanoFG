import vcf as pyvcf
import argparse
import sys
import requests
import json
import time
import sys
from EnsemblRestClient import EnsemblRestClient

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('-v', '--vcf', required=True, type=str, help='Path to fusion gene vcf')
parser.add_argument('-d', '--directory', required=True, type=str, help='Path to primer directory')
parser.add_argument('-o', '--offset', default=0, type=int, help='Offset [default: 0]')
parser.add_argument('-f', '--flank', default=200, type=int, help='Flank [default: 200]')
parser.add_argument('-m', '--mask',action='store_true')
args = parser.parse_args()
vcf = args.vcf
directory = args.directory
offset = args.offset
flank = args.flank
mask = args.mask

#############################################   Mask common variation locations in the sequence  #############################################
def mask_seq( chr, start, end, strand, seq ):
    ext = "/overlap/region/human/"+str(chr)+":"+str(start)+"-"+str(end)+":"+str(strand)+"?feature=variation"
    headers={ "Content-Type" : "application/json"}
    data=EnsemblRestClient.perform_rest_action(server, ext, headers)
    masked_seq = seq

    for snp in data:
        snp_pos= snp['start']-start
        if str(strand) == '-1':
            snp_pos = len(masked_seq)-snp_pos-1
        masked_seq = masked_seq[:snp_pos]+'n'+masked_seq[snp_pos+1:];
    return( masked_seq )

#############################################   Gather the reference sequence of a specific area in the genome  #############################################
def get_seq(chr, start, end, strand):
    ext = "/info/assembly/homo_sapiens/"+str(chr)
    headers={ "Content-Type" : "application/json"}
    data=EnsemblRestClient.perform_rest_action(server, ext, headers)
    chrlength = data['length']

    ext = "/sequence/region/human/"+str(chr)+":"+str(start)+"-"+str(end)+":"+str(strand)+""
    if start < 1:
        sys.stderr.write("\tFailed to get seq, because POS is too close to START of the chr\n")
        seq = False
    elif end > chrlength:
        sys.stderr.write("\tFailed to get seq, because ENDPOS is too close to END of the chr\n")
        seq = False
    else:
        data=EnsemblRestClient.perform_rest_action(server, ext, headers)
        seq = data['seq']
        if mask:
            seq = mask_seq(chr, start, end, strand, seq )
    return( seq )

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

#############################################   RUNNING CODE   #############################################
EnsemblRestClient=EnsemblRestClient()
server = "http://grch37.rest.ensembl.org"

with open(vcf, "r") as fusion_vcf:
    vcf_reader = pyvcf.Reader(fusion_vcf)

    ### DETERMINE IF VCF IS PRODUCED BY SNIFFLES OR NANOSV
    if "source" in vcf_reader.metadata:
        if vcf_reader.metadata["source"][0].lower()=="sniffles":
            vcf_type="Sniffles"
    elif "cmdline" in vcf_reader.metadata:
        if "nanosv" in vcf_reader.metadata["cmdline"][0].lower():
            vcf_type="NanoSV"

    for record in vcf_reader:
        if "FUSION" in record.INFO:
            ### GATHER THE SEQUENCE AROUND BREAKPOINTS THAT PRODUCE A VALID FUSION GENE
            if vcf_type=="NanoSV":
                pos1_orientation=record.ALT[0].orientation
                pos2_orientation=record.ALT[0].remoteOrientation

            #SNIFFLES DOES NOT SHOW A CORRECT BND STRUCTURE FOR ALL BREAKPOINTS. FOR THAT REASON, THE STRANDS VALUE IN THE INFO FIELD IS USED TO PRODUCE A CORRECT BND STRUCTURE
            elif vcf_type=="Sniffles":
                if record.INFO["STRANDS"][0][0]=="+":
                    pos1_orientation=False
                else:
                    pos1_orientation=True
                if record.INFO["STRANDS"][0][1]=="+":
                    pos2_orientation=False
                else:
                    pos2_orientation=True

            seq1 = ''
            seq2 = ''
            if pos1_orientation:
                seq1 = get_seq(record.CHROM, record.POS+offset, record.POS+offset+flank, -1)
            else:
                seq1 = get_seq(record.CHROM, record.POS-flank-offset, record.POS-offset, 1)
            if not seq1:
                continue
            if pos2_orientation:
                seq2 = get_seq(record.ALT[0].chr, record.ALT[0].pos+offset, record.ALT[0].pos+offset+flank, 1)
            else:
                seq2 = get_seq(record.ALT[0].chr, record.ALT[0].pos-offset-flank, record.ALT[0].pos-offset, -1)
            if not seq2:
                continue

            with open(directory+"/"+record.ID+"."+record.INFO["FUSION"][0].replace("-","_")+".fasta", "w") as primer_output:
                primer_output.write(">"+record.ID+"."+record.INFO["FUSION"][0].replace("-","_")+"\n")
                primer_output.write(seq1+"[]"+seq2+"\n")
