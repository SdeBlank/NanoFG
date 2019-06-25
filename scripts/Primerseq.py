import vcf as pyvcf
import argparse
import sys
import requests
import json
import time
import sys
from EnsemblRestClient import EnsemblRestClient

server = "http://grch37.rest.ensembl.org"

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('-v', '--vcf', required=True, type=str, help='Path to fusion gene vcf')
parser.add_argument('-o', '--output', required=True, type=str, help='Path to output')
parser.add_argument('-f', '--flank', default=200, type=int, help='Flank [default: 200]')
args = parser.parse_args()
vcf = args.vcf
output=args.output
flank = args.flank

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


########################################### MAIN CODE
EnsemblRestClient=EnsemblRestClient()

with open(vcf, "r") as fusion_vcf, open(output, "w") as primer_output:
    vcf_reader = pyvcf.Reader(vcf)

    if "source" in vcf_reader.metadata:
        if vcf_reader.metadata["source"][0].lower()=="sniffles":
            vcf_type="Sniffles"
    elif "cmdline" in vcf_reader.metadata:
        if "nanosv" in vcf_reader.metadata["cmdline"][0].lower():
            vcf_type="NanoSV"

    for record in vcf_reader:
        if "FUSION" in record.INFO:
            if isinstance(record.ALT[0], pyvcf.model._SV) or isinstance( record.ALT[0], pyvcf.model._Substitution ):
                record = alt_convert( record )
            if not record:
                continue

            if vcf_type=="NanoSV":
                pos1_orientation=record.ALT[0].orientation
                pos2_orientation=record.ALT[0].remoteOrientation
            elif vcf_type=="Sniffles":
                if record.INFO["STRANDS"][0][0]=="+":
                    pos1_orientation=False
                else:
                    pos1_orientation=True
                if record.INFO["STRANDS"][0][1]=="+":
                    pos2_orientation=False
                else:
                    pos2_orientation=True

            if pos1_orientation:
                region1=record.CHROM+":"str(record.POS)+"-"+str(record.POS+flank)
            else:
                region1=record.CHROM+":"str(record.POS)+"-"+str(record.POS-flank)

            if pos2_orientation:
                region2=record.ALT[0].chr+":"str(record.ALT[0].pos)+"-"+str(record.ALT[0].pos+flank)
            else:
                region2=record.ALT[0].chr+":"str(record.ALT[0].pos)+"-"+str(record.ALT[0].pos-flank)

            primer_output.write(record.ID+"\t"+region1+"\t"+region2+"\n")
