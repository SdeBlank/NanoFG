#!/usr/bin/python

import argparse
import datetime
from EnsemblRestClient import EnsemblRestClient
import sys
import re

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('-b', '--bed', type=str, help='Output bed file', required=True)
parser.add_argument('-r', '--region', type=str, help='List of genes or regions to select from the bam-file', required=True)

args = parser.parse_args()

bed=args.bed
selection=args.region
selection=selection.split(",")
EnsemblRestClient=EnsemblRestClient()

### Create a bed-file that contains all the regions that need to be selected from the full bam-file
with open(bed, "w") as bed:
    for item in selection:

        ### If a region is given in a chr:pos1-pos2 format, this format is directly used
        if re.match(".+:\d+-\d+", item):
            chrom=item.split(":")[0]
            pos1=item.split(":")[1].split("-")[0]
            pos2=item.split(":")[1].split("-")[1]

        else:
            server='http://grch37.rest.ensembl.org'

            ### If a desired region is given by a ensembl identifier (e.g. ENSG00000141510) or the gene name, the positions are taken from ensembl
            if re.match("^ENS", item):
                endpoint="/lookup/id/"+item
            else:
                endpoint="/lookup/symbol/homo_sapiens/"+item

            headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
            params={"expand": "1"}

            ensembl_request=EnsemblRestClient.perform_rest_action(server, endpoint, headers, params)
            chrom=str(ensembl_request["seq_region_name"])
            pos1=str(ensembl_request["start"])
            pos2=str(ensembl_request["end"])

        try:
            if pos1<pos2:
                bed.write("\t".join([chrom, pos1, pos2])+"\n")
            else:
                bed.write("\t".join([chrom, pos2, pos1])+"\n")
        except:
            ### If pos1 or pos2 does not exist, the given region is not correct or the gene cannot be found.
            sys.exit(item+" is not a valid selection. Use a valid ensembl id, gene name or region (e.g. 3:1234-2345)")
