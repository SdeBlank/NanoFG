import requests
import json
import vcf as pyvcf
import argparse
import datetime
print(1, datetime.datetime.now())
parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('vcf', help='VCF file')

args = parser.parse_args()
VCF=args.vcf

def parse_vcf(vcf):
    with open(VCF, "r") as vcf:
        VCF_READER=pyvcf.Reader(vcf)

        for record in VCF_READER:
            CHROM1=record.CHROM
            POS1=record.POS
            CHROM2=record.ALT[0].chr
            POS2=record.ALT[0].pos
            POS1_ORIENTATION=record.ALT[0].orientation
            POS2_ORIENTATION=record.ALT[0].remoteOrientation
            #if int(record.ID)==31116:
            breakend1_info=breakend_annotation(CHROM1, POS1, POS1_ORIENTATION, 1)
            breakend2_info=breakend_annotation(CHROM2, POS2, POS2_ORIENTATION, 2)
            print(record.ID)
            print(len(breakend1_info))
            print(len(breakend2_info))

def breakend_annotation(CHROM, POS, orientation, order):
    ### Request genes, transcripts and exons
    SERVER="https://GRCh37.rest.ensembl.org/overlap/region/"
    SPECIES="human"
    HEADERS={"Content-Type" : "application/json"}
    FEATURE="gene"
    request = requests.get(SERVER+SPECIES+"/"+str(CHROM)+":"+str(POS)+"-"+str(POS)+"?feature="+FEATURE, headers=HEADERS)
    response = request.text
    data=json.loads(response)
    GENES={}
    for hit in data:
        hit["transcripts"]={}
        GENES[hit["gene_id"]]=hit

    ### Shorten time to request gene and transcript together? Is gene always above transcript in request?
    FEATURE="transcript"
    request = requests.get(SERVER+SPECIES+"/"+str(CHROM)+":"+str(POS)+"-"+str(POS)+"?feature="+FEATURE, headers=HEADERS)
    response = request.text
    data=json.loads(response)
    for hit in data:
        for gene_id, gene in GENES.items():
            if hit["Parent"]==gene["gene_id"]:
                hit["exons"]={}
                hit["total_exons"]=0
                if "tag" not in hit:
                    hit["tag"]="none"
                GENES[gene_id]["transcripts"][hit["transcript_id"]]=hit

    FEATURE="exon"
    SERVER="https://GRCh37.rest.ensembl.org/overlap/id/"
    for gene_id, gene in GENES.items():
        if gene["biotype"]=="protein_coding":
            request = requests.get(SERVER+gene["id"]+"?feature="+FEATURE, headers=HEADERS)
            response = request.text
            data=json.loads(response)
            for transcript_id, transcript in GENES[gene_id]["transcripts"].items():
                for hit in data:
                    if hit["Parent"]==transcript["id"]:
                        GENES[gene_id]["transcripts"][transcript_id]["total_exons"]+=1
                        GENES[gene_id]["transcripts"][transcript_id]["exons"][hit["rank"]]=hit

    ### Select overlapping regions
    HITS=[]
    for gene_id, gene in GENES.items():
        TRANSCRIPTS=gene["transcripts"]
        if gene["biotype"] != "protein_coding":
            continue
        for transcript_id, transcript in TRANSCRIPTS.items():
            INFO={}
            EXONS=transcript["exons"]
            INFO["Gene_id"]=gene["gene_id"]
            INFO["Gene_name"]=gene["external_name"]
            INFO["Strand"]=gene["strand"]
            INFO["Transcript_id"]=transcript["transcript_id"]
            INFO["Biotype"]=transcript["biotype"]
            INFO["Total_exons"]=transcript["total_exons"]
            INFO["Flag"]=transcript["tag"]
            if INFO["Biotype"] != "protein_coding":
                continue
            for exon_id, exon in EXONS.items():
                MATCH=False
                if exon["ensembl_end_phase"] is not -1:
                    if exon["ensembl_phase"] is not -1:
                        INFO["CDS"]=True
                        if exon["strand"]==1:
                            if POS>=exon["start"] and POS<=exon["end"]:
                                INFO["Type"]="exon"
                                INFO["Rank"]=exon["rank"]
                                if orientation:
                                    INFO["Phase"]=(3-((exon["end"]-exon["ensembl_end_phase"]-POS+1)%3))%3
                                else:
                                    INFO["Phase"]=(POS-exon["start"]+(3-exon["ensembl_phase"])%3)%3
                                MATCH=True
                            elif POS<exon["start"]:
                                INFO["Type"]="intron"
                                if orientation:
                                    INFO["Rank"]=exon["rank"]
                                    INFO["Phase"]=exon["ensembl_phase"]
                                else:
                                    INFO["Rank"]=exon["rank"]-1
                                    INFO["Phase"]=EXONS[INFO["Rank"]]["ensembl_end_phase"]
                                MATCH=True
                        elif exon["strand"]==-1:
                            if POS>=exon["end"] and POS<=exon["start"]:
                                INFO["Type"]="exon"
                                INFO["Rank"]=exon["rank"]
                                if orientation:
                                    INFO["Phase"]=(exon["start"]-POS+(3-exon["ensembl_phase"])%3)%3
                                else:
                                    INFO["Phase"]=(3-((POS-exon["end"]+exon["ensembl_end_phase"])%3))%3
                                MATCH=True
                            elif POS>exon["start"]:
                                INFO["Type"]="intron"
                                if orientation:
                                    INFO["Rank"]=exon["rank"]-1
                                    INFO["Phase"]=EXONS[INFO["Rank"]]["ensembl_end_phase"]
                                else:
                                    INFO["Rank"]=exon["rank"]
                                    INFO["Phase"]=exon["ensembl_phase"]
                                MATCH=True
                else:
                    INFO["CDS"]=False
                if MATCH:
                    INFO["Exon_id"]=EXONS[INFO["Rank"]]["exon_id"]
                    INFO["Exon_start"]=EXONS[INFO["Rank"]]["start"]
                    INFO["Exon_end"]=EXONS[INFO["Rank"]]["end"]
                    break
            HITS.append(INFO)
    return(HITS)

# CHROM=7
# POS=34035524
# TYPE="TH"
#breakend_annotation(CHROM, POS)
parse_vcf(VCF)
#print(breakend_annotation(CHROM, POS))
print(2, datetime.datetime.now())
