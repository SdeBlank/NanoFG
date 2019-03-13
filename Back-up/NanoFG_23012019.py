import requests
import json
import vcf as pyvcf
import argparse
import datetime

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('-v', '--vcf', type=str, help='input VCF', required=True)

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
            if "SVLEN" in record.INFO:
                SVLEN=record.INFO["SVLEN"][0]
            else:
                SVLEN=0
            if (SVLEN<1000 and SVLEN!=0) or not isinstance(record.ALT[0], pyvcf.model._Breakend) or len(record.FILTER)>0:
                continue

            breakend1_info=breakend_annotation(CHROM1, POS1, POS1_ORIENTATION, 1)
            if len(breakend1_info)==0:
                continue
            breakend2_info=breakend_annotation(CHROM2, POS2, POS2_ORIENTATION, 2)
            if len(breakend2_info)==0:
                continue


            for gene1 in breakend1_info:
                for gene2 in breakend2_info:
                    if gene1["Phase"]==gene2["Phase"] and gene1["Strand"]==gene2["Strand"] and gene1["Gene_id"]!=gene2["Gene_id"]:
                        print(record.ID)
                        print(breakend1_info)
                        print(breakend2_info)
                        print("JAAAAA##########")

def breakend_annotation(CHROM, POS, orientation, order):
    ### Request genes, transcripts and exons
    SERVER="https://GRCh37.rest.ensembl.org/overlap/region/"
    SPECIES="human"
    HEADERS={"Content-Type" : "application/json"}
    FEATURE="gene"
    genes_request=requests.get(SERVER+SPECIES+"/"+str(CHROM)+":"+str(POS)+"-"+str(POS)+"?feature="+FEATURE, headers=HEADERS)
    genes_response = genes_request.text
    genes_data=json.loads(genes_response)

    HITS=[]

    for hit in genes_data:
        INFO={}
        GENE_ID=hit["gene_id"]
        if hit["biotype"]=="protein_coding":
            SERVER="https://GRCh37.rest.ensembl.org/lookup/id/"
            request = requests.get(SERVER+str(GENE_ID)+"?expand=1"+FEATURE, headers=HEADERS)
            response = request.text
            gene_info=json.loads(response)
            INFO["Gene_id"]=gene_info["id"]
            INFO["Gene_name"]=gene_info["display_name"]
            INFO["Strand"]=gene_info["strand"]

            for transcript in gene_info["Transcript"]:
                if transcript["is_canonical"]==1:
                    INFO["Transcript_id"]=transcript["id"]
                    INFO["Biotype"]=transcript["biotype"]
                    INFO["Total_exons"]=len(transcript["Exon"])
                    if POS>transcript["Translation"]["start"] and POS<transcript["Translation"]["end"]:
                        INFO["CDS"]="Within"
                    elif INFO["Strand"]==1:
                        if POS<transcript["Translation"]["start"]:
                            INFO["CDS"]="Before"
                            INFO["Distance_to_CDS"]=abs(transcript["Translation"]["start"]-POS)
                        elif POS>transcript["Translation"]["end"]:
                            INFO["CDS"]="After"
                            INFO["Distance_to_CDS"]=abs(POS-transcript["Translation"]["end"])
                    elif INFO["Strand"]==-1:
                        if POS<transcript["Translation"]["start"]:
                            INFO["CDS"]="After"
                            INFO["Distance_to_CDS"]=abs(transcript["Translation"]["start"]-POS)
                        elif POS>transcript["Translation"]["end"]:
                            INFO["CDS"]="Before"
                            INFO["Distance_to_CDS"]=abs(POS-transcript["Translation"]["end"])

            FEATURE="exon"
            SERVER="https://GRCh37.rest.ensembl.org/overlap/id/"
            request = requests.get(SERVER+INFO["Transcript_id"]+"?feature="+FEATURE, headers=HEADERS)
            response = request.text
            data=json.loads(response)
            EXONS={}
            for hit in data:
                if hit["Parent"]==INFO["Transcript_id"]:
                    EXONS[hit["rank"]]=hit

            for exon_id, exon in EXONS.items():
                MATCH=False
                if not exon["ensembl_end_phase"]==-1 and not exon["ensembl_phase"]==-1:
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
                            elif exon["rank"]>1:
                                INFO["Rank"]=exon["rank"]-1
                                INFO["Phase"]=EXONS[INFO["Rank"]]["ensembl_end_phase"]              ### What happens if breakend before CDS but in-frame?
                            else:
                                break
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
                        elif POS>exon["start"] and exon["rank"]>1:
                            INFO["Type"]="intron"
                            if orientation:
                                INFO["Rank"]=exon["rank"]-1
                                INFO["Phase"]=EXONS[INFO["Rank"]]["ensembl_end_phase"]
                            else:
                                INFO["Rank"]=exon["rank"]
                                INFO["Phase"]=exon["ensembl_phase"]
                            MATCH=True
                else:
                    if exon["ensembl_end_phase"]==-1 and exon["ensembl_phase"]==-1:
                        continue
                    elif exon["ensembl_end_phase"]==-1:
                        if INFO["CDS"]=="within":


                if MATCH:
                    INFO["Exon_id"]=EXONS[INFO["Rank"]]["exon_id"]
                    INFO["Exon_start"]=EXONS[INFO["Rank"]]["start"]
                    INFO["Exon_end"]=EXONS[INFO["Rank"]]["end"]
                    HITS.append(INFO)
                    break

    return(HITS)

# CHROM=7
# POS=34035524
# TYPE="TH"
#breakend_annotation(CHROM, POS)
print("Start:", datetime.datetime.now())
parse_vcf(VCF)
print("End:", datetime.datetime.now())
#print(breakend_annotation(CHROM, POS))
