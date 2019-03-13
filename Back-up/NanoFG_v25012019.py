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
            #print(record.ID)
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

            fusion_check(record, breakend1_info, breakend2_info, POS1_ORIENTATION, POS2_ORIENTATION)


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
                    INFO["CDS_start"]=transcript["Translation"]["start"]
                    INFO["CDS_end"]=transcript["Translation"]["end"]
                    if POS>transcript["Translation"]["start"] and POS<transcript["Translation"]["end"]:
                        INFO["CDS"]="Within"
                    elif INFO["Strand"]==1:
                        if POS<transcript["Translation"]["start"]:
                            INFO["CDS"]="Before"
                            INFO["Distance_to_CDS"]=abs(transcript["Translation"]["start"]-POS)
                        elif POS>transcript["Translation"]["end"]:
                            INFO["CDS"]="After"
                    elif INFO["Strand"]==-1:
                        if POS<transcript["Translation"]["start"]:
                            INFO["CDS"]="After"
                        elif POS>transcript["Translation"]["end"]:
                            INFO["CDS"]="Before"
                            INFO["Distance_to_CDS"]=abs(POS-transcript["Translation"]["end"])

            if INFO["CDS"]=="After":
                continue

            FEATURE="exon"
            SERVER="https://GRCh37.rest.ensembl.org/overlap/id/"
            request = requests.get(SERVER+INFO["Transcript_id"]+"?feature="+FEATURE, headers=HEADERS)
            response = request.text
            data=json.loads(response)
            EXONS={}
            for hit in data:
                if hit["Parent"]==INFO["Transcript_id"]:
                    if hit["strand"]==-1:
                        start=hit["start"]
                        hit["start"]=hit["end"]
                        hit["end"]=start
                    hit["Start_CDS"]=False
                    hit["End_CDS"]=False
                    if hit["ensembl_phase"]==-1 and hit["ensembl_end_phase"]!=-1:
                        hit["Start_CDS"]=True
                        START_CDS=hit["rank"]
                        hit["CDS_length"]=abs(hit["end"]-INFO["CDS_start"]+1)
                    elif hit["ensembl_phase"]!=-1 and hit["ensembl_phase"]!=-1:
                        hit["CDS_length"]=abs(hit["end"]-hit["start"]+1)
                    elif hit["ensembl_phase"]!=-1 and hit["ensembl_phase"]==-1:
                        hit["End_CDS"]=True
                        END_CDS=hit["rank"]
                        hit["CDS_length"]=abs(hit["start"]-INFO["CDS_end"]+1)
                    elif hit["ensembl_phase"]==-1 and hit["ensembl_phase"]==-1:
                        hit["CDS_length"]=0
                    EXONS[hit["rank"]]=hit

            for rank in range(INFO["Total_exons"]):
                RANK=rank+1
                exon=EXONS[RANK]
                MATCH=False
                #print(exon["start"], exon["end"], hit["strand"])
                if INFO["CDS"]=="Within":
                    if exon["strand"]==1:
                        if POS>=exon["start"] and POS<=exon["end"]:
                            INFO["Type"]="exon"
                            INFO["Rank"]=exon["rank"]
                            if orientation:
                                INFO["Phase"]=(3-((exon["end"]-exon["ensembl_end_phase"]-POS+1)%3))%3
                                INFO["CDS_length"]=abs(exon["end"]-POS)+1
                                for rank in range(INFO["Rank"]+1, INFO["Total_exons"]+1):
                                    INFO["CDS_length"]+=EXONS[rank]["CDS_length"]
                            else:
                                INFO["Phase"]=(POS-exon["start"]+(3-exon["ensembl_phase"])%3)%3
                                INFO["CDS_length"]=abs(POS-exon["start"])+1
                                if INFO["Start_CDS"]:
                                    INFO["CDS_length"]=abs(POS-exon["start"])+1
                                for rank in range(1, INFO["Rank"]):
                                    INFO["CDS_length"]+=EXONS[rank]["CDS_length"]
                            MATCH=True
                        elif POS<exon["start"]:
                            INFO["Type"]="intron"
                            if orientation:
                                INFO["Rank"]=exon["rank"]
                                INFO["Phase"]=exon["ensembl_phase"]
                                INFO["CDS_length"]=0
                                for rank in range(INFO["Rank"], INFO["Total_exons"]+1):
                                    INFO["CDS_length"]+=EXONS[rank]["CDS_length"]
                            elif exon["rank"]>1:
                                INFO["Rank"]=exon["rank"]-1
                                INFO["Phase"]=EXONS[INFO["Rank"]]["ensembl_end_phase"]              ### What happens if breakend before CDS but in-frame?
                                INFO["CDS_length"]=0
                                for rank in range(1, INFO["Rank"]+1):
                                    INFO["CDS_length"]+=EXONS[rank]["CDS_length"]
                            else:
                                break
                            MATCH=True
                    elif exon["strand"]==-1:
                        if POS>=exon["end"] and POS<=exon["start"]:
                            INFO["Type"]="exon"
                            # if INFO["CDS"]=="before":
                            #     type=UTR
                            #     Distance_to_CDS=sadadsa
                            INFO["Rank"]=exon["rank"]
                            if orientation:
                                INFO["Phase"]=(exon["start"]-POS+(3-exon["ensembl_phase"])%3)%3
                                INFO["CDS_length"]=abs(POS-exon["start"])+1
                                for rank in range(1, INFO["Rank"]):
                                    INFO["CDS_length"]+=EXONS[rank]["CDS_length"]
                            else:
                                INFO["Phase"]=(3-((POS-exon["end"]+exon["ensembl_end_phase"])%3))%3
                                INFO["CDS_length"]=abs(exon["end"]-POS)+1
                                for rank in range(INFO["Rank"]+1, INFO["Total_exons"]+1):
                                    INFO["CDS_length"]+=EXONS[rank]["CDS_length"]
                            MATCH=True
                        elif POS>exon["start"]:
                            INFO["Type"]="intron"
                            if orientation:
                                if exon["rank"]>1:
                                    INFO["Rank"]=exon["rank"]-1
                                    INFO["Phase"]=EXONS[INFO["Rank"]]["ensembl_end_phase"]
                                    INFO["CDS_length"]=0
                                    for rank in range(1, INFO["Rank"]+1):
                                        INFO["CDS_length"]+=EXONS[rank]["CDS_length"]
                            else:
                                INFO["Rank"]=exon["rank"]
                                INFO["Phase"]=exon["ensembl_phase"]
                                INFO["CDS_length"]=0
                                for rank in range(INFO["Rank"], INFO["Total_exons"]+1):
                                    INFO["CDS_length"]+=EXONS[rank]["CDS_length"]

                            MATCH=True

                elif INFO["CDS"]=="Before":
                    if exon["strand"]==1:
                        if orientation:
                            if POS>=exon["start"] and POS<=exon["end"]:
                                INFO["Type"]="exon"
                                INFO["Rank"]=exon["rank"]
                                INFO["Distance_to_CDS"]=abs(POS-exon["end"])+1
                                for rank in range(INFO["RANK"]+1, START_CDS):
                                    if rank==START_CDS
                                        INFO["Distance_to_CDS"]+=abs(EXONS[rank]["start"]-INFO["CDS_start"])+1
                                    else:
                                        INFO["Distance_to_CDS"]+=abs(EXONS[rank]["end"]-EXONS[rank]["start"])+1
                                INFO["Phase"]=(3-((exon["end"]-exon["ensembl_end_phase"]-POS+1)%3))%3
                                INFO["CDS_length"]=0
                                for rank in range(1, INFO["Total_exons"]+1):
                                    INFO["CDS_length"]+=EXONS[rank]["CDS_length"]
                            elif POS<exon["start"]:
                    elif exon["strand"]==-1:
                        if not orientation:
                            INFO["Type"]="5'UTR"
                            INFO["Rank"]=1




                if MATCH:
                    INFO["Exon_id"]=EXONS[INFO["Rank"]]["exon_id"]
                    INFO["Exon_start"]=EXONS[INFO["Rank"]]["start"]
                    INFO["Exon_end"]=EXONS[INFO["Rank"]]["end"]
                    HITS.append(INFO)
                    break

    return(HITS)



def fusion_check(Record, Breakend1, Breakend2, Orientation1, Orientation2):
    for annotation1 in Breakend1:
        for annotation2 in Breakend2:
            if annotation1["Gene_id"]!=annotation2["Gene_id"]:
                if Orientation1 and Orientation2:
                    if (annotation1["Strand"]!=annotation2["Strand"] and annotation1["Phase"]==annotation2["Phase"] and
                            annotation1["Type"]==annotation2["Type"] and annotation1["CDS"]==annotation2["CDS"]):
                        #print(Record.ID + "\t" + annotation1["Gene_name"] +"-"+annotation2["Gene_name"])
                        print(annotation1["Gene_name"], annotation2["Gene_name"])
                        print(annotation1)
                        print(annotation2)
                        #return True
                    #else:
                        #return False
                elif Orientation1 and not Orientation2:
                    if (annotation1["Strand"]==annotation2["Strand"] and annotation1["Phase"]==annotation2["Phase"] and
                            annotation1["Type"]==annotation2["Type"] and annotation1["CDS"]==annotation2["CDS"]):
                        #print(Record.ID + "\t" + annotation1["Gene_name"] +"-"+annotation2["Gene_name"])
                        print(annotation1["Gene_name"], annotation2["Gene_name"])
                        print(annotation1)
                        print(annotation2)
                        #return True
                    #else:
                        #return False
                if not Orientation1 and Orientation2:
                    if (annotation1["Strand"]==annotation2["Strand"] and annotation1["Phase"]==annotation2["Phase"] and
                            annotation1["Type"]==annotation2["Type"] and annotation1["CDS"]==annotation2["CDS"]):
                        #print(Record.ID + "\t" + annotation1["Gene_name"] +"-"+annotation2["Gene_name"])
                        print(annotation1["Gene_name"], annotation2["Gene_name"])
                        print(annotation1)
                        print(annotation2)
                        #return True
                    #else:
                        #return False
                if not Orientation1 and not Orientation2:
                    if (annotation1["Strand"]!=annotation2["Strand"] and annotation1["Phase"]==annotation2["Phase"] and
                            annotation1["Type"]==annotation2["Type"] and annotation1["CDS"]==annotation2["CDS"]):
                        #print(Record.ID + "\t" + annotation1["Gene_name"] +"-"+annotation2["Gene_name"])
                        print(annotation1["Gene_name"], annotation2["Gene_name"])
                        print(annotation1)
                        print(annotation2)
                        #return True
                    #else:
                        #return False


# CHROM=7
# POS=34035524
# TYPE="TH"
#breakend_annotation(CHROM, POS)
print("Start:", datetime.datetime.now())

parse_vcf(VCF)

print("End:", datetime.datetime.now())
#print(breakend_annotation(CHROM, POS))
