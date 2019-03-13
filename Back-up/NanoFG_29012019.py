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

            # Do not activate filter step yet, testing on SOMATIC set, so filter will contain BPI-SOMATIC and PCR-SOMATIC
            # if not isinstance(record.ALT[0], pyvcf.model._Breakend) or len(record.FILTER)>0:
            #     continue

            breakend1_annotation=overlap_annotation(CHROM1, POS1)
            if not breakend1_annotation:
                continue
            breakend2_annotation=overlap_annotation(CHROM2, POS2)
            if not breakend2_annotation:
                continue

            breakend1_info=breakend_annotation(CHROM1, POS1, POS1_ORIENTATION, breakend1_annotation)
            breakend2_info=breakend_annotation(CHROM2, POS2, POS2_ORIENTATION, breakend2_annotation)

            fusion_check(record, breakend1_info, breakend2_info, POS1_ORIENTATION, POS2_ORIENTATION)

def overlap_annotation(CHROM, POS):
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
        GENE_ID=hit["gene_id"]
        if hit["biotype"]=="protein_coding":
            INFO={}
            SERVER="https://GRCh37.rest.ensembl.org/lookup/id/"
            request = requests.get(SERVER+str(GENE_ID)+"?expand=1"+FEATURE, headers=HEADERS)
            response = request.text
            gene_info=json.loads(response)
            INFO["Gene_id"]=gene_info["id"]
            INFO["Gene_name"]=gene_info["display_name"]
            INFO["Strand"]=gene_info["strand"]
            INFO["Gene_start"]=gene_info["start"]
            INFO["Gene_end"]=gene_info["end"]

            for transcript in gene_info["Transcript"]:
                if transcript["is_canonical"]==1:
                    CDS=False
                    INFO["Transcript_id"]=transcript["id"]
                    INFO["Biotype"]=transcript["biotype"]
                    INFO["Transcript_start"]=transcript["start"]
                    INFO["Transcript_end"]=transcript["end"]
                    INFO["Total_exons"]=len(transcript["Exon"])
                    if INFO["Strand"]==1:
                        INFO["CDS_start"]=transcript["Translation"]["start"]
                        INFO["CDS_end"]=transcript["Translation"]["end"]
                    else:
                        INFO["CDS_start"]=transcript["Translation"]["end"]
                        INFO["CDS_end"]=transcript["Translation"]["start"]
                    INFO["Exons"]=[]

                    ### EXONS
                    for rank, exon in enumerate(transcript["Exon"]):
                        EXON_INFO={}
                        EXON_INFO["Rank"]=rank+1
                        EXON_INFO["Type"]="exon"
                        CHRON_START=exon["start"]
                        CHRON_END=exon["end"]
                        if INFO["Strand"]==1:
                            EXON_INFO["Start"]=exon["start"]
                            EXON_INFO["End"]=exon["end"]
                        else:
                            EXON_INFO["Start"]=exon["end"]
                            EXON_INFO["End"]=exon["start"]
                        EXON_INFO["Contains_start_CDS"]=False
                        EXON_INFO["Contains_end_CDS"]=False

                        if CDS:
                            if INFO["CDS_end"]>=CHRON_START and INFO["CDS_end"]<=CHRON_END:
                                EXON_INFO["Contains_end_CDS"]=True
                                EXON_INFO["CDS"]=True
                                EXON_INFO["Start_phase"]=PHASE
                                EXON_INFO["End_phase"]=-1
                                EXON_INFO["CDS_length"]=abs(INFO["CDS_end"]-EXON_INFO["Start"])+1
                                CDS=False
                            else:
                                EXON_INFO["CDS"]=True
                                EXON_INFO["Start_phase"]=PHASE
                                EXON_INFO["End_phase"]=(abs(EXON_INFO["End"]-EXON_INFO["Start"])+1+PHASE)%3
                                EXON_INFO["CDS_length"]=abs(EXON_INFO["End"]-EXON_INFO["Start"])+1
                        elif INFO["CDS_start"]>=CHRON_START and INFO["CDS_start"]<=CHRON_END:
                            EXON_INFO["Contains_start_CDS"]=True
                            EXON_INFO["CDS"]=True
                            EXON_INFO["Start_phase"]=-1
                            EXON_INFO["End_phase"]=(abs(EXON_INFO["End"]-INFO["CDS_start"])+1)%3
                            EXON_INFO["CDS_length"]=abs(EXON_INFO["End"]-INFO["CDS_start"])+1
                            CDS=True

                        else:
                            EXON_INFO["CDS"]=False
                            EXON_INFO["Start_phase"]=-1
                            EXON_INFO["End_phase"]=-1
                            EXON_INFO["CDS_length"]=0
                        PHASE=EXON_INFO["End_phase"]
                        INFO["Exons"].append(EXON_INFO)

                        ### INTRONS
                        if rank<len(transcript["Exon"])-1:
                            INTRON_INFO={}
                            INTRON_INFO["Type"]="intron"
                            INTRON_INFO["Rank"]=rank+1
                            INTRON_INFO["Start"]=transcript["Exon"][rank]["end"]
                            INTRON_INFO["End"]=transcript["Exon"][rank+1]["start"]
                            INTRON_INFO["CDS"]=CDS
                            INFO["Exons"].append(INTRON_INFO)
            HITS.append(INFO)
    return HITS

def breakend_annotation(CHROM, POS, orientation, Info):
    HITS=[]
    for gene in Info:
        BND_INFO={k:v for k,v in gene.items() if k!="Exons"}
        # What if BND lands outside transcript but inside the gene (Possibly only promoter)
        #print(gene)
        if gene["Strand"]==1:
            if POS<gene["Transcript_start"]:
                BND_INFO["Breakpoint_location"]="before_transcript"
            elif POS<gene["CDS_start"]:
                BND_INFO["Breakpoint_location"]="5'UTR"
            elif POS>gene["CDS_start"] and POS<gene["CDS_end"]:
                BND_INFO["Breakpoint_location"]="CDS"
            elif POS>gene["Transcript_end"]:
                BND_INFO["Breakpoint_location"]="after_transcript"
            elif POS>gene["CDS_end"]:
                BND_INFO["Breakpoint_location"]="3'UTR"
        else:
            orientation= not orientation
            if POS>gene["Transcript_start"]:
                BND_INFO["Breakpoint_location"]="before_transcript"
            elif POS>gene["CDS_start"]:
                BND_INFO["Breakpoint_location"]="5'UTR"
            elif POS<gene["CDS_start"] and POS>gene["CDS_end"]:
                BND_INFO["Breakpoint_location"]="CDS"
            elif POS<gene["Transcript_end"]:
                BND_INFO["Breakpoint_location"]="after_transcript"
            elif POS<gene["CDS_end"]:
                BND_INFO["Breakpoint_location"]="3'UTR"

        if orientation:
            BND_INFO["Order"]="3'"
        else:
            BND_INFO["Order"]="5'"

        for idx, piece in enumerate(gene["Exons"]):
            if gene["Strand"]==1:
                CHRON_START=piece["Start"]
                CHRON_END=piece["End"]
            else:
                CHRON_START=piece["End"]
                CHRON_END=piece["Start"]

            if not BND_INFO["Breakpoint_location"]=="after_transcript" and not BND_INFO["Breakpoint_location"]=="before_transcript":
                if POS>=CHRON_START and POS<=CHRON_END:
                    BND_INFO["Type"]=piece["Type"]
                    if BND_INFO["Type"]=="exon":
                        BND_INFO["Rank"]=piece["Rank"]
                        if orientation:
                            BND_INFO["Phase"]=(abs(piece["End"]-POS)+1+piece["End_phase"])%3
                            BND_INFO["CDS_length"]=abs(POS-piece["End"])+1
                            for rank in range(idx+1, len(gene["Exons"])):
                                if gene["Exons"][rank]["Type"]=="exon":
                                    BND_INFO["CDS_length"]+=gene["Exons"][rank]["CDS_length"]
                        else:
                            BND_INFO["Phase"]=(abs(POS-piece["Start"])+1+piece["Start_phase"])%3
                            BND_INFO["CDS_length"]=abs(POS-piece["Start"])+1
                            for rank in range(0, idx):
                                if gene["Exons"][rank]["Type"]=="exon":
                                    BND_INFO["CDS_length"]+=gene["Exons"][rank]["CDS_length"]
                    else:
                        BND_INFO["Rank"]=piece["Rank"]
                        if orientation:
                            BND_INFO["Phase"]=gene["Exons"][idx+1]["Start_phase"]
                            BND_INFO["CDS_length"]=0
                            for rank in range(idx, len(gene["Exons"])):
                                if gene["Exons"][rank]["Type"]=="exon":
                                    BND_INFO["CDS_length"]+=gene["Exons"][rank]["CDS_length"]
                        else:
                            BND_INFO["Phase"]=gene["Exons"][idx-1]["End_phase"]
                            BND_INFO["CDS_length"]=0
                            for rank in range(0, idx):
                                if gene["Exons"][rank]["Type"]=="exon":
                                    BND_INFO["CDS_length"]+=gene["Exons"][rank]["CDS_length"]
        if BND_INFO["Breakpoint_location"]=="CDS":
            HITS.append(BND_INFO)
    return HITS

def fusion_check(Record, Breakend1, Breakend2, Orientation1, Orientation2):
    for annotation1 in Breakend1:
        for annotation2 in Breakend2:
            if annotation1["Gene_id"]!=annotation2["Gene_id"] and annotation1["Gene_name"]!=annotation2["Gene_name"]:
                if Orientation1 and Orientation2:
                    if (annotation1["Strand"]!=annotation2["Strand"] and annotation1["Phase"]==annotation2["Phase"] and annotation1["Type"]==annotation2["Type"] and annotation1["Breakpoint_location"]==annotation2["Breakpoint_location"]):
                        print(Record.ID)
                        print(annotation1["Order"], annotation1["Gene_name"], annotation2["Gene_name"], annotation2["Order"])
                        print(annotation1)
                        print(annotation2)

                elif Orientation1 and not Orientation2:
                    if (annotation1["Strand"]==annotation2["Strand"] and annotation1["Phase"]==annotation2["Phase"] and annotation1["Type"]==annotation2["Type"] and annotation1["Breakpoint_location"]==annotation2["Breakpoint_location"]):
                        print(Record.ID)
                        print(annotation1["Order"], annotation1["Gene_name"], annotation2["Gene_name"], annotation2["Order"])
                        print(annotation1)
                        print(annotation2)

                if not Orientation1 and Orientation2:
                    if (annotation1["Strand"]==annotation2["Strand"] and annotation1["Phase"]==annotation2["Phase"] and annotation1["Type"]==annotation2["Type"] and annotation1["Breakpoint_location"]==annotation2["Breakpoint_location"]):
                        print(Record.ID)
                        print(annotation1["Order"], annotation1["Gene_name"], annotation2["Gene_name"], annotation2["Order"])
                        print(annotation1)
                        print(annotation2)

                if not Orientation1 and not Orientation2:
                    if (annotation1["Strand"]!=annotation2["Strand"] and annotation1["Phase"]==annotation2["Phase"] and annotation1["Type"]==annotation2["Type"] and annotation1["Breakpoint_location"]==annotation2["Breakpoint_location"]):
                        print(Record.ID)
                        print(annotation1["Order"], annotation1["Gene_name"], annotation2["Gene_name"], annotation2["Order"])
                        print(annotation1)
                        print(annotation2)


print("Start:", datetime.datetime.now())
parse_vcf(VCF)
print("End:", datetime.datetime.now())
