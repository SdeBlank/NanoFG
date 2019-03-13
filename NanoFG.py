import requests
import json
import vcf as pyvcf
import argparse
import datetime
import sys

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('-v', '--vcf', type=str, help='input VCF', required=True)

args = parser.parse_args()
VCF=args.vcf

def parse_vcf(vcf):
    with open(VCF, "r") as vcf:
        VCF_READER=pyvcf.Reader(vcf)

        for record in VCF_READER:
            # Do not activate filter step yet, testing on SOMATIC set, so filter will contain BPI-SOMATIC and PCR-SOMATIC
            if not isinstance(record.ALT[0], pyvcf.model._Breakend):# or len(record.FILTER)>0:
                continue
            CHROM1=record.CHROM
            POS1=record.POS
            CHROM2=record.ALT[0].chr
            POS2=record.ALT[0].pos
            POS1_ORIENTATION=record.ALT[0].orientation
            POS2_ORIENTATION=record.ALT[0].remoteOrientation

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
    TRY=1
    while TRY <= 10:
        try:
            genes_request=requests.get(SERVER+SPECIES+"/"+str(CHROM)+":"+str(POS)+"-"+str(POS)+"?feature="+FEATURE, headers=HEADERS)
            genes_response = genes_request.text
            genes_data=json.loads(genes_response)
            break
        except:
            if TRY==10:
                sys.exit("Error while requesting from ENSEMBL database after "+str(TRY)+" tries")
            TRY +=1
    HITS=[]

    for hit in genes_data:
        GENE_ID=hit["gene_id"]
        if hit["biotype"]=="protein_coding":
            if hit["description"]:
                if "readthrough" not in hit["description"]:
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
                                        #EXON_INFO["End_phase"]=-1
                                        EXON_INFO["End_phase"]=0
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
                                    #EXON_INFO["Start_phase"]=-1
                                    EXON_INFO["Start_phase"]=0
                                    if INFO["CDS_end"]>=CHRON_START and INFO["CDS_end"]<=CHRON_END:
                                        EXON_INFO["Contains_end_CDS"]=True
                                        EXON_INFO["End_phase"]=0
                                        EXON_INFO["CDS_length"]=abs(INFO["CDS_end"]-INFO["CDS_start"])+1
                                    else:
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
                                    INTRON_INFO["Phase"]=EXON_INFO["End_phase"]
                                    INTRON_INFO["Start"]=transcript["Exon"][rank]["end"]
                                    INTRON_INFO["End"]=transcript["Exon"][rank+1]["start"]
                                    INTRON_INFO["CDS"]=CDS
                                    INFO["Exons"].append(INTRON_INFO)
                    if INFO["Biotype"]=="protein_coding":
                        HITS.append(INFO)
    return HITS

def breakend_annotation(CHROM, POS, orientation, Info):
    HITS=[]
    for gene in Info:
        BND_INFO={k:v for k,v in gene.items() if k!="Exons"}                                                        #!!!!! ADD INFO AS WELL.. NEED INFORMATION ABOUT SURROUNDING EXONS FOR EXON_INTRON FUSIONS
        #BND_INFO={k:v for k,v in gene.items()}
        # What if BND lands outside transcript but inside the gene (Possibly only promoter)
        if gene["Strand"]==1:
            if POS<gene["Transcript_start"]:
                BND_INFO["Breakpoint_location"]="before_transcript"
                BND_INFO["Distance from transcript start"]=abs(POS-gene["Transcript_start"])
            elif POS<gene["CDS_start"]:
                BND_INFO["Breakpoint_location"]="5'UTR"
            elif POS>gene["CDS_start"] and POS<gene["CDS_end"]:
                if abs(POS-gene["CDS_end"])<3:
                    BND_INFO["Breakpoint_location"]="StopCodon"
                elif abs(POS-gene["CDS_end"])>3 and abs(POS-gene["CDS_end"])<10:                 #### ADDED pos inside stop codon==different and pos close to translation end gives 3'UTR fusion instead
                    BND_INFO["Breakpoint_location"]="3'UTR"
                else:
                    BND_INFO["Breakpoint_location"]="CDS"
            elif POS>gene["Transcript_end"]:
                BND_INFO["Breakpoint_location"]="after_transcript"
                BND_INFO["Distance from transcript start"]=abs(POS-gene["Transcript_end"])
            elif POS>gene["CDS_end"]:
                BND_INFO["Breakpoint_location"]="3'UTR"
        else:
            orientation= not orientation
            if POS<gene["Transcript_start"]:
                BND_INFO["Breakpoint_location"]="after_transcript"
                BND_INFO["Distance from transcript end"]=abs(POS-gene["Transcript_start"])
            elif POS<gene["CDS_start"]:
                BND_INFO["Breakpoint_location"]="3'UTR"
            elif POS>gene["CDS_start"] and POS<gene["CDS_end"]:
                if abs(POS-gene["CDS_start"])<3:
                    BND_INFO["Breakpoint_location"]="StopCodon"
                elif abs(POS-gene["CDS_start"])>3 and abs(POS-gene["CDS_start"])<10:
                    BND_INFO["Breakpoint_location"]="3'UTR"
                else:
                    BND_INFO["Breakpoint_location"]="CDS"
            elif POS>gene["Transcript_end"]:
                BND_INFO["Breakpoint_location"]="before_transcript"
                BND_INFO["Distance from transcript start"]=abs(POS-gene["Transcript_end"])
            elif POS>gene["CDS_end"]:
                BND_INFO["Breakpoint_location"]="5'UTR"

        if orientation:
            BND_INFO["Order"]="3'"
        else:
            BND_INFO["Order"]="5'"

        for idx, sequence in enumerate(gene["Exons"]):
            if gene["Strand"]==1:
                CHRON_START=sequence["Start"]
                CHRON_END=sequence["End"]
            else:
                CHRON_START=sequence["End"]
                CHRON_END=sequence["Start"]

            if not BND_INFO["Breakpoint_location"]=="after_transcript" and not BND_INFO["Breakpoint_location"]=="before_transcript":
                if POS>=CHRON_START and POS<=CHRON_END:
                    BND_INFO["Type"]=sequence["Type"]

                    if BND_INFO["Type"]=="exon":
                        BND_INFO["Rank"]=sequence["Rank"]

                        # if sequence["Contains_start_CDS"]==True:
                        #     BND_INFO["Phase"]=(abs(INFO["CDS_start"]-POS)+sequence["Start_phase"])%3
                        # else:
                        #     BND_INFO["Phase"]=(abs(sequence["Start"]-POS)+sequence["Start_phase"])%3
                        if BND_INFO["Breakpoint_location"]=="5'UTR" or BND_INFO["Breakpoint_location"]=="3'UTR":
                            BND_INFO["Phase"]=-1
                        else:
                            if orientation:
                                if sequence["Contains_end_CDS"]==True:
                                    #BND_INFO["Phase"]=(((abs(gene["CDS_end"]-POS)+1+(3-sequence["End_phase"]))%3)*2)%3
                                    BND_INFO["Phase"]=(abs(gene["CDS_start"]-POS)+sequence["Start_phase"])%3
                                    BND_INFO["CDS_length"]=abs(POS-gene["CDS_end"])+1
                                    print (gene["Gene_id"], POS)
                                else:
                                    #BND_INFO["Phase"]=(((abs(sequence["End"]-POS)+1+(3-sequence["End_phase"]))%3)*2)%3
                                    BND_INFO["Phase"]=(abs(sequence["Start"]-POS)+sequence["Start_phase"])%3
                                    BND_INFO["CDS_length"]=abs(POS-sequence["End"])+1
                                    print (gene["Gene_id"], POS)
                                    #BND_INFO["Phase"]=(abs(sequence["Start"]-POS)+sequence["Start_phase"])%3
                                #BND_INFO["Phase"]=(abs(sequence["Start"]-POS)+1+sequence["Start_phase"])%3

                                for rank in range(idx+1, len(gene["Exons"])):
                                    if gene["Exons"][rank]["Type"]=="exon":
                                        BND_INFO["CDS_length"]+=gene["Exons"][rank]["CDS_length"]
                            else:
                                if sequence["Contains_start_CDS"]==True:
                                    BND_INFO["Phase"]=(abs(POS-gene["CDS_start"])+1+sequence["Start_phase"])%3
                                    BND_INFO["CDS_length"]=abs(POS-gene["CDS_start"])+1
                                else:
                                    BND_INFO["Phase"]=(abs(POS-sequence["Start"])+1+sequence["Start_phase"])%3
                                    BND_INFO["CDS_length"]=abs(POS-sequence["Start"])+1
                                # BND_INFO["Phase"]=(abs(POS-sequence["Start"])+1+sequence["Start_phase"])%3

                                for rank in range(0, idx):
                                    if gene["Exons"][rank]["Type"]=="exon":
                                        BND_INFO["CDS_length"]+=gene["Exons"][rank]["CDS_length"]

                    else:
                        BND_INFO["Rank"]=sequence["Rank"]
                        BND_INFO["Phase"]=sequence["Phase"]
                        BND_INFO["CDS_length"]=0
                        if orientation:
                            for rank in range(idx, len(gene["Exons"])):
                                if gene["Exons"][rank]["Type"]=="exon":
                                    BND_INFO["CDS_length"]+=gene["Exons"][rank]["CDS_length"]
                        else:
                            for rank in range(0, idx):
                                if gene["Exons"][rank]["Type"]=="exon":
                                    BND_INFO["CDS_length"]+=gene["Exons"][rank]["CDS_length"]
                #else: UTR

        HITS.append(BND_INFO)
    return HITS

def fusion_check(Record, Breakend1, Breakend2, Orientation1, Orientation2):
    for annotation1 in Breakend1:
        for annotation2 in Breakend2:
            if annotation1["Gene_id"]!=annotation2["Gene_id"] and annotation1["Gene_name"]!=annotation2["Gene_name"]:
                if annotation1["Breakpoint_location"]==annotation2["Breakpoint_location"]:
                    print("\n")
                    print(Record.CHROM, annotation1["Breakpoint_location"], annotation2["Breakpoint_location"])
                    print(Record.POS, annotation1["Transcript_start"], annotation1["Transcript_end"])
                    print(Record.ALT[0].pos, annotation2["Transcript_start"], annotation2["Transcript_end"])
                    print(annotation1["Gene_name"], annotation2["Gene_name"])


                if Orientation1 and Orientation2:
                    if (annotation1["Strand"]!=annotation2["Strand"] and annotation1["Breakpoint_location"]==annotation2["Breakpoint_location"]):
                        if annotation1["Breakpoint_location"]=="CDS"  and annotation1["Type"]==annotation2["Type"]:
                            ##### Exon fusion with not the same phase can still be right, as calling might not be 100% accurate
                            if annotation1["Phase"]==annotation2["Phase"]:
                                print(Record.ID)
                                print("Fusion protein")
                                print(annotation1["Order"], annotation1["Gene_name"], annotation2["Gene_name"], annotation2["Order"])
                                print(annotation1)
                                print(annotation2)
                            elif annotation1["Type"]=="exon":
                                print(Record.ID)
                                print("Fusion protein (OUT OF PHASE)")
                                print(annotation1["Order"], annotation1["Gene_name"], annotation2["Gene_name"], annotation2["Order"])
                                print(annotation1)
                                print(annotation2)
                        elif annotation1["Breakpoint_location"]=="5'UTR":
                            print(Record.ID)
                            print("Promoter fusion")
                            print(annotation1["Order"], annotation1["Gene_name"], annotation2["Gene_name"], annotation2["Order"])
                            print(annotation1)
                            print(annotation2)
                        elif annotation1["Breakpoint_location"]=="3'UTR" and annotation1["Type"]==annotation2["Type"]:
                            print(Record.ID)
                            print("3'UTR fusion")
                            print(annotation1["Order"], annotation1["Gene_name"], annotation2["Gene_name"], annotation2["Order"])
                            print(annotation1)
                            print(annotation2)

                elif Orientation1 and not Orientation2:
                    if (annotation1["Strand"]==annotation2["Strand"] and annotation1["Breakpoint_location"]==annotation2["Breakpoint_location"]):
                        if annotation1["Breakpoint_location"]=="CDS" and annotation1["Type"]==annotation2["Type"]:
                            if annotation1["Phase"]==annotation2["Phase"]:
                                print(Record.ID)
                                print("Fusion protein")
                                print(annotation1["Order"], annotation1["Gene_name"], annotation2["Gene_name"], annotation2["Order"])
                                print(annotation1)
                                print(annotation2)
                            elif annotation1["Type"]=="exon":
                                print(Record.ID)
                                print("Fusion protein (OUT OF PHASE)")
                                print(annotation1["Order"], annotation1["Gene_name"], annotation2["Gene_name"], annotation2["Order"])
                                print(annotation1)
                                print(annotation2)
                        elif annotation1["Breakpoint_location"]=="5'UTR":
                            print(Record.ID)
                            print("Promoter fusion")
                            print(annotation1["Order"], annotation1["Gene_name"], annotation2["Gene_name"], annotation2["Order"])
                            print(annotation1)
                            print(annotation2)
                        elif annotation1["Breakpoint_location"]=="3'UTR" and annotation1["Type"]==annotation2["Type"]:
                            print(Record.ID)
                            print("3'UTR fusion")
                            print(annotation1["Order"], annotation1["Gene_name"], annotation2["Gene_name"], annotation2["Order"])
                            print(annotation1)
                            print(annotation2)

                if not Orientation1 and Orientation2:
                    if (annotation1["Strand"]==annotation2["Strand"] and annotation1["Breakpoint_location"]==annotation2["Breakpoint_location"]):
                        if annotation1["Breakpoint_location"]=="CDS" and annotation1["Type"]==annotation2["Type"]:
                            if annotation1["Phase"]==annotation2["Phase"]:
                                print(Record.ID)
                                print("Fusion protein")
                                print(annotation1["Order"], annotation1["Gene_name"], annotation2["Gene_name"], annotation2["Order"])
                                print(annotation1)
                                print(annotation2)
                            elif annotation1["Type"]=="exon":
                                print(Record.ID)
                                print("Fusion protein (OUT OF PHASE)")
                                print(annotation1["Order"], annotation1["Gene_name"], annotation2["Gene_name"], annotation2["Order"])
                                print(annotation1)
                                print(annotation2)
                        elif annotation1["Breakpoint_location"]=="5'UTR":
                            print(Record.ID)
                            print("Promoter fusion")
                            print(annotation1["Order"], annotation1["Gene_name"], annotation2["Gene_name"], annotation2["Order"])
                            print(annotation1)
                            print(annotation2)
                        elif annotation1["Breakpoint_location"]=="3'UTR" and annotation1["Type"]==annotation2["Type"]:
                            print(Record.ID)
                            print("3'UTR fusion")
                            print(annotation1["Order"], annotation1["Gene_name"], annotation2["Gene_name"], annotation2["Order"])
                            print(annotation1)
                            print(annotation2)

                if not Orientation1 and not Orientation2:
                    if (annotation1["Strand"]!=annotation2["Strand"] and annotation1["Breakpoint_location"]==annotation2["Breakpoint_location"]):
                        if annotation1["Breakpoint_location"]=="CDS" and annotation1["Type"]==annotation2["Type"]:
                            if annotation1["Phase"]==annotation2["Phase"]:
                                print(Record.ID)
                                print("Fusion protein")
                                print(annotation1["Order"], annotation1["Gene_name"], annotation2["Gene_name"], annotation2["Order"])
                                print(annotation1)
                                print(annotation2)
                            elif annotation1["Type"]=="exon":
                                print(Record.ID)
                                print("Fusion protein (OUT OF PHASE)")
                                print(annotation1["Order"], annotation1["Gene_name"], annotation2["Gene_name"], annotation2["Order"])
                                print(annotation1)
                                print(annotation2)
                        elif annotation1["Breakpoint_location"]=="5'UTR":
                            print(Record.ID)
                            print("Promoter fusion")
                            print(annotation1["Order"], annotation1["Gene_name"], annotation2["Gene_name"], annotation2["Order"])
                            print(annotation1)
                            print(annotation2)
                        elif annotation1["Breakpoint_location"]=="3'UTR" and annotation1["Type"]==annotation2["Type"]:
                            print(Record.ID)
                            print("3'UTR fusion")
                            print(annotation1["Order"], annotation1["Gene_name"], annotation2["Gene_name"], annotation2["Order"])
                            print(annotation1)
                            print(annotation2)


print("Start:", datetime.datetime.now())
parse_vcf(VCF)
print("End:", datetime.datetime.now())
