import vcf as pyvcf
import argparse
import datetime
import sys
from EnsemblRestClient import EnsemblRestClient

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('-v', '--vcf', type=str, help='input VCF', required=True)
parser.add_argument('-fo', '--fusion_output', type=str, help='Fusion gene output file (table)', required=True)
parser.add_argument('-o', '--output', type=str, help='Annotated VCF output', required=True)

args = parser.parse_args()

def ParseVcf(vcf, vcf_output, info_output):
    with open(vcf, "r") as vcf, open(vcf_output, "w") as vcf_output, open(info_output, "w") as fusion_output:
        VCF_READER=pyvcf.Reader(vcf)
        VCF_READER.infos['FUSION']=pyvcf.parser._Info('FUSION', ".", "String", "Gene names of the fused genes reported if present", "NanoSV", "X")
        VCF_WRITER=pyvcf.Writer(vcf_output, VCF_READER, lineterminator='\n')
        fusion_output.write("\t".join(["ID","Fusion_type", "Flags", "ENSEMBL IDS", "5'_gene", "5'_BND","5'_CDS length", "5' Original_CDS_length","3'_gene", "3'_BND","3'_CDS length", "3' Original_CDS_length"])+"\n")
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

            print(record.ID)
            #Gather all ENSEMBL information on genes that overlap with the BND
            breakend1_annotation=EnsemblAnnotation(CHROM1, POS1)
            if not breakend1_annotation:
                continue
            breakend2_annotation=EnsemblAnnotation(CHROM2, POS2)
            if not breakend2_annotation:
                continue
            #Use requested information to calculate BND specific features
            breakend1_info=BreakendAnnotation(CHROM1, POS1, POS1_ORIENTATION, breakend1_annotation)
            breakend2_info=BreakendAnnotation(CHROM2, POS2, POS2_ORIENTATION, breakend2_annotation)

            #Cross-compare all BND1 hits against all BND2 hits, determine correct fusions and produce output
            Fusions, Vcf_fusion_info=BreakpointAnnotation(record, breakend1_info, breakend2_info, POS1_ORIENTATION, POS2_ORIENTATION, info_output)

            #Produce output
            for fusion in Fusions:
                fusion_output.write(fusion)

            if len(Vcf_fusion_info)>0:
                record.INFO["FUSION"]=Vcf_fusion_info
            VCF_WRITER.write_record(record)

# Gather all information about protein coding genes that overlap with a break-end, independent of the location of the breakend in the gene
def EnsemblAnnotation(CHROM, POS):
    SERVER='http://grch37.rest.ensembl.org'
    ENDPOINT="/overlap/region/human/"+str(CHROM)+":"+str(POS)+"-"+str(POS)
    HEADERS={"Content-Type" : "application/json"}
    PARAMS={"feature": "transcript"}
    genes_data=EnsemblRestClient.perform_rest_action(SERVER, ENDPOINT, HEADERS, PARAMS)

    transcript_ccds={}
    UNIQUE_GENES=[]
    for hit in genes_data:
        if hit["biotype"]=="protein_coding":
            if "ccdsid" in hit:
                transcript_ccds[hit["id"]]=hit["ccdsid"]
            else:
                transcript_ccds[hit["id"]]=None

            if hit["Parent"] not in UNIQUE_GENES:
                UNIQUE_GENES.append(hit["Parent"])

    HITS=[]
    for GENE_ID in UNIQUE_GENES:
        #GENE_ID=hit
        ENDPOINT="/lookup/id/"+str(GENE_ID)
        PARAMS={"expand": "1"}
        gene_info=EnsemblRestClient.perform_rest_action(SERVER, ENDPOINT, HEADERS, PARAMS)
        if gene_info["biotype"]=="protein_coding":
            INFO={}
            INFO["Gene_id"]=gene_info["id"]
            INFO["Gene_name"]=gene_info["display_name"]
            INFO["Strand"]=gene_info["strand"]
            INFO["Gene_start"]=gene_info["start"]
            INFO["Gene_end"]=gene_info["end"]
            INFO["Chromosome"]=CHROM
            INFO["Flags"]=[]

            if "description" in gene_info:
                if "readthrough" in gene_info["description"]:
                    INFO["Flags"].append("Readthrough")

            ##### FLAG for CTC-...  and RP..... proteins (Often not well characterized or readthrough genes)


            for transcript in gene_info["Transcript"]:
                if transcript["is_canonical"]==1:
                    if transcript["id"] in transcript_ccds:
                        if transcript_ccds[transcript["id"]] is None:
                            INFO["Flags"].append("No CCDS")
                    LENGTH_CDS=0
                    CDS=False
                    INFO["Transcript_id"]=transcript["id"]
                    INFO["Biotype"]=transcript["biotype"]
                    INFO["Transcript_start"]=transcript["start"]
                    INFO["Transcript_end"]=transcript["end"]
                    INFO["Total_exons"]=len(transcript["Exon"])
                    INFO["Original_CDS_length"]=(transcript["Translation"]["length"]*3)+3
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
                            EXON_INFO["Start_phase"]="-1"                                         #ADD SOMETHING ELSE SO AN ERROR IS GIVEN IF USED
                            EXON_INFO["End_phase"]="-1"
                            EXON_INFO["CDS_length"]=0
                        PHASE=EXON_INFO["End_phase"]
                        LENGTH_CDS+=EXON_INFO["CDS_length"]
                        INFO["Exons"].append(EXON_INFO)

                        ### INTRONS
                        if rank<len(transcript["Exon"])-1:
                            INTRON_INFO={}
                            INTRON_INFO["Type"]="intron"
                            INTRON_INFO["Rank"]=rank+1
                            INTRON_INFO["Phase"]=EXON_INFO["End_phase"]
                            if INFO["Strand"]==1:
                                INTRON_INFO["Start"]=transcript["Exon"][rank]["end"]+1
                                INTRON_INFO["End"]=transcript["Exon"][rank+1]["start"]-1
                            else:
                                INTRON_INFO["Start"]=transcript["Exon"][rank]["start"]-1
                                INTRON_INFO["End"]=transcript["Exon"][rank+1]["end"]+1
                            # INTRON_INFO["Start"]=transcript["Exon"][rank]["end"]
                            # INTRON_INFO["End"]=transcript["Exon"][rank+1]["start"]
                            INTRON_INFO["CDS"]=CDS
                            INFO["Exons"].append(INTRON_INFO)
                    #print(INFO["Gene_id"],LENGTH_CDS, transcript["Translation"]["length"]*3)
                    if LENGTH_CDS-3!=transcript["Translation"]["length"]*3:
                        INFO["Flags"].append("Possible incomplete CDS")             #as currently I am unable to request the given ENSEMBL flags. Bias towards incomplete but bases%3=0
            HITS.append(INFO)
    return HITS

# Use the gene information gathered in EnsemblAnnotation and add specific breakend information (e.g. what exon, what frame is created etc.)
def BreakendAnnotation(CHROM, POS, orientation, Info):
    HITS=[]
    for gene in Info:
        BND_INFO={k:v for k,v in gene.items()}
        if gene["Strand"]==1:
            ORIENTATION=orientation
            if POS<gene["Transcript_start"]:
                BND_INFO["Breakpoint_location"]="before_transcript"
                BND_INFO["Distance from transcript start"]=abs(POS-gene["Transcript_start"])
                BND_INFO["Type"]=None
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
                BND_INFO["Type"]=None
            elif POS>gene["CDS_end"]:
                BND_INFO["Breakpoint_location"]="3'UTR"
        else:
            ORIENTATION= not orientation
            if POS<gene["Transcript_start"]:
                BND_INFO["Breakpoint_location"]="after_transcript"
                BND_INFO["Distance from transcript end"]=abs(POS-gene["Transcript_start"])
                BND_INFO["Type"]=None
            elif POS<gene["CDS_end"]:
                BND_INFO["Breakpoint_location"]="3'UTR"
            elif POS>gene["CDS_end"] and POS<gene["CDS_start"]:
                if abs(POS-gene["CDS_end"])<3:
                    BND_INFO["Breakpoint_location"]="StopCodon"
                elif abs(POS-gene["CDS_end"])>3 and abs(POS-gene["CDS_end"])<10:
                    BND_INFO["Breakpoint_location"]="3'UTR"
                else:
                    BND_INFO["Breakpoint_location"]="CDS"
            elif POS>gene["Transcript_end"]:
                BND_INFO["Breakpoint_location"]="before_transcript"
                BND_INFO["Distance from transcript start"]=abs(POS-gene["Transcript_end"])
                BND_INFO["Type"]=None
            elif POS>gene["CDS_start"]:
                BND_INFO["Breakpoint_location"]="5'UTR"

        if ORIENTATION:
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
                        BND_INFO["Exon_start_phase"]=sequence["Start_phase"]
                        BND_INFO["Exon_end_phase"]=sequence["End_phase"]

                        if BND_INFO["Breakpoint_location"]=="5'UTR" or BND_INFO["Breakpoint_location"]=="3'UTR":
                            BND_INFO["Phase"]=-1                                                               # Perhaps change to non-integer so it crashes if due to some reason calculation happens with this
                        else:
                            if ORIENTATION:
                                if sequence["Contains_start_CDS"]:
                                    BND_INFO["Phase"]=(abs(gene["CDS_start"]-POS)+sequence["Start_phase"])%3
                                    BND_INFO["CDS_length"]=abs(POS-gene["CDS_end"])+1
                                else:
                                    BND_INFO["Phase"]=(abs(sequence["Start"]-POS)+sequence["Start_phase"])%3
                                    BND_INFO["CDS_length"]=abs(POS-sequence["End"])+1
                                if sequence["Contains_end_CDS"]:
                                    BND_INFO["CDS_length"]=abs(POS-gene["CDS_end"])+1

                                for rank in range(idx+1, len(gene["Exons"])):
                                    if gene["Exons"][rank]["Type"]=="exon":
                                        BND_INFO["CDS_length"]+=gene["Exons"][rank]["CDS_length"]
                            else:
                                if sequence["Contains_start_CDS"]:
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
                        if BND_INFO["Order"]=="3'":
                            for rank in range(idx, len(gene["Exons"])):
                                if gene["Exons"][rank]["Type"]=="exon":
                                    BND_INFO["CDS_length"]+=gene["Exons"][rank]["CDS_length"]
                        elif BND_INFO["Order"]=="5'":
                            for rank in range(0, idx):
                                if gene["Exons"][rank]["Type"]=="exon":
                                    BND_INFO["CDS_length"]+=gene["Exons"][rank]["CDS_length"]
                #else: UTR

        HITS.append(BND_INFO)
    return HITS

#Takes the information from two breakend gathered by BreakendAnnotation and what type of fusion is present and if it is in phase
def FusionCheck(Annotation1, Annotation2):
    if Annotation1["Breakpoint_location"]=="CDS" and Annotation1["Type"]==Annotation2["Type"]:
        Annotation1_CDS_length=Annotation1["CDS_length"]
        Annotation2_CDS_length=Annotation2["CDS_length"]
        if Annotation1["Type"]=="exon":
            if Annotation1["Phase"]==Annotation2["Phase"]:
                Fusion_type="exon-exon"
            else:
                Fusion_type="exon-exon (OUT OF FRAME)"
        elif Annotation1["Type"]=="intron":
            if Annotation1["Phase"]==Annotation2["Phase"]:
                Fusion_type="intron-intron"
            else:
                Fusion_type="intron-intron (OUT OF FRAME)"
        else:
            sys.exit("Unknown error in fusion check - 1")
    elif Annotation1["Type"]=="exon" and Annotation2["Type"]=="intron":
        Annotation1_CDS_length=0
        Annotation2_CDS_length=Annotation2["CDS_length"]
        if Annotation1["Order"]=="5'":
            for rank in range(0, (Annotation1["Rank"]*2)-2):
                if Annotation1["Exons"][rank]["Type"]=="exon":
                    Annotation1_CDS_length+=Annotation1["Exons"][rank]["CDS_length"]
            if Annotation1["Exon_start_phase"]==Annotation2["Phase"]:
                Fusion_type="exon-intron"
            else:
                Fusion_type="exon-intron (OUT OF FRAME)"
        elif Annotation1["Order"]=="3'":
            for rank in range((Annotation1["Rank"]*2)-1, len(Annotation1["Exons"])):
                if Annotation1["Exons"][rank]["Type"]=="exon":
                    Annotation1_CDS_length+=Annotation1["Exons"][rank]["CDS_length"]

            if Annotation1["Exon_end_phase"]==Annotation2["Phase"]:
                Fusion_type="intron-exon"
            else:
                Fusion_type="intron-exon (OUT OF FRAME)"
        else:
            sys.exit("Unknown error in fusion check - 2")
    elif Annotation2["Type"]=="exon" and Annotation1["Type"]=="intron":
        Annotation1_CDS_length=Annotation1["CDS_length"]
        Annotation2_CDS_length=0
        if Annotation2["Order"]=="5'":
            for rank in range(0, (Annotation2["Rank"]*2)-2):
                if Annotation2["Exons"][rank]["Type"]=="exon":
                    Annotation2_CDS_length+=Annotation2["Exons"][rank]["CDS_length"]

            if Annotation2["Exon_start_phase"]==Annotation1["Phase"]:
                Fusion_type="exon-intron"
            else:
                Fusion_type="exon-intron (OUT OF FRAME)"

        elif Annotation2["Order"]=="3'":
            for rank in range((Annotation2["Rank"]*2)-1, len(Annotation2["Exons"])):
                if Annotation2["Exons"][rank]["Type"]=="exon":
                    Annotation2_CDS_length+=Annotation2["Exons"][rank]["CDS_length"]
            if Annotation2["Exon_end_phase"]==Annotation1["Phase"]:
                Fusion_type="intron-exon"
            else:
                Fusion_type="intron-exon (OUT OF FRAME)"
        else:
            sys.exit("Unknown error in fusion check - 3")

    elif Annotation1["Breakpoint_location"]=="5'UTR" and Annotation1["Type"]==Annotation2["Type"]:
        Fusion_type="Promoter fusion"
        if Annotation1["Order"]=="5'":
            Annotation1_CDS_length=0
            Annotation2_CDS_length=Annotation2["Original_CDS_length"]
        else:
            Annotation1_CDS_length=Annotation2["Original_CDS_length"]
            Annotation2_CDS_length=0

    elif Annotation1["Breakpoint_location"]=="3'UTR" and Annotation1["Type"]==Annotation2["Type"]:
        Fusion_type="3'UTR fusion"
        if Annotation1["Order"]=="5'":
            Annotation1_CDS_length=Annotation2["Original_CDS_length"]
            Annotation2_CDS_length=0
        else:
            Annotation1_CDS_length=0
            Annotation2_CDS_length=Annotation2["Original_CDS_length"]
    else:
        sys.exit("Unknown fusion in vcf")

    return (Fusion_type, Annotation1_CDS_length, Annotation2_CDS_length)

def BreakpointAnnotation(Record, Breakend1, Breakend2, Orientation1, Orientation2, Output):
    CHROM1=Record.CHROM
    POS1=Record.POS
    CHROM2=Record.ALT[0].chr
    POS2=Record.ALT[0].pos
    Fusion_output=[]
    Vcf_output=[]
    for annotation1 in Breakend1:
        annotation1["BND"]=str(CHROM1)+":"+str(POS1)
        for annotation2 in Breakend2:
            annotation2["BND"]=str(CHROM2)+":"+str(POS2)

            #Add an extra flag based on both fused genes and produce a list of flags
            FLAGS=annotation1["Flags"]+annotation2["Flags"]

            if ((annotation1["Gene_start"]>annotation2["Gene_start"] and annotation1["Gene_start"]<annotation2["Gene_end"] and
                annotation1["Gene_end"]>annotation2["Gene_start"] and annotation1["Gene_end"]<annotation2["Gene_end"]) or
                (annotation2["Gene_start"]>annotation1["Gene_start"] and annotation2["Gene_start"]<annotation1["Gene_end"] and
                annotation2["Gene_end"]>annotation1["Gene_start"] and annotation2["Gene_end"]<annotation1["Gene_end"])):
                    FLAGS.append("Gene-within-Gene")

            if len(FLAGS)==0:
                FLAGS=["None"]

            #Discard fusions of the same gene and discard fusions where fused genes lie on the same strand and both breakends are in both fusion partners
            if (annotation1["Gene_id"]!=annotation2["Gene_id"] and annotation1["Gene_name"]!=annotation2["Gene_name"] and not
                (POS1 > annotation1["Gene_start"] and POS1 < annotation1["Gene_end"] and POS2 > annotation1["Gene_start"] and POS2 < annotation1["Gene_end"] and
                POS1 > annotation2["Gene_start"] and POS1 < annotation2["Gene_end"] and POS2 > annotation2["Gene_start"] and POS2 < annotation2["Gene_end"] and
                annotation1["Strand"]==annotation2["Strand"])):
                if annotation1["Order"]=="5'" and annotation2["Order"]=="3'":
                    FIVE_PRIME_GENE=annotation1
                    THREE_PRIME_GENE=annotation2
                elif annotation2["Order"]=="5'" and annotation1["Order"]=="3'":
                    FIVE_PRIME_GENE=annotation2
                    THREE_PRIME_GENE=annotation1
                else:
                    continue

                if ((Orientation1 and Orientation2 and annotation1["Strand"]!=annotation2["Strand"] and annotation1["Breakpoint_location"]==annotation2["Breakpoint_location"]) or
                    (Orientation1 and not Orientation2 and annotation1["Strand"]==annotation2["Strand"] and annotation1["Breakpoint_location"]==annotation2["Breakpoint_location"]) or
                    (not Orientation1 and Orientation2 and annotation1["Strand"]==annotation2["Strand"] and annotation1["Breakpoint_location"]==annotation2["Breakpoint_location"]) or
                    (not Orientation1 and not Orientation2 and annotation1["Strand"]!=annotation2["Strand"] and annotation1["Breakpoint_location"]==annotation2["Breakpoint_location"])):
                        FUSION_TYPE, annotation1["CDS_length"], annotation2["CDS_length"] = FusionCheck(annotation1, annotation2)

                        Fusion_output.append("\t".join([str(Record.ID), FUSION_TYPE, ";".join(FLAGS), FIVE_PRIME_GENE["Gene_id"]+"-"+THREE_PRIME_GENE["Gene_id"] ,FIVE_PRIME_GENE["Gene_name"], FIVE_PRIME_GENE["BND"],str(FIVE_PRIME_GENE["CDS_length"]), str(FIVE_PRIME_GENE["Original_CDS_length"]),
                        THREE_PRIME_GENE["Gene_name"], THREE_PRIME_GENE["BND"], str(THREE_PRIME_GENE["CDS_length"]), str(THREE_PRIME_GENE["Original_CDS_length"])])+"\n")

                        Vcf_output.append(FIVE_PRIME_GENE["Gene_id"]+"-"+THREE_PRIME_GENE["Gene_id"])

    return (Fusion_output, Vcf_output)

print("Start:", datetime.datetime.now())
VCF_IN=args.vcf
VCF_OUTPUT=args.output
INFO_OUTPUT=args.fusion_output
EnsemblRestClient=EnsemblRestClient()
ParseVcf(VCF_IN, VCF_OUTPUT, INFO_OUTPUT)
print("End:", datetime.datetime.now())