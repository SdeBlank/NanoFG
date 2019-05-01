#!/usr/bin/python

import argparse
import datetime
import vcf as pyvcf
from EnsemblRestClient import EnsemblRestClient
import sys
import nltk
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib import gridspec
from matplotlib.backends.backend_pdf import PdfPages
from difflib import SequenceMatcher
import copy
import re

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('-v', '--vcf', type=str, help='Input NanoSV vcf file', required=True)
parser.add_argument('-fo', '--fusion_output', type=str, help='Fusion gene info output file', required=True)
parser.add_argument('-o', '--output', type=str, help='Fusion gene annotated vcf file', required=True)
parser.add_argument('-p', '--pdf', type=str, help='Fusion gene pdf file', required=True)

args = parser.parse_args()

def parse_vcf(vcf, vcf_output, info_output, pdf):
    with open(vcf, "r") as vcf, open(vcf_output, "w") as vcf_output, open(info_output, "w") as fusion_output, PdfPages(pdf) as output_pdf:
        vcf_reader=pyvcf.Reader(vcf)
        vcf_reader.infos['FUSION']=pyvcf.parser._Info('FUSION', ".", "String", "Gene names of the fused genes reported if present", "NanoSV", "X")
        vcf_writer=pyvcf.Writer(vcf_output, vcf_reader, lineterminator='\n')
        fusion_output.write("\t".join(["ID","Fusion_type", "Flags", "ENSEMBL_IDS", "5'_gene", "5'_Breakpoint_location" ,"5'_BND", "5'_CDS_length", "5'_Original_CDS_length","3'_gene", "3'_Breakpoint_location", "3'_BND","3'_CDS_length", "3'_Original_CDS_length"])+"\n")
        for record in vcf_reader:
            chrom1=record.CHROM
            pos1=record.POS
            chrom2=record.ALT[0].chr
            pos2=record.ALT[0].pos
            pos1_orientation=record.ALT[0].orientation
            pos2_orientation=record.ALT[0].remoteOrientation

            print(record.ID)
            #Gather all ENSEMBL information on genes that overlap with the BND
            breakend1_annotation=ensembl_annotation(chrom1, pos1)
            if not breakend1_annotation:
                vcf_writer.write_record(record)
                continue
            breakend2_annotation=ensembl_annotation(chrom2, pos2)
            if not breakend2_annotation:
                vcf_writer.write_record(record)
                continue

            #Use requested information to calculate BND specific features
            breakend1_info=breakend_annotation(chrom1, pos1, pos1_orientation, breakend1_annotation)
            breakend2_info=breakend_annotation(chrom2, pos2, pos2_orientation, breakend2_annotation)

            #Cross-compare all BND1 hits against all BND2 hits, determine correct fusions and produce output
            fusions, vcf_fusion_info=breakpoint_annotation(record, breakend1_info, breakend2_info, pos1_orientation, pos2_orientation, info_output)

            #Produce output
            for fusion in fusions:
                fusion_output.write("\t".join([str(record.ID), fusion["Fusion_type"], ";".join(fusion["Flags"]), fusion["5'"]["Gene_id"]+"-"+fusion["5'"]["Gene_id"] ,fusion["5'"]["Gene_name"],
                fusion["5'"]["Type"]+" "+str(fusion["5'"]["Rank"])+"-"+str(fusion["5'"]["Rank"]+1), fusion["5'"]["BND"], str(fusion["5'"]["CDS_length"]), str(fusion["5'"]["Original_CDS_length"]),
                fusion["3'"]["Gene_name"], fusion["3'"]["Type"]+" "+str(fusion["3'"]["Rank"])+"-"+str(fusion["3'"]["Rank"]+1), fusion["3'"]["BND"], str(fusion["3'"]["CDS_length"]),
                str(fusion["3'"]["Original_CDS_length"])])+"\n")
                visualisation(fusion, record, output_pdf)

            if len(vcf_fusion_info)>0:
                record.INFO["FUSION"]=vcf_fusion_info
            vcf_writer.write_record(record)

# Gather all information about protein coding genes that overlap with a break-end, independent of the location of the breakend in the gene
def ensembl_annotation(CHROM, POS):
    server='http://grch37.rest.ensembl.org'
    endpoint="/overlap/region/human/"+str(CHROM)+":"+str(POS)+"-"+str(POS)
    headers={"Content-Type" : "application/json"}
    params={"feature": "transcript"}
    genes_data=EnsemblRestClient.perform_rest_action(server, endpoint, headers, params)

    transcript_ccds={}
    unqiue_genes=[]
    for hit in genes_data:
        if hit["biotype"]=="protein_coding":
            if "ccdsid" in hit:
                transcript_ccds[hit["id"]]=hit["ccdsid"]
            else:
                transcript_ccds[hit["id"]]=None

            if hit["Parent"] not in unqiue_genes:
                unqiue_genes.append(hit["Parent"])

    HITS=[]
    for gene_id in unqiue_genes:
        #gene_id=hit
        endpoint="/lookup/id/"+str(gene_id)
        params={"expand": "1"}
        gene_info=EnsemblRestClient.perform_rest_action(server, endpoint, headers, params)
        if gene_info["biotype"]=="protein_coding":
            ensembl_info={}
            ensembl_info["Gene_id"]=gene_info["id"]
            ensembl_info["Gene_name"]=gene_info["display_name"]
            ensembl_info["Strand"]=gene_info["strand"]
            ensembl_info["Gene_start"]=gene_info["start"]
            ensembl_info["Gene_end"]=gene_info["end"]
            ensembl_info["Chromosome"]=CHROM
            ensembl_info["Flags"]=[]

            if "description" in gene_info:
                if "readthrough" in gene_info["description"]:
                    ensembl_info["Flags"].append("fusion-with-readthrough")

            ##### FLAG for CTC-...  and RP..... proteins (Often not well characterized or readthrough genes)


            for transcript in gene_info["Transcript"]:
                if transcript["is_canonical"]==1:
                    if transcript["id"] in transcript_ccds:
                        if transcript_ccds[transcript["id"]] is None:
                            ensembl_info["Flags"].append("No-CCDS")
                    cds_length=0
                    CDS=False
                    ensembl_info["Transcript_id"]=transcript["id"]
                    ensembl_info["Translation_id"]=transcript["Translation"]["id"]
                    ensembl_info["Biotype"]=transcript["biotype"]
                    ensembl_info["Transcript_start"]=transcript["start"]
                    ensembl_info["Transcript_end"]=transcript["end"]
                    ensembl_info["Total_exons"]=len(transcript["Exon"])
                    ensembl_info["Total_exon_length"]=0
                    ensembl_info["Original_CDS_length"]=(transcript["Translation"]["length"]*3)+3
                    if ensembl_info["Strand"]==1:
                        ensembl_info["CDS_start"]=transcript["Translation"]["start"]
                        ensembl_info["CDS_end"]=transcript["Translation"]["end"]
                    else:
                        ensembl_info["CDS_start"]=transcript["Translation"]["end"]
                        ensembl_info["CDS_end"]=transcript["Translation"]["start"]
                    ensembl_info["Exons"]=[]

                    ### EXONS
                    for rank, exon in enumerate(transcript["Exon"]):
                        exon_info={}
                        exon_info["Rank"]=rank+1
                        exon_info["Type"]="exon"
                        CHRON_START=exon["start"]
                        CHRON_END=exon["end"]
                        ensembl_info["Total_exon_length"]+=abs(exon["end"]-exon["start"])
                        if ensembl_info["Strand"]==1:
                            exon_info["Start"]=exon["start"]
                            exon_info["End"]=exon["end"]
                        else:
                            exon_info["Start"]=exon["end"]
                            exon_info["End"]=exon["start"]
                        exon_info["Contains_start_CDS"]=False
                        exon_info["Contains_end_CDS"]=False

                        if CDS:
                            if ensembl_info["CDS_end"]>=CHRON_START and ensembl_info["CDS_end"]<=CHRON_END:
                                exon_info["Contains_end_CDS"]=True
                                exon_info["CDS"]=True
                                exon_info["Start_phase"]=phase
                                #exon_info["End_phase"]=-1
                                exon_info["End_phase"]=0
                                exon_info["CDS_length"]=abs(ensembl_info["CDS_end"]-exon_info["Start"])+1
                                CDS=False
                            else:
                                exon_info["CDS"]=True
                                exon_info["Start_phase"]=phase
                                exon_info["End_phase"]=(abs(exon_info["End"]-exon_info["Start"])+1+phase)%3
                                exon_info["CDS_length"]=abs(exon_info["End"]-exon_info["Start"])+1
                        elif ensembl_info["CDS_start"]>=CHRON_START and ensembl_info["CDS_start"]<=CHRON_END:
                            exon_info["Contains_start_CDS"]=True
                            exon_info["CDS"]=True
                            #exon_info["Start_phase"]=-1
                            exon_info["Start_phase"]=0
                            if ensembl_info["CDS_end"]>=CHRON_START and ensembl_info["CDS_end"]<=CHRON_END:
                                exon_info["Contains_end_CDS"]=True
                                exon_info["End_phase"]=0
                                exon_info["CDS_length"]=abs(ensembl_info["CDS_end"]-ensembl_info["CDS_start"])+1
                            else:
                                exon_info["End_phase"]=(abs(exon_info["End"]-ensembl_info["CDS_start"])+1)%3
                                exon_info["CDS_length"]=abs(exon_info["End"]-ensembl_info["CDS_start"])+1
                            CDS=True

                        else:
                            exon_info["CDS"]=False
                            exon_info["Start_phase"]="-1"
                            exon_info["End_phase"]="-1"
                            exon_info["CDS_length"]=0
                        phase=exon_info["End_phase"]
                        cds_length+=exon_info["CDS_length"]
                        ensembl_info["Exons"].append(exon_info)

                        ### INTRONS
                        if rank<len(transcript["Exon"])-1:
                            intron_info={}
                            intron_info["Type"]="intron"
                            intron_info["Rank"]=rank+1
                            intron_info["Phase"]=exon_info["End_phase"]
                            if ensembl_info["Strand"]==1:
                                intron_info["Start"]=transcript["Exon"][rank]["end"]+1
                                intron_info["End"]=transcript["Exon"][rank+1]["start"]-1
                            else:
                                intron_info["Start"]=transcript["Exon"][rank]["start"]-1
                                intron_info["End"]=transcript["Exon"][rank+1]["end"]+1
                            intron_info["CDS"]=CDS
                            ensembl_info["Exons"].append(intron_info)
                    #print(ensembl_info["Gene_id"],cds_length, transcript["Translation"]["length"]*3)
                    if cds_length-3!=transcript["Translation"]["length"]*3:
                        ensembl_info["Flags"].append("Possible-incomplete-CDS")             #as currently I am unable to request the given ENSEMBL Flags. Bias towards incomplete but bases%3=0
            HITS.append(ensembl_info)
    return HITS

# Use the gene information gathered in ensembl_annotation and add specific breakend information (e.g. what exon, what frame is created etc.)
def breakend_annotation(CHROM, POS, orientation, Info):
    HITS=[]
    for gene in Info:
        BND_INFO={k:v for k,v in gene.items()}
        BND_INFO["Breakend_position"]=POS
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
                # elif abs(POS-gene["CDS_end"])<10:                 #### ADDED pos inside stop codon==different and pos close to translation end gives 3'UTR fusion instead
                #     BND_INFO["Breakpoint_location"]="3'UTR"
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

#Takes the information from two breakend gathered by breakend_annotation and what type of fusion is present and if it is in phase
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

def breakpoint_annotation(Record, Breakend1, Breakend2, Orientation1, Orientation2, Output):
    chrom1=Record.CHROM
    pos1=Record.POS
    chrom2=Record.ALT[0].chr
    pos2=Record.ALT[0].pos
    Fusion_output=[]
    Vcf_output=[]
    for annotation1 in Breakend1:
        annotation1["BND"]=str(chrom1)+":"+str(pos1)
        for annotation2 in Breakend2:
            annotation2["BND"]=str(chrom2)+":"+str(pos2)

            #Add an extra flag based on both fused genes and produce a list of flags
            FLAGS=annotation1["Flags"]+annotation2["Flags"]

            if (((annotation1["Gene_start"]>annotation2["Gene_start"] and annotation1["Gene_start"]<annotation2["Gene_end"] and
                annotation1["Gene_end"]>annotation2["Gene_start"] and annotation1["Gene_end"]<annotation2["Gene_end"]) or
                (annotation2["Gene_start"]>annotation1["Gene_start"] and annotation2["Gene_start"]<annotation1["Gene_end"] and
                annotation2["Gene_end"]>annotation1["Gene_start"] and annotation2["Gene_end"]<annotation1["Gene_end"])) and
                annotation1['Strand']==annotation2['Strand']):

                FLAGS.append("Gene-within-Gene")

            Gene_difference=nltk.edit_distance(annotation1["Gene_name"], annotation2["Gene_name"])
            if Gene_difference/len(annotation1["Gene_name"]) <= 0.4 and Gene_difference/len(annotation2["Gene_name"]) <= 0.4:
                FLAGS.append("Similar-genes")

            #Discard fusions of the same gene and discard fusions where fused genes lie on the same strand and both breakends are in both fusion partners
            if (annotation1["Gene_id"]==annotation2["Gene_id"] or annotation1["Gene_name"]==annotation2["Gene_name"] or
                (pos1 > annotation1["Gene_start"] and pos1 < annotation1["Gene_end"] and pos2 > annotation1["Gene_start"] and pos2 < annotation1["Gene_end"] and
                pos1 > annotation2["Gene_start"] and pos1 < annotation2["Gene_end"] and pos2 > annotation2["Gene_start"] and pos2 < annotation2["Gene_end"] and
                annotation1["Strand"]==annotation2["Strand"])):
                continue
            elif (annotation1["Gene_id"]==annotation2["Gene_id"] or annotation1["Gene_name"]==annotation2["Gene_name"] or
                (pos1 > annotation1["Gene_start"] and pos1 < annotation1["Gene_end"] and pos2 > annotation1["Gene_start"] and pos2 < annotation1["Gene_end"] and
                pos1 > annotation2["Gene_start"] and pos1 < annotation2["Gene_end"] and pos2 > annotation2["Gene_start"] and pos2 < annotation2["Gene_end"] and
                annotation1["Strand"]!=annotation2["Strand"])):
                continue
            else:
                if annotation1["Order"]=="5'" and annotation2["Order"]=="3'":
                    five_prime_gene=annotation1
                    three_prime_gene=annotation2
                elif annotation2["Order"]=="5'" and annotation1["Order"]=="3'":
                    five_prime_gene=annotation2
                    three_prime_gene=annotation1
                else:
                    continue

                if ((Orientation1 and Orientation2 and annotation1["Strand"]!=annotation2["Strand"] and annotation1["Breakpoint_location"]==annotation2["Breakpoint_location"]) or
                    (Orientation1 and not Orientation2 and annotation1["Strand"]==annotation2["Strand"] and annotation1["Breakpoint_location"]==annotation2["Breakpoint_location"]) or
                    (not Orientation1 and Orientation2 and annotation1["Strand"]==annotation2["Strand"] and annotation1["Breakpoint_location"]==annotation2["Breakpoint_location"]) or
                    (not Orientation1 and not Orientation2 and annotation1["Strand"]!=annotation2["Strand"] and annotation1["Breakpoint_location"]==annotation2["Breakpoint_location"])):
                        FUSION_TYPE, annotation1["CDS_length"], annotation2["CDS_length"] = FusionCheck(annotation1, annotation2)
                        if len(FLAGS)==0:
                            FLAGS=["None"]
                        Fusion_output.append({"5'": five_prime_gene, "3'": three_prime_gene, "Fusion_type": FUSION_TYPE, "Flags":FLAGS})
                        Vcf_output.append(five_prime_gene["Gene_id"]+"-"+three_prime_gene["Gene_id"])

    return (Fusion_output, Vcf_output)

def get_domains(protein_id, coding_exons, relative_length, intron_relative_length):
    server='http://grch37.rest.ensembl.org'
    endpoint="/overlap/translation/"+protein_id
    headers={"Content-Type" : "application/json"}
    ensembl_domains=EnsemblRestClient.perform_rest_action(server, endpoint, headers)

    domain_layers=[[] for n in range(12)]
    cds_start=coding_exons["rel_pos"][0][0]
    cds_end=coding_exons["rel_pos"][-1][1]
    used_domains={}
    to_be_used_domains={}
    for domain in ensembl_domains:
        next_domain=False
        if domain["description"]=="":
            continue

        for index, exon in enumerate(coding_exons["pos"]):
            nr_of_introns=index
            if domain["start"]*3>=exon[0] and domain["start"]*3<=exon[1]:
                pos_start=coding_exons["rel_pos"][index][0]+abs(domain["start"]*3-exon[0])*relative_length
            if domain["end"]*3>=exon[0] and domain["end"]*3<=exon[1]:
                pos_end=coding_exons["rel_pos"][index][0]+abs(domain["end"]*3-exon[0])*relative_length

        text_begin=pos_start+0.5*abs(pos_end-pos_start)-(0.5*len(domain["description"])/3)
        text_end=pos_start+0.5*abs(pos_end-pos_start)+(0.5*len(domain["description"])/3)

        #print(domain["description"])
        for layer, domains in enumerate(domain_layers):
            next_layer=False
            nr_of_items=len(domains)
            if nr_of_items==0:
                domain_layers[layer].append([pos_start, pos_end, domain["description"],0.075*layer, domain["start"], domain["end"], text_begin, text_end])
                break
            else:
                index=0
                while index < nr_of_items:
                    pos=domains[index]
                    if domain["description"]==pos[2] or SequenceMatcher(None, domain["description"], pos[2]).ratio()>0.75:
                        if abs(domain["start"]-pos[4])<=50 and abs(domain["end"]-pos[5])<=50 or (domain["start"]>pos[4] and domain["end"]<pos[5]):
                            if domain["start"]<pos[4]:
                                domain_layers[layer][index][0]=pos_start
                                domain_layers[layer][index][4]=domain["start"]
                            if domain["end"]>pos[5]:
                                domain_layers[layer][index][1]=pos_end
                                domain_layers[layer][index][5]=domain["end"]
                            next_layer=True
                            next_domain=True
                            break
                    if (abs(pos[6]-text_begin)<=20 or abs(pos[6]-text_end)<=20 or abs(pos[7]-text_begin)<=20 or abs(pos[7]-text_end)<=20 or
                    abs(pos[0]-pos_start)<=10 or abs(pos[0]-pos_end)<=10 or abs(pos[1]-pos_start)<=10 or abs(pos[1]-pos_end)<=10 or
                    (text_begin>pos[6] and text_begin<pos[7]) or (text_end>pos[6] and text_end<pos[7]) or
                    (pos_start>pos[0] and pos_start<pos[1]) or (pos_end>pos[0] and pos_end<pos[1]) or
                    (text_begin<pos[6] and text_end>pos[7])):
                        next_layer=True
                    index+=1
                if next_domain:
                    break
                if not next_layer:
                    domain_layers[layer].append([pos_start, pos_end, domain["description"],0.075*layer, domain["start"], domain["end"], text_begin, text_end])
                    break

    for index, domains in enumerate(domain_layers):
        for domain in domains:
            if domain[2] not in to_be_used_domains:
                to_be_used_domains[domain[2]]=[domain]
            else:
                to_be_used_domains[domain[2]].append(domain)
    return to_be_used_domains


def visualisation(annotated_breakpoints, Record, pdf):
    rowplots = 7
    colplots = 3
    xmin = 0
    xmax = 100
    yvalue = .5
    cds_height=0.2
    non_cds_height=0.1
    Intron_size=1

    def domains_plot( domains, n, color, order ):
        ax = fig.add_subplot(gs[n, 1])
        max_y=0
        for domain in domains:
            yline=0.1+domain[3]
            ytext=yline+0.015
            if order=="3'":
                ytext=yline-0.015
            x0 = domain[0]/float(xmax)
            x1 = domain[1]/float(xmax)
            ax.axhline(y=yline, xmin=x0, xmax=x1,color=color, linewidth=1)
            ax.text((domain[0]+domain[1])/2,ytext,domain[2],horizontalalignment='center',verticalalignment='bottom', size=6)
        plt.xlim( (xmin,xmax) )
        if order=="3'":
            matplotlib.axes.Axes.invert_yaxis(ax)
        plt.axis('off')

    def splicing_plot(x1, x2, n, color):
        y1 = yvalue
        y2 = yvalue+.5
        ax = fig.add_subplot(gs[n, 1])
        ax.plot([x1, x2], [y1, y2], color, linestyle="--", zorder=3, linewidth=1)
        plt.ylim( (yvalue-0.5, yvalue+0.5) )
        plt.xlim( (xmin,xmax) )
        plt.axis('off')

    def introns_plot( start, end, n, color ):
        x0 = start/float(xmax)
        x1 = end/float(xmax)
        ax = fig.add_subplot(gs[n, 1])
        ax.axhline(y=yvalue, xmin=x0, xmax=x1, color=color, zorder=1, linewidth=1)

    def exons_plot( exons, n ):
        ax = fig.add_subplot(gs[n, 1])
        for exon in exons:
            if exon[3]:
                y0 = yvalue-cds_height
                y1 = 2*cds_height
            else:
                y0 = yvalue-non_cds_height
                y1 = 2*non_cds_height

            if exon[4]=="5'":
                color="blue"
            elif exon[4]=="3'":
                color="red"

            skipped=exon[5]
            if skipped:
                ax.add_patch(patches.Rectangle((exon[0], y0),exon[1],y1,color=color, zorder=2, alpha=0.5) )
            else:
                ax.add_patch(patches.Rectangle((exon[0], y0),exon[1],y1,color=color, zorder=2) )
            if exon[2] is not None:
                ax.text(exon[0]+(exon[1]/2),yvalue,exon[2],horizontalalignment='center', verticalalignment='center', color='white', size=6)
        plt.xlim( (xmin,xmax) )
        plt.axis('off')

    def breakpoint_plot(breakpoint_x, n):
        ax = fig.add_subplot(gs[n, 1])
        ax.axvline(x=breakpoint_x, ymin=yvalue-0.5, ymax=yvalue+0.5, color="black")
        plt.xlim( (xmin,xmax) )
        plt.axis('off')

    def gene_visualisation(gene_info, order):
        exon_list=[]
        domain_list=[]

        breaks=gene_info["Total_exons"]-1
        part=(xmax-breaks)/gene_info["Total_exon_length"]

        cds_exons={"rel_pos": [], "pos":[], "introns":0}
        space=0
        exon_length=0
        exon_start=0

        for index, exon in enumerate(gene_info["Exons"]):
            gene_info["Exons"][index]["Origin"]=order
            origin=order
            if gene_info["Type"]==exon["Type"] and gene_info["Rank"]==exon["Rank"]:
                if gene_info["Type"]=="intron":
                    breakpoint=exon_start+exon_length+0.5*space
                elif gene_info["Type"]=="exon":
                    breakpoint=exon_start+exon_length+space+(abs(gene_info["Breakend_position"]-exon["Start"])+1)*part
            if exon["Type"]=="exon":
                exon_start=exon_start+exon_length+space
                exon_length=(abs(exon["End"]-exon["Start"])+1)*part
                label=exon["Rank"]
                if exon["CDS"]:
                    if exon["Contains_start_CDS"]:
                        cds_end=abs(exon["End"]-gene_info["CDS_start"])+1
                        cds_exons["rel_pos"].append([exon_start+(abs(gene_info["CDS_start"]-exon["Start"])+1)*part, exon_start+exon_length])
                        cds_exons["pos"].append([0, abs(exon["End"]-gene_info["CDS_start"])+1])
                        cds_exons["introns"]+=1
                        CDS_start=exon_start+abs(gene_info["CDS_start"]-exon["Start"])*part
                        coding_length=abs((exon["End"]-gene_info["CDS_start"])+1)*part
                        exon_list.append((exon_start, exon_length, None, False, origin, False))
                        exon_list.append((CDS_start, coding_length, label, True, origin, False))
                    elif exon["Contains_end_CDS"]:
                        cds_exons["rel_pos"].append([exon_start, exon_start+exon_length])
                        cds_exons["pos"].append([cds_end+1, abs(gene_info["CDS_end"]-exon["Start"])+1+cds_end])
                        coding_length=abs(exon["Start"]-gene_info["CDS_end"])*part
                        exon_list.append((exon_start, exon_length, None, False, origin, False))
                        exon_list.append((exon_start, coding_length, label, True, origin, False))
                    else:
                        cds_exons["rel_pos"].append([exon_start, exon_start+exon_length])
                        cds_exons["pos"].append([cds_end+1, abs(exon["End"]-exon["Start"])+1+cds_end])
                        cds_exons["introns"]+=1
                        cds_end=abs(exon["End"]-exon["Start"])+1+cds_end
                        exon_list.append((exon_start, exon_length, label, True, origin, False))
                else:
                    exon_list.append((exon_start, exon_length, label, False, origin, False))
                space=Intron_size

        ensembl_domains=get_domains(gene_info["Translation_id"], cds_exons, part, space)
        for info in ensembl_domains.values():
            for positions in info:
                domain_list.append(positions)

        return (exon_list, domain_list, breakpoint)

############################################################################# Plot Title and flags

    fig = plt.figure(figsize=(11.7,8.3), edgecolor="black")
    gs  = gridspec.GridSpec(rowplots, colplots, hspace=0,height_ratios=[0.5,2, 0.5 ,0.5, 0.5, 2,1], width_ratios=[1,99, 0.1])

    ax = fig.add_subplot(gs[0,:])
    if "OUT OF FRAME" in annotated_breakpoints["Fusion_type"]:
        ax.text(0.5*xmax, 0, Record.ID+". " + annotated_breakpoints["5'"]["Gene_name"]+"-"+annotated_breakpoints["3'"]["Gene_name"], horizontalalignment='center',verticalalignment='bottom', size=16, fontweight='bold')
        ax.text(0.5*xmax, 0.35, "(OUT OF FRAME)", horizontalalignment='center',verticalalignment='bottom', size=14, color="red")
        ax.text(0.5*xmax, 0.6, annotated_breakpoints["5'"]["Gene_id"]+"-"+annotated_breakpoints["3'"]["Gene_id"], horizontalalignment='center',verticalalignment='bottom', size=9)
        ax.axhline(y=0.7, xmin=0, xmax=100, color="black", linewidth=1)
    else:
        ax.text(0.5*xmax, 0, Record.ID+". " + annotated_breakpoints["5'"]["Gene_name"]+"-"+annotated_breakpoints["3'"]["Gene_name"], horizontalalignment='center',verticalalignment='bottom', size=16, fontweight='bold')
        ax.text(0.5*xmax, 0.3, annotated_breakpoints["5'"]["Gene_id"]+"-"+annotated_breakpoints["3'"]["Gene_id"], horizontalalignment='center',verticalalignment='bottom', size=9)
        ax.axhline(y=0.4, xmin=0, xmax=100, color="black", linewidth=1)

    if annotated_breakpoints["Flags"]!=["None"]:
        ax.text(0, 1, "; ".join(annotated_breakpoints["Flags"]), horizontalalignment='left',verticalalignment='center', size=12, color="red")

    ax.axhline(y=1.4, xmin=0, xmax=100, color="black", linewidth=1, linestyle="--")

    plt.xlim( (xmin,xmax) )
    plt.axis('off')
    matplotlib.axes.Axes.invert_yaxis(ax)


############################################################################# DONOR GENE VISUALISATION
    donor_exons, donor_domains, breakpoint = gene_visualisation(annotated_breakpoints["5'"], "5'")

    ax = fig.add_subplot(gs[2, 0])
    ax.text(0, yvalue, annotated_breakpoints["5'"]["Gene_name"], horizontalalignment='center',verticalalignment='bottom', size=9)
    ax.text(0, yvalue-0.1, "("+str(int(annotated_breakpoints["5'"]["Original_CDS_length"]/3))+" aa)", horizontalalignment='center',verticalalignment='center', size=7)
    plt.axis('off')

    ax = fig.add_subplot(gs[2, 2])
    if annotated_breakpoints["5'"]["Strand"]==1:
        ax.text(0, yvalue, "+", horizontalalignment='center',verticalalignment='center', size=15)
    else:
        ax.text(0, yvalue, "_", horizontalalignment='center',verticalalignment='bottom', size=15)
    plt.axis('off')

    domains_plot( donor_domains, 1 , "blue", "5'")
    exons_plot( donor_exons , 2 )
    introns_plot( 0, 100, 2, "blue" )

    breakpoint_plot(breakpoint, 2)

############################################################################# ACCEPTOR GENE VISUALISATION

    acceptor_exons, acceptor_domains, breakpoint = gene_visualisation(annotated_breakpoints["3'"], "3'")

    Protein_info=annotated_breakpoints["3'"]["Gene_name"]+"\n"+"("+str(int(annotated_breakpoints["3'"]["Original_CDS_length"]/3))+" aa)"
    ax = fig.add_subplot(gs[4, 0])
    ax.text(0, yvalue, annotated_breakpoints["3'"]["Gene_name"], horizontalalignment='center',verticalalignment='bottom', size=9)
    ax.text(0, yvalue-0.1, "("+str(int(annotated_breakpoints["3'"]["Original_CDS_length"]/3))+" aa)", horizontalalignment='center',verticalalignment='center', size=7)
    plt.axis('off')

    ax = fig.add_subplot(gs[4, 2])
    if annotated_breakpoints["3'"]["Strand"]==1:
        ax.text(0, yvalue, "+", horizontalalignment='center',verticalalignment='center', size=15)
    else:
        ax.text(0, yvalue, "_", horizontalalignment='center',verticalalignment='bottom', size=15)
    plt.axis('off')

    introns_plot( 0, 100, 4, "red" )
    exons_plot( acceptor_exons, 4 )
    breakpoint_plot(breakpoint, 4)
    domains_plot( acceptor_domains, 5 , "red", "3'")

############################################################################# FUSION GENE VISUALISATION

    fused_exons = []
    print(annotated_breakpoints["Fusion_type"])
    if "exon-exon" in annotated_breakpoints["Fusion_type"]:
        fused_protein=annotated_breakpoints["5'"]["Exons"][:annotated_breakpoints["5'"]["Rank"]*2-1]+annotated_breakpoints["3'"]["Exons"][annotated_breakpoints["3'"]["Rank"]*2-2:]
        fused_length=0
        fused_nr_of_exons=0
        for index, exon in enumerate(fused_protein):
            if exon["Type"]=="exon":
                fused_nr_of_exons+=1
                origin=exon["Origin"]
                if exon["Origin"]=="5'" and exon["Rank"]==annotated_breakpoints["5'"]["Rank"]:
                    fused_protein[index]["End"]=annotated_breakpoints["5'"]["Breakend_position"]
                    exon["End"]=annotated_breakpoints["5'"]["Breakend_position"]
                if exon["Origin"]=="3'" and exon["Rank"]==annotated_breakpoints["3'"]["Rank"]:
                    fused_protein[index]["Start"]=annotated_breakpoints["3'"]["Breakend_position"]
                    exon["Start"]=annotated_breakpoints["3'"]["Breakend_position"]
                fused_length+=abs(exon["End"]-exon["Start"])
        breaks=fused_nr_of_exons-2
        part=(xmax-breaks)/fused_length
        space=0
        exon_length=0
        exon_start=0
        origin="5'"
        for exon in fused_protein:
            if exon["Type"]=="exon":
                if exon["Origin"]!=origin:
                    BND=exon_start+exon_length
                    exon_start=exon_start+exon_length
                else:
                    exon_start=exon_start+exon_length+space
                origin=exon["Origin"]
                exon_length=abs(exon["End"]-exon["Start"])*part
                label=exon["Rank"]
                if exon["CDS"]:
                    if exon["Contains_start_CDS"] and exon["Origin"]=="5'":
                        fused_exons.append((exon_start, exon_length, None, False, origin, False))
                        CDS_start=abs(annotated_breakpoints["5'"]["CDS_start"]-exon["Start"])*part
                        coding_length=abs(exon["End"]-annotated_breakpoints["5'"]["CDS_start"])*part
                        fused_exons.append((CDS_start, coding_length, label, True,origin, False))
                    elif exon["Contains_end_CDS"] and exon["Origin"]=="3'":
                        fused_exons.append((exon_start, exon_length, None, False, origin, False))
                        coding_length=abs(exon["Start"]-annotated_breakpoints["3'"]["CDS_end"])*part
                    else:
                        fused_exons.append((exon_start, exon_length, label, True, origin, False))
                else:
                    fused_exons.append((exon_start, exon_length, label, False, origin, False))
                space=Intron_size
    elif "intron-intron" in annotated_breakpoints["Fusion_type"]:
        fused_protein=annotated_breakpoints["5'"]["Exons"][:annotated_breakpoints["5'"]["Rank"]*2-1]+annotated_breakpoints["3'"]["Exons"][annotated_breakpoints["3'"]["Rank"]*2-1:]
        fused_length=0
        fused_nr_of_exons=0
        for exon in fused_protein:
            if exon["Type"]=="exon":
                fused_length+=abs(exon["End"]-exon["Start"])
                fused_nr_of_exons+=1
        breaks=fused_nr_of_exons-1
        part=(xmax-breaks)/fused_length
        space=0
        exon_length=0
        exon_start=0
        origin="5'"
        for exon in fused_protein:
            if exon["Type"]=="exon":
                if exon["Origin"]!=origin:
                    BND=exon_start+exon_length+(space/2)
                origin=exon["Origin"]
                exon_start=exon_start+exon_length+space
                exon_length=abs(exon["End"]-exon["Start"])*part
                label=exon["Rank"]
                if exon["CDS"]:
                    if exon["Contains_start_CDS"]:
                        fused_exons.append((exon_start, exon_length, None, False, origin, False))
                        CDS_start=abs(annotated_breakpoints["5'"]["CDS_start"]-exon["Start"])*part
                        coding_length=abs(exon["End"]-annotated_breakpoints["5'"]["CDS_start"])*part
                        fused_exons.append((CDS_start, coding_length, label, True,origin, False))
                    elif exon["Contains_end_CDS"]:
                        fused_exons.append((exon_start, exon_length, None, False, origin, False))
                        coding_length=abs(exon["Start"]-annotated_breakpoints["3'"]["CDS_end"])*part
                        fused_exons.append((exon_start, coding_length, label, True, origin, False))
                    else:
                        fused_exons.append((exon_start, exon_length, label, True, origin, False))
                else:
                    fused_exons.append((exon_start, exon_length, label, False, origin, False))
                space=Intron_size

    elif "exon" in annotated_breakpoints["Fusion_type"] and "intron" in annotated_breakpoints["Fusion_type"]:
        if "exon-intron" in annotated_breakpoints["Fusion_type"]:
            fused_protein=copy.deepcopy(annotated_breakpoints["5'"]["Exons"][:annotated_breakpoints["5'"]["Rank"]*2]+annotated_breakpoints["3'"]["Exons"][annotated_breakpoints["3'"]["Rank"]*2-1:])
        elif "intron-exon" in annotated_breakpoints["Fusion_type"]:
            fused_protein=copy.deepcopy(annotated_breakpoints["5'"]["Exons"][:annotated_breakpoints["5'"]["Rank"]*2]+annotated_breakpoints["3'"]["Exons"][annotated_breakpoints["3'"]["Rank"]*2-2:])
        fused_length=0
        fused_nr_of_exons=0
        for index, exon in enumerate(fused_protein):
            if exon["Type"]=="exon":
                fused_nr_of_exons+=1
                origin=exon["Origin"]
                if exon["Origin"]=="5'" and exon["Rank"]==annotated_breakpoints["5'"]["Rank"] and "exon-intron" in annotated_breakpoints["Fusion_type"]:
                    fused_protein[index]["End"]=annotated_breakpoints["5'"]["Breakend_position"]
                    exon["End"]=annotated_breakpoints["5'"]["Breakend_position"]
                if exon["Origin"]=="3'" and exon["Rank"]==annotated_breakpoints["3'"]["Rank"] and "intron-exon" in annotated_breakpoints["Fusion_type"]:
                    fused_protein[index]["Start"]=annotated_breakpoints["3'"]["Breakend_position"]
                    exon["Start"]=annotated_breakpoints["3'"]["Breakend_position"]
                fused_length+=abs(exon["End"]-exon["Start"])
        breaks=fused_nr_of_exons-1
        part=(xmax-breaks)/fused_length
        space=0
        exon_length=0
        exon_start=0
        origin="5'"
        for index, exon in enumerate(fused_protein):
            if exon["Type"]=="exon":
                if exon["Origin"]!=origin:
                    if "exon-intron" in annotated_breakpoints["Fusion_type"]:
                        BND=exon_start+exon_length
                        splicing_plot(exon_start-space, exon_start+0.5*exon_length, 3, "black")
                        splicing_plot(exon_start+exon_length+space, exon_start+0.5*exon_length, 3, "black")
                    else:
                        BND=exon_start+exon_length+space
                        splicing_plot(exon_start+exon_length, exon_start+exon_length+space+0.5*abs(exon["End"]-exon["Start"])*part, 3, "black")
                        splicing_plot(exon_start+exon_length+2*space+abs(exon["End"]-exon["Start"])*part, exon_start+exon_length+space+0.5*abs(exon["End"]-exon["Start"])*part, 3, "black")
                origin=exon["Origin"]
                exon_start=exon_start+exon_length+space
                exon_length=abs(exon["End"]-exon["Start"])*part
                label=exon["Rank"]
                if exon["CDS"]:
                    if exon["Contains_start_CDS"] and exon["Origin"]=="5'":
                        fused_exons.append((exon_start, exon_length, None, False, origin, False))
                        CDS_start=exon_start+abs(annotated_breakpoints["5'"]["CDS_start"]-exon["Start"])*part
                        coding_length=abs(exon["End"]-annotated_breakpoints["5'"]["CDS_start"])*part
                        fused_exons.append((CDS_start, coding_length, label, True,origin, False))
                    elif exon["Contains_end_CDS"] and exon["Origin"]=="3'":
                        fused_exons.append((exon_start, exon_length, None, False, origin, False))
                        coding_length=abs(exon["Start"]-annotated_breakpoints["3'"]["CDS_end"])*part
                        fused_exons.append((exon_start, coding_length, label, True, origin, False))
                    else:
                        fused_exons.append((exon_start, exon_length, label, True, origin, False))
                else:
                    fused_exons.append((exon_start, exon_length, label, False, origin, False))
                space=Intron_size

    Protein_info=str(int((annotated_breakpoints["5'"]["CDS_length"]+annotated_breakpoints["3'"]["CDS_length"])/3))+" aa"
    ax = fig.add_subplot(gs[3, 0])
    ax.text(0, yvalue, annotated_breakpoints["5'"]["Gene_name"]+"-"+annotated_breakpoints["3'"]["Gene_name"], horizontalalignment='center',verticalalignment='bottom', size=9)
    ax.text(0, yvalue-0.1, "("+Protein_info+")", horizontalalignment='center',verticalalignment='center', size=7)
    plt.axis('off')

    exons_plot( fused_exons, 3 )
    introns_plot( 0, BND, 3, "blue")
    introns_plot( BND, 100, 3, "red" )

    ############################################################################# Additional info visualisation
    ax = fig.add_subplot(gs[6, :])
    ax.axhline(y=0, xmin=0, xmax=100, color="black", linewidth=1, linestyle="--")
    SV_ID=re.findall("^\d+", Record.INFO["ALT_READ_IDS"][0])[0]

    ax.text(0, 0.2, "Original SV_ID:", horizontalalignment='left',verticalalignment='top', size=10, fontweight='bold')
    ax.text(0, 0.5, SV_ID, horizontalalignment='left',verticalalignment='center', size=9)
    ax.text(0.2, 0.2, "Fusion type:", horizontalalignment='left',verticalalignment='top', size=10, fontweight='bold')
    ax.text(0.2, 0.5, annotated_breakpoints["Fusion_type"].split(" ")[0], horizontalalignment='left',verticalalignment='center', size=9)
    ax.text(0.4, 0.2, "5' Breakpoint:", horizontalalignment='left',verticalalignment='top', size=10, fontweight='bold')
    ax.text(0.4, 0.5, annotated_breakpoints["5'"]["BND"], horizontalalignment='left',verticalalignment='center', size=9)
    ax.text(0.6, 0.2, "3' Breakpoint:", horizontalalignment='left',verticalalignment='top', size=10, fontweight='bold')
    ax.text(0.6, 0.5, annotated_breakpoints["3'"]["BND"], horizontalalignment='left',verticalalignment='center', size=9)
    ax.text(0.8, 0.2, "Consensus sequence:", horizontalalignment='left',verticalalignment='top', size=10, fontweight='bold')
    ax.text(0.8, 0.5, SV_ID+"_wtdbg2.ctg.fa", horizontalalignment='left',verticalalignment='center', size=9)

    matplotlib.axes.Axes.invert_yaxis(ax)
    plt.axis('off')

    pdf.savefig()
    plt.close()

#####################################################

print("Start:", datetime.datetime.now())
VCF_IN=args.vcf
VCF_OUTPUT=args.output
INFO_OUTPUT=args.fusion_output
PDF=args.pdf
EnsemblRestClient=EnsemblRestClient()
parse_vcf(VCF_IN, VCF_OUTPUT, INFO_OUTPUT, PDF)
print("End:", datetime.datetime.now())
