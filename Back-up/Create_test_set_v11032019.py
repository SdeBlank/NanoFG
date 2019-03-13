import requests
import json
import argparse
import random
import sys

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Put here a description.')
parser.add_argument('-o', '--output', type=str, help='Output with breakpoint sequences', required=True)

args = parser.parse_args()
VCF=args.output

def select_exons(gene_info):
    for transcript in gene_info["Transcript"]:
        if transcript["is_canonical"]==1:
            EXONS=[]
            INTRONS=[]
            CDS=False

            if transcript["strand"]==1:
                CDS_start=transcript["Translation"]["start"]
                CDS_end=transcript["Translation"]["end"]
            else:
                CDS_start=transcript["Translation"]["end"]
                CDS_end=transcript["Translation"]["start"]

            ### EXONS
            for rank, exon in enumerate(transcript["Exon"]):
                EXON_INFO={}
                EXON_INFO["Type"]="exon"
                EXON_INFO["Rank"]=rank+1
                CHRON_START=exon["start"]
                CHRON_END=exon["end"]
                if gene_info["strand"]==1:
                    EXON_INFO["Start"]=exon["start"]
                    EXON_INFO["End"]=exon["end"]
                else:
                    EXON_INFO["Start"]=exon["end"]
                    EXON_INFO["End"]=exon["start"]

                EXON_INFO["Contains_start_CDS"]=False
                EXON_INFO["Contains_end_CDS"]=False

                if CDS:
                    if CDS_end>=CHRON_START and CDS_end<=CHRON_END:
                        EXON_INFO["Contains_end_CDS"]=True
                        EXON_INFO["CDS"]=True
                        EXON_INFO["Start_phase"]=PHASE
                        EXON_INFO["End_phase"]=0
                        EXON_INFO["CDS_length"]=abs(CDS_end-EXON_INFO["Start"])+1
                        CDS=False
                    else:
                        EXON_INFO["CDS"]=True
                        EXON_INFO["Start_phase"]=PHASE
                        EXON_INFO["End_phase"]=(abs(EXON_INFO["End"]-EXON_INFO["Start"])+1+PHASE)%3
                        EXON_INFO["CDS_length"]=abs(EXON_INFO["End"]-EXON_INFO["Start"])+1
                elif CDS_start>=CHRON_START and CDS_start<=CHRON_END:
                    EXON_INFO["Contains_start_CDS"]=True
                    EXON_INFO["CDS"]=True
                    EXON_INFO["Start_phase"]=0
                    if CDS_end>=CHRON_START and CDS_end<=CHRON_END:
                        EXON_INFO["Contains_end_CDS"]=True
                        EXON_INFO["End_phase"]=0
                        EXON_INFO["CDS_length"]=abs(CDS_end-CDS_start)+1
                    else:
                        EXON_INFO["End_phase"]=(abs(EXON_INFO["End"]-CDS_start)+1)%3
                        EXON_INFO["CDS_length"]=abs(EXON_INFO["End"]-CDS_start)+1
                    CDS=True

                else:
                    EXON_INFO["CDS"]=False
                    EXON_INFO["Start_phase"]=-1
                    EXON_INFO["End_phase"]=-1
                    EXON_INFO["CDS_length"]=0
                PHASE=EXON_INFO["End_phase"]
                EXONS.append(EXON_INFO)

                ### INTRONS
                if rank<len(transcript["Exon"])-1:
                    INTRON_INFO={}
                    INTRON_INFO["Type"]="intron"
                    INTRON_INFO["Rank"]=rank+1
                    INTRON_INFO["Phase"]=EXON_INFO["End_phase"]
                    INTRON_INFO["Start"]=transcript["Exon"][rank]["end"]
                    INTRON_INFO["End"]=transcript["Exon"][rank+1]["start"]
                    INTRON_INFO["CDS"]=CDS
                    INTRONS.append(INTRON_INFO)
    return EXONS, INTRONS

def select_random_fusion(TRANS_CHANCE, SET, FUSION_TYPE, chromosomes, chromosome_lengths):
    while True:
        TRY_AGAIN=False
        #print("request 1")
        while True:
            try:
                CHROM1=random.choice(chromosomes)
                POS=random.randrange(0, chromosome_lengths[CHROM1], 1)
                POS2=POS+4999999
                #print(CHROM, POS)

                SERVER="https://GRCh37.rest.ensembl.org/overlap/region/"
                SPECIES="human"
                HEADERS={"Content-Type" : "application/json"}
                FEATURE="gene"


                TRY=1
                while TRY <= 10:
                    try:
                        genes_request=requests.get(SERVER+SPECIES+"/"+str(CHROM1)+":"+str(POS)+"-"+str(POS2)+"?feature="+FEATURE, headers=HEADERS)
                        genes_response = genes_request.text
                        genes_data=json.loads(genes_response)
                        break
                    except:
                        if TRY==10:
                            sys.exit("Error while requesting from ENSEMBL database after "+str(TRY)+" tries")
                        TRY +=1

                protein_genes=[]

                for gene in genes_data:
                    if gene["biotype"]=="protein_coding":
                        protein_genes.append(gene)

                choice=random.randrange(0,len(protein_genes), 1)
                gene1=protein_genes[choice]

                SERVER="https://GRCh37.rest.ensembl.org/lookup/id/"

                TRY=1
                while TRY <= 10:
                    try:
                        request = requests.get(SERVER+str(gene1["gene_id"])+"?expand=1"+FEATURE, headers=HEADERS)
                        response = request.text
                        gene1_info=json.loads(response)
                        break
                    except:
                        if TRY==10:
                            sys.exit("Error while requesting from ENSEMBL database after "+str(TRY)+" tries")
                        TRY +=1

                RETRY=False
                for transcript in gene1_info["Transcript"]:
                    if transcript["is_canonical"]==1:
                        #print(1, transcript["biotype"])
                        if transcript["biotype"]!="protein_coding":
                            RETRY=True

                        gene1_info["Translation_start"]=transcript["Translation"]["start"]
                        gene1_info["Translation_end"]=transcript["Translation"]["end"]

                if RETRY:
                    continue

                break
            except:
                continue
        #print("request 2")
        while True:
            try:
                if random.random() < TRANS_CHANCE:
                    CHROM2=random.choice(chromosomes)
                    POS=random.randrange(0, chromosome_lengths[CHROM2], 1)
                else:
                    CHROM2=CHROM1
                    POS=random.randrange(0, chromosome_lengths[CHROM2], 1)
                POS2=POS+4999999

                SERVER="https://GRCh37.rest.ensembl.org/overlap/region/"
                #print("request 3")
                TRY=1
                while TRY <= 10:
                    try:
                        genes_request=requests.get(SERVER+SPECIES+"/"+str(CHROM2)+":"+str(POS)+"-"+str(POS2)+"?feature="+FEATURE, headers=HEADERS)
                        genes_response = genes_request.text
                        genes_data=json.loads(genes_response)
                        break
                    except:
                        if TRY==10:
                            sys.exit("Error while requesting from ENSEMBL database after "+str(TRY)+" tries")
                        TRY +=1

                protein_genes=[]
                for gene in genes_data:
                    if gene["biotype"]=="protein_coding":
                        protein_genes.append(gene)

                choice=random.randrange(0,len(protein_genes), 1)
                gene2=protein_genes[choice]

                SERVER="https://GRCh37.rest.ensembl.org/lookup/id/"
                #print("request 4")
                TRY=1
                while TRY <= 10:
                    try:
                        request = requests.get(SERVER+str(gene2["gene_id"])+"?expand=1"+FEATURE, headers=HEADERS)
                        response = request.text
                        gene2_info=json.loads(response)
                        break
                    except:
                        if TRY==10:
                            sys.exit("Error while requesting from ENSEMBL database after "+str(TRY)+" tries")
                        TRY +=1

                RETRY=False
                for transcript in gene2_info["Transcript"]:
                    if transcript["is_canonical"]==1:
                        #print(2, transcript["biotype"])
                        if transcript["biotype"]!="protein_coding":
                            RETRY=True
                        gene2_info["Translation_start"]=transcript["Translation"]["start"]
                        gene2_info["Translation_end"]=transcript["Translation"]["end"]

                if RETRY or gene1_info["display_name"]==gene2_info["display_name"] or gene1_info["id"]==gene2_info["id"]:
                    continue

                break
            except:
                continue

        gene1_exon_info, gene1_intron_info=select_exons(gene1_info)
        gene2_exon_info, gene2_intron_info=select_exons(gene2_info)

        if not gene1_exon_info or not gene1_intron_info or not gene2_exon_info or not gene2_intron_info:
            continue

        if FUSION_TYPE=="intron-intron":
            TRY=0
            while True:
                TRY+=1
                if TRY==10:
                    TRY_AGAIN=True
                    break

                intron1=random.choice(gene1_intron_info)
                intron2=random.choice(gene2_intron_info)

                if not intron1["CDS"] or not intron2["CDS"]:
                    continue

                if gene1["strand"]==1:
                    start1=intron1["Start"]
                    end1=intron1["End"]
                else:
                    start1=intron1["End"]
                    end1=intron1["Start"]

                if gene2["strand"]==1:
                    start2=intron2["Start"]
                    end2=intron2["End"]
                else:
                    start2=intron2["End"]
                    end2=intron2["Start"]

                pos1=random.randrange(start1, end1, 1)
                pos2=random.randrange(start2, end2, 1)

                if gene1["strand"]==gene2["strand"]:
                    read1=[CHROM1, gene1["strand"], pos1-50, pos1, gene1["gene_id"], intron1["Phase"], intron1]
                    read2=[CHROM2, gene2["strand"], pos2, pos2+50, gene2["gene_id"], intron2["Phase"], intron2]
                else:
                    read1=[CHROM1, gene1["strand"], pos1-50, pos1, gene1["gene_id"], intron1["Phase"], intron1]
                    read2=[CHROM2, gene2["strand"], pos2-50, pos2, gene2["gene_id"], intron2["Phase"], intron2]
                # print("intron1 phase:", intron1["Phase"])
                # print("intron2 phase:", intron2["Phase"])

                if intron1["Phase"]==intron2["Phase"] and SET:
                    break
                elif intron1["Phase"]==intron2["Phase"] and not SET:
                    continue
                elif intron1["Phase"]!=intron2["Phase"] and SET:
                    continue
                elif intron1["Phase"]!=intron2["Phase"] and not SET:
                    break
            if TRY_AGAIN==True:
                continue




        if FUSION_TYPE=="exon-exon":
            TRY=0
            while True:
                TRY+=1
                #print("TRY:", TRY)
                if TRY==10:
                    TRY_AGAIN=True
                    break

                exon1=random.choice(gene1_exon_info)
                exon2=random.choice(gene2_exon_info)

                if not exon1["CDS"] or not exon2["CDS"]:
                    continue

                if gene1["strand"]==1:
                    start1=exon1["Start"]
                    end1=exon1["End"]
                else:
                    start1=exon1["End"]
                    end1=exon1["Start"]

                if gene2["strand"]==1:
                    start2=exon2["Start"]
                    end2=exon2["End"]
                else:
                    start2=exon2["End"]
                    end2=exon2["Start"]

                pos1=random.randrange(start1, end1, 1)
                #print(start1, end1)
                #print(gene1_info["Translation_start"], gene1_info["Translation_end"])
                if gene1["strand"]==1:
                    while pos1<gene1_info["Translation_start"] or pos1>gene1_info["Translation_end"]-3:
                        pos1=random.randrange(start1, end1, 1)
                else:
                    while pos1<gene1_info["Translation_start"]+3 or pos1>gene1_info["Translation_end"]:
                        pos1=random.randrange(start1, end1, 1)

                pos2=random.randrange(start2, end2, 1)

                if gene2["strand"]==1:
                    while pos2<gene2_info["Translation_start"] or pos2>gene2_info["Translation_end"]-3:
                        pos2=random.randrange(start2, end2, 1)
                else:
                    while pos2<gene2_info["Translation_start"]+3 or pos2>gene2_info["Translation_end"]:
                        pos2=random.randrange(start2, end2, 1)

                #if exon2["Start_phase"]!=-1 and exon1["End_phase"]!=-1:
                if gene1["strand"]==gene2["strand"]:
                    if gene1["strand"]==1:
                        if exon1["Contains_start_CDS"]:
                            phase1=(abs(gene1_info["Translation_start"]-pos1)+1+exon1["Start_phase"])%3
                        else:
                            phase1=(abs(exon1["Start"]-pos1)+1+exon1["Start_phase"])%3
                    else:
                        if exon1["Contains_start_CDS"]:
                            phase1=(abs(gene1_info["Translation_start"]-pos1)+exon1["Start_phase"])%3
                        else:
                            phase1=(abs(exon1["Start"]-pos1)+exon1["Start_phase"])%3

                    if gene2["strand"]==1:
                        if exon1["Contains_start_CDS"]:
                            phase2=(abs(gene2_info["Translation_start"]-pos2)+exon2["Start_phase"])%3
                        else:
                            phase2=(abs(exon2["Start"]-pos2)+exon2["Start_phase"])%3
                    else:
                        if exon2["Contains_start_CDS"]:
                            phase2=(abs(gene2_info["Translation_start"]-pos2)+1+exon2["Start_phase"])%3
                        else:
                            phase2=(abs(exon2["Start"]-pos2)+1+exon2["Start_phase"])%3

                else:
                    if gene1["strand"]==1:
                        if exon1["Contains_start_CDS"]:
                            phase1=(abs(gene1_info["Translation_start"]-pos1)+1+exon1["Start_phase"])%3
                        else:
                            phase1=(abs(exon1["Start"]-pos1)+1+exon1["Start_phase"])%3

                        if exon2["Contains_start_CDS"]:
                            phase2=(abs(gene2_info["Translation_start"]-pos2)+exon2["Start_phase"])%3
                        else:
                            phase2=(abs(exon2["Start"]-pos2)+exon2["Start_phase"])%3

                    else:
                        if exon1["Contains_start_CDS"]:
                            phase1=(abs(gene1_info["Translation_start"]-pos1)+exon1["Start_phase"])%3
                        else:
                            phase1=(abs(exon1["Start"]-pos1)+exon1["Start_phase"])%3

                        if exon2["Contains_start_CDS"]:
                            phase2=(abs(gene2_info["Translation_start"]-pos2)+1+exon2["Start_phase"])%3
                        else:
                            phase2=(abs(exon2["Start"]-pos2)+1+exon2["Start_phase"])%3

                if gene1["strand"]==gene2["strand"]:
                    read1=[CHROM1, gene1["strand"], pos1-50, pos1, gene1["gene_id"], phase1, exon1]
                    read2=[CHROM2, gene2["strand"], pos2, pos2+50, gene2["gene_id"], phase2, exon2]
                else:
                    read1=[CHROM1, gene1["strand"], pos1-50, pos1, gene1["gene_id"], phase1, exon1]
                    read2=[CHROM2, gene2["strand"], pos2-50, pos2, gene2["gene_id"], phase2, exon2]
                #print("phase1:",phase1)
                #print("phase2:",phase2)
                if phase1==phase2 and SET:
                    print("POS1:", pos1)
                    print("POS2:", pos2)
                    break
                elif phase1==phase2 and not SET:
                    #print("phase1==phase2 and not SET")
                    continue
                elif phase1!=phase2 and SET:
                    #print("phase1!=phase2 and SET")
                    continue
                elif phase1!=phase2 and not SET:
                    print("POS1:", pos1)
                    print("POS2:", pos2)
                    break


        #print("TRY_AGAIN:",TRY_AGAIN)
            if TRY_AGAIN==True:
                continue

            # print("exon1 phase:", exon1["Start_phase"], exon1["End_phase"])
            # print("exon2 phase:", exon2["Start_phase"], exon2["End_phase"])
            # print("phase:",phase1)
            # print("phase:",phase2)

        if FUSION_TYPE=="intron-exon":
            TRY=0
            while True:
                TRY+=1
                #print("TRY:", TRY)
                if TRY==10:
                    TRY_AGAIN=True
                    break

                exon1=random.choice(gene1_exon_info)
                intron2=random.choice(gene2_intron_info)

                if not exon1["CDS"] or not intron2["CDS"] or exon1["Rank"]==1:
                    continue

                if gene1["strand"]==1:
                    if exon1["Rank"]==1:
                        continue
                    start1=exon1["Start"]
                    end1=exon1["End"]
                else:
                    if exon1["Rank"]==len(gene1_exon_info):
                        continue
                    start1=exon1["End"]
                    end1=exon1["Start"]

                if gene2["strand"]==1:
                    start2=intron2["Start"]
                    end2=intron2["End"]
                else:
                    start2=intron2["End"]
                    end2=intron2["Start"]

                pos1=random.randrange(start1, end1, 1)
                #print(start1, end1)
                #print(gene1_info["Translation_start"], gene1_info["Translation_end"])

                if gene1["strand"]==1:
                    while pos1<gene1_info["Translation_start"] or pos1>gene1_info["Translation_end"]-3:
                        pos1=random.randrange(start1, end1, 1)
                else:
                    while pos1<gene1_info["Translation_start"]+3 or pos1>gene1_info["Translation_end"]:
                        pos1=random.randrange(start1, end1, 1)

                pos2=random.randrange(start2, end2, 1)

                if gene1["strand"]==1:
                    if exon1["Contains_start_CDS"]:
                        continue
                    else:
                        #phase1=gene1_exon_info[exon1["Rank"]-2]["End_phase"]
                        phase1=exon1["Start_phase"]
                else:
                    #phase1=gene1_exon_info[exon1["Rank"]]["Start_phase"]
                    phase1=exon1["End_phase"]

                phase2=intron2["Phase"]

                #print("phase1:",phase1)
                #print("phase2:",phase2)
                if gene1["strand"]==gene2["strand"]:
                    read1=[CHROM1, gene1["strand"], pos1-50, pos1, gene1["gene_id"], phase1, exon1]
                    read2=[CHROM2, gene2["strand"], pos2, pos2+50, gene2["gene_id"], phase2, intron2]
                else:
                    read1=[CHROM1, gene1["strand"], pos1-50, pos1, gene1["gene_id"], phase1, exon1]
                    read2=[CHROM2, gene2["strand"], pos2-50, pos2, gene2["gene_id"], phase2, intron2]

                if phase1==phase2 and SET:
                    break
                elif phase1==phase2 and not SET:
                    #print("phase1==phase2 and not SET")
                    continue
                elif phase1!=phase2 and SET:
                    #print("phase1!=phase2 and SET")
                    continue
                elif phase1!=phase2 and not SET:
                    break
        #print("TRY_AGAIN:",TRY_AGAIN)
            if TRY_AGAIN==True:
                continue

        return read1,read2

def get_sequence(range):
    chrom=range[0]
    strand=range[1]
    pos1=range[2]
    pos2=range[3]

    SERVER = "https://grch37.rest.ensembl.org/sequence/region/human/"+str(chrom)+":"+str(pos1)+"-"+str(pos2)+":1"

    TRY=1
    while TRY <= 10:
        try:
            request = requests.get(SERVER, headers={ "Content-Type" : "application/json"})
            data = json.loads(request.text)
            break
        except:
            if TRY==10:
                sys.exit("Error while requesting from ENSEMBL database after "+str(TRY)+" tries")
            TRY +=1
    seq=data["seq"]
    return seq



CHROMOSOMES=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X"]
CHR_LENGTHS={}
for chr in CHROMOSOMES:
    SERVER = "https://GRCh37.rest.ensembl.org/info/assembly/homo_sapiens/"+str(chr)
    TRY=1
    while TRY <= 10:
        try:
            request = requests.get(SERVER, headers={ "Content-Type" : "application/json"})
            data = json.loads(request.text)
            break
        except:
            if TRY==10:
                sys.exit("Error while requesting from ENSEMBL database after "+str(TRY)+" tries")
        TRY +=1

    CHR_LENGTHS[chr]=int(data['length'])

NUMBER=10
TYPES=["exon-exon"]#, "5'UTR", "3'UTR"]
TRANS_PERC=0.2
TRUTH=True

for type in TYPES:
    for x in range(NUMBER):
        #print(str(x+1)+": Getting 2 random genes")
        GENE1, GENE2 = select_random_fusion(TRANS_PERC, TRUTH, type, CHROMOSOMES, CHR_LENGTHS)
        #print(str(x+1)+": Getting fusion sequence")
        SEQ1 = get_sequence(GENE1)
        SEQ2 = get_sequence(GENE2)
        if GENE1[1]!=GENE2[1]:
            complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
            SEQ2=''.join([complement[base] for base in SEQ2[::-1]])
        FULL_SEQ=SEQ1+SEQ2
        #print("FUSION " + str(x+1))
        print(">"+type+"_"+str(TRUTH)+"_"+str(x+1) + "\t" + type + "\t" + GENE1[4]+"-"+GENE2[4])
        print(GENE1[0:6])
        print(GENE1[6])
        print(GENE2[0:6])
        print(GENE2[6])
        print(SEQ1, SEQ2)
