import pysam
import vcf as pyvcf
import argparse
import re
import copy

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Input parameters for NanoFG.')
parser.add_argument('-v', '--vcf', type=str, help='Input NanoSV vcf file', required=True)
parser.add_argument('-b', '--bam', type=str, help='Input bam file', required=True)
parser.add_argument('-o', '--output', type=str, help='Output complex SV vcf file', required=True)
args = parser.parse_args()

with open(args.vcf, "r") as vcf, open(args.output, "w") as vcf_output:
    vcf_reader=pyvcf.Reader(vcf)
    vcf_reader.infos['FUSION']=pyvcf.parser._Info('FUSION', ".", "String", "Gene names of the fused genes reported if present", "NanoSV", "X")
    vcf_reader.infos['ORIGINAL_SVID']=pyvcf.parser._Info('ORIGINAL_SVID', "1", "Integer", "SVID in the vcf of the full vcf", "NanoSV", "X")
    vcf_writer=pyvcf.Writer(vcf_output, vcf_reader, lineterminator='\n')

    RECORDS={}

    ### DETERMINE IF VCF IS PRODUCED BY SNIFFLES OR NANOSV
    if "source" in vcf_reader.metadata:
        if vcf_reader.metadata["source"][0].lower()=="sniffles":
            vcf_type="Sniffles"
            vcf_reader.infos['RNAMES']=pyvcf.parser._Info('RNAMES', ".", "String", "Names of reads supporting SVs (comma separated)", "Sniffles", "X")
    elif "cmdline" in vcf_reader.metadata:
        if "nanosv" in vcf_reader.metadata["cmdline"][0].lower():
            vcf_type="NanoSV"
    else:
        sys.exit("Unknown VCF format or re-run sniffle with '-n -1' to get supporting reads")

    svs_per_read={}
    for record in vcf_reader:
        RECORDS[record.ID]=record
        if vcf_type=="NanoSV":
            compared_id=re.findall("^\d+", record.INFO["ALT_READ_IDS"][0])[0]
            pos1_orientation=record.ALT[0].orientation
            pos2_orientation=record.ALT[0].remoteOrientation
            ### Use set to remove reads that support the same breakpoint twice. This is most likely a result of WGA (whole genome amplification).
            ### Natural occurrences can be due to tandem duplications, but currently seem unlikely to be very important in fusions
            good_reads=[read for read in record.INFO["ALT_READ_IDS"] if record.INFO["ALT_READ_IDS"].count(read)==1]
            for supporting_read in good_reads:
                if supporting_read not in svs_per_read:
                    svs_per_read[supporting_read]=[record.ID]
                elif record.ID not in svs_per_read[supporting_read]:
                    svs_per_read[supporting_read].append(record.ID)

    ### Select reads that support multiple breakpoints
    complex_reads={}
    for key, value in svs_per_read.items():
        if not len(value)<2:
            complex_sv="-".join(value)
            if complex_sv not in complex_reads:
                complex_reads[complex_sv]=[key]
            else:
                complex_reads[complex_sv].append(key)

    alignments={}
    not_unique=[]

    ### Get information on all reads from the bam file
    all_reads={}
    bamfile = pysam.AlignmentFile(args.bam, "rb" )
    for read in bamfile.fetch(until_eof=True):
        if not read.seq == None and not read.is_unmapped:
            if read not in not_unique:
                not_unique.append(read)
                if read.is_reverse:
                    left_clipped=int(re.findall("\d+", read.cigarstring)[-1])
                    right_clipped=int(re.findall("\d+", read.cigarstring)[0])
                    start=read.reference_end
                    end=read.reference_start
                    strand="-"
                else:
                    left_clipped=int(re.findall("\d+", read.cigarstring)[0])
                    right_clipped=int(re.findall("\d+", read.cigarstring)[-1])
                    start=read.reference_start
                    end=read.reference_end
                    strand="+"

                alignment_info=[read.query_name, read.reference_name, int(start), int(end) , left_clipped, right_clipped, strand, abs(int(start)-int(end))]
                if read.query_name not in all_reads:
                    all_reads[read.query_name]=[alignment_info]
                else:
                    compare=copy.deepcopy(all_reads[read.query_name])
                    for index, alignment in enumerate(compare):
                        if left_clipped<alignment[4]:
                            all_reads[read.query_name].insert(index, alignment_info)
                        elif index+1==len(all_reads[read.query_name]):
                            all_reads[read.query_name].append(alignment_info)
    bamfile.close()
    ### Go through every complex SV, currently take the first reads that supports the complex SV and check orientation
    for complex_sv, reads in complex_reads.items():
        completed=False
        for selected_read in reads:
            alignments[selected_read]=all_reads[selected_read]

            piece_length=0
            for piece in alignments[selected_read][1:-1]:
                piece_length+=piece[7]

            bnd_info={}
            BND=1
            for loc in range(len(alignments[selected_read])-1):
                bnd_info[BND]=[selected_read , alignments[selected_read][loc][1], int(alignments[selected_read][loc][3]),
                                alignments[selected_read][loc+1][1], int(alignments[selected_read][loc+1][2])]
                BND+=1
            if len(bnd_info.keys())>len(complex_sv.split("-")):
                continue
            else:
                completed=True
                break
        if not completed:
            print("Number of aligments in the read do not overlap with the number of SVs detected by the used SV caller")
            continue
        Likeliness={}
        for key, value in bnd_info.items():
            Likeliness[key]={}
            for SV in complex_sv.split("-"):
                if (str(RECORDS[SV].CHROM)==value[1] and str(RECORDS[SV].ALT[0].chr)==value[3]) or (str(RECORDS[SV].CHROM)==value[3] and str(RECORDS[SV].ALT[0].chr)==value[1]):
                    Likeliness[key][SV]=abs(RECORDS[SV].POS-value[2])+abs(RECORDS[SV].ALT[0].pos-value[4])

        order=[]
        for rank, score in Likeliness.items():
            closest=100000000000
            for sv_id,distance in score.items():
                if distance<closest:
                    most_likely=sv_id
                    closest=distance
            order.append(most_likely)

        ## Creating a hardcopy produces an error in writing the final vcf. Currently no solution has been found to write both the separate SVs and the linked SV
        #first_SV=copy.deepcopy(RECORDS[order[0]])
        #last_SV=copy.deepcopy(RECORDS[order[-1]])
        first_SV=RECORDS[order[0]]
        last_SV=RECORDS[order[-1]]
        complex_record=first_SV
        orientation = first_SV.ALT[0].orientation
        remoteOrientation = last_SV.ALT[0].remoteOrientation
        complex_record.ID = complex_sv
        complex_record.ALT = [ pyvcf.model._Breakend( last_SV.ALT[0].chr, last_SV.ALT[0].pos, orientation, remoteOrientation, record.REF, True ) ]
        complex_record.FILTER.append("complex")

        vcf_writer.write_record(complex_record)
