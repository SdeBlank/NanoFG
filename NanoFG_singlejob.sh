#!/bin/bash

usage() {
echo "

NanoFG_singlejob.sh -b BAM [-v VCF] [-s SELECTED_GENES_OR_REGIONS] [-df]

Required parameters:
    -b|--bam                                                           Path to bam file

Optional parameters:

GENERAL
    -h|--help                                                          Shows help-v|--vcf
    -v|--vcf		                                                       Path to vcf file
    -t|--threads
    -e|--venv                                                          Path to virtual environment[${VENV}]

SELECTION AND FILTERING
    -s|--selection                                                     Select genes or areas to check for fusion genes.
                                                                       Insert a list of genes or areas, separated by a comma
                                                                       e.g. 'BRAF,TP53' or 'ENSG00000157764,ENSG00000141510' or '17:7565097-7590yy856'
    -dc|--dont_clean                                                   Don't clean up the intermediate files
    -df|--dont_filter                                                  Don't filter out all non-PASS SVs
    -cc|--consensus_calling                                            Create a consensus sequence of the fusion-supporting reads. Not recommended on low-coverage data.

OUTPUT
    -o|--outputdir                                                     Path to output directory
    -io|--info_output                                                  Path to the NanoFG output info file []
    -vo|--vcf_output                                                   Path to the NanoFG output vcf file []
    -p|--pdf                                                           Path to the NanoFG output pdf file []

REQUIRED TOOLS
    -sv|--sv_caller                                                    NanoSV or path to Sniffles [${SV_CALLER}]
    -sa|--samtools                                                     Path to sambamba|samtools [${SAMTOOLS}]
    -mm2|--minimap2                                                    Path to minimap2 [${MINIMAP2}]
    -l|--last_dir                                                      Path to LAST directory [${LAST_DIR}]
    -w|--wtdbg2_dir                                                    Path to wtdbg2 directory [${WTDBG2_DIR}]

SCRIPTS
    -fres|--fusion_read_extraction_script                              Path to the fusion_read_extraction.py script [${FUSION_READ_EXTRACTION_SCRIPT}]
    -fcs|--fusion_check_script                                         Path to vcf_primer_filter.py [$FUSION_CHECK_SCRIPT]

CONSENSUS CALLING
    -ccws|--consensus_calling_wtdbg2_settings                          Wtdbg2 settings [${CONSENSUS_CALLING_WTDBG2_SETTINGS}]

SV CALLING
    -nmc|--nanosv_minimap2_config                                               Path to config to use for nanosv [$NANOSV_MINIMAP2_CONFIG]
    -ss|--sniffles_settings                                            Settings to use for sniffles [$SNIFFLES_SETTINGS]

LAST MAPPING
    -lmr|--last_mapping_refgenome                                      Reference genome [${LAST_MAPPING_REFGENOME}]
    -lmrd|--last_mapping_refdict                                       Reference genome .dict file [${LAST_MAPPING_REFDICT}]
    -lms|--last_mapping_settings                                       LAST settings [${LAST_MAPPING_SETTINGS}]
    -lmt|--last_mapping_threads                                        Number of threads [${LAST_MAPPING_THREADS}]

PRIMER DESIGN
    -pdf|--primer_design_flank                                         Flanking distance around to breakpoint to extract [$PRIMER_DESIGN_FLANK]
    -pdhv|--primer_design_h_vmem                                       Primer design memory [$PRIMER_DESIGN_MEM]
    -pdhr|--primer_design_h_rt                                         Primer design time [$PRIMER_DESIGN_TIME]
    -pdd|--primer_design_dir                                           Path to primer3 directory [$PRIMER_DESIGN_DIR]
    -pdb|--primer_design_bindir                                        Path to primer3 bin dir [$PRIMER_DESIGN_BINDIR]
    -pdpt|--primer_design_pcr_type                                     PCR type [$PRIMER_DESIGN_PCR_TYPE]
    -pdtp|--primer_design_tilling_params                               Tilling params [$PRIMER_DESIGN_TILLING_PARAMS]
    -pdp|--primer_design_psr                                           PSR [$PRIMER_DESIGN_PSR]
    -pdgp|--primer_design_guix_profile                                 Path to guix profile [$PRIMER_DESIGN_GUIX_PROFILE]
    -pdpc|--primer_design_primer3_core                                 Path to primer3_core [$PRIMER_DESIGN_PRIMER3_CORE]
    -pdm|--primer_design_mispriming                                    Path to mispriming [$PRIMER_DESIGN_MISPRIMING]
"
}

POSITIONAL=()

#GENERAL DEFAULTS
NANOFG_DIR=$(realpath $(dirname ${BASH_SOURCE[0]}))

source $NANOFG_DIR/paths.ini

PIPELINE_DIR=$NANOFG_DIR/pipeline
FILES_DIR=$NANOFG_DIR/files
SCRIPT_DIR=$NANOFG_DIR/scripts
VENV=${NANOFG_DIR}/venv/bin/activate

THREADS=8
CONSENSUS_CALLING=false
DONT_CLEAN=false
DONT_FILTER=false

OUTPUTDIR=$(realpath ./)

#TOOL PATH DEFAULTS
SAMTOOLS=$PATH_SAMTOOLS
MINIMAP2=$PATH_MINIMAP2
LAST_DIR=$PATH_LAST_DIR
WTDBG2_DIR=$PATH_WTDBG2_DIR

SV_CALLER='/hpc/cog_bioinf/kloosterman/tools/NanoSV/nanosv/NanoSV.py'
NANOSV_MINIMAP2_CONFIG=$FILES_DIR/nanosv_minimap2_config.ini
NANOSV_LAST_CONFIG=$FILES_DIR/nanosv_last_config.ini
NANOSV_LAST_CONSENSUS_CONFIG=$FILES_DIR/nanosv_last_consensus_config.ini
SNIFFLES_SETTINGS='-s 2 -n -1 --genotype'
#SNIFFLES_SETTINGS='-s 2 -n -1 --genotype -d 1'     SETTINGS TO DETECT RECIPROCAL TRANSLOCATIONS

#REGION SELECTION
REGION_SELECTION_SCRIPT=$SCRIPT_DIR/RegionSelection.py
REGION_SELECTION_BED_OUTPUT=$OUTPUTDIR/regions.bed
REGION_SELECTION_BAM_OUTPUT=$OUTPUTDIR/regions.bam

#FUSION READ EXTRACTION DEFAULTS
FUSION_READ_EXTRACTION_SCRIPT=$SCRIPT_DIR/FusionReadExtraction.py

#CONSENSUS CALLING DEFAULTS
CONSENSUS_CALLING_WTDBG2_SETTINGS='-x ont -g 3g -q'

#LAST MAPPING DEFAULTS
LAST_MAPPING_REFFASTA=$PATH_HOMO_SAPIENS_REFFASTA
LAST_MAPPING_REFGENOME=$PATH_HOMO_SAPIENS_REFGENOME
LAST_MAPPING_REFDICT=$PATH_HOMO_SAPIENS_REFDICT
LAST_MAPPING_SETTINGS="-Q 0 -p ${LAST_DIR}/last_params"
LAST_MAPPING_THREADS=1

#FUSION CHECK DEFAULTS
FUSION_CHECK_SCRIPT=$SCRIPT_DIR/FusionCheck.py

#PRIMER DESIGN DEFAULTS
PRIMER_DESIGN_DIR=$PATH_PRIMER_DESIGN_DIR
PRIMER_DESIGN_BINDIR=$PRIMER_DESIGN_DIR/primers
PRIMER_DESIGN_GUIX_PROFILE=$PRIMER_DESIGN_DIR/emboss/.guix-profile
PRIMER_DESIGN_PRIMER3_CORE=$PRIMER_DESIGN_DIR/primer3/src/primer3_core
PRIMER_DESIGN_MISPRIMING=$PRIMER_DESIGN_DIR/repbase/current/empty.ref
PRIMER_DESIGN_PCR_TYPE='single'
PRIMER_DESIGN_TILLING_PARAMS=''
PRIMER_DESIGN_PSR='60-200'
PRIMER_DESIGN_FLANK=200


while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
    -h|--help)
    usage
    exit
    shift # past argument
    ;;
    -v|--vcf)
    VCF="$2"
    shift # past argument
    shift # past value
    ;;
    -b|--bam)
    BAM="$2"
    shift # past argument
    shift # past value
    ;;
    -s|--selection)
    SELECTION="$2"
    shift # past argument
    shift # past value
    ;;
    -io|--info_output)
    FUSION_CHECK_INFO_OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
    -vo|--vcf_output)
    FUSION_CHECK_VCF_OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
    -p|--pdf)
    FUSION_CHECK_PDF_OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--outputdir)
    OUTPUTDIR="$2"
    shift # past argument
    shift # past value
    ;;
    -d|--nanofg_dir)
    NANOFG_DIR="$2"
    shift # past argument
    shift # past value
    ;;
    -cc|--consensus_calling)
    CONSENSUS_CALLING=true
    shift # past argument
    ;;
    -dc|--dont_clean)
    DONT_CLEAN=true
    shift # past argument
    ;;
    -df|--dont_filter)
    DONT_FILTER=true
    shift # past argument
    ;;
    -e|--venv)
    VENV="$2"
    shift # past argument
    shift # past value
    ;;
    -sa|--samtools)
    SAMTOOLS="$2"
    shift # past argument
    shift # past value
    ;;
    -mm2|--minimap2)
    MINIMAP2="$2"
    shift # past argument
    shift # past value
    ;;
    -l|--last_dir)
    LAST_DIR="$2"
    shift # past argument
    shift # past value
    ;;
    -w|--wtdbg2_dir)
    WTDBG2_DIR="$2"
    shift # past argument
    shift # past value
    ;;
    -sv|--sv_caller)
    SV_CALLER="$2"
    shift # past argument
    shift # past value
    ;;
    -nmc|--nanosv_minimap2_config)                                           #### NOT USED
    NANOSV_MINIMAP2_CONFIG="$2"
    shift # past argument
    shift # past value
    ;;
    -ss|--sniffles_settings)
    SNIFFLES_SETTINGS="$2"
    shift # past argument
    shift # past value
    ;;
    -fres|--fusion_read_extraction_script)
    FUSION_READ_EXTRACTION_SCRIPT="$2"
    shift # past argument
    shift # past value
    ;;
    -ccws|--consensus_calling_wtdbg2_settings)
    CONSENSUS_CALLING_WTDBG2_SETTINGS="$2"
    shift # past argument
    shift # past value
    ;;
    -lmr|--last_mapping_refgenome)
    LAST_MAPPING_REFGENOME="$2"
    shift # past argument
    shift # past value
    ;;
    -lmrd|--last_mapping_refdict)
    LAST_MAPPING_REFDICT="$2"
    shift # past argument
    shift # past value
    ;;
    -lms|--last_mapping_settings)
    LAST_MAPPING_SETTINGS="$2"
    shift # past argument
    shift # past value
    ;;
    -fcs|--fusion_check_script)
    FUSION_CHECK_SCRIPT="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

if [ -z $BAM ]; then
    echo "Missing -b|--bam parameter"
    usage
    exit
fi

echo -e "`date` \t Running on `uname -n`"

SAMPLE=$(basename $BAM)
SAMPLE=${SAMPLE/.bam/}

if [ -z $FUSION_CHECK_VCF_OUTPUT ]; then
    FUSION_CHECK_VCF_OUTPUT=$OUTPUTDIR/${SAMPLE}_FusionGenes.vcf
fi

if [ -z $FUSION_CHECK_INFO_OUTPUT ]; then
    FUSION_CHECK_INFO_OUTPUT=$OUTPUTDIR/${SAMPLE}_FusionGenesInfo.txt
fi

if [ -z $FUSION_CHECK_PDF_OUTPUT ]; then
    FUSION_CHECK_PDF_OUTPUT=$OUTPUTDIR/${SAMPLE}_FusionGenes.pdf
fi

VCF_OUTPUT=$OUTPUTDIR/${SAMPLE}_FusionGenes.vcf
INFO_OUTPUT=$OUTPUTDIR/${SAMPLE}_FusionGenesInfo.txt

PIPELINE_DIR=$NANOFG_DIR/pipeline
SCRIPT_DIR=$NANOFG_DIR/scripts
CANDIDATE_DIR=$OUTPUTDIR/candidate_fusions
PRIMER_DIR=$OUTPUTDIR/primers

BAM_MERGE_OUT=$OUTPUTDIR/candidate_fusion_genes.bam
SV_CALLING_OUT=$OUTPUTDIR/candidate_fusion_genes.vcf

if [ -d $CANDIDATE_DIR ]; then
  rm $CANDIDATE_DIR/*
else
  mkdir -p $CANDIDATE_DIR
fi

if [ ! -d $CANDIDATE_DIR ]; then
    exit
fi

. $VENV

################################################## SELECTION OF REGIONS THAT ARE SPECIFIED BY THE USER (OPTIONAL)
if [ ! -z $VCF ];then
  REGION_SELECTION_BAM_OUTPUT=$BAM
  echo "### vcf (-v) already provided. Skipping selection and sv calling"
elif [ -z $SELECTION ];then
  REGION_SELECTION_BAM_OUTPUT=$BAM
  echo "### No selection parameter (-s) provided. Using all mapped reads"
else
  echo -e "`date` \t 'Selecting regions to check for fusion genes..."
  python $REGION_SELECTION_SCRIPT \
  -b $REGION_SELECTION_BED_OUTPUT \
  -r $SELECTION

  REGION_SELECTION_SAM_OUTPUT=${REGION_SELECTION_BAM_OUTPUT/.bam/.sam}
  $SAMTOOLS view -H $BAM > $REGION_SELECTION_SAM_OUTPUT
  $SAMTOOLS view -@ $THREADS -L $REGION_SELECTION_BED_OUTPUT $BAM | cut -f 1 | sort -k1n | uniq > reads.tmp
  $SAMTOOLS view  -@ $THREADS $BAM | grep -f reads.tmp >> $REGION_SELECTION_SAM_OUTPUT
  $SAMTOOLS view  -@ $THREADS -h -S -f bam $REGION_SELECTION_SAM_OUTPUT | $SAMTOOLS sort -o $REGION_SELECTION_BAM_OUTPUT /dev/stdin
  $SAMTOOLS index $REGION_SELECTION_BAM_OUTPUT
  rm reads.tmp
  rm $REGION_SELECTION_SAM_OUTPUT
fi

################################################## CALLING OF SVS ON GIVEN BAM FILE (OPTIONAL)
if [ -z $VCF ]; then
  echo -e "`date` \t SV calling..."
  VCF=$OUTPUTDIR/${SAMPLE}.vcf
  bash $PIPELINE_DIR/sv_calling.sh \
    -sv $SV_CALLER \
    -b $REGION_SELECTION_BAM_OUTPUT \
    -t $THREADS \
    -s $SAMTOOLS \
    -v $VENV \
    -c $NANOSV_MINIMAP2_CONFIG \
    -ss "$SNIFFLES_SETTINGS" \
    -o $VCF
fi

################################################## REMOVAL OF INSERTIONS (MANDATORY) AND SVS WITHOUT THE PASS FILTER (OPTIONAL)

VCF_FILTERED=${OUTPUTDIR}/$(basename $VCF)

if [ $DONT_FILTER = false ];then
  echo -e "`date` \t Removing insertions and all variants without a PASS filter..."
  VCF_FILTERED=${VCF_FILTERED/.vcf/_noINS_PASS.vcf}
  grep "^#" $VCF > $VCF_FILTERED
  grep -v "^#" $VCF | awk '$5!="<INS>"' | awk '$7=="PASS"' >> $VCF_FILTERED
else
  echo -e "`date` \t Removing insertions..."
  VCF_FILTERED=${VCF_FILTERED/.vcf/_noINS.vcf}
  grep "^#" $VCF > $VCF_FILTERED
  grep -v "^#" $VCF | awk '$5!="<INS>"' >> $VCF_FILTERED
fi

##################################################  EXTRACTION OF READS THAT SUPPORT POSSIBLE FUSION-GENERATING SVS, BASED ON ORIENTATION AND STRAND COMBINATION
echo -e "`date` \t Extracting reads that support candidate fusion genes..."

python $FUSION_READ_EXTRACTION_SCRIPT \
  -b $BAM \
  -v $VCF_FILTERED \
  -o $CANDIDATE_DIR

################################################## CREATION OF A CONSENSUS SEQUENCE OF THE READS THAT SUPPORT A SV (OPTIONAL, NOT RECOMMENDED IF LOW COVERAGE DATA)
if [ $CONSENSUS_CALLING = true ];then
  echo -e "`date` \t Producing consensus if possible..."
else
  echo -e "`date` \t Consensus calling not activated. Use '-cc' to turn consensus calling on"
fi

for FASTA in $CANDIDATE_DIR/*.fasta; do
  if [ $CONSENSUS_CALLING = true ];then
    bash $PIPELINE_DIR/consensus_calling.sh \
      -f $FASTA \
      -t $THREADS \
      -w $WTDBG2_DIR \
      -ws  "$CONSENSUS_CALLING_WTDBG2_SETTINGS"
  fi

  if ! [ -s ${FASTA/.fasta/_wtdbg2.ctg.fa} ] && [ $CONSENSUS_CALLING = false ];then
    ln -s $FASTA ${FASTA/.fasta/.no_ctg.fa}
  fi
done

################################################## MAPPING OF THE READS (OR CONSENSUS) USING LAST FOR NANOSV AND MINIMAP2 FOR SNIFFLES.
################################################## LAST IS PREFERRED FOR ACCURATE EXON-EXON FUSION DETECTION
echo -e "`date` \t Mapping candidate fusion genes..."

if [[ $SV_CALLER == *"nanosv"* ]] || [[ $SV_CALLER == *"NanoSV"* ]]; then
  MAPPING_ARGS="-t $LAST_MAPPING_THREADS -r $LAST_MAPPING_REFGENOME -rd $LAST_MAPPING_REFDICT -l $LAST_DIR -ls '$LAST_MAPPING_SETTINGS' -s $SAMTOOLS"
  for FA in $CANDIDATE_DIR/*.fa; do
    echo $FA;
  done | \
  xargs -I{} --max-procs $THREADS bash -c "bash $PIPELINE_DIR/last_mapping.sh -f {} $MAPPING_ARGS; exit 1;"

elif [[ $SV_CALLER == *"sniffles"* ]] || [[ $SV_CALLER == *"Sniffles"* ]]; then
  MAPPING_ARGS="-mm2 $MINIMAP2 -r $LAST_MAPPING_REFFASTA -t $THREADS -s $SAMTOOLS"
  for FA in $CANDIDATE_DIR/*.fa; do
    echo $FA;
  done | \
  xargs -I{} --max-procs $THREADS bash -c "bash $PIPELINE_DIR/minimap2_mapping.sh -f {} $MAPPING_ARGS; exit 1;"
fi

################################################## MERGING ALL SEPARATE BAM FILES OF ALL FUSION GENE CANDIDATES INTO A SINGLE BAM FILE
echo -e "`date` \t Merging bams..."

NUMBER_OF_BAMS=$(ls $CANDIDATE_DIR/*.last.sorted.bam | wc -l)

if [[ NUMBER_OF_BAMS -gt 1 ]];then
  $SAMTOOLS merge -f $BAM_MERGE_OUT $CANDIDATE_DIR/*.last.sorted.bam
  if [[ $SAMTOOLS == *"samtools"* ]] || [[ $SAMTOOLS == *"Samtools"* ]]; then
    $SAMTOOLS index $BAM_MERGE_OUT
  fi
elif [[ NUMBER_OF_BAMS -eq 1 ]];then
  cp $CANDIDATE_DIR/*.last.sorted.bam $BAM_MERGE_OUT
  cp $CANDIDATE_DIR/*.last.sorted.bam.bai ${BAM_MERGE_OUT}.bai
  echo "NO CANDIDATE FUSION GENES FOUND"
else
  exit
fi

if [[ $SV_CALLER == *"sniffles"* ]] || [[ $SV_CALLER == *"Sniffles"* ]]; then
  $SAMTOOLS calmd $BAM_MERGE_OUT $LAST_MAPPING_REFFASTA -b > temp.bam
  mv temp.bam $BAM_MERGE_OUT
  $SAMTOOLS index $BAM_MERGE_OUT
fi
################################################## CALLING SVS FOR THE MAPPED FUSION CANDIDATES
echo -e "`date` \t Calling SVs..."

if [ $CONSENSUS_CALLING = true ];then
  bash $PIPELINE_DIR/sv_calling.sh \
    -sv $SV_CALLER \
    -b $BAM_MERGE_OUT \
    -t $THREADS \
    -s $SAMTOOLS \
    -v $VENV \
    -c $NANOSV_LAST_CONSENSUS_CONFIG \
    -ss "$SNIFFLES_SETTINGS" \
    -o $SV_CALLING_OUT
else
  bash $PIPELINE_DIR/sv_calling.sh \
    -sv $SV_CALLER \
    -b $BAM_MERGE_OUT \
    -t $THREADS \
    -s $SAMTOOLS \
    -v $VENV \
    -c $NANOSV_LAST_CONFIG \
    -ss "$SNIFFLES_SETTINGS" \
    -o $SV_CALLING_OUT
fi

# if [ $CONSENSUS_CALLING = true ];then
#   SV_CALLING_SETTINGS="-sv $SV_CALLER -t 1 -s $SAMTOOLS -v $VENV -c $NANOSV_LAST_CONSENSUS_CONFIG -ss '$SNIFFLES_SETTINGS'"
# else
#   SV_CALLING_SETTINGS="-sv $SV_CALLER -t 1 -s $SAMTOOLS -v $VENV -c $NANOSV_LAST_CONFIG -ss '$SNIFFLES_SETTINGS'"
# fi
#
# for CANDIDATE_BAM in $CANDIDATE_DIR/*.last.sorted.bam; do
#   echo ${CANDIDATE_BAM/.bam/};
# done | \
# xargs -I{} --max-procs $THREADS bash -c "bash $PIPELINE_DIR/sv_calling.sh -b {}.bam -o {}.vcf $SV_CALLING_SETTINGS; exit 1;"
#
# grep "^#" $(ls $CANDIDATE_DIR/*.last.sorted.vcf | head -n 1) > $SV_CALLING_OUT
# for CANDIDATE_VCF in $CANDIDATE_DIR/*.last.sorted.vcf; do
#   grep -v "^#" $CANDIDATE_VCF >> $SV_CALLING_OUT
# done

################################################## REMOVAL OF INSERTIONS (MANDATORY) AND SVS WITHOUT THE PASS FILTER (OPTIONAL)

if [ $DONT_FILTER = false ];then
  SV_CALLING_OUT_FILTERED=${SV_CALLING_OUT/.vcf/_noINS_PASS.vcf}
  echo -e "`date` \t Removing insertions and all variants without a PASS filter..."
  grep "^#" $SV_CALLING_OUT > $SV_CALLING_OUT_FILTERED
  grep -v "^#" $SV_CALLING_OUT | awk '$5!="<INS>"' | awk '$7=="PASS"' >> $SV_CALLING_OUT_FILTERED
else
  SV_CALLING_OUT_FILTERED=${SV_CALLING_OUT/.vcf/_noINS.vcf}
  echo -e "`date` \t Removing insertions..."
  grep "^#" $SV_CALLING_OUT > $SV_CALLING_OUT_FILTERED
  grep -v "^#" $SV_CALLING_OUT | awk '$5!="<INS>"' >> $SV_CALLING_OUT_FILTERED
fi

################################################### CHECKING THE FUSION GENE FOR ADDITIONAL INFORMATION (FRAME, SIMILARITY, GENE OVERLAP, ETC.)
echo -e "`date` \t Checking fusion candidates..."

bash $PIPELINE_DIR/fusion_check.sh \
  -v $SV_CALLING_OUT_FILTERED \
  -ov $VCF \
  -o $FUSION_CHECK_VCF_OUTPUT \
  -fo $FUSION_CHECK_INFO_OUTPUT \
  -p $FUSION_CHECK_PDF_OUTPUT \
  -s $FUSION_CHECK_SCRIPT \
  -e $VENV

################################################### DESIGNING PRIMERS FOR THE DETECTED FUSIONS
echo -e "`date` \t Designing primers around fusion breakpoints..."
mkdir -p $PRIMER_DIR
if [ ! -d $PRIMER_DIR ]; then
  exit
fi

python $SCRIPT_DIR/Primerseq.py \
 -v $FUSION_CHECK_VCF_OUTPUT \
 -d $PRIMER_DIR \
 -f $PRIMER_DESIGN_FLANK

mkdir -p $PRIMER_DIR/tmp
if [ ! -d $PRIMER_DIR/tmp ]; then
 exit
fi

if [ -d $PRIMER_DESIGN_DIR ];then
  cd $PRIMER_DIR/tmp
  for FUSION_FASTA in $PRIMER_DIR/*.fasta; do
    if [ -z $PRIMER_DESIGN_TILLING_PARAMS ]; then
      bash $PIPELINE_DIR/primer_design.sh \
       -f $FUSION_FASTA \
       -o ${FUSION_FASTA/.fasta/.primers} \
       -pdb $PRIMER_DESIGN_BINDIR \
       -pdpt $PRIMER_DESIGN_PCR_TYPE \
       -pdp $PRIMER_DESIGN_PSR \
       -pdgp $PRIMER_DESIGN_GUIX_PROFILE \
       -pdpc $PRIMER_DESIGN_PRIMER3_CORE \
       -pdm $PRIMER_DESIGN_MISPRIMING
    else
      bash $PIPELINE_DIR/primer_design.sh \
       -f $FUSION_FASTA \
       -o ${FUSION_FASTA/.fasta/.primers} \
       -pdb $PRIMER_DESIGN_BINDIR \
       -pdpt $PRIMER_DESIGN_PCR_TYPE \
       -pdtp $PRIMER_DESIGN_TILLING_PARAMS \
       -pdp $PRIMER_DESIGN_PSR \
       -pdgp $PRIMER_DESIGN_GUIX_PROFILE \
       -pdpc $PRIMER_DESIGN_PRIMER3_CORE \
       -pdm $PRIMER_DESIGN_MISPRIMING
    fi
    mv $PRIMER_DIR/tmp/primer3.out ${FUSION_FASTA/.fasta/_primerinfo.txt}
    rm $PRIMER_DIR/tmp/*
  done
  rmdir $PRIMER_DIR/tmp/
  cat $PRIMER_DIR/*.primers > $OUTPUTDIR/${SAMPLE}_FusionGenes.primers
  cd $OUTPUTDIR
else
  echo "Path to primer3 does not exist. Giving only breakpoint sequences"
  cat $PRIMER_DIR/*.fasta > $OUTPUTDIR/${SAMPLE}_FusionGenesBNDseq.fasta
fi

################################################### IF -DC IS SPECIFIED, ALL FILES EXCEPT THE FINAL FILES ARE DELETED TO PROVIDE A CLEAN OUTPUT

if [ $DONT_CLEAN = false ];then
  rm $VCF_FILTERED
  rm $CANDIDATE_DIR/*
  rmdir $CANDIDATE_DIR
  rm $BAM_MERGE_OUT
  rm $BAM_MERGE_OUT.bai
  rm $SV_CALLING_OUT
  rm $SV_CALLING_OUT_FILTERED
  rm $PRIMER_DIR/*
  rmdir $PRIMER_DIR
fi

deactivate

echo -e "`date` \t Done"
