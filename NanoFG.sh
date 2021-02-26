#!/bin/bash

usage() {
echo "
bash NanoFG.sh -f </path/to/fastq> [-s SELECTED_GENES_OR_REGIONS] [-cc] [-df] [-dc]
OR
bash NanoFG.sh -b </path/to/bam> [-v </path/to/vcf>] [-s SELECTED_GENES_OR_REGIONS] [-cc] [-df] [-dc]


Required parameters:
    -f|--fastqdir                                                      Path to fastq directory
    OR
    -b|--bam                                                           Path to bam file

Optional parameters:

GENERAL
    -h|--help                                                          Shows help
    -n|--name                                                          Name of the sample used
    -v|--vcf		                                                       Path to vcf file
    -t|--threads                                                       Number of threads
    -e|--venv                                                          Path to virtual environment[${VENV}]
    -wl|--without_last                                                 Do not remap fusion candidates with LAST, but use minimap2

SELECTION AND FILTERING
    -nc|--non_coding                                                   Also include non-coding fusions in the results (Not fully tested yet)
    -s|--selection                                                     Select genes or areas to check for fusion genes.
                                                                       Insert a list of genes or areas, separated by a comma
                                                             e.g. 'BRAF,TP53' or 'ENSG00000157764,ENSG00000141510' or '17:7565097-7590856'
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

CONSENSUS CALLING
    -ccws|--consensus_calling_wtdbg2_settings                          Wtdbg2 settings [${CONSENSUS_CALLING_WTDBG2_SETTINGS}]

SV CALLING SETTINGS
    -nmc|--nanosv_minimap2_config                                      Path to config to use for nanosv [$NANOSV_MINIMAP2_CONFIG]
    -nlc|--nanosv_last_config                                          NanoSV config to detect SVs in the LAST mapped fusion candidates [${NANOSV_LAST_NOCONSENSUS_CONFIG}]
    -ss|--sniffles_settings                                            Settings to use for sniffles [$SNIFFLES_SETTINGS]

MAPPING
    -rf|--reffasta                                                     Reference genome fasta [${REFFASTA}]
    -rg|--refgenome                                                    Reference genome [${REFGENOME}]
    -rd|--refdict                                                      Reference genome .dict file [${REFDICT}]
    -mm2s|--minimap2_settings                                          Minimap2 settings [${MINIMAP2_SETTINGS}]
    -lms|--last_mapping_settings                                       LAST settings [${LAST_MAPPING_SETTINGS}]
    -lmt|--last_mapping_threads                                        Number of threads [${LAST_MAPPING_THREADS}]

PRIMER DESIGN
    -pdf|--primer_design_flank                                         Flanking distance around to breakpoint to extract [$PRIMER_DESIGN_FLANK]
    -pdd|--primer_design_dir                                           Path to primer3 directory [$PRIMER_DESIGN_DIR]
    -pdp|--primer_design_psr                                           PSR [$PRIMER_DESIGN_PSR]
    -pdmtf|--primer_Desing_minimal_target_flank                        Minimal bases flanking breakpoint that has to be included in primer product [$PRIMER_DESIGN_MINIMAL_TARGET_FLANK]
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

THREADS=1

NON_CODING=false
CONSENSUS_CALLING=false
DONT_CLEAN=false
DONT_FILTER=false
USE_LAST=true

OUTPUTDIR=$(realpath ./)

#TOOL PATH DEFAULTS
SAMTOOLS=$PATH_SAMTOOLS
MINIMAP2=$PATH_MINIMAP2
LAST_DIR=$PATH_LAST_DIR
WTDBG2_DIR=$PATH_WTDBG2_DIR
SV_CALLER=$PATH_SV_CALLER

#REFERENCE DEFAULTS
REFFASTA=$PATH_HOMO_SAPIENS_REFFASTA
REFGENOME=$PATH_HOMO_SAPIENS_REFGENOME
REFDICT=$PATH_HOMO_SAPIENS_REFDICT

#SV CALLING DEFAULTS
NANOSV_MINIMAP2_CONFIG=$FILES_DIR/nanosv_minimap2_config.ini
NANOSV_MINIMAP2_CONSENSUS_CONFIG=$FILES_DIR/nanosv_minimap2_consensus_config.ini
NANOSV_LAST_NOCONSENSUS_CONFIG=$FILES_DIR/nanosv_last_config.ini
NANOSV_LAST_CONSENSUS_CONFIG=$FILES_DIR/nanosv_last_consensus_config.ini

SNIFFLES_SETTINGS='-s 2 -n -1 -d 10 --genotype'

#REGION SELECTION DEFAULTS
REGION_SELECTION_SCRIPT=$SCRIPT_DIR/RegionSelection.py

#FUSION READ EXTRACTION DEFAULTS
FUSION_READ_EXTRACTION_SCRIPT=$SCRIPT_DIR/FusionReadExtraction.py
FUSION_CHECK_SCRIPT=$SCRIPT_DIR/FusionCheck.py

#CONSENSUS CALLING DEFAULTS
CONSENSUS_CALLING_WTDBG2_SETTINGS='-x ont -g 3g -q'

#MAPPING DEFAULTS
MINIMAP2_SETTINGS='-x map-ont -a --MD'

LAST_MAPPING_SETTINGS="-Q 0 -p ${LAST_DIR}/last_params"
LAST_MAPPING_THREADS=1

#PRIMER DESIGN DEFAULTS
PRIMER_DESIGN_GETSEQ_SCRIPT=$SCRIPT_DIR/PrimerFlankDesign.py
PRIMER_DESIGN_DIR=$PATH_PRIMER_DESIGN_DIR
PRIMER_DESIGN_PRIMER3_CORE=$PRIMER_DESIGN_DIR/primer3/src/primer3_core
PRIMER_DESIGN_PSR='100-200'
PRIMER_DESIGN_FLANK='200'
PRIMER_DESIGN_MINIMAL_TARGET_FLANK='10'


while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
    -h|--help)
    usage
    exit
    shift # past argument
    ;;
    -f|--fastqdir)
    FASTQDIR="$2"
    shift # past argument
    shift # past value
    ;;
    -b|--bam)
    BAM="$2"
    shift # past argument
    shift # past value
    ;;
    -n|--name)
    SAMPLE="$2"
    shift # past argument
    shift # past value
    ;;
    -v|--vcf)
    VCF="$2"
    shift # past argument
    shift # past value
    ;;
    -t|--threads)
    THREADS="$2"
    shift # past argument
    shift # past value
    ;;
    -nc|--non_coding)
    NON_CODING=true
    shift # past argument
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
    -wl|--without_last)
    USE_LAST=false
    shift # past argument
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
    -nmc|--nanosv_minimap2_config)
    NANOSV_MINIMAP2_CONFIG="$2"
    shift # past argument
    shift # past value
    ;;
    -nlc|--nanosv_last_config)
    NANOSV_LAST_CONFIG="$2"
    shift # past argument
    shift # past value
    ;;
    -ss|--sniffles_settings)
    SNIFFLES_SETTINGS="$2"
    shift # past argument
    shift # past value
    ;;
    -ccws|--consensus_calling_wtdbg2_settings)
    CONSENSUS_CALLING_WTDBG2_SETTINGS="$2"
    shift # past argument
    shift # past value
    ;;
    -rf|--reffasta)
    REFFASTA="$2"
    shift # past argument
    shift # past value
    ;;
    -rg|--refgenome)
    REFGENOME="$2"
    shift # past argument
    shift # past value
    ;;
    -rd|--refdict)
    REFDICT="$2"
    shift # past argument
    shift # past value
    ;;
    -mm2s|--minimap2_settings)
    MINIMAP2_SETTINGS="$2"
    shift # past argument
    shift # past value
    ;;
    -lms|--last_mapping_settings)
    LAST_MAPPING_SETTINGS="$2"
    shift # past argument
    shift # past value
    ;;
    -pdf|--primer_design_flank)
    PRIMER_DESIGN_FLANK="$2"
    shift # past argument
    shift # past value
    ;;
    -pdd|--primer_design_dir)
    PRIMER_DESIGN_DIR="$2"
    shift # past argument
    shift # past value
    ;;
    -pdp|--primer_design_psr)
    PRIMER_DESIGN_PSR="$2"
    shift # past argument
    shift # past value
    ;;
    -pdmtf|--primer_Desing_minimal_target_flank)
    PRIMER_DESIGN_MINIMAL_TARGET_FLANK="$2"
    PRIMER_DESIGN_PSR="$2"
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

if [ -z $BAM ] && [ -z $FASTQDIR ]; then
    echo "Missing -f|--fastqdir OR -b|--bam parameter"
    usage
    exit
fi

OUTPUTDIR=$(realpath ${OUTPUTDIR})

echo -e "`date` \t Running on `uname -n`"

if [ -z $SAMPLE ];then
  if [ -z $BAM ]; then
    SAMPLE=$(ls $FASTQDIR/*.fastq | head -n 1)
    SAMPLE=$(basename $SAMPLE)
    SAMPLE=${SAMPLE/.fastq/}
  else
    SAMPLE=$(basename $BAM)
    SAMPLE=${SAMPLE/.bam/}
  fi
fi

if [ -z $FUSION_CHECK_VCF_OUTPUT ]; then
    FUSION_CHECK_VCF_OUTPUT=$OUTPUTDIR/${SAMPLE}_FusionGenes.vcf
fi

if [ -z $FUSION_CHECK_INFO_OUTPUT ]; then
    FUSION_CHECK_INFO_OUTPUT=$OUTPUTDIR/${SAMPLE}_FusionGenesInfo.txt
fi

if [ -z $FUSION_CHECK_PDF_OUTPUT ]; then
    FUSION_CHECK_PDF_OUTPUT=$OUTPUTDIR/${SAMPLE}_FusionGenes.pdf
fi

REGION_SELECTION_BED_OUTPUT=$OUTPUTDIR/regions.bed
REGION_SELECTION_BAM_OUTPUT=$OUTPUTDIR/regions.bam
VCF_OUTPUT=$OUTPUTDIR/${SAMPLE}_FusionGenes.vcf
INFO_OUTPUT=$OUTPUTDIR/${SAMPLE}_FusionGenesInfo.txt

PIPELINE_DIR=$NANOFG_DIR/pipeline
SCRIPT_DIR=$NANOFG_DIR/scripts
CANDIDATE_DIR=$OUTPUTDIR/candidate_fusions
PRIMER_DIR=$OUTPUTDIR/primers

BAM_MERGE_OUT=$OUTPUTDIR/candidate_fusion_genes.bam
SV_CALLING_OUT=$OUTPUTDIR/candidate_fusion_genes.vcf


if [ -d "$CANDIDATE_DIR" ]; then
  rm -f $CANDIDATE_DIR/*
else
  mkdir -p $CANDIDATE_DIR
fi

if [ -d "$PRIMER_DIR" ]; then
  rm -f $PRIMER_DIR/*
else
  mkdir -p $PRIMER_DIR
fi

if [ ! -d $CANDIDATE_DIR ]; then
    exit
fi

. $VENV

################################################## MAPPING OF THE NANOPORE READS USING MINIMAP2 (OPTIONAL)
if [ ! -z $BAM ];then
  echo "### bam (-b) already provided. Skipping minimap2 mapping"
else
  echo -e "`date` \t Mapping all reads using minimap2..."

  FASTQ=${OUTPUTDIR}/${SAMPLE}.merged.fastq
  cat $FASTQDIR/*.fastq > $FASTQ
  BAM=${OUTPUTDIR}/${SAMPLE}.bam

  bash $PIPELINE_DIR/minimap2_mapping.sh \
    -f $FASTQ \
    -o $BAM \
    -mm2 $MINIMAP2 \
    -mm2s "$MINIMAP2_SETTINGS" \
    -r $REFFASTA \
    -t $THREADS \
    -s $SAMTOOLS
fi

################################################## SELECTION OF REGIONS THAT ARE SPECIFIED BY THE USER (OPTIONAL)
if [ ! -z $VCF ];then
  REGION_SELECTION_BAM_OUTPUT=$BAM
  echo "### vcf (-v) already provided. Skipping selection and sv calling"
elif [ -z $SELECTION ];then
  REGION_SELECTION_BAM_OUTPUT=$BAM
  echo "### No selection parameter (-s) provided. Using all mapped reads"
else
  echo -e "`date` \t Selecting regions to check for fusion genes..."
  python $REGION_SELECTION_SCRIPT \
  -b $REGION_SELECTION_BED_OUTPUT \
  -r $SELECTION

  if ! [ $? -eq 0 ]; then
    echo "!!! REGION SELECTION NOT CORRECTLY COMPLETED... exiting"
    exit
  fi

  REGION_SELECTION_SAM_OUTPUT=${REGION_SELECTION_BAM_OUTPUT/.bam/.sam}
  $SAMTOOLS view -H $BAM > $REGION_SELECTION_SAM_OUTPUT
  $SAMTOOLS view -@ $THREADS -L $REGION_SELECTION_BED_OUTPUT $BAM | cut -f 1 | sort -k1n | uniq > $OUTPUTDIR/reads.tmp
  $SAMTOOLS view  -@ $THREADS $BAM | grep -f reads.tmp >> $REGION_SELECTION_SAM_OUTPUT
  $SAMTOOLS view  -@ $THREADS -h -S -f bam $REGION_SELECTION_SAM_OUTPUT | $SAMTOOLS sort -o $REGION_SELECTION_BAM_OUTPUT /dev/stdin
  $SAMTOOLS index $REGION_SELECTION_BAM_OUTPUT
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
    if ! [ $? -eq 0 ]; then
      echo "!!! SV CALLING NOT CORRECTLY COMPLETED... exiting"
      exit
    fi
fi

################################################## REMOVAL OF INSERTIONS (ALWAYS) AND SVS WITHOUT THE PASS FILTER (OPTIONAL)

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

if [ $NON_CODING = true ];then
  python $FUSION_READ_EXTRACTION_SCRIPT \
    -b $BAM \
    -v $VCF_FILTERED \
    -o $CANDIDATE_DIR \
    -nc $NON_CODING
else
  python $FUSION_READ_EXTRACTION_SCRIPT \
    -b $BAM \
    -v $VCF_FILTERED \
    -o $CANDIDATE_DIR
fi

if ! [ $? -eq 0 ]; then
  echo "!!! FUSION READ EXTRACTION NOT CORRECTLY COMPLETED... exiting"
  exit
fi
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
    ln -s $(realpath $FASTA) ${FASTA/.fasta/.no_ctg.fa}
  fi
done

################################################## MAPPING OF THE READS (OR CONSENSUS) USING LAST FOR NANOSV AND MINIMAP2 FOR SNIFFLES.
################################################## LAST IS PREFERRED FOR ACCURATE EXON-EXON FUSION DETECTION
echo -e "`date` \t Mapping candidate fusion genes..."

if [[ $SV_CALLER == *"nanosv"* ]] || [[ $SV_CALLER == *"NanoSV"* ]]; then
  if [ $USE_LAST = true ]; then
    echo -e "`date` \t Using LAST for remapping"
    MAPPING_ARGS="-t $LAST_MAPPING_THREADS -r $REFGENOME -rd $REFDICT -l $LAST_DIR -ls '$LAST_MAPPING_SETTINGS' -s $SAMTOOLS"
    for FA in $CANDIDATE_DIR/*.fa; do
      echo $FA;
    done | \
    xargs -I{} --max-procs $THREADS bash -c "bash $PIPELINE_DIR/last_mapping.sh -f {} $MAPPING_ARGS; exit 1;"
  else
    echo -e "`date` \t Using minimap2 for remapping"
    MAPPING_ARGS="-mm2 $MINIMAP2 -oc -mm2s '$MINIMAP2_SETTINGS' -r $REFFASTA -t $LAST_MAPPING_THREADS -s $SAMTOOLS"
    for FA in $CANDIDATE_DIR/*.fa; do
      echo ${FA/.fa/};
    done | \
    xargs -I{} --max-procs $THREADS bash -c "bash $PIPELINE_DIR/minimap2_mapping.sh -f {}.fa -o {}.sorted.bam $MAPPING_ARGS; exit 1;"
  fi

elif [[ $SV_CALLER == *"sniffles"* ]] || [[ $SV_CALLER == *"Sniffles"* ]]; then
  echo -e "`date` \t Using minimap2 for remapping"
  MAPPING_ARGS="-mm2 $MINIMAP2 -oc -mm2s '$MINIMAP2_SETTINGS' -r $REFFASTA -t $LAST_MAPPING_THREADS -s $SAMTOOLS"
  for FA in $CANDIDATE_DIR/*.fa; do
    echo ${FA/.fa/};
  done | \
  xargs -I{} --max-procs $THREADS bash -c "bash $PIPELINE_DIR/minimap2_mapping.sh -f {}.fa -o {}.sorted.bam $MAPPING_ARGS; exit 1;"
fi
################################################## MERGING ALL SEPARATE BAM FILES OF ALL FUSION GENE CANDIDATES INTO A SINGLE BAM FILE
echo -e "`date` \t Merging bams..."

NUMBER_OF_BAMS=$(ls $CANDIDATE_DIR/*.sorted.bam | wc -l)

if [[ NUMBER_OF_BAMS -gt 1 ]];then
  $SAMTOOLS merge -f $BAM_MERGE_OUT $CANDIDATE_DIR/*.sorted.bam
  if [[ $SAMTOOLS == *"samtools"* ]] || [[ $SAMTOOLS == *"Samtools"* ]]; then
    $SAMTOOLS index $BAM_MERGE_OUT
  fi
elif [[ NUMBER_OF_BAMS -eq 1 ]];then
  cp $CANDIDATE_DIR/*.sorted.bam $BAM_MERGE_OUT
  cp $CANDIDATE_DIR/*.sorted.bam.bai ${BAM_MERGE_OUT}.bai
else
  echo "NO CANDIDATE FUSION GENES FOUND"
  exit
fi
################################################## CALLING SVS FOR THE MAPPED FUSION CANDIDATES
echo -e "`date` \t Calling SVs..."
if [ -z $NANOSV_LAST_CONFIG ]; then
  if [ $CONSENSUS_CALLING = true ];then
    if [ $USE_LAST = true ]; then
      NANOSV_CONFIG=$NANOSV_LAST_CONSENSUS_CONFIG
    else
      NANOSV_CONFIG=$NANOSV_MINIMAP2_CONSENSUS_CONFIG
    fi
    SNIFFLES_SETTINGS='-s 1 -n -1 --genotype'
  else
    if [ $USE_LAST = true ]; then
      NANOSV_CONFIG=$NANOSV_LAST_NOCONSENSUS_CONFIG
    else
      NANOSV_CONFIG=$NANOSV_MINIMAP2_CONFIG
    fi
    SNIFFLES_SETTINGS='-s 2 -n -1 --genotype'
  fi
fi

bash $PIPELINE_DIR/sv_calling.sh \
  -sv $SV_CALLER \
  -b $BAM_MERGE_OUT \
  -t $THREADS \
  -s $SAMTOOLS \
  -v $VENV \
  -c $NANOSV_CONFIG \
  -ss "$SNIFFLES_SETTINGS" \
  -o $SV_CALLING_OUT

if ! [ $? -eq 0 ]; then
  echo "!!! SV CALLING NOT CORRECTLY COMPLETED... exiting"
  exit
fi

################################################## REMOVAL OF INSERTIONS (ALWAYS) AND SVS WITHOUT THE PASS FILTER (OPTIONAL)

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

################################################### COMBINING SVS ON THE SAME READ FOR THE DETECTION OF COMPLEX FUSIONS
echo -e "`date` \t Linking and combining SVs for complex fusion detection"
VCF_COMPLEX=${OUTPUTDIR}/complex.vcf

python $SCRIPT_DIR/CombineSVs.py \
-v $SV_CALLING_OUT_FILTERED \
-b $BAM_MERGE_OUT \
-o $VCF_COMPLEX

grep -v "^#" $VCF_COMPLEX >> $SV_CALLING_OUT_FILTERED


################################################### CHECKING THE CANDIDATE FUSION GENES FOR ADDITIONAL INFORMATION (FRAME, SIMILARITY, GENE OVERLAP, ETC.)
echo -e "`date` \t Checking fusion candidates..."

if [ $NON_CODING = true ];then
  python $FUSION_CHECK_SCRIPT \
    -v $SV_CALLING_OUT_FILTERED \
    -ov $VCF \
    -o $FUSION_CHECK_VCF_OUTPUT \
    -fo $FUSION_CHECK_INFO_OUTPUT \
    -p $FUSION_CHECK_PDF_OUTPUT \
    -nc
else
  python $FUSION_CHECK_SCRIPT \
    -v $SV_CALLING_OUT_FILTERED \
    -ov $VCF \
    -o $FUSION_CHECK_VCF_OUTPUT \
    -fo $FUSION_CHECK_INFO_OUTPUT \
    -p $FUSION_CHECK_PDF_OUTPUT
fi

if ! [ $? -eq 0 ]; then
  echo "!!! FUSION CHECK NOT CORRECTLY COMPLETED... exiting"
  exit
fi

################################################### DESIGNING PRIMERS FOR THE DETECTED FUSIONS
echo -e "`date` \t Designing primers around fusion breakpoints..."

python $PRIMER_DESIGN_GETSEQ_SCRIPT \
 -v $FUSION_CHECK_VCF_OUTPUT \
 -d $PRIMER_DIR \
 -f $PRIMER_DESIGN_FLANK

if ! [ $? -eq 0 ]; then
 echo "!!! PRIMER FLANK DESIGN NOT CORRECTLY COMPLETED... exiting"
 exit
fi

if [ -d $PRIMER_DESIGN_DIR ];then
  echo -e "FUSION_ID\tFORWARD_PRIMER_ID\tFORWARD_PRIMER_SEQ\tREVERSE_PRIMER_ID\tREVERSE_PRIMER_SEQ\tPRODUCT_SIZE" > $OUTPUTDIR/${SAMPLE}_FusionGenes.primers
  for FUSION_FASTA in $PRIMER_DIR/*.fasta; do
    bash $PIPELINE_DIR/primer_design.sh \
     -f $FUSION_FASTA \
     -o $OUTPUTDIR/${SAMPLE}_FusionGenes.primers \
     -psr $PRIMER_DESIGN_PSR \
     -pdpc $PRIMER_DESIGN_PRIMER3_CORE \
     -mtf $PRIMER_DESIGN_MINIMAL_TARGET_FLANK
  done
  cd $OUTPUTDIR
else
  echo "Path to primer3 does not exist. Giving only breakpoint sequences"
  cat $PRIMER_DIR/*.fasta > $OUTPUTDIR/${SAMPLE}_FusionGenesBNDseq.fasta
fi

################################################### IF -DC IS NOT SPECIFIED, ALL FILES EXCEPT THE OUTPUT FILES ARE DELETED TO PROVIDE A CLEAN OUTPUT

if [ $DONT_CLEAN = false ];then
  if [ -f $OUTPUTDIR/reads.tmp ];then
    rm $OUTPUTDIR/reads.tmp
  fi
  rm -f $VCF_FILTERED
  rm -f $CANDIDATE_DIR/*
  rmdir -f $CANDIDATE_DIR
  rm -f $BAM_MERGE_OUT
  rm -f $BAM_MERGE_OUT.bai
  rm -f $SV_CALLING_OUT
  rm -f $SV_CALLING_OUT_FILTERED
  # rm $PRIMER_DIR/*
  # rmdir $PRIMER_DIR
fi

deactivate
echo -e "`date` \t Done"
