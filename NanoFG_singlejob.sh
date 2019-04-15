#!/bin/bash

usage() {
echo "
Required parameters:
    -v|--vcf		                                                                  Path to vcf file
    -b|--bam                                                                      Path to bam file

Optional parameters:

GENERAL
    -h|--help                                                                     Shows help
    -o|--outputdir                                                                Path to output directory
    -d|--nanofg_dir                                                               Directory that contains NanoFG [${NANOFG_DIR}]
    -e|--venv                                                                     Path to virtual environment[${VENV}]
    -m|--mail                                                                     Email adress
    -dc|--dont_clean                                                              Don't clean up the intermediate files

REQUIRED TOOLS
    -n|--nanosv                                                                   Path to NanoSV [${NANOSV}]
    -s|--sambamba                                                                 Path to sambamba|samtools [${SAMBAMBA}]
    -l|--last_dir                                                                 Path to LAST directory [${LAST_DIR}]
    -w|--wtdbg2_dir                                                               Path to wtdbg2 directory [${WTDBG2_DIR}]

VCF SPLIT
    -vsl|--vcf_split_lines                                                        Override the nr of lines to split the original vcf into
    -vst|vcf_split_threads                                                        Number of threads [${VCF_SPLIT_THREADS}]
    -vshv|vcf_split_h_vmem                                                        Vcf split memory [${VCF_SPLIT_MEMORY}]
    -vshr|vcf_split_h_rt                                                          Vcf split time [${VCF_SPLIT_TIME}]

FUSION READ EXTRACTION
    -fres|--fusion_read_extraction_script                                         Path to the fusion_read_extraction.py script [${FUSION_READ_EXTRACTION_SCRIPT}]
    -fret|--fusion_read_extraction_threads                                        Number of threads [${FUSION_READ_EXTRACTION_THREADS}]
    -frehv|--fusion_read_extraction_h_vmem                                        Fusion read extraction memory [${FUSION_READ_EXTRACTION_MEMORY}]
    -freht|--fusion_read_extraction_h_rt                                          Fusion read extraction time [${FUSION_READ_EXTRACTION_TIME}]

CONSENSUS MAPPING
    -cmr|--consensus_mapping_refgenome                                            Reference genome [${REF}]
    -cmrd|--consensus_mapping_refdict                                             Reference genome .dict file [${REF_DICT}]
    -cmws|--consensus_mapping_wtdbg2_settings                                     wtdbg2 settings [${WTDBG2_SETTINGS}]
    -cmls|--consensus_mapping_last_settings                                       LAST settings [${LAST_SETTINGS}]
    -cmt|--consensus_mapping_threads                                              Number of threads [${CONSENSUS_MAPPING_THREADS}]
    -cmhv|--consensus_mapping_h_vmem                                              Consensus mapping memory [${CONSENSUS_MAPPING_MEMORY}]
    -cmhr|--consensus_mapping_h_rt                                                Consensus mapping time [${CONSENSUS_MAPPING_TIME}]

BAM MERGE
    -bmt|vcf_split_threads                                                        Number of threads [${VCF_SPLIT_THREADS}]
    -bmhv|vcf_split_h_vmem                                                        Bam merge memory [${VCF_SPLIT_MEMORY}]
    -bmhr|vcf_split_h_rt                                                          Bam merge time [${VCF_SPLIT_TIME}]

SV CALLING
    -sct|sv_calling_threads                                                       Number of threads [${SV_CALLING_THREADS}]
    -schv|sv_calling_h_vmem                                                       SV calling memory [${SV_CALLING_MEMORY}]
    -schr|sv_calling_h_rt                                                         SV calling time [${SV_CALLING_TIME}]

FUSION CHECK
    -fcio|--fusion_check_info_output                                              Path to the NanoFG output info file []
    -fcvo|--fusion_check_vcf_output                                               Path to the NanoFG output vcf file []
    -fcs|--fusion_check_script                                                    Path to vcf_primer_filter.py [$FUSION_CHECK_SCRIPT]
    -fct|--fusion_check_threads                                                   VCF output file [$FUSION_CHECK_THREADS]
    -fchv|--fusion_check_h_vmem                                                   VCF output file [$FUSION_CHECK_MEMORY]
    -fcht|--fusion_check_h_rt                                                     VCF output file [$FUSION_CHECK_TIME]
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

OUTPUTDIR=$(realpath ./)
echo $OUTPUTDIR
DONT_CLEAN=false

#TOOL PATH DEFAULTS
#SAMBAMBA=/hpc/local/CentOS7/cog_bioinf/sambamba_v0.6.5/sambamba
SAMBAMBA=$PATH_SAMBAMBA
#LAST_DIR=/hpc/cog_bioinf/kloosterman/tools/last-921
LAST_DIR=$PATH_LAST_DIR
#WTDBG2_DIR=/hpc/cog_bioinf/kloosterman/tools/wtdbg2_v2.2
WTDBG2_DIR=$PATH_WTDBG2_DIR

echo $SAMBAMBA
echo $LAST_DIR
echo $WTDBG2_DIR
#VCF SPLIT DEFAULTS
VCF_SPLIT_THREADS=1
VCF_SPLIT_TIME=0:5:0
VCF_SPLIT_MEMORY=10G

#FUSION READ EXTRACTION DEFAULTS
FUSION_READ_EXTRACTION_SCRIPT=$SCRIPT_DIR/FusionReadExtraction.py
FUSION_READ_EXTRACTION_THREADS=1
FUSION_READ_EXTRACTION_TIME=0:10:0
FUSION_READ_EXTRACTION_MEMORY=10G

#CONSENSUS MAPPING DEFAULTS
# CONSENSUS_MAPPING_REFGENOME=/hpc/cog_bioinf/GENOMES/LAST/human_GATK_GRCh37
CONSENSUS_MAPPING_REFGENOME=$PATH_HOMO_SAPIENS_REFGENOME
# CONSENSUS_MAPPING_REFDICT=/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.dict
CONSENSUS_MAPPING_REFDICT=$PATH_HOMO_SAPIENS_REFDICT
echo $CONSENSUS_MAPPING_REFGENOME
echo $PATH_HOMO_SAPIENS_REFDICT
CONSENSUS_MAPPING_WTDBG2_SETTINGS='-x ont -g 3g'
CONSENSUS_MAPPING_LAST_SETTINGS="-Q 0 -p ${LAST_DIR}/last_params"
CONSENSUS_MAPPING_THREADS=8
CONSENSUS_MAPPING_TIME=0:15:0
CONSENSUS_MAPPING_MEMORY=20G

#BAM MERGE DEFAULTS
BAM_MERGE_THREADS=1
BAM_MERGE_TIME=0:10:0
BAM_MERGE_MEMORY=40G

#SV CALLING DEFAULTS
SV_CALLING_THREADS=1
SV_CALLING_TIME=0:10:0
SV_CALLING_MEMORY=10G
SV_CALLING_CONFIG=$FILES_DIR/nanosv_last_config.ini

#FUSION CHECK DEFAULTS
FUSION_CHECK_SCRIPT=$SCRIPT_DIR/FusionCheck.py
FUSION_CHECK_THREADS=1
FUSION_CHECK_TIME=0:30:0
FUSION_CHECK_MEMORY=10G

#CHECK NANOFG DEFAULTS
CHECK_NANOFG_THREADS=1
CHECK_NANOFG_TIME=0:5:0
CHECK_NANOFG_MEMORY=5G


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
    -e|--venv)
    VENV="$2"
    shift # past argument
    shift # past value
    ;;
    -m|--mail)
    MAIL="$2"
    shift # past argument
    shift # past value
    ;;
    -dc|--dont_clean)
    DONT_CLEAN=true
    shift # past argument
    ;;
    -s|--sambamba)
    SAMBAMBA="$2"
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
    -vsl|--vcf_split_lines)
    VCF_SPLIT_LINES=$2
    shift # past argument
    shift # past value
    ;;
    -vst|vcf_split_threads)
    VCF_SPLIT_THREADS="$2"
    shift # past argument
    shift # past value
    ;;
    -vshv|vcf_split_h_vmem)
    VCF_SPLIT_MEMORY="$2"
    shift # past argument
    shift # past value
    ;;
    -vshr|vcf_split_h_rt)
    VCF_SPLIT_TIME="$2"
    shift # past argument
    shift # past value
    ;;
    -fres|--fusion_read_extraction_script)
    FUSION_READ_EXTRACTION_SCRIPT="$2"
    shift # past argument
    shift # past value
    ;;
    -fret|--fusion_read_extraction_threads)
    FUSION_READ_EXTRACTION_THREADS="$2"
    shift # past argument
    shift # past value
    ;;
    -frehv|--fusion_read_extraction_h_vmem)
    FUSION_READ_EXTRACTION_MEMORY="$2"
    shift # past argument
    shift # past value
    ;;
    -freht|--fusion_read_extraction_h_rt)
    FUSION_READ_EXTRACTION_TIME="$2"
    shift # past argument
    shift # past value
    ;;
    -cmr|--consensus_mapping_refgenome)
    CONSENSUS_MAPPING_REFGENOME="$2"
    shift # past argument
    shift # past value
    ;;
    -cmrd|--consensus_mapping_refdict)
    CONSENSUS_MAPPING_REFDICT="$2"
    shift # past argument
    shift # past value
    ;;
    -cmws|--consensus_mapping_wtdbg2_settings)
    CONSENSUS_MAPPING_WTDBG2_SETTINGS="$2"
    shift # past argument
    shift # past value
    ;;
    -cmls|--consensus_mapping_last_settings)
    CONSENSUS_MAPPING_LAST_SETTINGS="$2"
    shift # past argument
    shift # past value
    ;;
    -cmt|--consensus_mapping_threads)
    CONSENSUS_MAPPING_THREADS="$2"
    shift # past argument
    shift # past value
    ;;
    -cmhv|--consensus_mapping_h_vmem)
    CONSENSUS_MAPPING_MEMORY="$2"
    shift # past argument
    shift # past value
    ;;
    -cmhr|--consensus_mapping_h_rt)
    CONSENSUS_MAPPING_TIME="$2"
    shift # past argument
    shift # past value
    ;;
    -bmt|vcf_split_threads)
    BAM_MERGE_THREADS="$2"
    shift # past argument
    shift # past value
    ;;
    -bmhv|vcf_split_h_vmem)
    BAM_MERGE_MEMORY="$2"
    shift # past argument
    shift # past value
    ;;
    -bmhr|vcf_split_h_rt)
    BAM_MERGE_TIME="$2"
    shift # past argument
    shift # past value
    ;;
    -sct|sv_calling_threads)
    SV_CALLING_THREADS="$2"
    shift # past argument
    shift # past value
    ;;
    -schv|sv_calling_h_vmem)
    SV_CALLING_MEMORY="$2"
    shift # past argument
    shift # past value
    ;;
    -schr|sv_calling_h_rt)
    SV_CALLING_TIME="$2"
    shift # past argument
    shift # past value
    ;;
    -fcio|--fusion_check_info_output)
    FUSION_CHECK_INFO_OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
    -fcvo|--fusion_check_vcf_output)
    FUSION_CHECK_VCF_OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
    -fcs|--fusion_check_script)
    FUSION_CHECK_SCRIPT="$2"
    shift # past argument
    shift # past value
    ;;
    -fct|--fusion_check_threads)
    FUSION_CHECK_THREADS="$2"
    shift # past argument
    shift # past value
    ;;
    -fchv|--fusion_check_h_vmem)
    FUSION_CHECK_MEMORY="$2"
    shift # past argument
    shift # past value
    ;;
    -fcht|--fusion_check_h_rt)
    FUSION_CHECK_TIME="$2"
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

if [ -z $VCF ]; then
    echo "Missing -v|--vcf parameter"
    usage
    exit
fi
if [ -z $BAM ]; then
    echo "Missing -b|--bam parameter"
    usage
    exit
fi

if [ -z $FUSION_CHECK_VCF_OUTPUT ]; then
    FUSION_CHECK_VCF_OUTPUT=./$(basename $VCF)
    FUSION_CHECK_VCF_OUTPUT=$OUTPUTDIR/${FUSION_CHECK_VCF_OUTPUT/.vcf/_FusionGenes.vcf}
fi

if [ -z $FUSION_CHECK_INFO_OUTPUT ]; then
    FUSION_CHECK_INFO_OUTPUT=./$(basename $VCF)
    FUSION_CHECK_INFO_OUTPUT=$OUTPUTDIR/${FUSION_CHECK_INFO_OUTPUT/.vcf/_FusionGenesInfo.txt}
fi

echo `date`: Running on `uname -n`

VCF_NAME=$(basename $VCF)
VCF_NAME=${VCF_NAME/.vcf/}
VCF_OUTPUT=$OUTPUTDIR/${VCF_NAME}_FusionGenes.vcf
INFO_OUTPUT=$OUTPUTDIR/${VCF_NAME}_FusionGenesInfo.txt

PIPELINE_DIR=$NANOFG_DIR/pipeline
SCRIPT_DIR=$NANOFG_DIR/scripts

SPLITDIR=$OUTPUTDIR/split_vcf

MERGE_BAMS_OUT=$OUTPUTDIR/consensus_last.sorted.bam

SV_CALLING_OUT=$OUTPUTDIR/consensus_nanosv.vcf

mkdir -p $SPLITDIR
if [ ! -d $SPLITDIR ]; then
    exit
fi

##################################################

VCF_NO_INS=${VCF/.vcf/_noINS.vcf}
VCF_NO_INS=${OUTPUTDIR}/$(basename $VCF_NO_INS)

grep "^#" $VCF > $VCF_NO_INS
grep -v "^#" $VCF | awk '$5!="<INS>"' >> $VCF_NO_INS

##################################################

. $VENV
python $FUSION_READ_EXTRACTION_SCRIPT \
  -b $BAM \
  -v $VCF_NO_INS \
  -o $SPLITDIR

##################################################

for FASTA in $SPLITDIR/*.fasta; do
  bash $PIPELINE_DIR/consensus_mapping.sh \
  -f $FASTA \
  -t $CONSENSUS_MAPPING_THREADS \
  -r $CONSENSUS_MAPPING_REFGENOME \
  -rd $CONSENSUS_MAPPING_REFDICT \
  -w $WTDBG2_DIR \
  -ws '$CONSENSUS_MAPPING_WTDBG2_SETTINGS' \
  -l $LAST_DIR \
  -ls '$CONSENSUS_MAPPING_LAST_SETTINGS' \
  -s $SAMBAMBA

##################################################

bash $PIPELINE_DIR/bam_merge.sh \
  -d $SPLITDIR \
  -s $SAMBAMBA \
  -o $MERGE_BAMS_OUT

##################################################

bash $PIPELINE_DIR/sv_calling.sh \
  -b $MERGE_BAMS_OUT \
  -n $NANOSV \
  -t $SV_CALLING_THREADS \
  -s $SAMBAMBA \
  -v $VENV \
  -c $SV_CALLING_CONFIG \
  -o $SV_CALLING_OUT

###################################################

bash $PIPELINE_DIR/fusion_check.sh \
  -v $SV_CALLING_OUT \
  -o $FUSION_CHECK_VCF_OUTPUT \
  -fo $FUSION_CHECK_INFO_OUTPUT \
  -s $FUSION_CHECK_SCRIPT \
  -e $VENV

###################################################
