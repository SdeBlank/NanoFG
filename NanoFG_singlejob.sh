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
    -io|--info_output                                                             Path to the NanoFG output info file []
    -vo|--vcf_output                                                              Path to the NanoFG output vcf file []
    -d|--nanofg_dir                                                               Directory that contains NanoFG [${NANOFG_DIR}]
    -t|--threads
    -e|--venv                                                                     Path to virtual environment[${VENV}]

REQUIRED TOOLS
    -n|--nanosv                                                                   Path to NanoSV [${NANOSV}]
    -s|--sambamba                                                                 Path to sambamba|samtools [${SAMBAMBA}]
    -l|--last_dir                                                                 Path to LAST directory [${LAST_DIR}]
    -w|--wtdbg2_dir                                                               Path to wtdbg2 directory [${WTDBG2_DIR}]

SCRIPTS
    -fres|--fusion_read_extraction_script                                         Path to the fusion_read_extraction.py script [${FUSION_READ_EXTRACTION_SCRIPT}]
    -fcs|--fusion_check_script                                                    Path to vcf_primer_filter.py [$FUSION_CHECK_SCRIPT]

CONSENSUS MAPPING
    -r|--consensus_marefgenome                                                    Reference genome [${REF}]
    -rd|-consensus_mapping_refdict                                                Reference genome .dict file [${REF_DICT}]
    -ws|--wtdbg2_settings                                                         wtdbg2 settings [${WTDBG2_SETTINGS}]
    -ls|--last_settings                                                           LAST settings [${LAST_SETTINGS}]
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

#TOOL PATH DEFAULTS
SAMBAMBA=$PATH_SAMBAMBA
LAST_DIR=$PATH_LAST_DIR
WTDBG2_DIR=$PATH_WTDBG2_DIR

#FUSION READ EXTRACTION DEFAULTS
FUSION_READ_EXTRACTION_SCRIPT=$SCRIPT_DIR/FusionReadExtraction.py

#CONSENSUS MAPPING DEFAULTS
CONSENSUS_MAPPING_REFGENOME=$PATH_HOMO_SAPIENS_REFGENOME
CONSENSUS_MAPPING_REFDICT=$PATH_HOMO_SAPIENS_REFDICT
CONSENSUS_MAPPING_WTDBG2_SETTINGS='-x ont -g 3g'
CONSENSUS_MAPPING_LAST_SETTINGS="-Q 0 -p ${LAST_DIR}/last_params"
CONSENSUS_MAPPING_THREADS=1
CONSENSUS_MAPPING_TIME=0:15:0
CONSENSUS_MAPPING_MEMORY=20G

SV_CALLING_CONFIG=$FILES_DIR/nanosv_last_config.ini

#FUSION CHECK DEFAULTS
FUSION_CHECK_SCRIPT=$SCRIPT_DIR/FusionCheck.py



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
    -fres|--fusion_read_extraction_script)
    FUSION_READ_EXTRACTION_SCRIPT="$2"
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
done

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
