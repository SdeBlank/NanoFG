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
    -p|--pdf                                                                      Path to the NanoFG output pdf file []
    -d|--nanofg_dir                                                               Directory that contains NanoFG [${NANOFG_DIR}]
    -t|--threads
    -cc|--consensus_calling                                                       Perform consensus sequencing on the reads to decrease runtime
    -df|--dont_filter                                                             Don't filter out all non-PASS SVs
    -e|--venv                                                                     Path to virtual environment[${VENV}]

REQUIRED TOOLS
    -n|--nanosv                                                                   Path to NanoSV [${NANOSV}]
    -s|--sambamba                                                                 Path to sambamba|samtools [${SAMBAMBA}]
    -l|--last_dir                                                                 Path to LAST directory [${LAST_DIR}]
    -w|--wtdbg2_dir                                                               Path to wtdbg2 directory [${WTDBG2_DIR}]

SCRIPTS
    -fres|--fusion_read_extraction_script                                         Path to the fusion_read_extraction.py script [${FUSION_READ_EXTRACTION_SCRIPT}]
    -fcs|--fusion_check_script                                                    Path to vcf_primer_filter.py [$FUSION_CHECK_SCRIPT]

CONSENSUS CALLING
    -ccws|--consensus_calling_wtdbg2_settings                                     Wtdbg2 settings [${CONSENSUS_CALLING_WTDBG2_SETTINGS}]

LAST MAPPING
    -lmr|--last_mapping_refgenome                                                 Reference genome [${LAST_MAPPING_REFGENOME}]
    -lmrd|--last_mapping_refdict                                                  Reference genome .dict file [${LAST_MAPPING_REFDICT}]
    -lms|--last_mapping_settings                                                  LAST settings [${LAST_MAPPING_SETTINGS}]
    -lmt|--last_mapping_threads                                                   Number of threads [${LAST_MAPPING_THREADS}]
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
DONT_FILTER=false

OUTPUTDIR=$(realpath ./)

#TOOL PATH DEFAULTS
SAMBAMBA=$PATH_SAMBAMBA
LAST_DIR=$PATH_LAST_DIR
WTDBG2_DIR=$PATH_WTDBG2_DIR

#FUSION READ EXTRACTION DEFAULTS
FUSION_READ_EXTRACTION_SCRIPT=$SCRIPT_DIR/FusionReadExtraction.py

#CONSENSUS CALLING DEFAULTS
CONSENSUS_CALLING_WTDBG2_SETTINGS='-x ont -g 3g'

#LAST MAPPING DEFAULTS
LAST_MAPPING_REFGENOME=$PATH_HOMO_SAPIENS_REFGENOME
LAST_MAPPING_REFDICT=$PATH_HOMO_SAPIENS_REFDICT
LAST_MAPPING_SETTINGS="-Q 0 -p ${LAST_DIR}/last_params"
LAST_MAPPING_THREADS=1

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
    -df|--dont_filter)
    DONT_FILTER=true
    shift # past argument
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
    FUSION_CHECK_VCF_OUTPUT=$(basename $VCF)
    FUSION_CHECK_VCF_OUTPUT=$OUTPUTDIR/${FUSION_CHECK_VCF_OUTPUT/.vcf/_FusionGenes.vcf}
fi

if [ -z $FUSION_CHECK_INFO_OUTPUT ]; then
    FUSION_CHECK_INFO_OUTPUT=$(basename $VCF)
    FUSION_CHECK_INFO_OUTPUT=$OUTPUTDIR/${FUSION_CHECK_INFO_OUTPUT/.vcf/_FusionGenesInfo.txt}
fi

if [ -z $FUSION_CHECK_PDF_OUTPUT ]; then
    FUSION_CHECK_PDF_OUTPUT=$(basename $VCF)
    FUSION_CHECK_PDF_OUTPUT=$OUTPUTDIR/${FUSION_CHECK_PDF_OUTPUT/.vcf/_FusionGenes.pdf}
fi

echo `date`: Running on `uname -n`

VCF_NAME=$(basename $VCF)
VCF_NAME=${VCF_NAME/.vcf/}
VCF_OUTPUT=$OUTPUTDIR/${VCF_NAME}_FusionGenes.vcf
INFO_OUTPUT=$OUTPUTDIR/${VCF_NAME}_FusionGenesInfo.txt

PIPELINE_DIR=$NANOFG_DIR/pipeline
SCRIPT_DIR=$NANOFG_DIR/scripts
CANDIDATE_DIR=$OUTPUTDIR/candidate_fusions

MERGE_BAMS_OUT=$OUTPUTDIR/consensus_last.sorted.bam
SV_CALLING_OUT=$OUTPUTDIR/consensus_nanosv.vcf

mkdir -p $CANDIDATE_DIR
if [ ! -d $CANDIDATE_DIR ]; then
    exit
fi

##################################################
echo 'Step 1. Insertion removal and possible filtering'

VCF_NO_INS=${VCF/.vcf/_noINS.vcf}
VCF_NO_INS=${OUTPUTDIR}/$(basename $VCF_NO_INS)

if [ $DONT_FILTER = false ];then
  grep "^#" $VCF > $VCF_NO_INS
  grep -v "^#" $VCF | awk '$5!="<INS>"' | awk '$7=="PASS"' >> $VCF_NO_INS
else
  grep "^#" $VCF > $VCF_NO_INS
  grep -v "^#" $VCF | awk '$5!="<INS>"' >> $VCF_NO_INS
fi

##################################################
echo 'Step 2. Candidate fusion read extraction'

. $VENV

python $FUSION_READ_EXTRACTION_SCRIPT \
  -b $BAM \
  -v $VCF_NO_INS \
  -o $CANDIDATE_DIR

##################################################
echo 'Step 3. Producing consensus if possible'

for FASTA in $CANDIDATE_DIR/*.fasta; do
  if [ $CONSENSUS_CALLING = true ]
    bash $PIPELINE_DIR/consensus_calling.sh \
      -f $FASTA \
      -t $THREADS \
      -w $WTDBG2_DIR \
      -ws  "$CONSENSUS_CALLING_WTDBG2_SETTINGS"
  fi

  if ! [ -s ${FASTA/.fasta/_wtdbg2.ctg.fa} ];then
    ln -s $FASTA ${FASTA/.fasta/.no_ctg.fa}
  fi
done

##################################################
echo 'Step 4. Mapping consensus sequence (all reads if no consensus is reached)'

LAST_MAPPING_ARGS="-t $LAST_MAPPING_THREADS -r $LAST_MAPPING_REFGENOME -rd $LAST_MAPPING_REFDICT -l $LAST_DIR -ls '$LAST_MAPPING_SETTINGS' -s $SAMBAMBA"

for FA in $CANDIDATE_DIR/*.fa; do
  echo $FA;
done | \
xargs -I{} --max-procs $THREADS bash -c "echo 'Start' {}; bash $PIPELINE_DIR/last_mapping.sh -f {} $LAST_MAPPING_ARGS; echo 'Done' {}; exit 1;"

##################################################
echo 'Step 5. Merging bams'
$SAMBAMBA merge $MERGE_BAMS_OUT $CANDIDATE_DIR/*.last.sorted.bam

##################################################
echo 'Step 6. Calling SVs'

bash $PIPELINE_DIR/sv_calling.sh \
  -b $MERGE_BAMS_OUT \
  -n $NANOSV \
  -t $THREADS \
  -s $SAMBAMBA \
  -v $VENV \
  -c $SV_CALLING_CONFIG \
  -o $SV_CALLING_OUT

###################################################
echo 'Step 7. Checking fusion candidates'

bash $PIPELINE_DIR/fusion_check.sh \
  -v $SV_CALLING_OUT \
  -ov $VCF \
  -o $FUSION_CHECK_VCF_OUTPUT \
  -fo $FUSION_CHECK_INFO_OUTPUT \
  -p $FUSION_CHECK_PDF_OUTPUT \
  -s $FUSION_CHECK_SCRIPT \
  -e $VENV

###################################################

deactivate
