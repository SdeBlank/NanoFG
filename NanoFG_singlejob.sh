#!/bin/bash

usage() {
echo "

NanoFG_singlejob.sh -b BAM [-v VCF] [-s SELECTED_GENES_OR_REGIONS] [-df]

Required parameters:
    -b|--bam                                                                      Path to bam file

Optional parameters:

GENERAL
    -h|--help                                                                     Shows help-v|--vcf
    -v|--vcf		                                                                  Path to vcf file		                                                                  Path to vcf file
    -t|--threads
    -e|--venv                                                                     Path to virtual environment[${VENV}]

SELECTION AND FILTERING
    -s|--selection                                                     Select genes or areas to check for fusion genes.
                                                                       Insert a list of genes or areas, separated by a comma
                                                                       e.g. 'BRAF,TP53' or 'ENSG00000157764,ENSG00000141510' or '17:7565097-7590yy856'
    -df|--dont_filter                                                  Don't filter out all non-PASS SVs
    -cc|--consensus_calling                                            Perform consensus sequencing on the reads to decrease runtime

OUTPUT
    -o|--outputdir                                                                Path to output directory
    -io|--info_output                                                             Path to the NanoFG output info file []
    -vo|--vcf_output                                                              Path to the NanoFG output vcf file []
    -p|--pdf                                                                      Path to the NanoFG output pdf file []

REQUIRED TOOLS
    -sv|--sv_caller                                                               NanoSV or path to Sniffles [${SV_CALLER}]
    -sa|--sambamba                                                                Path to sambamba|samtools [${SAMBAMBA}]
    -l|--last_dir                                                                 Path to LAST directory [${LAST_DIR}]
    -w|--wtdbg2_dir                                                               Path to wtdbg2 directory [${WTDBG2_DIR}]

SCRIPTS
    -fres|--fusion_read_extraction_script                                         Path to the fusion_read_extraction.py script [${FUSION_READ_EXTRACTION_SCRIPT}]
    -fcs|--fusion_check_script                                                    Path to vcf_primer_filter.py [$FUSION_CHECK_SCRIPT]

CONSENSUS CALLING
    -ccws|--consensus_calling_wtdbg2_settings                                     Wtdbg2 settings [${CONSENSUS_CALLING_WTDBG2_SETTINGS}]

SV CALLING
    -nsc|--nanosv_config                                               Path to config to use for nanosv [$NANOSV_CONFIG]
    -ss|--sniffles_settings                                            Settings to use for sniffles [$SNIFFLES_SETTINGS]

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

SELECTION=""
CONSENSUS_CALLING=false
DONT_FILTER=false

OUTPUTDIR=$(realpath ./)

#TOOL PATH DEFAULTS
SAMBAMBA=$PATH_SAMBAMBA
LAST_DIR=$PATH_LAST_DIR
WTDBG2_DIR=$PATH_WTDBG2_DIR

SV_CALLER='/hpc/cog_bioinf/kloosterman/tools/NanoSV/nanosv/NanoSV.py'
NANOSV_MINIMAP2_CONFIG=$FILES_DIR/nanosv_minimap2_config.ini
NANOSV_LAST_CONFIG=$FILES_DIR/nanosv_last_config.ini
SNIFFLES_SETTINGS='-s 2 -n -1 --genotype'

#REGION SELECTION
REGION_SELECTION_SCRIPT=$SCRIPT_DIR/RegionSelection.py
REGION_SELECTION_BED_OUTPUT=$OUTPUTDIR/regions.bed
REGION_SELECTION_BAM_OUTPUT=$OUTPUTDIR/regions.bam

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
    -df|--dont_filter)
    DONT_FILTER=true
    shift # past argument
    ;;
    -e|--venv)
    VENV="$2"
    shift # past argument
    shift # past value
    ;;
    -sa|--sambamba)
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
    -sv|--sv_caller)
    SV_CALLER="$2"
    shift # past argument
    shift # past value
    ;;
    -nsc|--nanosv_config)
    NANOSV_CONFIG="$2"
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

. $VENV

##################################################
echo 'Step 1. (Optional) Selecting regions to check for fusion genes'

if [ $SELECTION = "" ] && ! [ -z $VCF ];then
  REGION_SELECTION_BAM_OUTPUT=$BAM
  echo "No selection parameter given or vcf input given. Using all mapped reads or given vcf"
else
  python $REGION_SELECTION_SCRIPT \
  -b $REGION_SELECTION_BED_OUTPUT \
  -r $SELECTION

  $SAMBAMBA view -h -f bam -L $REGION_SELECTION_BED_OUTPUT -o $REGION_SELECTION_BAM_OUTPUT $BAM
fi

##################################################
echo 'Step 2. (Optional) SV calling'
if [ -z $VCF ]; then
  VCF=${BAM/.bam/.vcf}
  bash $PIPELINE_DIR/sv_calling.sh \
    -sv $SV_CALLER \
    -b $REGION_SELECTION_BAM_OUTPUT \
    -t $THREADS \
    -s $SAMBAMBA \
    -v $VENV \
    -c $NANOSV_CONFIG \
    -ss $SNIFFLES_SETTINGS \
    -o $VCF
else
  echo "vcf already given. skipping step"
fi

##################################################
echo 'Step 3. Insertion removal and possible filtering'

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
echo 'Step 4. Candidate fusion read extraction'

python $FUSION_READ_EXTRACTION_SCRIPT \
  -b $BAM \
  -v $VCF_NO_INS \
  -o $CANDIDATE_DIR

##################################################
echo 'Step 5. Producing consensus if possible'

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
echo 'Step 6. Mapping consensus sequence (all reads if no consensus is reached)'

LAST_MAPPING_ARGS="-t $LAST_MAPPING_THREADS -r $LAST_MAPPING_REFGENOME -rd $LAST_MAPPING_REFDICT -l $LAST_DIR -ls '$LAST_MAPPING_SETTINGS' -s $SAMBAMBA"

for FA in $CANDIDATE_DIR/*.fa; do
  echo $FA;
done | \
xargs -I{} --max-procs $THREADS bash -c "echo 'Start' {}; bash $PIPELINE_DIR/last_mapping.sh -f {} $LAST_MAPPING_ARGS; echo 'Done' {}; exit 1;"

##################################################
echo 'Step 7. Merging bams'
$SAMBAMBA merge $MERGE_BAMS_OUT $CANDIDATE_DIR/*.last.sorted.bam

##################################################
echo 'Step 8. Calling SVs'

bash $PIPELINE_DIR/sv_calling.sh \
  -sv $SV_CALLER \
  -b $BAM_MERGE_OUT \
  -t $THREADS \
  -s $SAMBAMBA \
  -v $VENV \
  -c $NANOSV_MINIMAP2_CONFIG \
  -ss $SNIFFLES_SETTINGS \
  -o $SV_CALLING_OUT

###################################################
echo 'Step 9. Checking fusion candidates'

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
