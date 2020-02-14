#!/bin/bash

# Create usage message
usage(){
echo "
Required parameters:
-f|--fasta                Input fasta file [${FASTA}]

Optional parameters:
-h|--help                 Shows help
-t|--threads              Number of threads [${THREADS}]
-r|--refgenome            Reference genome [${REF}]
-rd|--refdict             Reference genome .dict file [${REF_DICT}]
-l|--last_dir             Path to LAST directory [${LAST_DIR}]
-ls|--last_settings       LAST settings [${LAST_SETTINGS}]
-s|--sambamba             Path to sambamba|samtools [${SAMTOOLS}]
"
exit
}

POSITIONAL=()
NANOFG_DIR=$(realpath $(dirname $(dirname ${BASH_SOURCE[0]})))
source ${NANOFG_DIR}/paths.ini
# DEFAULT SETTINGS
THREADS=1
REF=${PATH_HOMO_SAPIENS_REFGENOME}
REF_DICT=${PATH_HOMO_SAPIENS_REFDICT}
LAST_DIR=${PATH_LAST_DIR}
LASTAL=${LAST_DIR}/src/lastal
LAST_SPLIT=${LAST_DIR}/src/last-split
LAST_PARAMS=${LAST_DIR}/last_params
MAF_CONVERT=${LAST_DIR}/scripts/maf-convert
LAST_SETTINGS="-Q 0 -p ${LAST_PARAMS}"
MAF_CONVERT=${LAST_DIR}/scripts/maf-convert
SAMTOOLS=${PATH_SAMTOOLS}

while [[ $# -gt 0 ]]; do
  KEY="$1"
  case ${KEY} in
    -h|--help)
    usage
    shift # past argument
    ;;
    -f|--fasta)
    FASTA="$2"
    shift # past argument
    shift # past value
    ;;
    -t|--threads)
    THREADS="$2"
    shift # past argument
    shift # past value
    ;;
    -r|--refgenome)
    REF="$2"
    shift # past argument
    shift # past value
    ;;
    -rd|--refdict)
    REF_DICT="$2"
    shift # past argument
    shift # past value
    ;;
    -l|--last_dir)
    LAST_DIR="$2"
    shift # past argument
    shift # past value
    ;;
    -ls|--last_settings)
    LAST_SETTINGS="$2"
    LAST_SETTINGS_OVERRIDE=True
    shift # past argument
    shift # past value
    ;;
    -s|--sambamba)
    SAMTOOLS="$2"
    shift # past argument
    shift # past value
    ;;
    *) # Unknown option
    POSITIONAL+=("$1") # Store unknown option
    shift
    ;;
  esac
done

set -- "${POSITIONAL[@]}" # Restore postional parameters

if [ -z ${FASTA} ]; then
  echo "Missing -f|--fasta parameter"
  usage
fi

LASTAL=${LAST_DIR}/src/lastal
LAST_SPLIT=${LAST_DIR}/src/last-split
LAST_PARAMS=${LAST_DIR}/last_params
MAF_CONVERT=${LAST_DIR}/scripts/maf-convert
LAST_SETTINGS="-Q 0 -p ${LAST_PARAMS}"
PREFIX=${FASTA/.fa/}

if [ -z $LAST_SETTINGS_OVERRIDE ];then
  LAST_SETTINGS=$(echo $LAST_SETTINGS | sed -e "s/-p [^ ]\+/-p ${LAST_PARAMS}/")
fi

LAST_COMMAND="${LASTAL} ${LAST_SETTINGS} ${REF} ${PREFIX}.fa | ${LAST_SPLIT} | ${MAF_CONVERT} -f ${REF_DICT} sam /dev/stdin | ${SAMTOOLS} view /dev/stdin -b | \
${SAMTOOLS} sort /dev/stdin -o ${PREFIX}.last.sorted.bam"
SAMTOOLS_INDEX_COMMAND="${SAMTOOLS} index ${PREFIX}.last.sorted.bam"

#echo ${LAST_COMMAND}
eval ${LAST_COMMAND}
#echo ${SAMTOOLS_INDEX_COMMAND}
eval ${SAMTOOLS_INDEX_COMMAND}
