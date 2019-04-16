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
-s|--sambamba             Path to sambamba|samtools [${SAMBAMBA}]
"
exit
}

POSITIONAL=()

# DEFAULT SETTINGS
THREADS=1
REF=/hpc/cog_bioinf/GENOMES/LAST/human_GATK_GRCh37
REF_DICT=/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.dict
LAST_DIR=/hpc/cog_bioinf/kloosterman/tools/last-921
LAST_SETTINGS="-Q 0 -p ${LAST_PARAMS}"
MAF_CONVERT=${LAST_DIR}/scripts/maf-convert
SAMBAMBA=/hpc/local/CentOS7/cog_bioinf/sambamba_v0.6.5/sambamba

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
    -rd|--refdict)
    REF_DICT="$2"
    shift # past argument
    shift # past value
    -l|--last_dir)
    LAST_DIR="$2"
    shift # past argument
    shift # past value
    ;;
    -ls|--last_settings)
    LAST_SETTINGS="$2"
    LAST_SETTINGS_OVERRIDE=True                 ##### FIND BETTER WAY???
    shift # past argument
    shift # past value
    ;;
    -s|--sambamba)
    SAMBAMBA="$2"
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

echo `date`: Running on `uname -n`

LASTAL=${LAST_DIR}/src/lastal
LAST_SPLIT=${LAST_DIR}/src/last-split
LAST_PARAMS=${LAST_DIR}/last_params
MAF_CONVERT=${LAST_DIR}/scripts/maf-convert

if [ -z $LAST_SETTINGS_OVERRIDE ];then
  LAST_SETTINGS=$(echo $LAST_SETTINGS | sed -e "s/-p [^ ]\+/-p ${LAST_PARAMS}/")
fi

LAST_COMMAND="${LASTAL} ${LAST_SETTINGS} ${REF} ${PREFIX}.ctg.fa | ${LAST_SPLIT} | ${MAF_CONVERT} -f ${REF_DICT} sam /dev/stdin | ${SAMBAMBA} view -S -f bam /dev/stdin | \
${SAMBAMBA} sort /dev/stdin -o ${PREFIX}.ctg.last.sorted.bam"
SAMBAMBA_INDEX_COMMAND="${SAMBAMBA} index ${PREFIX}.ctg.last.sorted.bam"

#echo ${LAST_COMMAND}
eval ${LAST_COMMAND}
#echo ${SAMBAMBA_INDEX_COMMAND}
eval ${SAMBAMBA_INDEX_COMMAND}
