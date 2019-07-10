#!/bin/bash

# Create usage message
usage(){
echo "
Required parameters:
-f|--fastq                Input fastq file
-o|--output               Name of bam output file

Optional parameters:
-h|--help                 Shows help
-mm2|--minimap2           Path to minimap2 [$MINIMAP2]
-mm2s|--minimap2          Minimap2 parameters [$MINIMAP2_SETTINGS]
-t|--threads              Number of threads [${THREADS}]
-r|--reffasta             Reference genome [${REFFASTA}]
-s|--sambamba             Path to sambamba|samtools [${SAMTOOLS}]
"
exit
}

POSITIONAL=()

NANOFG_DIR=$(realpath $(dirname ${BASH_SOURCE[0]}))/../
source $NANOFG_DIR/paths.ini

# DEFAULT SETTINGS
THREADS=1
REF=$PATH_HOMO_SAPIENS_REFFASTA
SAMTOOLS=$PATH_SAMTOOLS
MINIMAP2=$PATH_MINIMAP2
MINIMAP2_SETTINGS='-ax map-ont'

while [[ $# -gt 0 ]]; do
  KEY="$1"
  case ${KEY} in
    -h|--help)
    usage
    shift # past argument
    ;;
    -f|--fastq)
    FASTQ="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output)
    OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
    -mm2|--minimap2)
    MINIMAP2="$2"
    shift # past argument
    shift # past value
    ;;
    -mm2s|--minimap2)
    MINIMAP2_SETTINGS="$2"
    shift # past argument
    shift # past value
    ;;
    -t|--threads)
    THREADS="$2"
    shift # past argument
    shift # past value
    ;;
    -r|--reffasta)
    REFFASTA="$2"
    shift # past argument
    shift # past value
    ;;
    -s|--samtools)
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

if [ -z ${FASTQ} ]; then
  echo "Missing -f|--fastq parameter"
  usage
fi

if [ -z ${OUTPUT} ]; then
  echo "Missing -o|--output parameter"
  usage
fi

#echo `date`: Running on `uname -n`

$MINIMAP2 -t $THREADS $MINIMAP2_SETTINGS $REFFASTA $FASTQ\
| $SAMTOOLS view -h -S -b -@ $THREADS /dev/stdin \
| $SAMTOOLS sort /dev/stdin -o ${OUTPUT} -@ $THREADS

$SAMTOOLS index ${OUTPUT} -@ $THREADS

#echo `date`: Done
