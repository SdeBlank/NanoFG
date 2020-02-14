#!/bin/bash

usage() {
echo "
Required parameters:

Optional parameters:
    -h|--help                   Shows help
    -d|--split_directory        Directory of the split vcf files [${SPLITDIR}]
    -s|--sambamba_path          Path to sambamba|samtools [${SAMBAMBA}]
    -o|--output                 Path to merged output bam [$OUTPUT]
"
}

NANOFG_DIR=$(realpath $(dirname $(dirname ${BASH_SOURCE[0]})))
source ${NANOFG_DIR}/paths.ini

POSITIONAL=()
SAMBAMBA=${PATH_SAMTOOLS}
OUTPUT='/dev/stdout'
SPLITDIR=./split_vcf

while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
    -h|--help)
    usage
    exit
    shift # past argument
    ;;
    -d|--split_directory)
    SPLITDIR="$2"
    shift # past argument
    shift # past value
    ;;
    -s|--sambamba_path)
    SAMBAMBA="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output)
    OUTPUT="$2"
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

#echo `date`: Running on `uname -n`

#echo "$SAMBAMBA merge $OUTPUT $SPLITDIR/*.last.sorted.bam"
$SAMBAMBA merge $OUTPUT $SPLITDIR/*.last.sorted.bam

#echo `date`: Done
