#!/bin/bash

usage() {
echo "
Required parameters:
    -v|--vcf		    Path to vcf file
    -b|--bam        Path to bam file

Optional parameters:
    -h|--help
    -s|--script
    -d|--read_directory
    -e|--venv
"
}

POSITIONAL=()

#DEFAULTS
READDIR=./
NANOFG_DIR=$(realpath $(dirname $(dirname {BASH_SOURCE[0]})))
VENV=${NANOFG_DIR}/venv/bin/activate
SCRIPT=${NANOFG_DIR}/scripts/fusion_gene_read_extraction.py

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
    -d|--read_directory)
    READDIR="$2"
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

echo `date`: Running on `uname -n`

. $VENV

python $SCRIPT -v $VCF -b $BAM -o $READDIR

echo `date`: Done
