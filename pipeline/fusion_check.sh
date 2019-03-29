#!/bin/bash

usage() {
echo "
Required parameters:
    -v|--vcf		    Path to vcf file


Optional parameters:
    -h|--help       Shows help
    -o|--output     Vcf output file [_FusionGenes.vcf]
    -fo|--fusion_output    Txt output file [_FusionGenesInfo.txt]
    -s|--script     Path to vcf_primer_filter.py [$SCRIPT]
    -e|--venv       Path to virtual environment[$VENV]
"
}

POSITIONAL=()

NANOFG_DIR=$(realpath $(dirname $(dirname ${BASH_SOURCE[0]})))
echo $NANOFG_DIR
SCRIPT=$NANOFG_DIR/scripts/NanoFG.py
VENV=$NANOFG_DIR/venv/bin/activate

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
    -o|--output)
    OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
    -fo|--fusion_output)
    FUSION_OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
    -s|--script)
    SCRIPT="$2"
    shift # past argument
    shift # past value
    ;;
    -e|--venv)
    VENV="$2"
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

if [ -z $OUTPUT ]; then
    OUTPUT=./$(basename $VCF)
    OUTPUT=${OUTPUT/.vcf/_FusionGenes.vcf}
fi

if [ -z $FUSION_OUTPUT ]; then
    FUSION_OUTPUT=./$(basename $VCF)
    FUSION_OUTPUT=${FUSION_OUTPUT/.vcf/_FusionGenesInfo.txt}
fi

echo `date`: Running on `uname -n`

. $VENV

python $SCRIPT -v $VCF -o $OUTPUT -fo $FUSION_OUTPUT

echo `date`: Done
