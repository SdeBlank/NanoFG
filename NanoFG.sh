#!/bin/bash

usage() {
echo "
Required parameters:
    -v|--vcf		    Path to vcf file

Optional parameters:
    -h|--help       Shows help
    -s|--script     Path to vcf_primer_filter.py [$SCRIPT]
    -l|--LINES      Number of lines to put in each split vcf file [$LINES]
    -o|--output     VCF output file [$OUTPUT]
    -vo|--vcf_output  VCF output file
    -e|--venv       Path to virtual environment[$VENV]
"
}

POSITIONAL=()

SCRIPT="/home/cog/sdeblank/Documents/github/NanoFG/NanoFG.py"
VENV="/data/sharc/venv/bin/activate"
LINES=100


while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
    -h|--help)
    usage
    shift # past argument
    ;;
    -v|--vcf)
    VCF="$2"
    shift # past argument
    shift # past value
    ;;
    -s|--script)
    SCRIPT="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output)
    OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
    -vo|--vcf_output)
    VCF_OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
    -l|--lines)
    LINES="$2"
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
    OUTPUT=${OUTPUT/.vcf/_FusionGenes.txt}
fi

if [ -z $VCF_OUTPUT ]; then
    VCF_OUTPUT=./$(basename $VCF)
    VCF_OUTPUT=${VCF_OUTPUT/.vcf/_FusionGenes.vcf}
fi

SPLITDIR="$(dirname $OUTPUT)/split_vcf"
JOBDIR="$(dirname $OUTPUT)/jobs"

echo `date`: Running on `uname -n`

. $VENV

if [ $LINES != None ]; then
  if [ ! -d $SPLITDIR ]; then
      mkdir $SPLITDIR
  fi

  if [ ! -d $JOBDIR ]; then
      mkdir $JOBDIR
  fi

  HEADER=$(grep "^#" $VCF)
  AWK="grep -v \"^#\" $VCF | awk -v HEADER=\"\$HEADER\" 'NR%$LINES==1 { file = \"$SPLITDIR/\" int(NR/$LINES)+1 \".vcf\"; print HEADER > file } { print > file }'"
  eval $AWK

  NUMBER_OF_LINES_VCF_1=$(grep -v "^#" $VCF | wc -l | grep -oP "(^\d+)")
  NUMBER_OF_LINES_VCF_2=$(cat $SPLITDIR/*.vcf | grep -v "^#" | wc -l | grep -oP "(^\d+)")

  if [ $NUMBER_OF_LINES_VCF_1 == $NUMBER_OF_LINES_VCF_2 ]; then
    for SPLIT_VCF in $SPLITDIR/*.vcf; do
      SPLIT_OUTPUT=${SPLIT_VCF/.vcf/_FusionGenes.txt}
      SPLIT_VCF_OUTPUT=${SPLIT_VCF/.vcf/_FusionGenes.vcf}
      python $SCRIPT -v $SPLIT_VCF -fo $SPLIT_OUTPUT -o $SPLIT_VCF_OUTPUT
    done
  fi

  NUMBER_SPLIT_VCF=$(ls -l $SPLITDIR/* | grep -cv "FusionGenes" | grep -oP "(^\d+)")
  NUMBER_SPLIT_OUTPUT=$(ls -l $SPLITDIR/* | grep -c "FusionGenes.txt" | grep -oP "(^\d+)")

  if [ $NUMBER_SPLIT_VCF == $NUMBER_SPLIT_OUTPUT ]; then
    head -n 1 $SPLITDIR/1_FusionGenes.txt > $OUTPUT
    for fusion_output in $SPLITDIR/*FusionGenes.txt; do
      tail -n +2 $fusion_output >> $OUTPUT
    done
  fi
else
  python $SCRIPT -v $VCF -fo $OUTPUT -o $VCF_OUTPUT
fi

echo `date`: Done
