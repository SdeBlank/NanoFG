#!/bin/bash

usage() {
echo "
Required parameters:
    -v|--vcf		    Path to vcf file

Optional parameters:
    -h|--help               Shows help
    -o|--outputdir                                                                Path to output directory
    -d|--split_directory        directory that contains NanoFG [$NANOFG_DIR]
    -l|--lines     Number of lines to put in each spit vcf file [Devides vcf in 50 files]
    -n|--number_of_files
"
}

POSITIONAL=()

#DEFAULTS
OUTPUTDIR=$(realpath ./)
SPLITDIR=${OUTPUTDIR}/split_vcf
NUMBER_OF_FILES=15

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
    -o|--vcf)
    VCF="$2"
    shift # past argument
    shift # past value
    ;;
    -d|--split_directory)
    SPLITDIR="$2"
    shift # past argument
    shift # past value
    ;;
    -l|--lines)
    LINES="$2"
    shift # past argument
    shift # past value
    ;;
    -n|--number_of_files)
    NUMBER_OF_FILES="$2"
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

if [ -z $LINES ]; then
  NUMBER_OF_SVS=$(grep -vc "^#" $VCF | grep -oP "(^\d+)")
  ((LINES = ($NUMBER_OF_SVS + $NUMBER_OF_FILES - 1) / $NUMBER_OF_FILES))
  if [ $LINES -lt 100 ]; then
    LINES=100
  fi
fi

echo `date`: Running on `uname -n`

HEADER=$(grep "^#" $VCF)
AWK="grep -v \"^#\" $VCF | awk -v HEADER=\"\$HEADER\" 'NR%$LINES==1 { file = \"$SPLITDIR/\" int(NR/$LINES)+1 \".vcf\"; print HEADER > file } { print > file }'"

eval $AWK

echo `date`: Done
