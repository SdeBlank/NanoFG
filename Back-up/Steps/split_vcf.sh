#!/bin/bash

usage() {
echo "
Required parameters:
    -v|--vcf		    Path to vcf file

Optional parameters:
    -h|--help       Shows help
    -l|--LINES      Number of lines to put in each spit vcf file [$LINES]
    -o|--outputdir     VCF output directory [$OUTPUT]
"
}

POSITIONAL=()

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
    -o|--outputdir)
    SPLITDIR="$2"
    shift # past argument
    shift # past value
    ;;
    -l|--lines)
    LINES="$2"
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

if [ -z $SPLITDIR ]; then
    SPLITDIR="./split_vcf"
fi

echo `date`: Running on `uname -n`

if [ ! -d $SPLITDIR ]; then
    mkdir $SPLITDIR
fi

HEADER=$(grep "^#" $VCF)
AWK="grep -v \"^#\" $VCF | awk -v HEADER=\"\$HEADER\" 'NR%$LINES==1 { file = \"$SPLITDIR/\" int(NR/$LINES)+1 \".vcf\"; print HEADER > file } { print > file }'"
eval $AWK

echo `date`: Done
