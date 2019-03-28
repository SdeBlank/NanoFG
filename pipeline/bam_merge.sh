#!/bin/bash

usage() {
echo "
Required parameters:
    -d|--split_directory
    -s|--sambamba_path
    -o|--output

Optional parameters:
    -h|--help       Shows help
"
}

POSITIONAL=()

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
    -s|--sambamba_path)
    SAMBAMBA="$2"
    shift # past argument
    shift # past value
    -o|--output)
    OUTPUT="$2"
    shift # past argument
    shift # past value
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
    esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

if [ -z $SPLITDIR ]; then
    echo "Missing -d|--split_directory"
    usage
    exit
fi

if [ -z $SAMBAMBA ]; then
    echo "Missing -s|--sambamba_path"
    usage
    exit
fi

if [ -z $OUTPUT ]; then
    echo "Missing -o|--output"
    usage
    exit
fi

echo `date`: Running on `uname -n`

for SPLIT_BAM in $SPLITDIR/*.bam; do
  if [ -z $BAM_NAMES ]; then
    BAM_NAMES=$SPLIT_BAM
  else
    BAM_NAMES=$BAM_NAMES' '$SPLIT_BAM
done

$SAMBAMBA merge $OUTPUT $BAM_NAMES

echo `date`: Done
