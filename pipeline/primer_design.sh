#!/bin/bash

usage() {
echo "
Required parameters:
    -f|--fasta		Path to fasta file
    -o|--output

Optional parameters:
    -h|--help     Shows help
    -psr|--psr  PSR [$PSR]
    -pdpc|--primer3_core   Primer3 core [$PRIMER3_CORE]
    -mtf|--minimal_target_flank Minimal flank around breakpoint that has to be included in PCR product [$MIN_TARGET_FLANK]
"
}

POSITIONAL=()

NANOFG_DIR=$(realpath $(dirname $(dirname ${BASH_SOURCE[0]})))
source ${NANOFG_DIR}/paths.ini

PRIMER_DESIGN_DIR=${PATH_PRIMER_DESIGN_DIR}

# DEFAULTS
PSR='60-200'
MIN_TARGET_FLANK=10
PRIMER3_CORE=$PRIMER_DESIGN_DIR/src/primer3_core

while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
    -h|--help)
    usage
    shift # past argument
    ;;
    -f|--fasta)
    FASTA="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output)
    OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
    -psr|--psr)
    PSR="$2"
    shift # past argument
    shift # past value
    ;;
    -mtf|--minimal_target_flank)
    MIN_TARGET_FLANK="$2"
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
if [ -z $FASTA ]; then
    echo "Missing -f|--fasta parameter"
    usage
    exit
fi

#echo `date`: Running on `uname -n`

FASTA_ID=$(grep ">" $FASTA)
FASTA_ID=${FASTA_ID//\>}
FASTA_SEQUENCE=$(grep -v ">" $FASTA)

bp_index() {
  x="${1%%[]*}"
  [[ "$x" = "$1" ]] && echo -1 || echo "${#x}"
}

SEQUENCE_TARGET=$(($(bp_index $FASTA_SEQUENCE)-$MIN_TARGET_FLANK+1))
SEQUENCE_TARGET_LENGTH=$(($MIN_TARGET_FLANK*2))
FASTA_SEQUENCE=${FASTA_SEQUENCE//[]}

PRIMER_OUTPUT=${FASTA//.fasta}.primers

$PRIMER3_CORE <<EOF > $PRIMER_OUTPUT
SEQUENCE_ID=$FASTA_ID
SEQUENCE_TEMPLATE=$FASTA_SEQUENCE
PRIMER_PRODUCT_SIZE_RANGE=$PSR
SEQUENCE_TARGET=$SEQUENCE_TARGET,$SEQUENCE_TARGET_LENGTH
PRIMER_NUM_RETURN=1
=
EOF

FORWARD_PRIMER_ID="$FASTA_ID-$(sed -n -e 's/^PRIMER_LEFT_0=//p' $PRIMER_OUTPUT)"
FORWARD_PRIMER_ID=${FORWARD_PRIMER_ID//,*}F
FORWARD_PRIMER_SEQ=$(sed -n -e 's/^PRIMER_LEFT_0_SEQUENCE=//p' $PRIMER_OUTPUT)
REVERSE_PRIMER_ID="$FASTA_ID-$(sed -n -e 's/^PRIMER_RIGHT_0=//p' $PRIMER_OUTPUT)"
REVERSE_PRIMER_ID=${REVERSE_PRIMER_ID//,*}R
REVERSE_PRIMER_SEQ=$(sed -n -e 's/^PRIMER_RIGHT_0_SEQUENCE=//p' $PRIMER_OUTPUT)
PRODUCT_SIZE=$(sed -n -e 's/^PRIMER_PAIR_0_PRODUCT_SIZE=//p' $PRIMER_OUTPUT)
echo -e "$FASTA_ID\t$FORWARD_PRIMER_ID\t$FORWARD_PRIMER_SEQ\t$REVERSE_PRIMER_ID\t$REVERSE_PRIMER_SEQ\t$PRODUCT_SIZE" >> $OUTPUT

#echo `date`: Done
