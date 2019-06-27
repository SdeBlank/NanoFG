#!/bin/bash

usage() {
echo "
Required parameters:
    -f|--fasta		Path to fasta file
    -o|--output

Optional parameters:
    -h|--help     Shows help
    -pdb|--bindir		Bindir [$BINDIR]
    -pdpt|--pcr_type   PCR type [$PCR_TYPE]
    -pdtp|--tilling_params Tilling parameters [$TILLING_PARAMS]
    -psr|--psr  PSR [$PSR]
    -pdgp|--guix_profile   GUIX profile [$GUIX_PROFILE]
    -pdpc|--primer3_core   Primer3 core [$PRIMER3_CORE]
    -pdm|--mispriming     Mispriming [$MISPRIMING]
    -pde|--primer_design_error
"
}

POSITIONAL=()

PRIMER_DESIGN_DIR="/hpc/cog_bioinf/cuppen/personal_data/jvalleinclan/tools_kloosterman/primer3/"
# DEFAULTS
BINDIR=$PRIMER_DESIGN_DIR/primers
PCR_TYPE='single'
TILLING_PARAMS=' '
PSR='60-200'
GUIX_PROFILE=$PRIMER_DESIGN_DIR/emboss/.guix-profile
PRIMER3_CORE=$PRIMER_DESIGN_DIR/primer3/src/primer3_core
MISPRIMING=$PRIMER_DESIGN_DIR/repbase/current/empty.ref
#MISPRIMING='/hpc/cog_bioinf/cuppen/personal_data/jvalleinclan/tools_kloosterman/primer3repbase/current/humrep.ref'

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
    -pdb|--bindir)
    BINDIR="$2"
    shift # past argument
    shift # past value
    ;;
    -pdpt|--pcr_type)
    PCR_TYPE="$2"
    shift # past argument
    shift # past value
    ;;
    -pdtp|--tilling_params)
    TILLING_PARAMS="$2"
    shift # past argument
    shift # past value
    ;;
    -psr|--psr)
    PSR="$2"
    shift # past argument
    shift # past value
    ;;
    -pdgp|--guix_profile)
    GUIX_PROFILE="$2"
    shift # past argument
    shift # past value
    ;;
    -pdpc|--primer3_core)
    PRIMER3_CORE="$2"
    shift # past argument
    shift # past value
    ;;
    -pdm|--mispriming)
    MISPRIMING="$2"
    shift # past argument
    shift # past value
    ;;
    -pde|--primer_design_error)
    PRIMER_DESIGN_ERR="$2"
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

guixr load-profile $GUIX_PROFILE --<<EOF
export EMBOSS_PRIMER3_CORE=$PRIMER3_CORE
$BINDIR/primerBATCH1 $MISPRIMING $PCR_TYPE $PSR $TILLING_PARAMS <$FASTA
EOF

grep -v "FAILED" ./primers.txt > ${OUTPUT}.tmp
paste <(cat ${OUTPUT}.tmp) <(grep "PRODUCT SIZE" ./primer3.out | grep -oP "\d+$") > $OUTPUT
rm ${OUTPUT}.tmp

#echo `date`: Done
