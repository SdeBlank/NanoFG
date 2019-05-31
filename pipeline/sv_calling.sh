#!/bin/bash

usage() {
echo "
Required parameters:
    -b|--bam		     Path to sorted bam file

Optional parameters:
    -h|--help		             Shows help
    -t|--threads	           Number of threads [$THREADS]
    -sv|--sv_caller          'NanoSV' or the path to Sniffles [$SV_CALLER]
    -s|--sambamba	           Path to sambamba [$SAMBAMBA]
    -v|--venv		             Path to virtual env of NanoSV [$VENV]
    -c|--config		           Path to config file [$CONFIG]
    -ss|--sniffles_settings  Settings for sniffles sv calling [$SNIFFLES_SETTINGS]
    -o|--output		           Path to vcf output file [$OUTPUT]
"
}

POSITIONAL=()

# DEFAULTS
NANOFG_DIR=$(realpath $(dirname $(dirname ${BASH_SOURCE[0]})))
FILES_DIR=$NANOFG_DIR/files

THREADS=1
SV_CALLER='/hpc/cog_bioinf/kloosterman/tools/NanoSV/nanosv/NanoSV.py'
SAMBAMBA='/hpc/local/CentOS7/cog_bioinf/sambamba_v0.6.5/sambamba'
OUTPUT='/dev/stdout'
VENV=$NANOFG_DIR/venv/bin/activate
CONFIG=$FILES_DIR/nanosv_last_config.ini
SNIFFLES_SETTINGS='-s 2 -n -1 --genotype'


while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
    -h|--help)
    usage
    shift # past argument
    ;;
    -b|--bam)
    BAM="$2"
    shift # past argument
    shift # past value
    ;;
    -t|--threads)
    THREADS="$2"
    shift # past argument
    shift # past value
    ;;
    -s|--sambamba)
    SAMBAMBA="$2"
    shift # past argument
    shift # past value
    ;;
    -sv|--sv_caller)
    SV_CALLER="$2"
    shift # past argument
    shift # past value
    ;;
    -v|--venv)
    VENV="$2"
    shift # past argument
    shift # past value
    ;;
    -c|--config)
    CONFIG="$2"
    shift # past argument
    shift # past value
    ;;
    -ss|--sniffles_settings)
    SNIFFLES_SETTINGS="$2"
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
if [ -z $BAM ]; then
    echo "Missing -b|--bam parameter"
    usage
    exit
fi

echo `date`: Running on `uname -n`

. $VENV

if [[ $SV_CALLER == *"nanosv"* ]] || [[ $SV_CALLER == *"NanoSV"* ]]; then
  #$SV_CALLER  \
  python $SV_CALLER \
  -s $SAMBAMBA \
  -c $CONFIG \
  -t $THREADS \
  -o $OUTPUT \
  $BAM
fi

if [[ $SV_CALLER == *"sniffles"* ]] || [[ $SV_CALLER == *"Sniffles"* ]]; then
  $SV_CALLER  \
  -v $OUTPUT \
  -m $BAM \
  -t $THREADS \
  $SNIFFLES_SETTINGS
fi

deactivate

echo `date`: Done
