#!/bin/bash

usage() {
echo "
Required parameters:
    -b|--bam		Path to sorted bam file

Optional parameters:
    -h|--help		Shows help
    -t|--threads	Number of threads [$THREADS]
    -s|--sambamba	Path to sambamba [$SAMBAMBA]
    -v|--venv		Path to virtual env of NanoSV [$VENV]
    -c|--config		Path to config file [$CONFIG]
    -o|--output		Path to vcf output file [$OUTPUT]
"
}

POSITIONAL=()

# DEFAULTS
NANOFG_DIR=$(realpath $(dirname $(dirname ${BASH_SOURCE[0]})))
FILES_DIR=$NANOFG_DIR/files

THREADS=1
SAMBAMBA='/hpc/local/CentOS7/cog_bioinf/sambamba_v0.6.5/sambamba'
OUTPUT='/dev/stdout'
VENV=$NANOFG_DIR/venv/bin/activate
CONFIG=$FILES_DIR/nanosv_last_config.ini


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

NanoSV  \
-s $SAMBAMBA \
-c $CONFIG \
-t $THREADS \
-o $OUTPUT \
$BAM

deactivate

echo `date`: Done
