#!/bin/bash

usage() {
echo "
Required parameters:
    -d|--split_directory

Optional parameters:
    -h|--help       Shows help
    -s|--script     Path to vcf_primer_filter.py [$SCRIPT]
    -e|--venv       Path to virtual environment[$VENV]
    -t|--threads
    -hv|--h_vmem
    -hr|--h_runtime
"
}

POSITIONAL=()

#DEFAULTS
SPLIT_FUSION_READ_EXTRACTION_THREADS=8
SPLIT_FUSION_READ_EXTRACTION_MEMORY=5G
SPLIT_FUSION_READ_EXTRACTION_TIME=0:30:0
NANOFG_DIR=$(realpath $(dirname $(dirname {BASH_SOURCE[0]})))
SCRIPT=$NANOFG_DIR/scripts/fusion_gene_read_extraction.py
VENV=NANOFG_DIR/venv/bin/activate

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
    -t|--threads)
    SPLIT_FUSION_READ_EXTRACTION_THREADS="$2"
    shift # past argument
    shift # past value
    ;;
    -hv|--h_vmem)
    SPLIT_FUSION_READ_EXTRACTION_MEMORY="$2"
    shift # past argument
    shift # past value
    ;;
    -hr|--h_rtvmem)
    SPLIT_FUSION_READ_EXTRACTION_TIME="$2"
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

if [ -z $SPLITDIR ]; then
    echo "Missing -d|--split_directory"
    usage
    exit
fi

for SPLIT_VCF in $SPLITDIR/*.vcf; do
  SPLIT_OUTPUT=${SPLIT_VCF/.vcf/_FusionGenes.txt}
  SPLIT_VCF_OUTPUT=${SPLIT_VCF/.vcf/_FusionGenes.vcf}

  SPLIT_FUSION_READ_EXTRACTION_JOBNAME=$(basename $SPLIT_VCF)
  SPLIT_FUSION_READ_EXTRACTION_JOBNAME=$SPLIT_FUSION_READ_EXTRACTION_JOBNAME_FUSION_READ_EXTRACTION
  SPLIT_FUSION_READ_EXTRACTION_JOB=$JOBDIR/$SPLIT_FUSION_READ_EXTRACTION_JOBNAME.sh
  SPLIT_FUSION_READ_EXTRACTION_ERR=$LOGDIR/$SPLIT_FUSION_READ_EXTRACTION_JOBNAME.err
  SPLIT_FUSION_READ_EXTRACTION_LOG=$LOGDIR/$SPLIT_FUSION_READ_EXTRACTION_JOBNAME.log

  if [ -z $NANOFG_SPLIT_JOBNAMES ]; then
    NANOFG_SPLIT_JOBNAMES=$NANOFG_SPLIT_JOBNAME
  else
    NANOFG_SPLIT_JOBNAMES=$NANOFG_SPLIT_JOBNAMES','$NANOFG_SPLIT_JOBNAME
  fi

  qsub << EOF
#!/bin/bash

#$ -N $SPLIT_FUSION_READ_EXTRACTION_JOBNAME
#$ -cwd
#$ -pe threaded $SPLIT_FUSION_READ_EXTRACTION_THREADS
#$ -l h_vmem=$SPLIT_FUSION_READ_EXTRACTION_MEMORY
#$ -l h_rt=$SPLIT_FUSION_READ_EXTRACTION_TIME
#$ -e $SPLIT_FUSION_READ_EXTRACTION_ERR
#$ -o $SPLIT_FUSION_READ_EXTRACTION_LOG

echo \`date\`: Running on \`uname -n\`

if [ ! -e ${LOGDIR}/${SPLIT_FUSION_READ_EXTRACTION_JOBNAME}.done ]
  . $VENV
  python $SCRIPT -v $SPLIT_VCF -fo $SPLIT_OUTPUT -o $SPLIT_VCF_OUTPUT

  FINISHED="\$(tail -n 1 \$LOGFILE | grep -o End)"
    if [ -z \$FINISHED ]; then
      echo "$SPLIT_FUSION_READ_EXTRACTION_JOB did not fully complete; Increase NANOFG_SPLIT_MEM or NANOFG_SPLIT_TIME" >&2
      exit
    else
      touch ${LOGDIR}/${SPLIT_FUSION_READ_EXTRACTION_JOBNAME}.done
    fi
else
  echo "$SPLIT_FUSION_READ_EXTRACTION_JOBNAME has been completed already"
fi

echo \`date\`: Done

EOF
done
