#!/bin/bash

usage() {
echo "
Required parameters:
    -sd|--split_directory
    -jd|--job_directory
    -ld|--directory
Optional parameters:
    -h|--help       Shows help
    -s|--script     Path to vcf_primer_filter.py [$SCRIPT]
    -o|--output     Fusion gene output file
    -vo|--vcf_output  VCF output file
"
}

POSITIONAL=()

SCRIPT="/home/cog/sdeblank/Documents/github/NanoFG/NanoFG.py"
VENV="/data/sharc/venv/bin/activate"

NANOFG_SPLIT_MEM=5G
NANOFG_SPLIT_TIME=0:10:0

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
    -sd|--split_directory)
    SPLITDIR="$2"
    shift # past argument
    shift # past value
    ;;
    -jd|--job_directory)
    JOBDIR="$2"
    shift # past argument
    shift # past value
    ;;
    -ld|--log_directory)
    JOBDIR="$2"
    shift # past argument
    shift # past value
    ;;
    -s|--script)
    SCRIPT="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output)
    SPLIT_OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
    -vo|--vcf_output)
    SPLIT_VCF_OUTPUT="$2"
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

if [ -z $SPLITDIR ]; then
    echo "Missing -sd|--split_directory parameter"
    usage
    exit
fi

if [ -z $JOBDIR ]; then
    echo "Missing -jd|--job_directory parameter"
    usage
    exit
fi

if [ -z $JOBDIR ]; then
    echo "Missing -ld|--log_directory parameter"
    usage
    exit
fi

echo `date`: Running on `uname -n`


for SPLIT_VCF in $SPLITDIR/*.vcf; do
  SPLIT_OUTPUT=${SPLIT_VCF/.vcf/_FusionGenes.txt}
  SPLIT_VCF_OUTPUT=${SPLIT_VCF/.vcf/_FusionGenes.vcf}

  NANOFG_SPLIT_JOBNAME=${SPLIT_VCF/.vcf/_NanoFG}
  NANOFG_SPLIT_JOB=$JOBDIR/$NANOFG_SPLIT_JOBNAME.sh
  NANOFG_SPLIT_ERR=$LOGDIR/$NANOFG_SPLIT_JOBNAME.log
  NANOFG_SPLIT_LOG=$LOGDIR/$NANOFG_SPLIT_JOBNAME.err

  NANOFG_SPLIT_JOBNAMES=$NANOFG_SPLIT_JOBNAMES','$NANOFG_SPLIT_JOBNAME

  cat << EOF > $NANOFG_SPLIT_JOB
  #!/bin/bash

  #$ -N $NANOFG_SPLIT_JOBNAME
  #$ -cwd
  #$ -l h_vmem=$NANOFG_SPLIT_MEM
  #$ -l h_rt=$NANOFG_SPLIT_TIME
  #$ -e $NANOFG_SPLIT_ERR
  #$ -o $NANOFG_SPLIT_LOG
  #$ -hold_jid $NANOFG_SPLIT_JOBNAMES

  . $VENV
  python $SCRIPT -v \$SPLIT_VCF -fo \$SPLIT_OUTPUT -o \$SPLIT_VCF_OUTPUT
EOF
  qsub $NANOFG_SPLIT_JOB
done

MERGE_OUTPUT_JOB << EOF
#!/bin/bash

#$ -N $MERGE_OUTPUT_JOBNAME
#$ -cwd
#$ -l h_vmem=$MERGE_OUTPUT_MEM
#$ -l h_rt=$MERGE_OUTPUT_TIME
#$ -e $MERGE_OUTPUT_ERR
#$ -o $MERGE_OUTPUT_LOG
#$ -hold_jid $NANOFG_SPLIT_JOBNAMES

NUMBER_SPLIT_VCF=$(ls -l $SPLITDIR/* | grep -cv "FusionGenes" | grep -oP "(^\d+)")
NUMBER_SPLIT_OUTPUT=$(ls -l $SPLITDIR/* | grep -c "FusionGenes.txt" | grep -oP "(^\d+)")

if [ $NUMBER_SPLIT_VCF == $NUMBER_SPLIT_OUTPUT ]; then
  head -n 1 $SPLITDIR/1_FusionGenes.txt > $OUTPUT
  for fusion_output in $SPLITDIR/*FusionGenes.txt; do
    tail -n +2 $fusion_output >> $OUTPUT
  done
fi
EOF

qsub $MERGE_OUTPUT_JOB

echo `date`: Done
