#!/bin/bash

usage() {
echo "
Required parameters:
    -v|--vcf		    Path to vcf file

Optional parameters:
    -h|--help       Shows help
    -d|--NanoFG_directory Directory that contains NanoFG [$NANOFG_DIR]
    -s|--script     Path to vcf_primer_filter.py [$SCRIPT]
    -l|--LINES      Number of lines to put in each spit vcf file [$LINES]
    -o|--output     VCF output file [$OUTPUT]
    -e|--venv       Path to virtual environment[$VENV]
"
}

POSITIONAL=()

SCRIPT="/home/cog/sdeblank/Documents/github/NanoFG/NanoFG.py"
VENV="/data/sharc/venv/bin/activate"
NANOFG_DIR=$(dirname $SCRIPT)
LINES=100

SPLIT_VCF_MEM=5G
SPLIT_VCF_TIME=0:5:0

NANOFG_MEM=5G
NANOFG_TIME=0:5:0

NANOFG_SPLIT_MEM=5G
NANOFG_SPLIT_TIME=0:15:0

MERGE_OUTPUT_MEM=5G
MERGE_OUTPUT_TIME=0:5:0

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

echo `date`: Running on `uname -n`

if [ -z $OUTPUT ]; then
    OUTPUT=./$(basename $VCF)
    OUTPUT=${OUTPUT/.vcf/_FusionGenes.txt}
fi

SPLITDIR="$(dirname $OUTPUT)/split_vcf"
JOBDIR="$(dirname $OUTPUT)/jobs"
LOGDIR="$(dirname $OUTPUT)/logs"

if [ ! -d $SPLITDIR ]; then
    mkdir $SPLITDIR
fi

if [ ! -d $JOBDIR ]; then
    mkdir $JOBDIR
fi

if [ ! -d $LOGDIR ]; then
    mkdir $LOGDIR
fi

VCF_NAME=$(basename $VCF)
# SPLIT_VCF_JOBNAME=${VCF_NAME/.vcf/_split_vcf}
# SPLIT_VCF_JOB=$JOBDIR/$SPLIT_VCF_JOBNAME.sh
# SPLIT_VCF_ERR=$LOGDIR/$SPLIT_VCF_JOBNAME.log
# SPLIT_VCF_LOG=$LOGDIR/$SPLIT_VCF_JOBNAME.err
#
# NANOFG_JOBNAME=${VCF_NAME/.vcf/_NanoFG}
# NANOFG_SPLIT_JOBNAMES=$NANOFG_JOBNAME
# NANOFG_JOB=$JOBDIR/$NANOFG_JOBNAME.sh
# NANOFG_ERR=$LOGDIR/$NANOFG_JOBNAME.log
# NANOFG_LOG=$LOGDIR/$NANOFG_JOBNAME.err

MERGE_OUTPUT_JOBNAME=${VCF_NAME/.vcf/_merge_output}
MERGE_OUTPUT_JOB=$JOBDIR/$MERGE_OUTPUT_JOBNAME.sh
MERGE_OUTPUT_ERR=$LOGDIR/$MERGE_OUTPUT_JOBNAME.err
MERGE_OUTPUT_LOG=$LOGDIR/$MERGE_OUTPUT_JOBNAME.log

HEADER=$(grep "^#" $VCF)
AWK="grep -v \"^#\" $VCF | awk -v HEADER=\"\$HEADER\" 'NR%$LINES==1 { file = \"$SPLITDIR/\" int(NR/$LINES)+1 \".vcf\"; print HEADER > file } { print > file }'"
eval $AWK

NUMBER_OF_LINES_VCF_1=$(grep -v "^#" $VCF | wc -l | grep -oP "(^\d+)")
NUMBER_OF_LINES_VCF_2=$(cat $SPLITDIR/*.vcf | grep -v "^#" | wc -l | grep -oP "(^\d+)")

if [ $NUMBER_OF_LINES_VCF_1 != $NUMBER_OF_LINES_VCF_2 ]; then
  echo "Number of lines in all split VCF files is not equal to the number of lines in the original VCF file"
  exit
fi

for SPLIT_VCF in $SPLITDIR/*.vcf; do
  SPLIT_OUTPUT=${SPLIT_VCF/.vcf/_FusionGenes.txt}
  SPLIT_VCF_OUTPUT=${SPLIT_VCF/.vcf/_FusionGenes.vcf}

  NANOFG_SPLIT_JOBNAME=$(basename $SPLIT_VCF)
  NANOFG_SPLIT_JOBNAME=NanoFG_$NANOFG_SPLIT_JOBNAME
  NANOFG_SPLIT_JOB=$JOBDIR/$NANOFG_SPLIT_JOBNAME.sh
  NANOFG_SPLIT_ERR=$LOGDIR/$NANOFG_SPLIT_JOBNAME.err
  NANOFG_SPLIT_LOG=$LOGDIR/$NANOFG_SPLIT_JOBNAME.log

  NANOFG_SPLIT_JOBNAMES=$NANOFG_SPLIT_JOBNAMES','$NANOFG_SPLIT_JOBNAME

  qsub << EOF
#!/bin/bash

#$ -N $NANOFG_SPLIT_JOBNAME
#$ -cwd
#$ -l h_vmem=$NANOFG_SPLIT_MEM
#$ -l h_rt=$NANOFG_SPLIT_TIME
#$ -e $NANOFG_SPLIT_ERR
#$ -o $NANOFG_SPLIT_LOG

. $VENV
python $SCRIPT -v $SPLIT_VCF -fo $SPLIT_OUTPUT -o $SPLIT_VCF_OUTPUT
EOF
done

qsub << EOF
#!/bin/bash

#$ -N $MERGE_OUTPUT_JOBNAME
#$ -cwd
#$ -l h_vmem=$MERGE_OUTPUT_MEM
#$ -l h_rt=$MERGE_OUTPUT_TIME
#$ -e $MERGE_OUTPUT_ERR
#$ -o $MERGE_OUTPUT_LOG
#$ -hold_jid $NANOFG_SPLIT_JOBNAMES

NUMBER_SPLIT_VCF=\$(ls -l $SPLITDIR/* | grep -cv "FusionGenes" | grep -oP "(^\d+)")
NUMBER_SPLIT_OUTPUT=\$(ls -l $SPLITDIR/* | grep -c "FusionGenes.txt" | grep -oP "(^\d+)")

for LOGFILE in $LOGDIR/*.log; do
  FINISHED=\$(tail -n 1 \$LOGFILE | grep End)
  if [ -z \$FINISHED ]; then
    echo "One or more of the NanoFG did not complete; Increase NANOFG_SPLIT_MEM or NANOFG_SPLIT_TIME"
    exit
  fi
done

if [ $NUMBER_SPLIT_VCF == $NUMBER_SPLIT_OUTPUT ]; then
  head -n 1 $SPLITDIR/1_FusionGenes.txt > $OUTPUT
  for fusion_output in $SPLITDIR/*FusionGenes.txt; do
    tail -n +2 \$fusion_output >> $OUTPUT
  done
else
  echo "Number of _FusionGenes.text files are not equal to the number of split vcfs"
  exit
fi
EOF

echo `date`: Done
