#!/bin/bash

usage() {
echo "
Required parameters:
    -v|--vcf		    Path to vcf file

Optional parameters:
    -h|--help       Shows help
    -s|--script     Path to vcf_primer_filter.py [$SCRIPT]
    -l|--LINES      Number of lines to put in each spit vcf file [$LINES]
    -o|--output     VCF output file [$OUTPUT]
    -sd|--split_directory Directory to place split vcfs in [./]
    -e|--venv       Path to virtual environment[$VENV]
"
}

POSITIONAL=()

SCRIPT="/home/cog/sdeblank/Documents/github/NanoFG/NanoFG.py"
VENV="/data/sharc/venv/bin/activate"
OUTPUT="/dev/stdout"
OUTPUTDIR="./split_vcf"
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
    -sd|--split_directory)
    OUTPUTDIR="$2"
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

. $VENV

if [ ! -d $OUTPUTDIR ]; then
    mkdir $OUTPUTDIR
fi

cat << EOF > $VCF_SPLIT_SH                  # ADD JOBS FOLDER ETC.
#!/bin/bash

#$ -N $VCF_SPLIT_JOBNAME
#$ -cwd
#$ -l h_vmem=$VCF_SPLIT_MEM
#$ -l h_rt=$VCF_SPLIT_TIME
#$ -e $VCF_SPLIT_ERR
#$ -o $VCF_SPLIT_LOG

HEADER=$(grep "^#" $VCF)
AWK="grep -v \"^#\" $VCF | awk -v HEADER=\"\$HEADER\" 'NR%$LINES==1 { file = \"$OUTPUTDIR/\" int(NR/$LINES)+1 \".vcf\"; print HEADER > file } { print > file }'"
eval $AWK
EOF

NUMBER_OF_LINES_VCF_1=$(grep -v "^#" $VCF | wc -l | grep -oP "(^\d+)")
NUMBER_OF_LINES_VCF_2=$(cat $OUTPUTDIR/*.vcf | grep -v "^#" | wc -l | grep -oP "(^\d+)")

if [ $NUMBER_OF_LINES_VCF_1 == $NUMBER_OF_LINES_VCF_2 ]; then
  for SPLIT_VCF in $OUTPUTDIR/*.vcf; do
    SPLIT_OUTPUT=${SPLIT_VCF/.vcf/_FusionGenes.txt}
    SPLIT_VCF_OUTPUT=${SPLIT_VCF/.vcf/_FusionGenes.vcf}
    qsub
    python $SCRIPT -v $SPLIT_VCF -fo $SPLIT_OUTPUT -o $SPLIT_VCF_OUTPUT
  done
fi

NUMBER_SPLIT_VCF=$(ls -l $OUTPUTDIR/* | grep -cv "FusionGenes" | grep -oP "(^\d+)")
NUMBER_SPLIT_OUTPUT=$(ls -l $OUTPUTDIR/* | grep -c "FusionGenes.txt" | grep -oP "(^\d+)")

if [ $NUMBER_SPLIT_VCF == $NUMBER_SPLIT_OUTPUT ]; then
  head -n 1 $OUTPUTDIR/1_FusionGenes.txt > $OUTPUT
  for fusion_output in $OUTPUTDIR/*FusionGenes.txt; do
    tail -n +2 $fusion_output >> $OUTPUT
  done
fi

echo `date`: Done
