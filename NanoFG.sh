#!/bin/bash

usage() {
echo "
Required parameters:
    -v|--vcf		    Path to vcf file

Optional parameters:
    -h|--help       Shows help
    -s|--script     Path to vcf_primer_filter.py [$SCRIPT]
    -o|--output     VCF output file [$OUTPUT]
    -e|--venv       Path to virtual environment[$VENV]
"
}

POSITIONAL=()

SCRIPT="/home/cog/sdeblank/Documents/github/All_scripts/NanoFG/NanoFG.py"
VENV="/data/sharc/venv/bin/activate"
OUTPUT="/dev/stdout"

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
    SCRIPT="$2"
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

python $SCRIPT -v $VCF


# LINES=20
# OUTPUTDIR="/home/cog/sdeblank/Documents/github/All_scripts/NanoFG/split_vcf"
# HEADER=$(grep "^#" $VCF)
# AWK="grep -v \"^#\" $VCF | awk -v HEADER=\"\$HEADER\" 'NR%$LINES==1 { file = \"$OUTPUTDIR/\" \"${VCF/.vcf/}.split\"int(NR/$LINES)+1 \".vcf\"; print HEADER > file } { print > file }'"
# eval $AWK

# for SPLIT_VCF in $OUTPUTDIR/*.vcf; do
#   BEDNAME=($(basename $BED | tr '.' ' '))
#   BEDNAME=${BEDNAME[0]}
#   BED_JOBNAME=${JOBNAMES[i]}
#   BED_SH=$JOBDIR/$BED_JOBNAME.sh
#   BED_ERR="$LOGDIR/${BED_JOBNAME}_\$TASK_ID.err"
#   BED_LOG="$LOGDIR/${BED_JOBNAME}_\$TASK_ID.log"
#   BED_IN="$VCF_SPLIT_OUTPUTDIR/\$ID.vcf"
#   BED_OUT="$VCF_SPLIT_OUTPUTDIR/\$ID.$BEDNAME.vcf"
#
# cat << EOF > $BED_SH
# #!/bin/bash
#
# #$ -N $BED_JOBNAME
# #$ -cwd
# #$ -t 1-$NUMBER_OF_SPLIT_FILES:1
# #$ -l h_vmem=$BED_MEM
# #$ -l h_rt=$BED_TIME
# #$ -e $BED_ERR
# #$ -o $BED_LOG
#
#
# ID=\$SGE_TASK_ID
#
# echo \`date\`: Running on \`uname -n\`
#
# if [ -e $BED_IN ]; then
#     if [ ! -e $BED_OUT.done ]; then
# 	     bash $STEPSDIR/bed_annotation.sh \\
# 	      -v $BED_IN \\
#         -b $BED \\
#         -s $BED_SCRIPT \\
#         -e $VENV \\
#         -o $BED_OUT
#     else
# 	     echo $BED_OUT already exists
#     fi
# else
#     echo $BED_IN does not exists
# fi
#
# echo \`date\`: Done
# EOF
#   qsub $BED_SH
#   ((i=i+1))
# done


#
# ### SPLIT VCF IN 10 SMALLER VCFS AND RUN PARALLEL
# for SPLIT_VCF in $OUTPUTDIR/*.feature.bed; do


echo `date`: Done
