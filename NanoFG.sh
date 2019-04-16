#!/bin/bash

usage() {
echo "
Required parameters:
    -v|--vcf		                                                                  Path to vcf file
    -b|--bam                                                                      Path to bam file

Optional parameters:

GENERAL
    -h|--help                                                                     Shows help
    -o|--outputdir                                                                Path to output directory
    -d|--nanofg_dir                                                               Directory that contains NanoFG [${NANOFG_DIR}]
    -e|--venv                                                                     Path to virtual environment[${VENV}]
    -m|--mail                                                                     Email adress
    -df|--dont_filter                                                             Don't filter out all non-PASS SVs
    -dc|--dont_clean                                                              Don't clean up the intermediate files

REQUIRED TOOLS
    -n|--nanosv                                                                   Path to NanoSV [${NANOSV}]
    -s|--sambamba                                                                 Path to sambamba|samtools [${SAMBAMBA}]
    -l|--last_dir                                                                 Path to LAST directory [${LAST_DIR}]
    -w|--wtdbg2_dir                                                               Path to wtdbg2 directory [${WTDBG2_DIR}]

FUSION READ EXTRACTION
    -fres|--fusion_read_extraction_script                                         Path to the fusion_read_extraction.py script [${FUSION_READ_EXTRACTION_SCRIPT}]
    -fret|--fusion_read_extraction_threads                                        Number of threads [${FUSION_READ_EXTRACTION_THREADS}]
    -frehv|--fusion_read_extraction_h_vmem                                        Fusion read extraction memory [${FUSION_READ_EXTRACTION_MEMORY}]
    -freht|--fusion_read_extraction_h_rt                                          Fusion read extraction time [${FUSION_READ_EXTRACTION_TIME}]

CONSENSUS MAPPING

    -ccws|--consensus_calling_wtdbg2_settings                                     wtdbg2 settings [${WTDBG2_SETTINGS}]
    -cct|--consensus_calling_threads                                              Number of threads [${CONSENSUS_CALLING_THREADS}]
    -cchv|--consensus_calling_h_vmem                                              Consensus mapping memory [${CONSENSUS_CALLING_MEMORY}]
    -cchr|--consensus_calling_h_rt                                                Consensus mapping time [${CONSENSUS_CALLING_TIME}]

LAST_MAPPING
    -lmr|--consensus_calling_refgenome                                            Reference genome [${REF}]
    -lmrd|--consensus_calling_refdict                                             Reference genome .dict file [${REF_DICT}]
    -lms|--consensus_calling_last_settings                                        LAST settings [${LAST_SETTINGS}]
    -lmt|--last_mapping_threads                                                   Number of threads[$LAST_MAPPING_THREADS]
    -lmhv|--last_mapping_h_vmem                                                   Last mapping memory[$LAST_MAPPING_MEMORY]
    -lmht|--last_mapping_h_rt                                                     Last mapping runtime[$LAST_MAPPING_TIME]

BAM MERGE
    -bmt|vcf_split_threads                                                        Number of threads [${VCF_SPLIT_THREADS}]
    -bmhv|vcf_split_h_vmem                                                        Bam merge memory [${VCF_SPLIT_MEMORY}]
    -bmhr|vcf_split_h_rt                                                          Bam merge time [${VCF_SPLIT_TIME}]

SV CALLING
    -sct|sv_calling_threads                                                       Number of threads [${SV_CALLING_THREADS}]
    -schv|sv_calling_h_vmem                                                       SV calling memory [${SV_CALLING_MEMORY}]
    -schr|sv_calling_h_rt                                                         SV calling time [${SV_CALLING_TIME}]

FUSION CHECK
    -fcio|--fusion_check_info_output                                              Path to the NanoFG output info file []
    -fcvo|--fusion_check_vcf_output                                               Path to the NanoFG output vcf file []
    -fcs|--fusion_check_script                                                    Path to vcf_primer_filter.py [$FUSION_CHECK_SCRIPT]
    -fct|--fusion_check_threads                                                   VCF output file [$FUSION_CHECK_THREADS]
    -fchv|--fusion_check_h_vmem                                                   VCF output file [$FUSION_CHECK_MEMORY]
    -fcht|--fusion_check_h_rt                                                     VCF output file [$FUSION_CHECK_TIME]
"
}

POSITIONAL=()

#GENERAL DEFAULTS
NANOFG_DIR=$(realpath $(dirname ${BASH_SOURCE[0]}))
source $NANOFG_DIR/paths.ini
PIPELINE_DIR=$NANOFG_DIR/pipeline
FILES_DIR=$NANOFG_DIR/files
SCRIPT_DIR=$NANOFG_DIR/scripts
VENV=${NANOFG_DIR}/venv/bin/activate

OUTPUTDIR=$(realpath ./)
echo $OUTPUTDIR
DONT_CLEAN=false

#TOOL PATH DEFAULTS
SAMBAMBA=$PATH_SAMBAMBA
LAST_DIR=$PATH_LAST_DIR
WTDBG2_DIR=$PATH_WTDBG2_DIR

#FUSION READ EXTRACTION DEFAULTS
FUSION_READ_EXTRACTION_SCRIPT=$SCRIPT_DIR/FusionReadExtraction.py
FUSION_READ_EXTRACTION_THREADS=1
FUSION_READ_EXTRACTION_TIME=0:10:0
FUSION_READ_EXTRACTION_MEMORY=10G

#CONSENSUS CALLING DEFAULTS
CONSENSUS_CALLING_WTDBG2_SETTINGS='-x ont -g 3g'
CONSENSUS_CALLING_THREADS=8
CONSENSUS_CALLING_TIME=0:15:0
CONSENSUS_CALLING_MEMORY=20G

#LAST MAPPING DEFAULTS
LAST_MAPPING_REFGENOME=$PATH_HOMO_SAPIENS_REFGENOME
LAST_MAPPING_REFDICT=$PATH_HOMO_SAPIENS_REFDICT
LAST_MAPPING_SETTINGS="-Q 0 -p ${LAST_DIR}/last_params"

#BAM MERGE DEFAULTS
BAM_MERGE_THREADS=1
BAM_MERGE_TIME=0:10:0
BAM_MERGE_MEMORY=40G

#SV CALLING DEFAULTS
SV_CALLING_THREADS=1
SV_CALLING_TIME=0:10:0
SV_CALLING_MEMORY=10G
SV_CALLING_CONFIG=$FILES_DIR/nanosv_last_config.ini

#FUSION CHECK DEFAULTS
FUSION_CHECK_SCRIPT=$SCRIPT_DIR/FusionCheck.py
FUSION_CHECK_THREADS=1
FUSION_CHECK_TIME=0:30:0
FUSION_CHECK_MEMORY=10G

#CHECK NANOFG DEFAULTS
CHECK_NANOFG_THREADS=1
CHECK_NANOFG_TIME=0:5:0
CHECK_NANOFG_MEMORY=5G


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
    -b|--bam)
    BAM="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--outputdir)
    OUTPUTDIR="$2"
    shift # past argument
    shift # past value
    ;;
    -d|--nanofg_dir)
    NANOFG_DIR="$2"
    shift # past argument
    shift # past value
    ;;
    -e|--venv)
    VENV="$2"
    shift # past argument
    shift # past value
    ;;
    -m|--mail)
    MAIL="$2"
    shift # past argument
    shift # past value
    ;;
    -df|--dont_filter)
    DONT_FILTER=true
    shift # past argument
    ;;
    -dc|--dont_clean)
    DONT_CLEAN=true
    shift # past argument
    ;;
    -s|--sambamba)
    SAMBAMBA="$2"
    shift # past argument
    shift # past value
    ;;
    -l|--last_dir)
    LAST_DIR="$2"
    shift # past argument
    shift # past value
    ;;
    -w|--wtdbg2_dir)
    WTDBG2_DIR="$2"
    shift # past argument
    shift # past value
    ;;
    -fres|--fusion_read_extraction_script)
    FUSION_READ_EXTRACTION_SCRIPT="$2"
    shift # past argument
    shift # past value
    ;;
    -fret|--fusion_read_extraction_threads)
    FUSION_READ_EXTRACTION_THREADS="$2"
    shift # past argument
    shift # past value
    ;;
    -frehv|--fusion_read_extraction_h_vmem)
    FUSION_READ_EXTRACTION_MEMORY="$2"
    shift # past argument
    shift # past value
    ;;
    -freht|--fusion_read_extraction_h_rt)
    FUSION_READ_EXTRACTION_TIME="$2"
    shift # past argument
    shift # past value
    ;;
    -ccws|--consensus_calling_wtdbg2_settings)
    CONSENSUS_CALLING_WTDBG2_SETTINGS="$2"
    shift # past argument
    shift # past value
    ;;
    -cct|--consensus_calling_threads)
    CONSENSUS_CALLING_THREADS="$2"
    shift # past argument
    shift # past value
    ;;
    -cchv|--consensus_calling_h_vmem)
    CONSENSUS_CALLING_MEMORY="$2"
    shift # past argument
    shift # past value
    ;;
    -cchr|--consensus_calling_h_rt)
    CONSENSUS_CALLING_TIME="$2"
    shift # past argument
    shift # past value
    ;;
    -lmt|--last_mapping_threads)
    LAST_MAPPING_THREADS="$2"
    shift # past argument
    shift # past value
    ;;
    -lmr|-LAST_MAPPING_REFGENOME)
    LAST_MAPPING_REFGENOME="$2"
    shift # past argument
    shift # past value
    ;;
    -lmrd|--last_mapping_refdict)
    LAST_MAPPING_REFDICT="$2"
    shift # past argument
    shift # past value
    ;;
    -lms|--last_mapping_settings)
    LAST_MAPPING_SETTINGS="$2"
    shift # past argument
    shift # past value
    ;;
    -bmt|vcf_split_threads)
    BAM_MERGE_THREADS="$2"
    shift # past argument
    shift # past value
    ;;
    -bmhv|vcf_split_h_vmem)
    BAM_MERGE_MEMORY="$2"
    shift # past argument
    shift # past value
    ;;
    -bmhr|vcf_split_h_rt)
    BAM_MERGE_TIME="$2"
    shift # past argument
    shift # past value
    ;;
    -sct|sv_calling_threads)
    SV_CALLING_THREADS="$2"
    shift # past argument
    shift # past value
    ;;
    -schv|sv_calling_h_vmem)
    SV_CALLING_MEMORY="$2"
    shift # past argument
    shift # past value
    ;;
    -schr|sv_calling_h_rt)
    SV_CALLING_TIME="$2"
    shift # past argument
    shift # past value
    ;;
    -fcio|--fusion_check_info_output)
    FUSION_CHECK_INFO_OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
    -fcvo|--fusion_check_vcf_output)
    FUSION_CHECK_VCF_OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
    -fcs|--fusion_check_script)
    FUSION_CHECK_SCRIPT="$2"
    shift # past argument
    shift # past value
    ;;
    -fct|--fusion_check_threads)
    FUSION_CHECK_THREADS="$2"
    shift # past argument
    shift # past value
    ;;
    -fchv|--fusion_check_h_vmem)
    FUSION_CHECK_MEMORY="$2"
    shift # past argument
    shift # past value
    ;;
    -fcht|--fusion_check_h_rt)
    FUSION_CHECK_TIME="$2"
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
if [ -z $BAM ]; then
    echo "Missing -b|--bam parameter"
    usage
    exit
fi

if [ -z $FUSION_CHECK_VCF_OUTPUT ]; then
    FUSION_CHECK_VCF_OUTPUT=./$(basename $VCF)
    FUSION_CHECK_VCF_OUTPUT=$OUTPUTDIR/${FUSION_CHECK_VCF_OUTPUT/.vcf/_FusionGenes.vcf}
fi

if [ -z $FUSION_CHECK_INFO_OUTPUT ]; then
    FUSION_CHECK_INFO_OUTPUT=./$(basename $VCF)
    FUSION_CHECK_INFO_OUTPUT=$OUTPUTDIR/${FUSION_CHECK_INFO_OUTPUT/.vcf/_FusionGenesInfo.txt}
fi

echo `date`: Running on `uname -n`

VCF_NAME=$(basename $VCF)
VCF_NAME=${VCF_NAME/.vcf/}
VCF_OUTPUT=$OUTPUTDIR/${VCF_NAME}_FusionGenes.vcf
INFO_OUTPUT=$OUTPUTDIR/${VCF_NAME}_FusionGenesInfo.txt

PIPELINE_DIR=$NANOFG_DIR/pipeline
SCRIPT_DIR=$NANOFG_DIR/scripts
JOB_DIR=$OUTPUTDIR/jobs
LOG_DIR=$OUTPUTDIR/logs

FUSION_READ_EXTRACTION_JOBNAME=${VCF_NAME}_FUSION_READ_EXTRACTION
FUSION_READ_EXTRACTION_SH=$JOB_DIR/$FUSION_READ_EXTRACTION_JOBNAME.sh
FUSION_READ_EXTRACTION_ERR=$LOG_DIR/$FUSION_READ_EXTRACTION_JOBNAME.err
FUSION_READ_EXTRACTION_LOG=$LOG_DIR/$FUSION_READ_EXTRACTION_JOBNAME.log

CONSENSUS_CALLING_JOBNAME=${VCF_NAME}_CONSENSUS_CALLING
CONSENSUS_CALLING_SH=$JOB_DIR/$CONSENSUS_CALLING_JOBNAME.sh
CONSENSUS_CALLING_ERR=$LOG_DIR/$CONSENSUS_CALLING_JOBNAME.err
CONSENSUS_CALLING_LOG=$LOG_DIR/$CONSENSUS_CALLING_JOBNAME.log

LAST_MAPPING_JOBNAME=${VCF_NAME}_LAST_MAPPING
LAST_MAPPING_SH=$JOB_DIR/$LAST_MAPPING_JOBNAME.sh
LAST_MAPPING_ERR=$LOG_DIR/$LAST_MAPPING_JOBNAME.err
LAST_MAPPING_LOG=$LOG_DIR/$LAST_MAPPING_JOBNAME.log

MERGE_BAMS_JOBNAME=${VCF_NAME}_MERGE_BAMS
MERGE_BAMS_SH=$JOB_DIR/$MERGE_BAMS_JOBNAME.sh
MERGE_BAMS_ERR=$LOG_DIR/$MERGE_BAMS_JOBNAME.err
MERGE_BAMS_LOG=$LOG_DIR/$MERGE_BAMS_JOBNAME.log
MERGE_BAMS_OUT=$OUTPUTDIR/consensus_last.sorted.bam

SV_CALLING_JOBNAME=${VCF_NAME}_SV_CALLING
SV_CALLING_SH=$JOB_DIR/$SV_CALLING_JOBNAME.sh
SV_CALLING_ERR=$LOG_DIR/$SV_CALLING_JOBNAME.err
SV_CALLING_LOG=$LOG_DIR/$SV_CALLING_JOBNAME.log
SV_CALLING_OUT=$OUTPUTDIR/consensus_nanosv.vcf

FUSION_CHECK_JOBNAME=${VCF_NAME}_FUSION_CHECK
FUSION_CHECK_SH=$JOB_DIR/$FUSION_CHECK_JOBNAME.sh
FUSION_CHECK_ERR=$LOG_DIR/$FUSION_CHECK_JOBNAME.err
FUSION_CHECK_LOG=$LOG_DIR/$FUSION_CHECK_JOBNAME.log

CHECK_NANOFG_JOBNAME=${VCF_NAME}_CHECKSHARC
CHECK_NANOFG_SH=$JOB_DIR/$CHECK_NANOFG_JOBNAME.sh
CHECK_NANOFG_ERR=$LOG_DIR/$CHECK_NANOFG_JOBNAME.err
CHECK_NANOFG_LOG=$LOG_DIR/$CHECK_NANOFG_JOBNAME.log
CHECK_NANOFG_OUT=$OUTPUTDIR/$VCF_NAME'.check'

mkdir -p $CANDIDATE_DIR
if [ ! -d $CANDIDATE_DIR ]; then
    exit
fi

mkdir -p $JOB_DIR
if [ ! -d $JOB_DIR ]; then
    exit
fi

mkdir -p $LOG_DIR
if [ ! -d $LOG_DIR ]; then
    exit
fi

fusion_read_extraction(){

cat << EOF > $FUSION_READ_EXTRACTION_SH
#!/bin/bash

#$ -N $FUSION_READ_EXTRACTION_JOBNAME
#$ -cwd
#$ -pe threaded $FUSION_READ_EXTRACTION_THREADS
#$ -l h_vmem=$FUSION_READ_EXTRACTION_MEMORY
#$ -l h_rt=$FUSION_READ_EXTRACTION_TIME
#$ -e $FUSION_READ_EXTRACTION_ERR
#$ -o $FUSION_READ_EXTRACTION_LOG

echo \`date\`: Running on \`uname -n\`

VCF_NO_INS=${VCF/.vcf/_noINS.vcf}
VCF_NO_INS=${OUTPUTDIR}/\$(basename \$VCF_NO_INS)

if [ -z $DONT_FILTER ];then
  grep "^#" $VCF > \$VCF_NO_INS
  grep -v "^#" $VCF | awk '\$5!="<INS>"' | awk '\$7=="PASS"' >> \$VCF_NO_INS
else
  grep "^#" $VCF > \$VCF_NO_INS
  grep -v "^#" $VCF | awk '\$5!="<INS>"' >> \$VCF_NO_INS
fi

if [ ! -e ${FUSION_READ_EXTRACTION_JOBNAME}.done ]; then
  . $VENV
  python $FUSION_READ_EXTRACTION_SCRIPT \
  -b $BAM \
  -v \$VCF_NO_INS \
  -o $CANDIDATE_DIR/

  FINISHED="\$(tail -n 2 ${FUSION_READ_EXTRACTION_ERR} | grep -o "End\|Done" | wc -l | grep -oP "(^\d+)")"

  if [ \$FINISHED==2 ]; then
    touch $LOG_DIR/${FUSION_READ_EXTRACTION_JOBNAME}.done
  else
    echo "Fusion read extraction did not complete; Increase FUSION_READ_EXTRACTION_MEMORY or FUSION_READ_EXTRACTION_TIME" >&2
    exit
  fi
fi

echo \`date\`: Done
EOF
qsub $FUSION_READ_EXTRACTION_SH
}


#PUT OUTPUTFILE FROM THIS SCRIPT TO MAKE OUTPUT, OTHERWISE CHECKING IF FILES EXIST IS HARDER
consensus_calling(){

cat << EOF > $CONSENSUS_CALLING_SH
#!/bin/bash

#$ -N $CONSENSUS_CALLING_JOBNAME
#$ -cwd
#$ -pe threaded $CONSENSUS_CALLING_THREADS
#$ -l h_vmem=$CONSENSUS_CALLING_MEMORY
#$ -l h_rt=$CONSENSUS_CALLING_TIME
#$ -e $CONSENSUS_CALLING_ERR
#$ -o $CONSENSUS_CALLING_LOG
#$ -hold_jid $FUSION_READ_EXTRACTION_JOBNAME

echo \`date\`: Running on \`uname -n\`

if [ -e $LOG_DIR/${FUSION_READ_EXTRACTION_JOBNAME}.done ]; then
  if [ ! -e $LOG_DIR/${CONSENSUS_CALLING_JOBNAME}.done ]; then
    for FASTA in $CANDIDATE_DIR/\$SGE_TASK_ID/*.fasta; do
      bash $PIPELINE_DIR/consensus_calling.sh \
      -f \$FASTA \
      -t $CONSENSUS_CALLING_THREADS \
      -w $WTDBG2_DIR \
      -ws '$CONSENSUS_CALLING_WTDBG2_SETTINGS'
    done

    NUMBER_FASTA_INPUT=\$(ls $CANDIDATE_DIR/*.fasta | wc -l | grep -oP "(^\d+)")
    NUMBER_CONTIG_OUTPUT=\$(ls $CANDIDATE_DIR/*.ctg.fa | wc -l | grep -oP "(^\d+)")

    if [ \$NUMBER_FASTA_INPUT==\$NUMBER_CONTIG_OUTPUT ];then
      touch $LOG_DIR/${CONSENSUS_CALLING_JOBNAME}.done
    else
      echo "Number of made contig files is not equal to the nubmer of input fasta files " >&2
    fi

  fi
fi

echo \`date\`: Done
EOF
qsub $CONSENSUS_CALLING_SH
}

last_mapping() {

cat << EOF > $LAST_MAPPING_SH
#!/bin/bash

#$ -N $LAST_MAPPING_JOBNAME
#$ -cwd
#$ -pe threaded $LAST_MAPPING_THREADS
#$ -l h_vmem=$LAST_MAPPING_MEMORY
#$ -l h_rt=$LAST_MAPPING_TIME
#$ -e $LAST_MAPPING_ERR
#$ -o $LAST_MAPPING_LOG
#$ -hold_jid $CONSENSUS_CALLING_JOBNAME

echo \`date\`: Running on \`uname -n\`

if [ -e $LOG_DIR/${CONSENSUS_CALLING_JOBNAME}.done ]; then
  if [ ! -e $LOG_DIR/${LAST_MAPPING_JOBNAME}.done ]; then
    LAST_MAPPING_ARGS="-t $LAST_MAPPING_THREADS -r $LAST_MAPPING_REFGENOME -rd $LAST_MAPPING_REFDICT -l $LAST_DIR -ls '$LAST_MAPPING_SETTINGS' -s $SAMBAMBA"

    for FA in $CANDIDATE_DIR/*.fa; do
        echo \$FA;
    done | \
    xargs -I{} --max-procs $THREADS bash -c "echo 'Start' {}; bash $PIPELINE_DIR/last_mapping.sh -f \$FA \$LAST_MAPPING_ARGS {}; echo 'Done' {}; exit 1;"

    NUMBER_FA_INPUT=\$(ls $CANDIDATE_DIR/*.fa | wc -l | grep -oP "(^\d+)")
    NUMBER_BAM_OUTPUT=\$(ls $CANDIDATE_DIR/*.ctg.last.sorted.bam | wc -l | grep -oP "(^\d+)")

    if [ \$NUMBER_FA_INPUT==\$NUMBER_BAM_OUTPUT ];then
      touch $LOG_DIR/${LAST_MAPPING_JOBNAME}.done
    else
      echo "Number of output bam files is not equal to the number of input contig fasta files" >&2
    fi
  fi
fi
EOF
qsub $LAST_MAPPING_SH
}

bam_merge() {

cat << EOF > $MERGE_BAMS_SH
#!/bin/bash

#$ -N $MERGE_BAMS_JOBNAME
#$ -cwd
#$ -pe threaded $BAM_MERGE_THREADS
#$ -l h_vmem=$BAM_MERGE_MEMORY
#$ -l h_rt=$BAM_MERGE_TIME
#$ -e $MERGE_BAMS_ERR
#$ -o $MERGE_BAMS_LOG
#$ -hold_jid $LAST_MAPPING_JOBNAME

echo \`date\`: Running on \`uname -n\`

if [ -e $LOG_DIR/$LAST_MAPPING_JOBNAME.done ]; then
  if [ ! -e $LOG_DIR/$MERGE_BAMS_JOBNAME.done ]; then
    bash $PIPELINE_DIR/bam_merge.sh \
      -d $CANDIDATE_DIR \
      -s $SAMBAMBA \
      -o $MERGE_BAMS_OUT

    if [ -e $MERGE_BAMS_OUT ];then
      touch $LOG_DIR/$MERGE_BAMS_JOBNAME.done
    fi
  fi
fi
echo \`date\`: Done
EOF
qsub $MERGE_BAMS_SH
}

sv_calling() {
cat << EOF > $SV_CALLING_SH
#!/bin/bash
#$ -N $SV_CALLING_JOBNAME
#$ -cwd
#$ -pe threaded $SV_CALLING_THREADS
#$ -l h_vmem=$SV_CALLING_MEMORY
#$ -l h_rt=$SV_CALLING_TIME
#$ -e $SV_CALLING_ERR
#$ -o $SV_CALLING_LOG
EOF

if [ ! -z $MERGE_BAMS_JOBNAME ]; then
cat << EOF >> $SV_CALLING_SH
#$ -hold_jid $MERGE_BAMS_JOBNAME
EOF
fi

cat << EOF >> $SV_CALLING_SH
echo \`date\`: Running on \`uname -n\`
if [ -e $LOG_DIR/$MERGE_BAMS_JOBNAME.done ]; then
    bash $PIPELINE_DIR/sv_calling.sh \
      -b $MERGE_BAMS_OUT \
      -n $NANOSV \
      -t $SV_CALLING_THREADS \
      -s $SAMBAMBA \
      -v $VENV \
      -c $SV_CALLING_CONFIG \
      -o $SV_CALLING_OUT
    NUMBER_OF_LINES_VCF=\$(grep -v "^#" $SV_CALLING_OUT | wc -l | grep -oP "(^\d+)")
    if [ \$NUMBER_OF_LINES_VCF != 0 ]; then
      touch $LOG_DIR/$SV_CALLING_JOBNAME.done
    fi
fi
echo \`date\`: Done
EOF
qsub $SV_CALLING_SH
}

fusion_check() {
cat << EOF > $FUSION_CHECK_SH
#!/bin/bash
#$ -N $FUSION_CHECK_JOBNAME
#$ -cwd
#$ -pe threaded $FUSION_CHECK_THREADS
#$ -l h_vmem=$FUSION_CHECK_MEMORY
#$ -l h_rt=$FUSION_CHECK_TIME
#$ -e $FUSION_CHECK_ERR
#$ -o $FUSION_CHECK_LOG
#$ -hold_jid $SV_CALLING_JOBNAME

echo \`date\`: Running on \`uname -n\`
if [ -e $LOG_DIR/$SV_CALLING_JOBNAME.done ];then
  bash $PIPELINE_DIR/fusion_check.sh \
    -v $SV_CALLING_OUT \
    -o $FUSION_CHECK_VCF_OUTPUT \
    -fo $FUSION_CHECK_INFO_OUTPUT \
    -s $FUSION_CHECK_SCRIPT \
    -e $VENV

  NUMBER_VCF_INPUT=\$(grep -v "^#" $SV_CALLING_OUT | wc -l | grep -oP "(^\d+)")
  NUMBER_VCF_OUTPUT=\$(grep -v "^#" $FUSION_CHECK_VCF_OUTPUT | wc -l | grep -oP "(^\d+)")

  FINISHED="\$(tail -n 2 $FUSION_CHECK_LOG | grep -o "End\|Done" | wc -l | grep -oP "(^\d+)")"

  if [ \$FINISHED==2 ]; then
    if [ $NUMBER_OF_LINES_VCF_1 != $NUMBER_OF_LINES_VCF_2 ]; then
      echo "Number of lines in all split VCF files is not equal to the number of lines in the original VCF file"
      exit
    else
      touch $LOG_DIR/$FUSION_CHECK_JOBNAME.done
    fi
  else
    echo "Fusion check did not complete; Increase FUSION_CHECK_MEMORY or FUSION_CHECK_TIME" >&2
    exit
  fi
fi
echo \`date\`: Done
EOF
qsub $FUSION_CHECK_SH
}

check_NanoFG() {
cat << EOF > $CHECK_NANOFG_SH
#!/bin/bash
#$ -N $CHECK_NANOFG_JOBNAME
#$ -cwd
#$ -pe threaded $CHECK_NANOFG_THREADS
#$ -l h_vmem=$CHECK_NANOFG_MEMORY
#$ -l h_rt=$CHECK_NANOFG_TIME
#$ -e $CHECK_NANOFG_ERR
#$ -o $CHECK_NANOFG_LOG
#$ -hold_jid $FUSION_CHECK_JOBNAME

echo \`date\`: Running on \`uname -n\`
CHECK_BOOL=true
echo "------------------------------------------------" >> $CHECK_NANOFG_OUT
echo "\`date\`" >> $CHECK_NANOFG_OUT
echo "Sample name: $VCF_NAME" >> $CHECK_NANOFG_OUT

if [ -e $LOG_DIR/$VCF_SPLIT_JOBNAME.done ]; then
    echo "Vcf split: Done" >> $CHECK_NANOFG_OUT
else
    echo "Vcf split: Fail" >> $CHECK_NANOFG_OUT
    CHECK_BOOL=false
fi

NUMBER_OF_FUSION_READ_EXTRACTION_JOBS=\$(ls $LOG_DIR/$FUSION_READ_EXTRACTION_JOBNAME_*.done | wc -l | grep -oP "(^\d+)")
if [ \$NUMBER_OF_FUSION_READ_EXTRACTION_JOBS==$NUMBER_OF_FILES ] && [ -e $LOG_DIR/$CONSENSUS_CALLING_JOBNAME.done ]; then
    echo "Fusion read extraction: Done" >> $CHECK_NANOFG_OUT
else
  echo "Fusion read extraction: Fail" >> $CHECK_NANOFG_OUT
  CHECK_BOOL=false
fi

NUMBER_OF_CONSENSUS_CALLING_JOBS=\$(ls $LOG_DIR/$CONSENSUS_CALLING_JOBNAME_*.done | wc -l | grep -oP "(^\d+)")
if [ \$NUMBER_OF_CONSENSUS_CALLING_JOBS==$NUMBER_OF_FILES ] && [ -e $LOG_DIR/$CONSENSUS_CALLING_JOBNAME.done ]; then
    echo "Consensus mapping: Done" >> $CHECK_NANOFG_OUT
else
    echo "Consensus mapping: Fail" >> $CHECK_NANOFG_OUT
    CHECK_BOOL=false
fi

if [ -e $LOG_DIR/$MERGE_BAMS_JOBNAME.done ]; then
    echo "Bam merge: Done" >> $CHECK_NANOFG_OUT
else
    echo "Bam merge: Fail" >> $CHECK_NANOFG_OUT
    CHECK_BOOL=false
fi

if [ -e $LOG_DIR/$SV_CALLING_JOBNAME.done ]; then
    echo "SV calling: Done" >> $CHECK_NANOFG_OUT
else
    echo "SV calling: Fail" >> $CHECK_NANOFG_OUT
    CHECK_BOOL=false
fi

if [ -e $LOG_DIR/$FUSION_CHECK_JOBNAME.done ]; then
    echo "Fusion check: Done" >> $CHECK_NANOFG_OUT
else
    echo "Fusion check: Fail" >> $CHECK_NANOFG_OUT
    CHECK_BOOL=false
fi

if [ \$CHECK_BOOL = true ]; then
  echo "Fusion check: Done" >> $CHECK_NANOFG_OUT
    touch $LOG_DIR/$CHECK_NANOFG_JOBNAME.done
    if [ $DONT_CLEAN = false ]; then
      rm -rf $SPLITDIR
      rm $MERGE_BAMS_OUT
      rm $SV_CALLING_OUT
    fi
else
  echo "NanoFG check: Fail" >> $CHECK_NANOFG_OUT
fi

if [ -z $MAIL ]; then
  tac $CHECK_NANOFG_OUT | sed '/^Qsub/q' | tac | mail -s 'NANOFG_${VCF_NAME}' $MAIL
fi

echo \`date\`: Done

sleep 20
EOF
qsub $CHECK_NANOFG_SH
}


if [ ! -e $LOG_DIR/$VCF_SPLIT_JOBNAME.done ]; then
    vcf_split
fi
if [ ! -e $LOG_DIR/$FUSION_READ_EXTRACTION_JOBNAME.done ]; then
    fusion_read_extraction
fi
if [ ! -e $LOG_DIR/$CONSENSUS_CALLING_JOBNAME.done ]; then
    consensus_calling
fi
if [ ! -e $LOG_DIR/$MERGE_BAMS_JOBNAME.done ]; then
    bam_merge
fi
if [ ! -e $LOG_DIR/$SV_CALLING_JOBNAME.done ]; then
    sv_calling
fi
if [ ! -e $LOG_DIR/$FUSION_CHECK_JOBNAME.done ]; then
    fusion_check
fi
if [ ! -e $LOG_DIR/$CHECK_NANOFG_JOBNAME.done ]; then
    check_NanoFG
fi

echo `date`: Done
