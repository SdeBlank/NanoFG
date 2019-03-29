#!/bin/bash

usage() {
echo "
Required parameters:
    -v|--vcf		                                                                Path to vcf file
    -b|--bam                                                                    Path to bam file

Optional parameters:

GENERAL
    -h|--help                                                                     Shows help
    -o|--outputdir                                                                Path to output directory
    -d|--nanofg_dir                                                         Directory that contains NanoFG [$NANOFG_DIR]
    -e|--venv                                                                     Path to virtual environment[$VENV]
    -dc|--dont_clean                                                              Don't clean up the intermediate files []

REQUIRED TOOLS
    -s|--sambamba                                                                 Path to sambamba|samtools [${SAMBAMBA}]
    -l|--last_dir                                                                 Path to LAST directory [${LAST_DIR}]
    -w|--wtdbg2_dir                                                               Path to wtdbg2 directory [${WTDBG2_DIR}]

VCF SPLIT
    -vsl|--vcf_split_lines                                                         Number of files to split the original vcf into
    -vst|vcf_split_threads                                                        Number of threads
    -vshv|vcf_split_h_vmem                                                        Vcf split memory
    -vshr|vcf_split_h_rt                                                          Vcf split time

FUSION READ EXTRACTION
    -fres|--fusion_read_extraction_script
    -fret|--fusion_read_extraction_threads                                        Number of threads
    -frehv|--fusion_read_extraction_h_vmem                                        Fusion read extraction memory
    -freht|--fusion_read_extraction_h_rt                                          Fusion read extraction time

CONSENSUS MAPPING
    -cmr|--consensus_mapping_refgenome                                            Reference genome [${REF}]
    -cmrd|--consensus_mapping_refdict                                             Reference genome .dict file [${REF_DICT}]
    -cmws|--consensus_mapping_wtdbg2_settings                                     wtdbg2 settings [${WTDBG2_SETTINGS}]
    -cmls|--consensus_mapping_last_settings                                       LAST settings [${LAST_SETTINGS}]
    -cmt|--consensus_mapping_threads                                              Number of threads [${THREADS}]
    -cmhv|--consensus_mapping_h_vmem                                              Number of threads [${THREADS}]
    -cmhr|--consensus_mapping_h_rt                                                Number of threads [${THREADS}]

BAM MERGE
    -bmt|vcf_split_threads                                                        Number of threads
    -bmhv|vcf_split_h_vmem                                                        Bam merge memory
    -bmhr|vcf_split_h_rt                                                          Bam merge time

SV CALLING
    -sct|sv_calling_threads                                                        Number of threads
    -schv|sv_calling_h_vmem                                                        SV calling memory
    -schr|sv_calling_h_rt                                                          SV calling time

FUSION CHECK
    -fcs|--fusion_check_script                                                    Path to vcf_primer_filter.py [$SCRIPT]
    -fct|--fusion_check_threads                                                   VCF output file [$OUTPUT]
    -fchv|--fusion_check_h_vmem                                                   VCF output file [$OUTPUT]
    -fcht|--fusion_check_h_rt                                                     VCF output file [$OUTPUT]
"
}

POSITIONAL=()

#GENERAL DEFAULTS
NANOFG_DIR=$(realpath $(dirname ${BASH_SOURCE[0]}))
PIPELINE_DIR=$NANOFG_DIR/pipeline
FILES_DIR=$NANOFG_DIR/files
SCRIPT_DIR=$NANOFG_DIR/scripts
VENV=${NANOFG_DIR}/venv/bin/activate

OUTPUTDIR='./'
DONT_CLEAN=false

#TOOL PATH DEFAULTS
SAMBAMBA=/hpc/local/CentOS7/cog_bioinf/sambamba_v0.6.5/sambamba
LAST_DIR=/hpc/cog_bioinf/kloosterman/tools/last-921
WTDBG2_DIR=/hpc/cog_bioinf/kloosterman/tools/wtdbg2_v2.2

#VCF SPLIT DEFAULTS
VCF_SPLIT_THREADS=1
VCF_SPLIT_TIME=0:30:0
VCF_SPLIT_MEMORY=10G

#CREATE FUSION READ EXTRACTION JOBS DEFAULTS
CREATE_FUSION_READ_EXTRACTION_JOBS_SCRIPT=$NANOFG_DIR/pipeline/fusion_gene_read_extraction.py
CREATE_FUSION_READ_EXTRACTION_JOBS_THREADS=8
CREATE_FUSION_READ_EXTRACTION_JOBS_TIME=0:10:0
CREATE_FUSION_READ_EXTRACTION_JOBS_MEMORY=5G

#FUSION READ EXTRACTION DEFAULTS
FUSION_READ_EXTRACTION_SCRIPT=$NANOFG_DIR/pipeline/fusion_gene_read_extraction.py
FUSION_READ_EXTRACTION_THREADS=1
FUSION_READ_EXTRACTION_TIME=0:10:0
FUSION_READ_EXTRACTION_MEMORY=5G

#CONSENSUS MAPPING DEFAULTS
CONSENSUS_MAPPING_REFGENOME=/hpc/cog_bioinf/GENOMES/LAST/human_GATK_GRCh37
CONSENSUS_MAPPING_REFDICT=/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.dict
CONSENSUS_MAPPING_WTDBG2_SETTINGS='-x ont -g 3g'
CONSENSUS_MAPPING_LAST_SETTINGS="-Q 0 -p ${LAST_DIR}/last_params"
CONSENSUS_MAPPING_THREADS=8
CONSENSUS_MAPPING_TIME=0:15:0
CONSENSUS_MAPPING_MEMORY=20G

#BAM MERGE DEFAULTS
BAM_MERGE_THREADS=8
BAM_MERGE_TIME=0:10:0
BAM_MERGE_MEMORY=10G

#SV CALLING DEFAULTS
SV_CALLING_THREADS=8
SV_CALLING_TIME=0:10:0
SV_CALLING_MEMORY=10G
SV_CALLING_CONFIG=$FILES_DIR/nanosv_last_config.ini

#FUSION CHECK DEFAULTS
FUSION_CHECK_SCRIPT=$NANOFG_DIR/pipeline/NanoFG.py
FUSION_CHECK_THREADS=8
FUSION_CHECK_TIME=0:20:0
FUSION_CHECK_MEMORY=10G



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
    -vsl|--vcf_split_lines)
    VCF_SPLIT_LINES=$2
    shift # past argument
    shift # past value
    ;;
    -vst|vcf_split_threads)
    VCF_SPLIT_THREADS="$2"
    shift # past argument
    shift # past value
    ;;
    -vshv|vcf_split_h_vmem)
    VCF_SPLIT_MEMORY="$2"
    shift # past argument
    shift # past value
    ;;
    -vshr|vcf_split_h_rt)
    VCF_SPLIT_TIME="$2"
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
    -cmr|--consensus_mapping_refgenome)
    CONSENSUS_MAPPING_REFGENOME="$2"
    shift # past argument
    shift # past value
    ;;
    -cmrd|--consensus_mapping_refdict)
    CONSENSUS_MAPPING_REFDICT="$2"
    shift # past argument
    shift # past value
    ;;
    -cmws|--consensus_mapping_wtdbg2_settings)
    CONSENSUS_MAPPING_WTDBG2_SETTINGS="$2"
    shift # past argument
    shift # past value
    ;;
    -cmls|--consensus_mapping_last_settings)
    CONSENSUS_MAPPING_LAST_SETTINGS="$2"
    shift # past argument
    shift # past value
    ;;
    -cmt|--consensus_mapping_threads)
    CONSENSUS_MAPPING_THREADS="$2"
    shift # past argument
    shift # past value
    ;;
    -cmhv|--consensus_mapping_h_vmem)
    CONSENSUS_MAPPING_MEMORY="$2"
    shift # past argument
    shift # past value
    ;;
    -cmhr|--consensus_mapping_h_rt)
    CONSENSUS_MAPPING_TIME="$2"
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
NUMBER_OF_SVS=$(grep -vc "^#" $VCF | grep -oP "(^\d+)")

if [ -z $VCF_SPLIT_LINES ]; then
  VCF_SPLIT_LINES=$(expr $NUMBER_OF_SVS / 100 + 1)
  NUMBER_OF_FILES=100
  if [ $VCF_SPLIT_LINES -lt 100 ]; then
    VCF_SPLIT_LINES=100
    NUMBER_OF_FILES=$(expr $NUMBER_OF_SVS / $VCF_SPLIT_LINES + 1)
  fi
else
  echo $NUMBER_OF_SVS
  echo $VCF_SPLIT_LINES
  NUMBER_OF_FILES=$(expr $NUMBER_OF_SVS / $VCF_SPLIT_LINES + 1)
fi

echo `date`: Running on `uname -n`

VCF_NAME=$(basename $VCF)
VCF_NAME=${VCF_NAME/.vcf/}
VCF_OUTPUT=$OUTPUTDIR/${VCF_NAME}_FusionGenes.vcf
INFO_OUTPUT=$OUTPUTDIR/${VCF_NAME}_FusionGenesInfo.txt

PIPELINE_DIR=$NANOFG_DIR/pipeline
SCRIPT_DIR=$NANOFG_DIR/scripts

JOBDIR=$OUTPUTDIR/jobs
LOGDIR=$OUTPUTDIR/logs

VCF_SPLIT_OUTDIR=$OUTPUTDIR/split_vcf
VCF_SPLIT_JOBNAME=${VCF_NAME}_SPLIT_VCF
VCF_SPLIT_SH=$JOBDIR/$VCF_SPLIT_JOBNAME.sh
VCF_SPLIT_ERR=$LOGDIR/$VCF_SPLIT_JOBNAME.err
VCF_SPLIT_LOG=$LOGDIR/$VCF_SPLIT_JOBNAME.log

CREATE_FUSION_READ_EXTRACTION_JOBS_JOBNAME=${VCF_NAME}_CREATE_FUSION_READ_EXTRACTION_JOBS
CREATE_FUSION_READ_EXTRACTION_JOBS_SH=$JOBDIR/$CREATE_FUSION_READ_EXTRACTION_JOBS_JOBNAME.sh
CREATE_FUSION_READ_EXTRACTION_JOBS_ERR=$LOGDIR/$CREATE_FUSION_READ_EXTRACTION_JOBS_JOBNAME.err
CREATE_FUSION_READ_EXTRACTION_JOBS_LOG=$LOGDIR/$CREATE_FUSION_READ_EXTRACTION_JOBS_JOBNAME.log

FUSION_READ_EXTRACTION_JOBNAME=${VCF_NAME}_FUSION_READ_EXTRACTION
FUSION_READ_EXTRACTION_SH=$JOBDIR/$FUSION_READ_EXTRACTION_JOBNAME.sh
FUSION_READ_EXTRACTION_ERR=$LOGDIR/$FUSION_READ_EXTRACTION_JOBNAME.err
FUSION_READ_EXTRACTION_LOG=$LOGDIR/$FUSION_READ_EXTRACTION_JOBNAME.log

CONSENSUS_MAPPING_JOBNAME=${VCF_NAME}_CONSENSUS_MAPPING
CONSENSUS_MAPPING_SH=$JOBDIR/$CONSENSUS_MAPPING_JOBNAME.sh
CONSENSUS_MAPPING_ERR=$LOGDIR/$CONSENSUS_MAPPING_JOBNAME.err
CONSENSUS_MAPPING_LOG=$LOGDIR/$CONSENSUS_MAPPING_JOBNAME.log

MERGE_BAMS_JOBNAME=${VCF_NAME}_MERGE_BAMS
MERGE_BAMS_SH=$JOBDIR/$MERGE_BAMS_JOBNAME.sh
MERGE_BAMS_ERR=$LOGDIR/$MERGE_BAMS_JOBNAME.err
MERGE_BAMS_LOG=$LOGDIR/$MERGE_BAMS_JOBNAME.log
MERGE_BAMS_OUT=$OUTPUTDIR/consensus_last.sorted.bam

SV_CALLING_JOBNAME=${VCF_NAME}_SV_CALLING
SV_CALLING_SH=$JOBDIR/$SV_CALLING_JOBNAME.sh
SV_CALLING_ERR=$LOGDIR/$SV_CALLING_JOBNAME.err
SV_CALLING_LOG=$LOGDIR/$SV_CALLING_JOBNAME.log
SV_CALLING_OUT=$OUTPUTDIR/consensus_nanosv.vcf

FUSION_CHECK_JOBNAME=${VCF_NAME}_FUSION_CHECK
FUSION_CHECK_SH=$JOBDIR/$FUSION_CHECK_JOBNAME.sh
FUSION_CHECK_ERR=$LOGDIR/$FUSION_CHECK_JOBNAME.err
FUSION_CHECK_LOG=$LOGDIR/$FUSION_CHECK_JOBNAME.log

mkdir -p $VCF_SPLIT_OUTDIR
if [ ! -d $VCF_SPLIT_OUTDIR ]; then
    exit
fi

mkdir -p $JOBDIR
if [ ! -d $JOBDIR ]; then
    exit
fi

mkdir -p $LOGDIR
if [ ! -d $LOGDIR ]; then
    exit
fi

vcf_split(){
cat << EOF > $VCF_SPLIT_SH
#!/bin/bash

#$ -N $VCF_SPLIT_JOBNAME
#$ -cwd
#$ -pe threaded $VCF_SPLIT_THREADS
#$ -l h_vmem=$VCF_SPLIT_MEMORY
#$ -l h_rt=$VCF_SPLIT_TIME
#$ -e $VCF_SPLIT_ERR
#$ -o $VCF_SPLIT_LOG

echo \`date\`: Running on \`uname -n\`

if [ ! -e $LOGDIR/$VCF_SPLIT_JOBNAME.done ];then
  if [ -e $VCF ];then
    bash $PIPELINE_DIR/vcf_split.sh \\
    -v $VCF \\
    -d $VCF_SPLIT_OUTDIR \\
    -l $VCF_SPLIT_LINES \\
  else
    echo "VCF file ($VCF) does not exist"
    exit
  fi

  NUMBER_OF_LINES_VCF=\$(grep -v "^#" $VCF | wc -l | grep -oP "(^\d+)")
  NUMBER_OF_LINES_SPLIT_VCFS=\$(cat $VCF_SPLIT_OUTDIR/*.vcf | grep -v "^#" | wc -l | grep -oP "(^\d+)")

  if [ \$NUMBER_OF_LINES_VCF == \$NUMBER_OF_LINES_SPLIT_VCFS ]; then
    touch $LOGDIR/$VCF_SPLIT_JOBNAME.done
  else
    echo "The total number of SVs in the split vcf files is different from the number of SVs in the original vcf" >&2
  fi
else
  echo "$VCF_SPLIT_JOBNAME has been completed already"
fi

echo \`date\`: Done
EOF
qsub $VCF_SPLIT_SH
}

create_fusion_read_extraction_jobs(){

cat << EOF > $FUSION_READ_EXTRACTION_SH
#!/bin/bash

#$ -N $FUSION_READ_EXTRACTION_JOBNAME
#$ -cwd
#$ -t 1:$NUMBER_OF_FILES:1
#$ -pe threaded $FUSION_READ_EXTRACTION_THREADS
#$ -l h_vmem=$FUSION_READ_EXTRACTION_MEMORY
#$ -l h_rt=$FUSION_READ_EXTRACTION_TIME
#$ -hold_jid $VCF_SPLIT_JOBNAME

echo \`date\`: Running on \`uname -n\`

if [ -e $LOGDIR/$VCF_SPLIT_JOBNAME.done ]; then
  if [ ! -e $FUSION_READ_EXTRACTION_JOBNAME.done ]; then
    mkdir -p $VCF_SPLIT_OUTDIR/\$SGE_TASK_ID
    python $SCRIPT_DIR/fusion_gene_read_extraction.py \
    -b $BAM \
    -v $VCF_SPLIT_OUTDIR/\$SGE_TASK_ID.vcf \
    -o $VCF_SPLIT_OUTDIR/\$SGE_TASK_ID

    touch $LOGDIR/$FUSION_READ_EXTRACTION_JOBNAME.done
  fi
fi
EOF
qsub $FUSION_READ_EXTRACTION_SH
}

consensus_mapping(){

cat << EOF > $CONSENSUS_MAPPING_SH
#!/bin/bash

#$ -N $CONSENSUS_MAPPING_JOBNAME
#$ -cwd
#$ -t 1:$NUMBER_OF_FILES:1
#$ -pe threaded $CONSENSUS_MAPPING_THREADS
#$ -l h_vmem=$CONSENSUS_MAPPING_MEMORY
#$ -l h_rt=$CONSENSUS_MAPPING_TIME
#$ -hold_jid $FUSION_READ_EXTRACTION_JOBNAME

echo \`date\`: Running on \`uname -n\`

if [ -e $LOGDIR/$FUSION_READ_EXTRACTION_JOBNAME.done ]; then
  if [ ! -e $LOGDIR/$CONSENSUS_MAPPING_JOBNAME.done ]; then
    for FASTA in $VCF_SPLIT_OUTDIR/\$SGE_TASK_ID/*.fasta; do
      bash $PIPELINE_DIR/Consensus_building_and_mapping.sh \
      -f \$FASTA \
      -t $CONSENSUS_MAPPING_THREADS \
      -r $CONSENSUS_MAPPING_REFGENOME \
      -rd $CONSENSUS_MAPPING_REFDICT \
      -w $WTDBG2_DIR \
      -ws '$CONSENSUS_MAPPING_WTDBG2_SETTINGS' \
      -l $LAST_DIR \
      -ls '$CONSENSUS_MAPPING_LAST_SETTINGS' \
      -s $SAMBAMBA
    done

    touch $LOGDIR/$CONSENSUS_MAPPING_JOBNAME.done

  fi
fi
EOF
qsub $CONSENSUS_MAPPING_SH
}

bam_merge(){

cat << EOF > $MERGE_BAMS_SH
#!/bin/bash

#$ -N $MERGE_BAMS_JOBNAME
#$ -cwd
#$ -pe threaded $BAM_MERGE_THREADS
#$ -l h_vmem=$BAM_MERGE_MEMORY
#$ -l h_rt=$BAM_MERGE_TIME
#$ -e $MERGE_BAMS_ERR
#$ -o $MERGE_BAMS_LOG
#$ -hold_jid $CONSENSUS_MAPPING_JOBNAME

echo \`date\`: Running on \`uname -n\`

if [ -e $LOGDIR/$CONSENSUS_MAPPING_JOBNAME.done ]; then
  if [ ! -e $LOGDIR/$MERGE_BAMS_JOBNAME.done ]; then
    bash $PIPELINE_DIR/bam_merge.sh \
      -d $VCF_SPLIT_OUTDIR \
      -s $SAMBAMBA \
      -o $MERGE_BAMS_OUT

    touch $LOGDIR/$MERGE_BAMS_JOBNAME.done
  fi
fi
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
if [ -e $MERGE_BAMS_JOBNAME.done ]; then
    bash $PIPELINE_DIR/sv_calling.sh -b $MERGE_BAMS_OUT -t $SV_CALLING_THREADS -s $SAMBAMBA -v $VENV -c $SV_CALLING_CONFIG -o $SV_CALLING_OUT
    NUMBER_OF_LINES_VCF=\$(grep -v "^#" $SV_CALLING_OUT | wc -l | grep -oP "(^\d+)")
    if [ \$NUMBER_OF_LINES_VCF != 0 ]; then
      touch $SV_CALLING_OUT.done
    fi
fi
echo \`date\`: Done
EOF
qsub $SV_CALLING_SH
}


#vcf_split

#create_fusion_read_extraction_jobs

#consensus_mapping

#bam_merge

sv_calling
