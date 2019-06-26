#!/bin/bash

# Create usage message
usage(){
echo "
Required parameters:
-f|--fasta                Input fasta file [${FASTA}]

Optional parameters:
-h|--help                 Shows help
-mm2|--minimap2             Path to minimap2
-t|--threads              Number of threads [${THREADS}]
-r|--refgenome            Reference genome [${REF}]
-rd|--refdict             Reference genome .dict file [${REF_DICT}]
-l|--last_dir             Path to LAST directory [${LAST_DIR}]
-ls|--last_settings       LAST settings [${LAST_SETTINGS}]
-s|--sambamba             Path to sambamba|samtools [${SAMTOOLS}]
"
exit
}

POSITIONAL=()

# DEFAULT SETTINGS
THREADS=1
REF=/hpc/cog_bioinf/GENOMES/LAST/human_GATK_GRCh37
SAMTOOLS=/hpc/local/CentOS7/cog_bioinf/sambamba_v0.6.5/sambamba
MINIMAP2=

while [[ $# -gt 0 ]]; do
  KEY="$1"
  case ${KEY} in
    -h|--help)
    usage
    shift # past argument
    ;;
    -f|--fasta)
    FASTA="$2"
    shift # past argument
    shift # past value
    ;;
    -mm2|--minimap2)
    MINIMAP2="$2"
    shift # past argument
    shift # past value
    ;;
    -t|--threads)
    THREADS="$2"
    shift # past argument
    shift # past value
    ;;
    -r|--refgenome)
    REF="$2"
    shift # past argument
    shift # past value
    -s|--samtools)
    SAMTOOLS="$2"
    shift # past argument
    shift # past value
    ;;
    *) # Unknown option
    POSITIONAL+=("$1") # Store unknown option
    shift
    ;;
  esac
done

set -- "${POSITIONAL[@]}" # Restore postional parameters

if [ -z ${FASTA} ]; then
  echo "Missing -f|--fasta parameter"
  usage
fi


#echo `date`: Running on `uname -n`

$MINIMAP2 -t $THREADS -ax map-ont $REF $FASTA\
| $SAMTOOLS view -h -S -b -t $THREADS /dev/stdin \
| $SAMTOOLS sort /dev/stdin -o ${bam/.bam/.sorted.bam} -t $THREADS

$SAMTOOLS index ${bam/.bam/.sorted.bam} -t $THREADS

#echo `date`: Done
