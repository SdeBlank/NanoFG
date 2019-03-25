#!/bin/bash

# Create usage message
usage(){
echo "
Optional parameters:

-h|--help                 Shows help
-f|--fasta                Input fasta file [${FASTA}]
-t|--threads              Number of threads [${THREADS}]
-r|--refgenome            Reference genome [${REF}]
-w|--wtdbg2_dir           Path to wtdbg2 directory [${WTDBG2_DIR}]
-ws|--wtdbg2_settings     wtdbg2 settings [${WTDBG2_SETTINGS}]
-m|--minimap2             Path to minimap2 [${MINIMAP2}]
-ms|--minimap2_settings   minimap2 settings [${MINIMAP2_SETTINGS}]
-s|--sambamba             Path to sambamba|samtools [${SAMBAMBA}]
"
exit
}

POSITIONAL=()

# DEFAULT SETTINGS
THREADS=1
REF=/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fa
WTDBG2_DIR=/hpc/cog_bioinf/kloosterman/tools/wtdbg2_v2.2
WTDBG2=${WTDBG2_DIR}/wtdbg2
WTPOA_CNS=${WTDBG2_DIR}/wtpoa-cns
WTDBG2_SETTINGS='-x ont -g 3g'
MINIMAP2=/hpc/cog_bioinf/kloosterman/tools/minimap2_v2.12/minimap2
MINIMAP2_SETTINGS='-x map-ont'
SAMBAMBA=/hpc/local/CentOS7/cog_bioinf/sambamba_v0.6.5/sambamba

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
    -t|--threads)
    THREADS="$2"
    shift # past argument
    shift # past value
    ;;
    -r|--refgenome)
    REF="$2"
    shift # past argument
    shift # past value
    ;;
    -w|--wtdbg2_dir)
    WTDBG2_DIR="$2"
    shift # past argument
    shift # past value
    ;;
    -ws|--wtdbg2_settings)
    WTDBG2_SETTINGS="$2"
    shift # past argument
    shift # past value
    ;;
    -m|--minimap2)
    MINIMAP2="$2"
    shift # past argument
    shift # past value
    ;;
    -ms|--minimap2_settings)
    MINIMAP2_SETTINGS="$2"
    shift # past argument
    shift # past value
    ;;
    -s|--sambamba)
    SAMBAMBA="$2"
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

PREFIX=${FASTA/.fasta/_wtdbg2}

WTDBG2_COMMAND="${WTDBG2} ${WTDBG2_SETTINGS} -i ${FASTA} -t ${THREADS} -o ${PREFIX}"
WTPOA_CNS_COMMAND="${WTPOA_CNS} -t ${THREADS} -i ${PREFIX}.ctg.lay.gz -o ${PREFIX}.ctg.fa"
MINIMAP2_COMMAND="${MINIMAP2} ${MINIMAP2_SETTING} -t ${THREADS} -a ${REF} ${PREFIX}.ctg.fa | ${SAMBAMBA} view -S -f bam /dev/stdin > ${PREFIX}.ctg.map.bam"
SAMBAMBA_INDEX_COMMAND="${SAMBAMBA} index ${PREFIX}.ctg.map.bam"

echo ${WTDBG2_COMMAND}
eval ${WTDBG2_COMMAND}
echo ${WTPOA_CNS_COMMAND}
eval ${WTPOA_CNS_COMMAND}
echo ${MINIMAP2_COMMAND}
eval ${MINIMAP2_COMMAND}
echo ${SAMBAMBA_INDEX_COMMAND}
eval ${SAMBAMBA_INDEX_COMMAND}
