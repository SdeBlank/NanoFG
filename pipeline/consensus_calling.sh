#!/bin/bash

# Create usage message
usage(){
echo "
Required parameters:
-f|--fasta                Input fasta file [${FASTA}]

Optional parameters:
-h|--help                 Shows help
-t|--threads              Number of threads [${THREADS}]
-w|--wtdbg2_dir           Path to wtdbg2 directory [${WTDBG2_DIR}]
-ws|--wtdbg2_settings     wtdbg2 settings [${WTDBG2_SETTINGS}]
"
exit
}

POSITIONAL=()

# DEFAULT SETTINGS
THREADS=1
WTDBG2_DIR=/hpc/cog_bioinf/kloosterman/tools/wtdbg2_v2.2
WTDBG2=${WTDBG2_DIR}/wtdbg2
WTPOA_CNS=${WTDBG2_DIR}/wtpoa-cns
WTDBG2_SETTINGS='-x ont -g 3g'

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

WTDBG2=${WTDBG2_DIR}/wtdbg2
WTPOA_CNS=${WTDBG2_DIR}/wtpoa-cns

PREFIX=${FASTA/.fasta/_wtdbg2}
SVID=`echo $(basename $FASTA) | cut -f 1 -d '_'`

WTDBG2_COMMAND="${WTDBG2} ${WTDBG2_SETTINGS} -i ${FASTA} -t ${THREADS} -fo ${PREFIX}"
WTPOA_CNS_COMMAND="${WTPOA_CNS} -t ${THREADS} -i ${PREFIX}.ctg.lay.gz -fo ${PREFIX}.ctg.fa"
SED_COMMAND="sed -i \"s/>ctg/>${SVID}_ctg/g\" ${PREFIX}.ctg.fa"

#echo ${WTDBG2_COMMAND}
eval ${WTDBG2_COMMAND}
#echo ${WTPOA_CNS_COMMAND}
eval ${WTPOA_CNS_COMMAND}
#echo ${SED_COMMAND}
eval ${SED_COMMAND}

#echo `date`: Done
