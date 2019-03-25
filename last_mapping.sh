PATH_LASTAL=/hpc/cog_bioinf/kloosterman/tools/last-921/src/lastal
PATH_LASTSPLIT=/hpc/cog_bioinf/kloosterman/tools/last-921/src/last-split
LASTAL_GENOME='/hpc/cog_bioinf/GENOMES/LAST/human_GATK_GRCh37'
LASTAL_SETTINGS="-Q 0 -p /hpc/cog_bioinf/kloosterman/tools/last-921/last_params $LASTAL_GENOME"

FASTA=$1
MAF="${FASTA/.fasta/.maf}"

$PATH_LASTAL $LASTAL_SETTINGS $FASTA | $PATH_LASTSPLIT > $MAF

PATH_MAFCONVERT=/hpc/cog_bioinf/kloosterman/tools/last-921/scripts/maf-convert
PATH_SAMBAMBA=/hpc/local/CentOS7/cog_bioinf/sambamba_v0.6.5/sambamba
CONVERSION_SETTINGS='-S --format=bam'

BAM="${MAF/.maf/.bam}"

$PATH_MAFCONVERT -d sam $MAF | $PATH_SAMBAMBA view $CONVERSION_SETTINGS /dev/stdin > $BAM

