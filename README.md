# NanoFG

## *--WORK IN PROGRESS--*
Third-generation fusion gene detection

STEP 1  DETECTION OF CANDIDATE FUSION GENES

Types of fusion gene detection to implement:

	Exon-exon *Done*

	Exon-intron *Done*

	Intron-intron *Done*

	3' UTR fusion

	5'UTR/Promoter fusion

	Issues:

	-Possible fusion if BND outside a gene and other inside exons due to skipping of the first or last exon. Might depend on the distance of the BND from the gene.

	-BND inside of start codon or stop codon

	-BND only few bp from start or end codon gives a protein fusion, but maybe more like 5' or 3' UTR fusion

STEP2  SCORING

Score on for example percentage of identical positions of confidence interval etc.
