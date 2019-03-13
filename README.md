# NanoFG
Third-generation fusion gene detection

Types of fusion gene detection to implement:

Exon-exon *Done*\n
Exon-intron *Done*
Intron-intron *Done*
3' UTR fusion
5'UTR/Promoter fusion


Issues:
-Possible fusion if BND outside a gene and other inside exons due to skipping of the first or last exon. Might depend on the distance of the BND from the gene.
-BND inside of start codon or stop codon
-BND only few bp from start or end codon gives a protein fusion, but maybe more like 5' or 3' UTR fusion
