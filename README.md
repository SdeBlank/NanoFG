# NanoFG

##### NanoFG integrates existing genomic tools and additonal code to provide an easy-to-use pipeline for the detection of fusion genes from Nanopore sequencing data

## INSTALL
```
git clone https://github.com/SdeBlank/NanoFG.git NanoFG
cd NanoFG

virtualenv venv -p python3
. venv/bin/activate
pip install -r requirements.txt
```
Adjust all paths in the paths.ini file to installed tools

## How to run

bash NanoFG.sh -f </path/to/fastqdir>  [-n SAMPLE_NAME ] [-s SELECTION] [-cf] [-cc] [-df] [-dc]
```
OR
```
bash NanoFG.sh -b </path/to/bam> [-v </path/to/vcf>] [-n SAMPLE_NAME ] [-s SELECTION] [-cf] [-cc] [-df] [-dc]

```
For more information, see the wiki:
https://github.com/SdeBlank/NanoFG/wiki
```

## Citation
Please see and cite: https://www.biorxiv.org/content/10.1101/807545v2

#### Additionally, NanoFG integrates different existing tools:
###### Samtools (1.7) - http://samtools.sourceforge.net/
_Li, H. et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 25, 2078–2079 (2009)._
###### Minimap2 (2.6) - https://github.com/lh3/minimap2
_Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100._

###### LAST (921) - http://last.cbrc.jp/doc/last.html
_Kielbasa, S. M., Wan, R., Sato, K., Horton, P. & Frith, M. C. Adaptive seeds tame genomic sequence comparison. Genome Research 21, 487–493 (2011)._

###### NanoSV (1.2.4) - https://github.com/mroosmalen/nanosv
_Cretu Stancu, M. et al. Mapping and phasing of structural variation in patient genomes using nanopore sequencing. Nat. Commun. 8, 1326 (2017)._

###### Wtdbg2 (2.2) - https://github.com/ruanjue/wtdbg2 
_Ruan, J. and Li, H. (2019) Fast and accurate long-read assembly with wtdbg2. Nat Methods_
