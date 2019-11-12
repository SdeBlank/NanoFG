# NanoFG Demo

This demonstration of NanoFG can be run on the simulation demo set demo.fastq. 

This data constitutes 3524 simulated nanopore reads containing fusion spanning reads for the gene fusion NEDD4L – CEP112.

## Dependencies
NanoFG has been tested on the following operating system:
- Operating system:       CentOS Linux release 7.7.1908 (Core)
- LSB Version:            :core-4.1-amd64:core-4.1-noarch

NanoFG requires installation of the following tools:

- Samtools (1.7) - http://samtools.sourceforge.net/
- Minimap2 (2.6) - https://github.com/lh3/minimap2 
- LAST (921) - http://last.cbrc.jp/doc/last.html
- NanoSV (1.2.4) - https://github.com/mroosmalen/nanosv

## Installation
```
git clone https://github.com/SdeBlank/NanoFG.git NanoFG
cd NanoFG

virtualenv venv -p python3
. venv/bin/activate
pip install -r requirements.txt
```
Adjust all paths in the paths.ini file to installed tools




## How to run

To run NanoFG on the demo set, enter this command in the bash command line:

```
NanoFG.sh -f NanoFG/demo/demo_fastqdir
```

## Expected results

The expected output of NanoFG on the demo.fastq should contain:
- demo_FusionGenesInfo.txt
- demo_FusionGenes.vcf
- demo_FusionGenes.primers
- demo_FusionGenes.pdf

NanoFG should detect 17 fusion spanning reads for the gene fusion NEDD4L – CEP112.

Estimated runtime: 5 minutes.
