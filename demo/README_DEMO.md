This demonstration of NanoFG can be run on the simulation demo set demo.fastq. 
This data constitutes simulated nanopore reads containing XX fusion spanning reads for the gene fusion XX – XX.

To run NanoFG on the demo set, enter this command in the bash command line:

```
NanoFG.sh -f </path/to/demo_fastqdir>
```

NanoFG has been tested on the following operating system:
- Operating system:       CentOS Linux release 7.7.1908 (Core)
- LSB Version:            :core-4.1-amd64:core-4.1-noarch

The expected output of NanoFG on the demo.fastq should contain:
- demo_FusionGeneInfo.txt
- demo_FusionGenes.vcf
- demo_FusionGenes.primers
- demo.FusionGenes.pdf

NanoFG should detect XX fusion spanning reads for the gene fusion XX – XX.

Estimated runtime: X minutes.
