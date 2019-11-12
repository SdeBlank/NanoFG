This demonstration of NanoFG can be run on the simulation demo set demo.fastq. 
This data constitutes 3524 simulated nanopore reads containing 17 fusion spanning reads for the gene fusion NEDD4L – CEP112.

To run NanoFG on the demo set, enter this command in the bash command line:

```
NanoFG.sh -f NanoFG/demo/demo_fastqdir
```

NanoFG has been tested on the following operating system:
- Operating system:       CentOS Linux release 7.7.1908 (Core)
- LSB Version:            :core-4.1-amd64:core-4.1-noarch

The expected output of NanoFG on the demo.fastq should contain:
- demo_FusionGenesInfo.txt
- demo_FusionGenes.vcf
- demo_FusionGenes.primers
- demo_FusionGenes.pdf

NanoFG should detect 17 fusion spanning reads for the gene fusion NEDD4L – CEP112.

Estimated runtime: 5 minutes.
