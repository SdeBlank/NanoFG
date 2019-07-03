# NanoFG

## *--WORK IN PROGRESS--*

## INSTALL
```
virtualenv venv -p </path/to/python>
. venv/bin/activate
pip install -r requirements.txt
```
## How to run
```
bash NanoFG.sh -b BAM [-v VCF] [-s SELECTION] [-cc] [-df] [-dc]
```

## INPUT

##### MANDATORY

**-b|--bam**

Path to bamfile

##### MAIN OPTIONAL 

**-v|--vcf**

Path to vcf

---

**-s|--selection  Regions to select from the bamfiles (separated by ',')**

Accepted formats:
- Direct region (e.g. 17:7565097-7590856 )
- Ensembl identifier (e.g. ENSG00000141510)
- Common gene name (e.g. TP53)

---

**-cc|--consensus_calling**

Creates a consensus of all supporting reads for a breakpoint before calling fusions. Increases the accuracy of breakpoint detection, which is especially important for exon-exon fusions. Only activate if there is sufficient coverage to create a consensus.

---

**-df|--dont_filter**

When activated, NanoFG does not filter breakpoints before and during its steps.

---

**-dc|--dont_clean**

When activated, NanoFG does not remove any intermediate files created during its process. Important if you want to keep the consensus sequences of the fusion gene after running NanoFG.




