# NanoFG

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

bash NanoFG.sh -f </path/to/fastq>  [-n SAMPLE_NAME ] [-s SELECTION] [-cf] [-cc] [-df] [-dc]
```
OR
```
bash NanoFG.sh -b </path/to/bam> [-v </path/to/vcf>] [-n SAMPLE_NAME ] [-s SELECTION] [-cf] [-cc] [-df] [-dc]

```
For more information, see the wiki:
https://github.com/SdeBlank/NanoFG/wiki
