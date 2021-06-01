# The tissue name will be included in the output:
TISSUE=NAcc
# All files ending in .parquet in this directory will be used:
PARQUET_DIR=example/NAcc/
# This VCF file should contain all variants in the tensorQTL results:
VCF_INFILE=example/geno.vcf.gz
TORUS_PATH=torus/src/torus
# This will be created as an intermediate file:
MATEQTL_FILE=example/mateqtl.$TISSUE.txt.gz
# This will be created as an intermediate file:
PIP_FILE=example/pip.$TISSUE.txt
# Omit the .gz extension, though actual output file will end in .gz:
VCF_OUTFILE=example/fastenloc.eqtl.anno.$TISSUE.vcf

set -e
set -o pipefail

echo 'Converting tensorQTL to Matrix eQTL format...'
python3 tensorqtl_to_matrixeqtl.py $PARQUET_DIR $MATEQTL_FILE
echo 'Running TORUS to get PIP...'
$TORUS_PATH -d $MATEQTL_FILE -dump_pip $PIP_FILE
echo 'Creating VCF for fastENLOC...'
python3 torus_to_fastenloc.py $PIP_FILE $VCF_INFILE $TISSUE $VCF_OUTFILE
echo 'Compressing the VCF file...'
bgzip $VCF_OUTFILE
