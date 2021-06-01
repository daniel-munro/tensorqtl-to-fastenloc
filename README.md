# tensorqtl-to-fastenloc

Convert [tensorQTL](https://github.com/broadinstitute/tensorqtl) output to eQTL annotation input for [fastENLOC](https://github.com/xqwen/fastenloc). This takes three steps:

1. A Python script converts tensorQTL `cis_nominal` mode output to Matrix eQTL output format.
2. TORUS produces posterior inclusion probabilities (PIP).
3. Another Python script prepares the results in the VCF format used by fastENLOC.

## Requirements

- [pandas](https://pypi.org/project/pandas/)
- [fastparquet](https://pypi.org/project/fastparquet/)
- [TORUS](https://github.com/xqwen/torus)
- bgzip

## Usage

After cloning this repo and installing TORUS, edit the first lines of [run.sh](run.sh) with your desired paths and other parameters. Then run `bash run.sh` on the command line from this directory.
