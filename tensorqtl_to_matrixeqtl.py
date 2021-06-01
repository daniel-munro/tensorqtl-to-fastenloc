import argparse
from fastparquet import ParquetFile
from glob import glob
import pandas as pd

p = argparse.ArgumentParser(description="Convert gene-variant pairs from tensorQTL to Matrix eQTL output format.")
p.add_argument("parquets", help="Directory containing parquet files with all tested gene-variant pairs")
p.add_argument("outfile", help="Output file name (*.txt.gz)")
p.add_argument("--perm", required=False, help="tensorQTL top association per gene-variant pair file (*.cis_qtl.txt.gz). Used to subset data to eGenes (optional).")
args = p.parse_args()

if args.perm is not None:
    perm = pd.read_csv(args.perm, sep="\t")
    perm = perm.loc[perm["qval"] < 0.05, :]
    genes = set(perm["phenotype_id"])

df = []
for f in glob(f"{args.parquets}/*.parquet"):
    d = ParquetFile(f).to_pandas()
    if args.perm is not None:
        d = d.loc[d["phenotype_id"].isin(genes)]
    # t-statistic is coefficient divided by its standard error (https://dss.princeton.edu/online_help/analysis/interpreting_regression.htm)
    d["t-stat"] = d["slope"] / d["slope_se"]
    d = d.rename(columns={"variant_id": "SNP", "phenotype_id": "gene", "slope": "beta", "pval_nominal": "p-value"})
    d = d[["SNP", "gene", "beta", "t-stat", "p-value"]]
    df.append(d)
df = pd.concat(df)
df.to_csv(args.outfile, sep="\t", index=False)
