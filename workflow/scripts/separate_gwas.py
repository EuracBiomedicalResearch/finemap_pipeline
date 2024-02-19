import pandas as pd
import argparse
import os
import subprocess as sb
import gzip
from tempfile import NamedTemporaryFile


def extract_chr_from_smstat(inputfile, chrom=None):
    # Extract header from sm stat
    myhead = pd.read_csv(inputfile, sep='\t', header=None, nrows=1)
    myhead = myhead.iloc[0, :].tolist()
    with NamedTemporaryFile(mode="w+b") as tmpfile:
        p1 = sb.run(["tabix", f"{inputfile}", f"{chrom}"], stdout=tmpfile)
        # p2 = sb.run(["gzip", "-c"], stdin=p1.stdout, stdout=tmpfile)
        try:
            p1.check_returncode()
            df = pd.read_csv(tmpfile.name, sep="\t", header=None)
            df.columns = myhead
        except sb.CalledProcessError:
            df = pd.DataFrame()
    return df


def main(inputfile, outdir, pval_thr=5e-8, pvalcol="LOG10P", pheno="",
         phenofile="", chrom=None):

    # Columns to write for the sumstat
    cols_to_write = ["CHROM", "GENPOS", "ID", pvalcol, "BETA", "SE", "N",
                     "CHISQ", "ALLELE0", "ALLELE1"]

    print(f"Read gwas: {inputfile}")
    mygwas = extract_chr_from_smstat(inputfile, chrom=chrom)

    # Compute p_values from LOGP
    print("Compute pvalues")
    mygwas['P'] = 10 ** -mygwas[pvalcol]
    ix = mygwas['P'] < pval_thr
    mygwas_filt = mygwas[ix]

    if mygwas_filt.shape[0] > 0:
        # Filter out NA values
        ix_nona = ~mygwas[pvalcol].isna()
        mygwas.loc[ix_nona, cols_to_write].to_csv(outdir, sep='\t',
                                                  index=False, header=True,
                                                  na_rep="NA")
    else:
        # Touch file
        f = open(outdir, "w")
        f.close()


if __name__ == '__main__':
    main(inputfile=snakemake.input.sumstat,
         outdir=snakemake.output[0],
         pval_thr=snakemake.params.pval_thr,
         pvalcol=snakemake.params.pval_col,
         pheno=snakemake.wildcards.pheno,
         chrom=snakemake.wildcards.chrom)

# parser = argparse.ArgumentParser()
# parser.add_argument('-i', '--infile', default=None, type=str)
# parser.add_argument('-o', '--outdir', default=".")
# parser.add_argument('--pcol', default="LOG10P")
# parser.add_argument('--pthr', default=5e-8, type=float)
# parser.add_argument('--pheno', default='')
# args = parser.parse_args()

# if not os.path.exists(args.outdir):
#     os.makedirs(args.outdir)

# main(inputfile=args.infile,
#      outdir=args.outdir,
#      pval_thr=args.pthr,
#      pvalcol=args.pcol,
#      pheno=args.pheno)
