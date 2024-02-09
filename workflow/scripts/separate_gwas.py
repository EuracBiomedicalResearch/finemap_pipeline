import pandas as pd
import argparse
import os


def main(inputfile, outdir, pval_thr=5e-8, pvalcol="LOG10P", pheno="",
         phenofile=""):

    # Check for directory existence, if not create it
    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

    print(f"Read gwas: {inputfile}")
    mygwas = pd.read_csv(inputfile, sep='\t', header=0,
                         usecols=["CHROM", "GENPOS", "ID", pvalcol,
                                  "BETA", "SE", "N", "CHISQ",
                                  "ALLELE0", "ALLELE1"])

    # Compute p_values from LOGP
    print("Compute pvalues")
    mygwas['P'] = 10 ** -mygwas[pvalcol]
    ix = mygwas['P'] < pval_thr
    mygwas_filt = mygwas[ix]
    chroms = mygwas_filt['CHROM'].unique()
    mygwas.set_index("CHROM", inplace=True)
    for cc in chroms:
        print(f"Running on chromosome: {cc}")
        g = mygwas.loc[cc]

        # Filter out NA values
        ix_nona = ~g[pvalcol].isna()
        g[ix_nona].to_csv(os.path.join(outdir, f'chr{cc}_sumstat.csv'),
                          sep='\t', index=True, index_label="CHROM",
                          header=True)



if __name__ == '__main__':
    main(inputfile=snakemake.input.sumstat,
         outdir=snakemake.output[0],
         pval_thr=snakemake.params.pval_thr,
         pvalcol=snakemake.params.pval_col,
         pheno=snakemake.wildcards.pheno)

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
