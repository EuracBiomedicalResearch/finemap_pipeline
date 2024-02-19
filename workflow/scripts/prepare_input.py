import pandas as pd
import numpy as np
import os
import glob


def prepare_inputs(tophits_dir, pval_thr=5e-8, outfile=None):
    try:
        flist = glob.glob(os.path.join(tophits_dir, "*.regenie.filtered.gz"))
    except FileNotFoundError as e:
        print(f"Cannot find tophits directory: {tophits_dir}", e)

    reslist = []
    for i, f in enumerate(flist):
        pheno = os.path.basename(f)
        pheno = pheno.replace(".regenie.filtered.gz", "")

        print(f"Running pheno: {pheno}, {i} / {len(flist)}", end="\r")
        tpdf = pd.read_csv(f, header=0, sep="\t", usecols=["CHROM", "LOG10P"])

        if tpdf.shape[0] > 0:
            tpdf["P"] = 10 ** -tpdf["LOG10P"]

            ix = tpdf['P'] < pval_thr
            tpdf_filt = tpdf[ix]
            if tpdf_filt.shape[0] > 0:
                chroms = tpdf_filt['CHROM'].unique()
                phenos = [pheno for j in range(len(chroms))]
                mydf = pd.DataFrame.from_dict({
                    "chrom": chroms,
                    "pheno": phenos
                }
                )
                reslist.append(mydf)

    if reslist:
        resdf = pd.concat(reslist, ignore_index=True)
    else:
        resdf = pd.DataFrame()

    if outfile:
        print(f"\nWriting file: {outfile}")
        resdf.to_csv(outfile, sep="\t", index=False, header=True)

    return resdf


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("--topdir", required=True)
    parser.add_argument("--pval_thr", default=5e-8)
    parser.add_argument("--outfile", default=None)

    args = parser.parse_args()

    # resdir = "/scratch/mfilosi/pQTL_somalogic/results_somalogic_HRC13K_resid/results"
    # tphits = os.path.join(resdir, "tophits")

    # pval_thr = 6.867189e-12

    prepare_inputs(args.tophits_dir, pval_thr=args.pval_thr, outfile=args.outfile)
