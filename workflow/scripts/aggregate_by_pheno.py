import pandas as pd
import os
from pathlib import Path


def main(summary_files, outfile, bypheno=True):
    # Open output file
    fo = open(outfile, "w")

    # Initialize variable to detect if it's the first time writing
    # in the output file
    firstw = True
    wmode = "w"

    # Cycle over input files
    for i, ff in enumerate(summary_files):
        # Exctract phenotype from path
        absp = Path(ff)
        pheno = absp.parts[-2]
        # Handle empty file with exception from pandas
        try:
            # Read input file
            print(f"Computing phenotype: {pheno}")
            try:
                sdf = pd.read_csv(ff, sep="\t")
            except pd.errors.EmptyDataError:
                sdf = pd.DataFrame()

            # Handling an empty data.frame
            if sdf.shape[0] > 0:
                if bypheno:
                    tmpdf = summarize_df(sdf)
                    tmpdf["pheno"] = pheno
                    tmpdf["sumfile"] = ff
                else:
                    tmpdf = sdf

                # Write results to output file
                tmpdf.to_csv(fo, sep="\t", index=False, header=firstw,
                             mode=wmode)

                # After the first writing set header to false and append mode
                firstw = False
                wmode = "a"
        except pd.errors.EmptyDataError:
            continue

    # Close file
    fo.close()


def summarize_df(df, cols=[]):
    """This function handle the summary report for credible sets
    produced by the runSusieR.R script.

    It extract the variable posterior probability and get the
    best hit for each credible set.

    Parameters
    ----------
    df : pandas.DataFrame
        contains the summary statistic filtered by credible sets computed
        with the susieR algorithm
    cols : list
        a list of characters containing a subset of the columns in the
        original sumstat to report in the output

    Returns
    -------
    dfres : pandas.DataFrame
        a summary dataframe. One line for each credible set.
    """
    gg = df.groupby(["CHROM", "cs"])
    idx = gg.agg({"variable_prob": "idxmax", "N": "size"})
    tmp = df.loc[idx["variable_prob"]]
    idx = idx.set_index("variable_prob")
    tmp["N_VARIANTS"] = idx["N"]
    if cols:
        try:
            dfres = tmp[cols]
        except KeyError:
            print("Some columns are not in the data.frame")
            dfres = tmp
    else:
        dfres = tmp
    return dfres


if __name__ == "__main__":

    main(snakemake.input, snakemake.output[0], snakemake.params.bypheno)
