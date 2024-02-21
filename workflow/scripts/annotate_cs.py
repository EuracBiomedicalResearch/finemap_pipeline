import pandas as pd
import os


def main(pheno, tophits_dir, cs_file):

    try:
        cs = pd.read_csv(cs_file, header=0, delimiter="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()

    tpfile = os.path.join(tophits_dir, f"{pheno}.regenie.filtered.gz")
    tpdf = pd.read_csv(tpfile, header=0, delimiter="\t",
                       usecols=["ID", "GENE_NAME"])
    cs.set_index("ID", inplace=True)
    tpdf.set_index("ID", inplace=True)

    cs_annot = cs.merge(tpdf, how="left", on="ID")

    return cs_annot


if __name__ == "__main__":

    cs_annot = main(pheno=snakemake.wildcards.pheno,
                    tophits_dir=snakemake.params.tophitsdir,
                    cs_file=snakemake.input[0])
    if cs_annot.shape[0] > 0:
        cs_annot.to_csv(snakemake.output[0], index_label="ID", sep="\t")
    else:
        f = open(snakemake.output[0], "w")
        f.close()
