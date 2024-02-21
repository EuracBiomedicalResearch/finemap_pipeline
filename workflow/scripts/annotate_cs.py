import pandas as pd


def main(pheno, tophits_dir, cs_file):

    cs = pd.read_csv(cs_file, header=0, delimiter="\t")
    tpfile = os.path(tophits_dir, f"{pheno}.regenie.filtered.gz")
    tpdf = pd.read_csv(tpfile, header=0, delimiter="\t",
                       usecols=["ID", "GENE_NAME"])
    cs.set_index("ID", inplace=True)
    tpdf.set_index("ID", inplace=True)
    cs_annot = pd.concat([cs, tpdf], ignore_index=False, axis=0)

    return cs_annot
if __name__ == "__main__":

    cs_annot = main(pheno=snakemake.wildcards.pheno,
                    tophits_dir=snakemake.params.tophitsdir,
                    cs_file=snakemake.input)
    if cs_annot.shape[0]:
        cs_annot.to_csv(snakemake.output, index_label="ID", sep="\t")
    else:
        f = open(snakemake.output, "w")
        f.close()