import pandas as pd


def main(phenofile, pheno, outfile):
    print(outfile)
    phenoval = pd.read_csv(phenofile, sep="\t", header=0,
                           usecols=["FID", "IID", pheno])
    phenoval.rename(columns={pheno: "Y"}, inplace=True)
    phenoval.to_csv(outfile, index=False, sep="\t")


if __name__ == "__main__":
    main(phenofile=snakemake.input.phenofile,
         pheno=snakemake.wildcards.pheno,
         outfile=snakemake.output[0])
    # import argparse

    # parser = argparse.ArgumentParser()
    # parser.add_argument("--phenofile", required=True)
    # parser.add_argument("--outfile", required=True)
    # parser.add_argument("--pheno", required=True)

    # # fake_args = ["--phenofile", snakemake.input.sumstat]
    # args = parser.parse_args()

    # main(phenofile=args.phenofile, pheno=args.pheno, outfile=args.outfile)
