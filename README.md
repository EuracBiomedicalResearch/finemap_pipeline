# Finemapping pipeline

This is meant to collect a set of script needed to run finemapping on genetic association studies.
We use `snakemake` workflow managment system based on  *python* language. For additional details, refers to the homepage 
of [`Snakemake`](https://snakemake.github.io)

In particular, we use the recently developed algorithm `susie` with its R implementation [`susieR`](https://stephenslab.github.io/susieR/index.html)

## Installation

### Requirements

```
- snakemake
- git
```

> **Optional** if already installed by the system administrator or already available in a conda environment.

See [**Install snakemake**](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for further 
information and specific parameters.

```bash
conda create -n snakemake bioconda::snakemake
```

### Pipeline installation

1. Now clone this repo into your working directory.

```bash
git clone https://github.com/EuracBiomedicalResearch/finemap_pipeline
cd finemap_pipeline
```

2. Write a configuration file

All the available parameters are defined through a configuration file written in `YAML` language.
Take the file [`config/config.yaml`](config/config.yaml) as an example and modify it according to your needs.


# References
 
> Wang, G., Sarkar, A., Carbonetto, P. & Stephens, M. (2020). A simple new approach to variable selection in regression, with application to genetic fine mapping. Journal of the Royal Statistical Society, Series B 82, 1273–1300. [https://doi.org/10.1111/rssb.12388](https://doi.org/10.1111/rssb.12388)
