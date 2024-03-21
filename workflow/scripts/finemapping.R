#!/shared/bioinf/R/bin/Rscript-4.3-BioC3.17

#---- Libraries ----
library(dplyr)
library(stringr)
library(glue)
library(susieR, lib.loc = "/home/mfilosi/R/rocker-rstudio/4.3")
library(Matrix)
library(rlog)
# library(Rfast, lib.loc = "/home/mfilosi/R/rocker-rstudio/4.3")
library(data.table)

#---- Functions ----
# Not needed any more
create_ldmat <- function(x, snplist){
  #---- Create LD matrix for all SNPs
  #---- Check if it's plink1.9 or plink2 input -----
  nn <- names(x)
  if (any(grepl("ID_", nn))){
    colix <- grep("ID_", nn)
  } else if (any(grepl("SNP_", nn))){
    colix <- grep("SNP_", nn)
  } else {
    stop("Not conformable ld file")
  }
  
  # Get Column for R value
  rcol <- grep("^(UNPHASED_R|R)$", nn, perl=TRUE)
  
  # snps <- unique(c(x$ID_A, x$ID_B))
  snps <- unique(unlist(x[, colix]))
  snps <- snps[!is.na(match(snps, snplist))]
  nsnps <- length(snps)
  snps <- data.frame(ID=snps) %>% mutate(i=1:n())
  
  j1 <- j2 <- c("ID")
  names(j1) <- nn[colix[1]]
  names(j2) <- nn[colix[2]]
  
  ldmat2 <- x %>% 
    left_join(snps, by=j1) %>% 
    left_join(snps, by=j2)
  
  cormat <- matrix(0, nrow=nsnps, ncol=nsnps)
  for (i in 1:nrow(ldmat2)){
    ix <- ldmat2[i, "i.x"]
    iy <- ldmat2[i, "i.y"]
    cormat[ix, iy] <- cormat[iy, ix] <- ldmat2[i, rcol]
  }
  diag(cormat) <- 1
  
  dimnames(cormat) <- list(snps$ID, snps$ID)
  # attr(cormat, "d") <- rep(nsnps, 2)
  
  return(cormat)
}

err_handling <- function(e){
  rlog::log_error(glue("{e}"))
  rlog::log_warn(glue("Creating empty output for clump: {clid}, chrom: {mychrom}"))
  # Return an empty list and set class to "susie"
  aa <- list()
  aa$sets <- list()
  class(aa) <- "susie"
  return(aa)
  }

#---- Setup ----
# cwd <- snakemake@params["cwd"]
ld_file <- snakemake@input[["ldfile"]]
sumstat_file <- snakemake@input[["smstat"]]
pheno_file <- snakemake@input[["phenofile"]]

# Load parameters for data handling
compute_ld <- snakemake@params[["use_ld"]]
CHRIS <- snakemake@params[["chris_id"]]

# Load parameters for susieR model
susie_min_abs_cor <- snakemake@params[["min_abs_corr"]]
susie_iter <- snakemake@params[["iter"]]

cs_smstat_file <- snakemake@output[["cs_smstat"]]
cs_report_file <- snakemake@output[["cs_report"]]
cs_rds_file <- snakemake@output[["cs_rds"]]
cs_log <- snakemake@output[["cs_log"]]

# RUN ONLY LOCALLY
# proj_path <- "/scratch/mfilosi/test_gwas_regenie/metaboGWAS"
# ld_file <- glue(proj_path, "/ld/meta_6/chr4.ld")
# geno_file <- glue(proj_path, "/ld/meta_6/chr4.ld_0.raw")
# snp_file <- glue(proj_path, "/ld/meta_6/chr4.ld_0.snplist")
# clump_file <- glue(proj_path, "/clumping/meta_6/chr4.clumps")
# finemap_file <- glue(proj_path, "/results/meta_6.finemapped")
# sumstat_file <- glue(proj_path, "/tmp-smstat/meta_6/chr4_sumstat.csv")
# pheno_file <- glue(proj_path, "/tmp-smstat/meta_6/phenotype.csv") 

rlog::log_info(glue("Running finemapping on file: {ld_file}"))

ld_file_size <- file.size(ld_file)
if (ld_file_size > 0) {
  rlog::log_debug("Start preprocessing for susieR")
  geno_info <- read.table(ld_file, header = TRUE, sep = "\t")

  #---- Read phenotype file ----
  rlog::log_debug(glue("Reading phenotype file: {pheno_file}"))
  pheno <- read.table(pheno_file, header = TRUE, sep = "\t")
  y <- pheno["Y"]
  sid <- pheno["IID"]
  if (CHRIS) {
    rlog::log_debug("Padding IDs...")
    sid <- sid %>% 
      mutate(IID = str_pad(IID, side = "left", width = 10, pad = 0))
  }

  # Remove NA from sample list
  nonasamp <- which(!is.na(y))
  sid <- sid[nonasamp, ]
  y <- y[nonasamp, ]

  if (length(nonasamp) == 0){
    rlog::log_error("No samples left for analysis...")
    rlog::log_error("Please check the overlap between sample ID in
    genotype files and phenotype file.")
    quit(save="no", status=1)
  } else {
    #---- read summary stat ----
    smstat <- fread(sumstat_file, header = TRUE, sep = "\t")

    #---- Prepare the result list ----
    cs_smstat_l <- list()
    cs_report_l <- list()
    cs_rssfit_l <- list()

    for (i in 1:nrow(geno_info)){
      #---- Get clump id and chromosome ----
      clid <- geno_info[i, "CLUMPID"]
      mychrom <- geno_info[i, "CHROM"]

      #---- Logging ----
      rlog::log_info(glue("Processing clump id {clid} on chrom: {mychrom}"))

      #---- Read clump genotype file ----
      geno_file <- geno_info[i, "GENOFILE"]
      # geno_file <- file.path(proj_path, geno_file)
      genos <- fread(geno_file)
      snps <- names(genos)[7:ncol(genos)]
      snps <- gsub("_.+", "", snps)
      rlog::log_info(glue("Analyzing: {length(snps)} snps!"))

      #---- Subset summary stat ----
      smstattmp <- smstat[ID %in% snps, ]
      ix <- match(snps, smstattmp$ID)
      smstattmp <- smstattmp[ix, ]

      #---- Match phenotype with genotypes ----
      gsix <- match(genos$IID, sid)
      genos_nona <- genos[which(!is.na(gsix))]
      ytmp <- y[gsix[!is.na(gsix)]]
      sidtmp <- sid[gsix[!is.na(gsix)]]
      X <- as.matrix(genos_nona[,7:ncol(genos_nona)]) * 1.0

      rlog::log_debug("Compute susieR model...\n")
      if (compute_ld){
        betas <- smstattmp$BETA
        sebetas <- smstattmp$SE
        n <- min(smstattmp$N, na.rm=TRUE)

        rlog::log_debug("Compute correlation matrix...")
        R <- cor(X)
        rlog::log_debug("Run susie model on LD.")
        rss_ress <- tryCatch(
          susie_rss(bhat = betas, shat = sebetas, n = n, R = R, max_iter = susie_iter, min_abs_corr = susie_min_abs_cor)
        , error = err_handling
          )
      } else {
        rlog::log_debug("Run susie model on genotype.")
        rss_ress <- tryCatch(
          susie(X, y = ytmp, max_iter = susie_iter,
                min_abs_corr = susie_min_abs_cor),
          error = err_handling
        )
      }

      cs <- rss_ress$sets
      if (length(cs$cs) > 0) {
        rlog::log_info("Found some credible sets!")
        vars <- summary(rss_ress)
        allsnps <- vars$vars
        allsnps <- allsnps %>% 
          mutate(
            ID = snps[variable],
            is_95_credset=variable_prob >= 0.95,
            is_99_credset=variable_prob >= 0.99)

        smstatcs <- smstattmp %>%
          left_join(allsnps, by = "ID") %>%
          filter(cs > 0) %>%
          mutate(csid = paste(CHROM, cs, sep = "_"))
        csrep <- vars$cs %>%
          mutate(CLUMPID = clid,
                 CHROM = mychrom,
                 csid = paste(CHROM, cs, sep = "_"))

        cs_smstat_l[[i]] <- data.table(smstatcs)
        cs_report_l[[i]] <- data.table(csrep)
        cs_rssfit_l[[i]] <- rss_ress
      } else {
        rlog::log_info("No credible set found!")
      }
    }

    cs_smstat <- data.table::rbindlist(cs_smstat_l)
    cs_report <- data.table::rbindlist(cs_report_l)

    #---- Saving output ----
    rlog::log_info("Saving output...")
    write.table(cs_report, file = cs_report_file, sep = "\t",
                col.names = TRUE, row.names = FALSE, quote = FALSE)
    write.table(cs_smstat, file = cs_smstat_file, sep = "\t",
                col.names = TRUE, row.names = FALSE, quote = FALSE)
    saveRDS(cs_rssfit_l, file = cs_rds_file)
  }
} else {
  rlog::log_debug("Touch output files...")
  cat(file = cs_report_file)
  cat(file = cs_smstat_file)
  cat(file = cs_rds_file)
}


rlog::log_info(cat("\n",
    "----------", "\n",
    "Done!!!", "\n",
    "----------", "\n"))

# susie_plot(rss_res,y = "PIP")
# susie_plot(zscores,y = "z")
# 
# lambda <- estimate_s_rss(zscores, R = ldmat, n=n)





