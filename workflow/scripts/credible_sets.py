import pandas as pd
import numpy as np
from scipy.stats import norm
import gcta


def abf_scipy(beta, se, priorsd=1, log=False):
    """Alternative implementation of the ABF
    """
    beta = float(beta)
    se = float(se)
    se2 = se**2
    prior2 = priorsd**2
    abf = norm.logpdf(beta, loc=0, scale=np.sqrt(prior2 + se2)) - \
    norm.logpdf(beta, loc=0, scale=se)
    if log:
        return abf
    else:
        return np.exp(abf)

def abf(beta, se, priorsd=1, log=False):
    """This is much faster then the abf function using `logpdf` function
    from scipy.stats.

    Parameters
    ----------
    beta: float
        beta estimation from a GWAS summary statistic
    se: float
        standard error for the beta estimations
    priorsd: float
        prior standard deviation of the phenotype
    logs: bool
        return the ABF factor (default) or the logABF if `log=True`

    Returns
    -------
    float
        depending on the `log` argument, return the ABF or the logABF

    """
    beta = float(beta)
    se = float(se)
    priorsd = float(priorsd)

    abf = abf_core(z=beta / se, se=se,  priorsd=priorsd, log=log)

    return abf

def abf_core(z, se, priorsd=1, log=False):
    zscore2 = z**2
    se2 = se**2
    prior2 = priorsd**2

    t1 = np.sqrt(se2 / (se2 + prior2))
    t2 = (zscore2 / 2) * (prior2 / (se2 + prior2))
    if log:
        abf = np.log(t1) + t2
    else:
        abf = t1 * np.exp(t2)
    return abf

def abf_from_pval(pval, maf, n, priorsd=1, **kwargs):
    zscore = norm.ppf(pval / 2.0)
    se = np.sqrt(get_var(maf, n, **kwargs))
    try:
        log = kwargs.pop("log")
    except KeyError:
        log = False
    abf = abf_core(zscore, se=se, priorsd=priorsd, log=log)
    return abf

def get_var(maf, N, prop_cases=None, **kwargs):
    """Estimate variance of the phenotype based on the MAF.
    Computation depends on the type of GWAS. Use this function
    when `beta` and `se` are not available.
    In particular if ABF is computed for case-control study it's better to
    use pvalue and and standard error estimated on the MAF.

    Parameters
    ----------
    maf: float
        minor allele frequency for the specific variant
    N: float, int
        total number of partecipants in the study
    prop_cases: None, float, int
        if provided the computation expect to be a case-control study

    """
    myden = 2 * N * maf * (1-maf)
    if prop_cases:
        myden *= prop_cases * (1 - prop_cases)
    return 1 / myden

def estimate_prior(sdY, type="quant", prop_cases=None):
    if type == "bin":
        return sdY * 0.2
    else:
        return sdY * 0.15

def log_sum(l):
    """
    This function comes from:
        https://github.com/opentargets/genetics-finemapping/blob/master/finemapping/credible_set.py#L183
    Calculates the log of the sum of the exponentiated logs taking out the
        max, i.e. insuring that the sum is not Inf
    Args:
        l (pandas Series)
    Returns:
        Sum of exponentiated logs
    """
    l_max = l.max()
    l_logsum = l_max + np.log(np.sum(np.exp(l - l_max)))
    return l_logsum

def normalize_abf(abf, log=False):
    """Inspiration for this function comes from the `R gtx` package:
        https://github.com/tobyjohnson/gtx/blob/master/R/abf.R
    """
    if log:
        x = abf - abf.max(numeric_only=True)
        x = np.exp(x)
    else:
        x = abf / abf.max(numeric_only=True)
    return x / np.sum(x)

def compute_credible_set(sumstat_window, prior_sd=1, cs_prob=0.95, **kwargs):
    mydata = sumstat_window
    mydata['logABF'] = mydata.apply(lambda x: abf(beta=x["BETA_COND"],
                                                     se=x["SE_COND"],
                                                     priorsd=prior_sd,
                                                     log=True),
                                    axis=1)
    mydata['normlogABF'] = normalize_abf(mydata['logABF'], log=True)
    mydata = mydata.sort_values(by="normlogABF", ascending=False)
    mydata['postProb'] = mydata['normlogABF'].cumsum()
    mydata['is_95_cred'] = mydata['postProb'].transform(lambda x: x>=0.95)
    mydata['is_99_cred'] = mydata['postProb'].transform(lambda x: x>=0.99)
    return mydata

def get_key(key, argdict, default=None):
    """Get argument from kwargs
    """
    try:
        myarg = argdict[key]
    except KeyError:
        myarg = default
    return myarg

def credible_set(index_var, sumstat, toploci, plinkfile, prior_sd=1,
                 cs_prob=0.95, **kwargs):
    """Function to compute credible set for a set of loci.

    Parameters
    ----------
    sumstat: pd.DataFrame
        a DataFrame containing with the summary statistic for the GWAS
        analysis. It has to contain conditional BETA, SE and PVAL_COND.
    toploci: pd.DataFrame
        dataFrame containing the list of toploci for the summary statistic
    prior_sd: float
        the prior standard deviation of the phenotype (default=1.0).
    cs_prob: float
        posterior probability for the credible set (default=0.95).
    **kwargs:
        a dictionary of named arguments. Expected arguments that will be used
        `cojo_wind`: region window span in kb
        `sdY`: variance of the phenotype for the specific GWAS
        `stype`: type of GWAS run, 'quant' or 'binary'
        other arguments will be ignored.

    Return
    ------
    pandas.DataFrame
        containing all the credible sets for the input toploci list
    """
    # Get parameters
    cojo_wind = get_key('cojo_wind', kwargs, default=500)

    # Initialize conditional analysis
    cc = gcta.ConditionalAnalysis(**kwargs)
    # for vv in toploci.to_dict(orient='records'):
    # Cycle over the index variants (if any)...
    print(f"computing: {index_var}")
    pos = index_var['GENPOS']

    # Get region around the index snp
    # -------------------------------
    myreg = ((sumstat["CHROM"] == index_var["CHROM"]) &
            (sumstat["GENPOS"] >= pos - (cojo_wind * 1000)) &
            (sumstat["GENPOS"] <= pos + (cojo_wind * 1000))
            )
    sumstat_wind = sumstat.loc[myreg, :]

    # Adjust p-value by condition beta and pvalue
    # -------------------------------------------
    cond_list = gcta.get_variant_to_cond(index_var=index_var["ID"],
                                         window_var=sumstat_wind['ID'],
                                         top_loci=toploci['ID'])
    df = cc.adjust_by_condition(sumstat_wind, plinkfile,
                                cond_list=cond_list)

    # Run estimation of credible set
    # ------------------------------
    stype = get_key('stype', kwargs, default="quant")
    priorsd = estimate_prior(prior_sd, type=stype)
    dfcs = compute_credible_set(df, prior_sd=priorsd)
    print(f"Conditioning snps: {len(cond_list)}")
    print(f"Credible sets at 95%: {dfcs['is_95_cred'].sum()}")
    print(f"Credible sets at 99%: {dfcs['is_99_cred'].sum()}")
    print(f"Tot. SNPs: {dfcs.shape[0]}")

    return dfcs


def susier_wrapper(index_var, sumstat, toploci, plinkfile, **kwargs):

    pth = os.path.basename(__file__)
