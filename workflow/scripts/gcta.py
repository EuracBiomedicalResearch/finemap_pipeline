import subprocess as sb
import numpy as np
import json
import os
from tempfile import mkdtemp
import logging
import shutil
import pandas as pd


class ConditionalAnalysis(object):
    def __init__(self, **kwargs):
        self.default_params = {'maf': {'flag': '--maf', 'value': '0.01'},
                               'cojo_collinear': {
                               'flag': '--cojo-collinear',
                               'value': '0.9'},
                               'cojo_p': {
                               'flag': '--cojo-p',
                               'value': str(5e-8)},
                               'cojo_wind': {
                               'flag': '--cojo-wind',
                               'value': str(10000)},
                               'extract': {
                               'flag': '--extract',
                               'value': None},
                               'cojo_top_SNPs': {
                               'flag': '--cojo-top-SNPs',
                               'value': 10}
                               }
        # Set arguments at init of the class
        self._gctaparams = {}
        # self.set_args(**kwargs)
        try:
            self.cmdbin = kwargs.pop("gctabin")
        except KeyError:
            self.cmdbin = 'gcta64'

        try:
            tmpdir = kwargs.pop("tmpdir")
        except KeyError:
            tmpdir = None

        try:
            os.path.exists(tmpdir)
            if os.path.exists(tmpdir):
                self.tmpdir = tmpdir
            else:
                self.tmpdir = mkdtemp(dir=tmpdir)
        except TypeError:
            self.tmpdir = mkdtemp()

        self.smstat = os.path.join(self.tmpdir, "gcta_sumstat.csv")
        self.snplist = os.path.join(self.tmpdir, "gcta_snplist.csv")

        logging.info(f"Create temporary directory to store files {self.tmpdir}")


    def set_args(self, **kwargs):
        for p, o in self.default_params.items():
            if p in kwargs:
                if p == 'extract':
                    if isinstance(kwargs[p], str):
                        # attlist += [o, kwargs[p]]
                        self._gctaparams[o['flag']] = kwargs[p]
                    elif kwargs[p] is None and self.snplist:
                        # attlist += [o, self.snplist]
                        self._gctaparams[o['flag']] = self.snplist
                else:
                    # attlist += [o, str(kwargs[p])]
                    self._gctaparams[o['flag']] = kwargs[p]

        if 'cojo_cond' in kwargs:
            if not isinstance(kwargs['cojo_cond'], str):
                raise NotImplementedError('Please provide a file name string. '
                                          'The file should contain a list of SNPs'
                                          'to condition on.')

            self._gctaparams['--cojo-cond'] = kwargs['cojo_cond']

            # Remove uncompatible arguments (if present)
            for k in ['--cojo-slct', '--cojo-top-SNPs']:
                try:
                    pp = self._gctaparams.pop(k)
                except KeyError:
                    pass
        else:
            self._gctaparams['--cojo-slct'] = None
            try:
                pp = self._gctaparams.pop('--cojo-cond')
            except KeyError:
                pass

    def set_default(self):
        """Set default arguments for GCTA algorithm
        """
        logging.info("Set defaul parameter for conditional analysis")
        paramdict = {k: v['value'] for k, v in self.default_params.items()}
        self.set_args(**paramdict)

    def regenie_to_gcta(self, smstat, fout=None, snplistfile=None,
                        lpval_thr=None):
        """ Read summary statistic from regenie

        This function read the sumstats in *.regenie.gz and create a summary
        stat format to run GCTA conditional analysis.

        Parameters
        ----------
        lpval_thr: float
            pvalue threshold to search for loci to report when translating from
             regenie to gcta

        Returns
        -------
        0
        """
        if not fout:
            fout = self.smstat

        if not snplistfile:
            snplistfile = self.snplist

        sumstat = smstat.copy()
        if lpval_thr and lpval_thr > 0:
            sumstat = sumstat.loc[sumstat['LOG10P'] > lpval_thr, :]

        sumstat['P'] = 10 ** (-sumstat['LOG10P'])

        sumstat = sumstat.rename(
            columns={
                'ID': 'SNP',
                'ALLELE1': 'A1',
                'ALLELE0': 'A2',
                'A1FREQ': 'freq',
                'BETA': 'b',
                'SE': 'se',
                'P': 'p'
            }
        )
        sumstat = sumstat.loc[:, ['SNP', 'A1', 'A2', 'freq', 'b',
                                  'se', 'p', 'N']]

        # Write Sumstat
        sumstat.to_csv(fout, sep='\t', index=None)

        # Write SNP list
        sumstat.SNP.to_csv(snplistfile, index=None, header=False)

        return 0

    def get_top_loci(self, sumstat, plinkfile, cojo_wind=500,
                     cojo_collinear=0.9, cojo_p=5e-8, chrom=None, **kwargs):

        logging.info("Search toploci using conditional analysis")

        # Empty parameters for conditional analysis
        self.reset_params()
        # Support providing sumstat as string or pandas dataframe
        # Check whether there are entries otherwise exit
        smstat = ck_smstat(sumstat)
        if smstat.shape[0] == 0:
            logging.warning("Empty summary stat, exit...")
            return 0

        # Set parameters for top loci
        paramdict = {'cojo_collinear': cojo_collinear,
                     'cojo_wind': cojo_wind,
                     'cojo_p': cojo_p}

        self.set_args(**paramdict)
        self._gctaparams['--bfile'] = plinkfile

        try:
            outfile = kwargs.pop('outfile')
        except KeyError:
            outfile = None
            pass

        if not outfile:
            outfile = os.path.join(self.tmpdir, "gcta_toploci_out")

        retcode = self.regenie_to_gcta(smstat, lpval_thr=-np.log10(cojo_p))

        # Set parameter for top loci identification
        self._gctaparams['--cojo-file'] = self.smstat
        self._gctaparams['--extract'] = self.snplist
        self._gctaparams['--out'] = outfile

        if chrom:
            self._gctaparams["--chr"] = chrom

        # Run the command
        self.run_cmd(**kwargs)
        try:
            jmadf = self.read_file(outfile, ftype="jma")

            # Get toploci SNP
            sel_snps = jmadf["SNP"]
            toploci = smstat.loc[smstat['ID'].isin(sel_snps), :]
        except FileNotFoundError:
            toploci = pd.DataFrame()

        # Remove temporary file
        logging.debug("Removing temporary file")
        os.remove(self.smstat)
        os.remove(self.snplist)

        return toploci

    def adjust_by_condition(self, sumstat_wind, plinkfile, cond_list=[],
                            chrom=None, **kwargs):
        """
        Parameters
        ----------
        sumstat_wind: pandas.DataFrame
            the summary statistic of a window around an index variant

        plinkfile: str
            a string specifying the plink file to use as LD reference
        cond_list: list
            a list of variants ID to condition on
        """
        try:
            outfile = kwargs.pop('outfile')
        except KeyError:
            outfile = None
            pass

        logging.info("Running adjusting")
        # Check summary statistic, if has to be read or is alread a pandasDF
        sumstat_wind = ck_smstat(sumstat_wind)

        # If no snps to condition on don't do nothing andjust create a correct
        # df
        if cond_list is None or len(cond_list) == 0:
            sm_cond = sumstat_wind.copy()
            sm_cond["BETA_COND"] = sm_cond['BETA']
            # sm_cond["PVAL_COND"] = 10**(-sm_cond["LOG10P"])
            sm_cond["PVAL_COND"] = sm_cond["P"]
            sm_cond["SE_COND"] = sm_cond["SE"]
        else:
            self.reset_params()

            # Write condition list
            gcta_cond = os.path.join(self.tmpdir, "gcta_cond_list.csv")
            with open(gcta_cond, "w") as gc:
                for s in cond_list:
                    gc.write(s + "\n")

            # Handle outfile
            if not outfile:
                outfile = os.path.join(self.tmpdir, "gcta_adjust")

            # Write summary stat window to file
            retcode = self.regenie_to_gcta(sumstat_wind)
            params = {"extract": self.snplist,
                      "cojo_cond": gcta_cond}

            self.set_args(**params)
            # self._gctaparams['--extract'] = self.snplist
            # self._gctaparams['--cojo-cond'] = gcta_cond
            self._gctaparams['--bfile'] = plinkfile
            self._gctaparams['--cojo-file'] = self.smstat
            self._gctaparams['--out'] = outfile

            # Run gcta command
            self.run_cmd(**kwargs)

            try:
                cmadf = self.read_file(outfile, ftype="cma")
                cmadf = cmadf.rename(columns={'SNP': 'ID'})
                sm_cond = pd.merge(sumstat_wind,
                                   cmadf.loc[:, ["ID", "bC", "bC_se", "pC"]])
                sm_cond = sm_cond.rename(columns={
                    'bC': 'BETA_COND', 'pC': 'PVAL_COND', 'bC_se': 'SE_COND'
                })
            except FileNotFoundError:
                print("Cannot find the file cma.cojo")
                sm_cond = sumstat_wind.copy()
                sm_cond["BETA_COND"] = sm_cond['BETA']
                sm_cond["PVAL_COND"] = sm_cond["P"]
                sm_cond["SE_COND"] = sm_cond["SE"]

        return sm_cond

    def reset_params(self):
        for k in self._gctaparams.keys():
            self._gctaparams.pop(k)

    def read_file(self, resfile, ftype="jma"):
        try:
            df = pd.read_csv(f'{resfile}.{ftype}.cojo', header=0, sep='\t')
        except FileNotFoundError as e:
            print("Problem with the file...")
            raise e
        return df

    def run_cmd(self, dry_run=False):
        cmd = self.create_cmd()

        if dry_run:
            print("---- This is a dry run ----")
            print(cmd)
        else:
            cp = sb.run(cmd)

            if cp.returncode == 0:
                print("Ok")
            else:
                print("Error in running gcta")
                print(cp.stderr)

    def create_cmd(self):
        """ This function creates a list of strings containing the command to
        run. We need this in order to standardize the process, since some
        datasets have multiple chromosomes and the conditional analysis should
        run on all the chromosomes.

        Returns
        -------
        none:
            it does not return nothing, just set the internal attribute
            `cmdlist` used to run the conditional analysis.
        """

        # Create command
        mycmd = [self.cmdbin]
        for k, v in self._gctaparams.items():
            if v:
                mycmd.extend([k, str(v)])
            else:
                mycmd.append(k)

        return mycmd

    def cleanup(self):
        shutil.rmtree(self.tmpdir)


def ck_smstat(sumstat):
    if isinstance(sumstat, str):
        try:
            smdf = pd.read_csv(sumstat, header=0, sep='\t')
        except FileNotFoundError as e:
            print(f"Could not find file {sumstat}")
            smdf = pd.DataFrame()
        return smdf
    elif isinstance(sumstat, pd.DataFrame):
        return sumstat


def get_variant_to_cond(index_var, window_var, top_loci):
    winvar = set(window_var)
    toploci = set(top_loci)
    condvar = winvar.intersection(toploci)
    condvar = list(condvar - set([index_var]))
    return condvar
