import os
import csv
import pandas as pd
import numpy as np
import subprocess as sb
from tempfile import NamedTemporaryFile


class Clump(object):
    def __init__(self, lead="", proxies=[], begin=0, end=0, **kwargs):
        self.lead = lead
        self.proxies = proxies

        if begin > end:
            self.begin = end
            self.end = begin
            self.strand = "-"
        else:
            self.begin = begin
            self.end = end
            self.strand = "+"

        self._span = self.end - self.begin

        try:
            pval = kwargs["pvalue"]
        except KeyError:
            pval = None
        self.pvalue = pval

    @classmethod
    def from_clumpfileline(cls, row):
        lead = row[2]
        proxies = get_snps_fromclump(row)
        try:
            pvalue = float(row[3])
        except ValueError:
            pvalue = 1

        return cls(lead=lead, proxies=proxies, pvalue=pvalue)

    @property
    def pvalue(self):
        return self._pvalue

    @property
    def lpvalue(self):
        return -np.log10(self._pvalue)

    @pvalue.setter
    def pvalue(self, new):
        try:
            if new >= 0 and new <= 1:
                self._pvalue = new
            else:
                raise ValueError
        except ValueError:
            print(f"Incorrect range for pvalue {new}")

    @property
    def span(self):
        return self._span


    def enlarge_clump(self, totsize=1000000):

        if self.span == 0:
            cllimits = self.get_span_from_clump()
            clspan = cllimits[1] - cllimits[0]
        else:
            cllimits = self.get_interval()
            clspan = self.span

        if clspan < totsize:
            newlimits = widen_span(cllimits, tot=totsize)
            # newspan = np.array(cllimits)
            # self.span = newspan
        else:
            newlimits = cllimits
        self.set_interval(newlimits)

        return self

    def get_span_from_clump(self):
        try:
            pos = [ll.split(":")[1] for ll in self.proxies + [self.lead]
                   if ll != "."]
            pos = [int(p) for p in pos]
            # print(pos)
            # pos.append(self.begin)
            # print(pos)
        except IndexError:
            # The list of proxies is empty
            print("No proxies provided...")
            pos = [0, 0]
        except ValueError:
            # The ID are not formatted as expected (chr:pos)
            print("Cannot retrieve position from proxies")
            pos = [0, 0]

        return [min(pos), max(pos)]

    def to_list(self):
        return [self.lead] + self.proxies

    def get_interval(self):
        return [self.begin, self.end]

    def set_interval(self, new):
        if isinstance(new, list):
            tmp = np.array(new)
            tmp.sort()
            try:
                tmp.astype("int64")
            except ValueError as e:
                print(f"Please provide a list of numerical values to set \
                      span: \n {e}")

        if isinstance(new, np.ndarray):
            tmp = new.copy()
            tmp.sort()

        if isinstance(new, dict):
            tmp = self._span
            try:
                tmp[0] = np.int64(new['begin'])
                tmp[1] = np.int64(new['end'])
            except KeyError:
                print("Please provide a correctly formed dictionary with \
                      keys 'begin' and 'end'")
            except ValueError:
                print("Some values in the dictionary are not numeric")

        self._span = tmp[-1] - tmp[0]
        self.begin = tmp[0]
        self.end = tmp[-1]

    def __add__(self, other):
        if isinstance(other, Clump):
            px = self.proxies + other.proxies
            ix = np.argsort(np.array([self.pvalue, other.pvalue]))

            if ix[0] == 0:
                ll = self.lead
                px += [other.lead]
                pval = self.pvalue
            else:
                ll = other.lead
                px += [self.lead]
                pval = other.pvalue

            # Initialize new clump
            res = Clump(lead=ll, proxies=px, pvalue=pval)

            # Update interval
            spanarra = np.array([self.get_interval(), other.get_interval()])
            st = spanarra.min(axis=0)[0]
            en = spanarra.max(axis=0)[-1]
            res.set_interval([st, en])

            return res
        else:
            raise NotImplementedError(f"Cannot add Clump class to \
                                      type {other.__class__}")


def get_snps_fromclump(line):
    if line[-1] == 'NONE':
        snps = []
    else:
        snps = line[-1].replace('(1)', '')
        snps = snps.split(',')

    return snps


def read_clump_file(file, enlarge=False, merge=False, **kwargs):
    """
    Parameters
    ----------
    file: str, FileBuffer
        read a clump file with extension `.clumped` resulted from
        plink clumping and merging over chromosomes
    enlarge: bool
        if set to True try to enlarge the the clumping size to a minimum of
        `totsize`. Parameter `totsize` can be specified through `kwargs`
    merge: bool
        if set to True search in the clumplist for overlapping regions and
        merge them if it finds any

    kwargs: Additional arguments can be passed to the function. For now only
    `totsize` a numerical value specifing the minimum span of a clumping region
    in bp (i.e. 1Mb: 1e6, 1kb: 1000). Default = 1Mb (1e6).

    Returns
    -------
    list of Clump instances
        a list of Clump() instances containing the clumping
        provided by the plink software
    """
    mylist = []

    # Handling empty files and do nothing...

    with open(file, mode='r') as fc:
        fcr = csv.reader(fc, delimiter='\t', skipinitialspace=True)
        next(fcr)

        for ll in fcr:
            if len(ll) > 0:
                # snps = get_snps_fromclump(ll)
                # cc = Clump(
                #     lead=ll[2],
                #     proxies=snps,
                #     pvalue=float(ll[3])
                # )
                cc = Clump.from_clumpfileline(ll)
                mylist.append(cc)

    # Enlarge clump to a minimum size of `totsize`
    # TODO: check if totsize is > 0...and other checks on totsize

    if enlarge:
        try:
            totsize = kwargs['totsize']
        except KeyError:
            totsize = 1e6
        # mylist = [enlarge_clump(cl, totsize=totsize) for cl in mylist]
        mylist = [cl.enlarge_clump(totsize=totsize) for cl in mylist]

    # Merge overlapping clumps
    if merge:
        mylist = sort_and_merge(mylist)

    return mylist


def compute_ld(snps, clid, plinkfile, dryrun=False, prefix="ld_clump",
               **kwargs):
    # Handle kwargs
    try:
        mem = kwargs["memory"]
        mem = int(mem)
    except KeyError:
        mem = 8192

    # plinkfile = os.path.join(genopath, f"chr{mychr}.forgwas.nofid")

    # Define file names
    ofile = f'{prefix}_{clid}'
    genofile = ofile + ".raw"
    snpfile = ofile + ".snplist"

    # Write snp list
    ff = open(snpfile, "w")
    # with NamedTemporaryFile(mode='w') as ftmp:
    for s in snps:
        ff.write(s)
        ff.write("\n")
        # ftmp.write(s)
        # ftmp.write("\n")
    ff.close()

    # DO NOT COMPUTE LD matrix with plink, sometimes it's not semi-positive
    # defined and can contain negative eigenvalues.
    # Thus extract the genotype matrix and compute correlation in python/R
    # See Issue #91: https://github.com/stephenslab/susieR/issues/91
    # ldcommand = ['plink', '--bfile', plinkfile, '--keep-allele-order',
    #                 # '--ld-window-r2', '0.1',
    #                 '--r', 'square',   '--extract', #    '--ld-snp-list'
    #                 snpfile, '--out', ofile, '--memory', str(mem)]
    ldcommand = ['plink', '--bfile', plinkfile, '--keep-allele-order',
                 '--extract', snpfile,
                 '--recode', 'A', '--out', ofile, '--memory', str(mem)]
    print(" ".join(ldcommand))

    if dryrun:
        print('---- This is a dry run ----')
        print('---- Running the following command...----')
        print(' '.join(ldcommand))
    else:
        if len(snps) > 1:
            myproc = sb.run(ldcommand)

            if myproc.returncode > 0:
                print("OOOpppps some issue with this process")

    return genofile, snpfile


def process_ld_file(ldfile):
    # Read the ld file and adjust column names
    # df = pd.read_csv(ldfile, sep="\t", header=0)
    print(ldfile)
    df = pd.read_csv(ldfile, sep="\s+", header=0)
    cc = df.columns
    df.columns = [c.replace("#", "") for c in cc]

    return df


def widen_span(span, tot):
    span2 = np.array(span)
    span2.sort()
    spandist = np.diff(span2)[0]
    totdiff = tot - spandist

    if totdiff > 0:
        span2[0] -= totdiff / 2.0
        span2[1] += totdiff / 2.0

    return span2


def intersect(cl1, cl2):
    cl1.sort()
    cl2.sort()

    # Handle possible cases

    if (cl1[0] > cl2[1]) or (cl1[1] < cl2[0]):
        # cl1 and cl2 are dop not overlap
        #  |-------|              cl1
        #             |-------|   cl2
        ck = False
    else:
        if (cl1[1] <= cl2[1]) and (cl1[1] >= cl2[0]):
            # cl1 end is within cl2 range
            #     |-------|        cl1
            #        |-------|     cl2
            ck = True

        if (cl1[0] >= cl2[0]) and (cl1[0] < cl2[1]):
            # cl1 end is within cl2 range
            #     |-------|        cl1
            #  |-------|           cl2
            ck = True

    return ck


def main_tmp(clumpfile, gendata, outfile="", **kwargs):
    # Get info on genetic data
    # gg = GenDataList.from_json("../../metaboGWAS/genetic_data.json")
    # gg = GenDataList.from_json(gendata)
    # ggg = gg.gendata["HRC13K"]

    genopath = {os.path.dirname(f) for f in gendata.get_plinkfiles()}
    genopath = genopath.pop()

    # Read clump file
    cldf = read_cl_file(clumpfile)
    clumplist = read_clump_file(clumpfile)
    dflist = []

    # Run LD computation for each clumping region

    for clump in cldf.iterrows():
        ldfile = compute_ld(clump, genopath=genopath,
                            prefix=outfile, **kwargs)
        # Get list of betas and sd for betas
        try:
            dflist.append(process_ld_file(ldfile, clump[0]))
        except FileNotFoundError:
            print(f"Cannot find file {ldfile}")
    dfall = pd.concat(dflist, ignore_index=True)

    print(f"Export csv file to {outfile}")
    dfall.to_csv(outfile, sep="\t", index=False, header=True)


def sort_and_merge(clumplist):
    nclumps = len(clumplist)
    res = np.zeros(nclumps) * np.nan
    clumplist.sort(key=lambda x: x.begin)
    newclumplist = []
    clustid = 0
    res[0] = clustid
    clnew = clumplist[0]

    for i in range(1, nclumps):
        if clumplist[i].begin < clumplist[i - 1].end:
            # Found an overlap, then merge clumps
            res[i] = clustid
            clnew += clumplist[i]
        else:
            # Assign merged clumps into the results
            newclumplist.append(clnew)

            # Set new pointer to new cluster
            clnew = clumplist[i]

            # Update cluster id
            clustid += 1
            res[i] = clustid

    # Add the last merge of the cycle
    newclumplist.append(clnew)

    return newclumplist


def main(clumpfile, plinkfile, chr, outfile="",  **kwargs):

    try:
        totsize = int(kwargs["totsize"])
    except KeyError:
        totsize = 1e6
    except ValueError as e:
        print("Invalid totsize value. Set it to default: 1e6", e)
        totsize = 1e6

    try:
        # Get clumping file size
        fs = os.path.getsize(clumpfile)
    except OSError:
        # If the file does not exists, then set size to 0
        fs = 0

    # Prepare the file for output
    fout = open(outfile, "w")
    # If file size is not 0...
    if fs > 0:
        clumplist = read_clump_file(clumpfile, enlarge=True, merge=True,
                                    totsize=1e6)
        fw = csv.writer(fout, delimiter="\t")

        # Output header
        myheader = ["CHROM", "CLUMPID", "GENOFILE", "SNPLIST"]
        firstw = True
        for i, cl in enumerate(clumplist):
            snps = cl.to_list()
            genofile, snpfile = compute_ld(snps, clid=i,
                                        plinkfile=plinkfile,
                                        prefix=outfile, **kwargs)

            # Check if the written file exists and is not empty
            try:
                if os.path.getsize(genofile) > 0:
                    # Write header if it's the first time writing
                    if firstw:
                        fw.writerow(myheader)
                        firstw = False
                    fw.writerow([chr, i, genofile, snpfile])
            except FileNotFoundError:
                pass

    fout.close()

    return outfile


if __name__ == "__main__":
    clumplist = main(clumpfile=snakemake.input[0],
                     plinkfile=snakemake.params.plinkfile,
                     outfile=snakemake.output[0],
                     chr=snakemake.wildcards.chrom,
                     memory=snakemake.resources.mem_mb,
                     totsize=snakemake.params.totsize)
    # import argparse
    # parser = argparse.ArgumentParser()

    # parser.add_argument("--clumpfile", default="", required=True)
    # parser.add_argument("--outfile", default="")  # , required=True)
    # parser.add_argument("--dryrun", action='store_true')
    # parser.add_argument("--plinkfile", required=True)
    # parser.add_argument("--memory", default=8192)
    # parser.add_argument("--chr", type=int)

    # args = parser.parse_args()

    # clumplist = main(args.clumpfile, plinkfile=args.plinkfile,
    #                  outfile=args.outfile, chr=args.chr)
