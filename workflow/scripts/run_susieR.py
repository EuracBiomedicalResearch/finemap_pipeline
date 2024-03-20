import os
import csv
import pandas as pd
import numpy as np
import subprocess as sb
from tempfile import NamedTemporaryFile, TemporaryDirectory
import clumps
import gcta
import credible_sets as cs


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
                cc = clumps.Clump.from_clumpfileline(ll)
                mylist.append(cc)

    # Enlarge clump to a minimum size of `totsize`
    # TODO: check if totsize is > 0...and other checks on totsize

    if enlarge:
        try:
            totsize = kwargs['totsize']
            totsize = np.abs(float(totsize))
        except KeyError:
            totsize = 1e6
        mylist = [cl.enlarge_clump(totsize=totsize) for cl in mylist]

    # Merge overlapping clumps
    if merge:
        mylist = sort_and_merge(mylist)

    for i, cl in enumerate(mylist):
        ii = cl.get_span_from_clump()
        cl.set_interval(ii)
        mylist[i] = cl

    return mylist


def plink_wrapper(sumstat, plinkfile, logp1=7.3, logp2=1.3, r2=0.1,
                  kb=10000, **kwargs):

    cmd = ['plink2',
           '--bfile', plinkfile,
           '--clump', sumstat,
           '--clump-log10', "'input-only'"]
    if logp1 > 0:
        cmd.extend(['--clump-log10-p1', str(logp1)])
    if logp2 > 0:
        cmd.extend(['--clump-log10-p2', str(logp2)])
    if r2 > 0 and r2 < 1:
        cmd.extend(['--clump-r2', str(r2)])
    if kb:
        cmd.extend(['--clump-kb', str(np.abs(kb))])

    try:
        pvalcol = kwargs['pvalcol']
        cmd.extend(['--clump-field', pvalcol])
    except KeyError:
        print('Use default p-value column in summary stats')

    try:
        snpsid = kwargs['ID']
        cmd.extend(['--clump-snp-field', 'ID'])
    except KeyError:
        pass

    try:
        snplist = kwargs['snplist']
    except KeyError:
        snplist = []

    try:
        memory = kwargs['memory']
        cmd.extend(['--memory', str(memory)])
    except KeyError:
        pass


    with TemporaryDirectory() as ftmp:
        if snplist:
            print(f"Select {len(snplist)} snps...")
            snpfile = os.path.join(ftmp, "snplist")
            with open(snpfile, 'w') as sp:
                for s in snplist:
                    sp.write(s + "\n")
            cmd.extend(['--extract', snpfile])

        # Add output file parameter
        ofile = os.path.join(ftmp, "clump")
        cmd.extend(['--out', ofile])

        # Run the command and read file
        cp = sb.run(' '.join(cmd), shell=True)

        clumplist = []
        if cp.returncode == 0:
            # Read clumping file
            clumplist = read_clump_file(ofile + ".clumps", **kwargs)

    return clumplist


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

def run_cred_set_susier(clumplist, sumstatdf, plinkfile, chrom,
                        clid=0, **kwargs):
    basename = os.path.dirname(__file__)
    reslist = []
    for i, cl in enumerate(clumplist):
        print(f"Processing clumpid: {clid} - {i}")
        span = cl.get_interval()
        lead = sumstatdf.loc[sumstatdf['ID'] == cl.lead, :]

        # Extract window of the sumstat
        sm_wind = ((sumstatdf['CHROM'] == lead.iloc[0]['CHROM']) &
                    (sumstatdf['GENPOS'] >= span[0]) &
                    (sumstatdf['GENPOS'] <= span[1]))
        snplist = sumstatdf.loc[sm_wind, 'ID']

        with TemporaryDirectory() as tmpd:
            myprefix = os.path.join(tmpd, "ld_clump")
            ld, snplistfile = compute_ld(snplist.to_list(), plinkfile=plinkfile,
                            clid=clid,
                            prefix=myprefix, **kwargs)
            smfile = os.path.join(tmpd, "smfile.csv")
            sumstatdf.to_csv(smfile, index=False, sep='\t')

            outfile = os.path.join(tmpd, "cs_out.csv")
            cmd = ['Rscript', os.path.join(basename, 'susier.R'),
                   ld, smfile,
                   outfile, str(clid)]

            print(f"File {ld} {os.path.exists(ld)}")
            # Run susier.py on the window
            cp = sb.run(' '.join(cmd), shell=True)
            if cp.returncode != 0:
                print("Azzzz....got an error!")
                resdf = None
            else:
                resdf = pd.read_csv(outfile, header=0, sep="\t")
                reslist.append(resdf)
    resdf = pd.concat(reslist)
    return resdf


def main_clumping(sumstat, plinkfile, outfile, chrom, memory=16000, **kwargs):

    # try:
    #     totsize = int(kwargs["totsize"])
    # except KeyError:
    #     totsize = 1e6
    # except ValueError as e:
    #     print("Invalid totsize value. Set it to default: 1e6", e)
    #     totsize = 1e6

    analysis_conf = kwargs["analysis_conf"]
    sumstatdf = pd.read_csv(sumstat, header=0, sep='\t')

    print("Clumping Step1")
    clumplist = plink_wrapper(sumstat, plinkfile, memory=memory,
                              **analysis_conf['clumping'], **kwargs,
                              enlarge=True, merge=True, )
    print("Clumping Step2")
    for i, cl1 in enumerate(clumplist):
        clumplist2 = plink_wrapper(sumstat, plinkfile, memory=memory,
                                logp1=-1, logp2=-1, kb=500, r2=0.5,
                                snplist=cl1.to_list(), **kwargs,
                                ID=analysis_conf["clumping"]["ID"],
                                pvalcol=analysis_conf["clumping"]["pvalcol"]
                                )
        resdf = run_cred_set_susier(clumplist2, sumstatdf, plinkfile, chrom=chrom,
                            clid=i)
        print(resdf)

    # print(spanres)
    # print(np.diff(spanres))
    # print(f"Clump1 span: {np.diff([mi1, ma1])}")

def main_conditional(sumstat, plinkfile, outfile="", chrom=None, memory=16000,
                     **kwargs):
    analysis_conf = kwargs.pop("analysis_conf")
    # analysis_conf = kwargs["analysis_config"]
    sumstatdf = pd.read_csv(sumstat, header=0, sep='\t')

    try:
        tmpdir = kwargs.pop("tmpdir")
    except KeyError:
        tmpdir = None
    # Get top loci based on conditional analysis
    cc = gcta.ConditionalAnalysis(tmpdir=tmpdir)
    toploci = cc.get_top_loci(sumstatdf, plinkfile, chrom=chrom, **kwargs)

    # Initialize the results
    cred_res = []
    for i, index_var in enumerate(toploci.to_dict(orient='records')):
        cred_set = cs.credible_set(index_var=index_var, sumstat=sumstatdf,
                                   toploci=toploci, plinkfile=plinkfile, prior_sd=1.0,
                                   cs_prob=0.95, stype="quant", tmpdir=tmpdir
                                   )
        if cred_set.shape[0] > 0:
            cred_set['index_var'] = index_var['ID']
            cred_set['csid'] = f"{chrom}_{i}"
            cred_res.append(cred_set)

    # Concatenate results
    cred_res_df = pd.concat(cred_res, axis=0)

    return cred_res_df


if __name__ == "__main__":
    finemap_config = snakemake.config

    clumplist = main(clumpfile=snakemake.input[0],
                     plinkfile=snakemake.params.plinkfile,
                     outfile=snakemake.output[0],
                     chr=snakemake.wildcards.chrom,
                     memory=snakemake.resources.mem_mb,
                     totsize=snakemake.params.totsize,
                     analysis_config=finemap_config)
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
