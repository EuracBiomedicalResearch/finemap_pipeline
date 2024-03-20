import pandas as pd
import numpy as np
import subprocess as sb


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


def widen_span(span, tot):
    span2 = np.array(span)
    span2.sort()
    spandist = np.diff(span2)[0]
    totdiff = tot - spandist

    if totdiff > 0:
        span2[0] -= totdiff / 2.0
        span2[1] += totdiff / 2.0

    return span2

# def run_clumping(sumstat, plinkfiles, **kwargs):

# plink2 --bfile {params.infile} --clump {input.smstat} \
# --clump-log10 'input-only' --clump-field {pvalcol} \
# --clump-log10-p1 {params.clump_logp1} --clump-log10-p2 {params.clump_logp2} \
# --clump-r2 {params.clump_r2} --clump-kb {params.clump_kb} \
# --clump-snp-field ID  --out {params.ofile} \
# --memory {resources.mem_mb} \
# --keep {params.sampfile}

