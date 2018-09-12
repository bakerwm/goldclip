#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import os
import io
import json
import tempfile


class Aligner_log(object):
    """Wrapper log file of aligner, bowtie, bowtie2, STAR
    report: total reads, unique mapped reads, multiple mapped reads

    Bowtie2:

    10000 reads; of these:
      10000 (100.00%) were unpaired; of these:
        166 (1.66%) aligned 0 times
        2815 (28.15%) aligned exactly 1 time
        7019 (70.19%) aligned >1 times
    98.34% overall alignment rate


    Bowtie:

    # reads processed: 10000
    # reads with at least one reported alignment: 3332 (33.32%)
    # reads that failed to align: 457 (4.57%)
    # reads with alignments suppressed due to -m: 6211 (62.11%)

    or:

    # reads processed: 10000
    # reads with at least one reported alignment: 9543 (95.43%)
    # reads that failed to align: 457 (4.57%)


    STAR:
    *final.Log.out

                                 Started job on |       Sep 12 11:08:57
                             Started mapping on |       Sep 12 11:11:27
                                    Finished on |       Sep 12 11:11:29
       Mapping speed, Million of reads per hour |       18.00

                          Number of input reads |       10000
                      Average input read length |       73
                                    UNIQUE READS:
                   Uniquely mapped reads number |       47
                        Uniquely mapped reads % |       0.47%
                          Average mapped length |       51.66
                       Number of splices: Total |       5
            Number of splices: Annotated (sjdb) |       0
                       Number of splices: GT/AG |       3
                       Number of splices: GC/AG |       0
                       Number of splices: AT/AC |       0
               Number of splices: Non-canonical |       2
                      Mismatch rate per base, % |       2.14%
                         Deletion rate per base |       0.04%
                        Deletion average length |       1.00
                        Insertion rate per base |       0.00%
                       Insertion average length |       0.00
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |       83
             % of reads mapped to multiple loci |       0.83%
        Number of reads mapped to too many loci |       19
             % of reads mapped to too many loci |       0.19%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |       0.02%
                 % of reads unmapped: too short |       98.31%
                     % of reads unmapped: other |       0.18%
                                  CHIMERIC READS:
                       Number of chimeric reads |       0
                            % of chimeric reads |       0.00%
    """


    def __init__(self, log):
        self.log = log
        # stat
        if isinstance(log, Aligner_log):
            self.stat = log.stat
        elif isinstance(log, dict):
            self.stat = log
        elif isinstance(log, io.TextIOWrapper):
            self.stat = self._log_parser()
        elif os.path.isfile(log):
            self.stat = self._log_parser()
        else:
            raise ValueError('not supported file')



    def _is_file(self):
        """Check the log file is exists, not empty
        """
        if os.path.isfile(self.log):
            return True
        else:
            return False



    def _is_non_empty(self):
        """Check if log file is empty"""
        if os.path.getsize(self.log) > 0:
            return True
        else:
            return False



    def _bowtie_log(self):
        """Wrapper bowtie log"""
        dd = {}
        with open(self.log, 'rt') as ff:
            for line in ff:
                if not ':' in line:
                    continue
                num = line.strip().split(':')[1]
                value = num.strip().split(' ')[0]
                value = int(value)
                if 'reads processed' in line:
                    dd['total'] = value
                elif 'at least one reported alignment' in line:
                    dd['map'] = value
                elif 'failed to align' in line:
                    dd['unmap'] = value
                elif 'alignments suppressed due to -m' in line:
                    dd['multiple'] = value
                else:
                    pass
        # unique and multiple
        if 'multiple' in dd:
            dd['unique'] = dd['map']
        else:
            dd['unique'] = dd['map'] # error
            dd['multiple'] = 0
        return dd



    def _bowtie2_log(self):
        """Wrapper bowtie2 log"""
        dd = {}
        with open(self.log, 'rt') as ff:
            for line in ff:
                value = line.strip().split(' ')[0]
                if '%' in value:
                    continue
                value = int(value)
                if 'reads; of these' in line:
                    dd['total'] = value
                elif 'aligned 0 times' in line:
                    dd['unmap'] = value
                elif 'aligned exactly 1 time' in line:
                    dd['unique'] = value
                elif 'aligned >1 times' in line:
                    dd['multiple'] = value
                else:
                    pass
        # unique and multiple
        dd['map'] = dd['unique'] + dd['multiple']
        return dd



    def _star_log(self):
        """Wrapper STAR *Final.log"""
        dd = {}
        with open(self.log, 'rt') as ff:
            for line in ff:
                value = line.strip().split('|')
                if not len(value) == 2:
                    continue
                value = value[1].strip()
                if 'Number of input reads' in line:
                    dd['total'] = int(value)
                elif 'Uniquely mapped reads number' in line:
                    dd['unique'] = int(value)
                elif 'Number of reads mapped to multiple loci' in line:
                    dd['multiple'] = int(value)
                else:
                    pass
        # unique and multiple
        dd['map'] = dd['unique'] + dd['multiple']
        dd['unmap'] = dd['total'] - dd['map']
        return dd


    def _log_parser(self):
        """Read log file as dict
        delimiter:
        bowtie:  ":"
        bowtie2:  ":"
        STAR:  "|"

        extra trimming
        1. trim "(10.00%)" 
        2. trim "blank" at both tails
        """
        log = self.log
        if isinstance(log, dict):
            return Aligner_log(log)
        elif isinstance(log, io.TextIOWrapper):
            line = next(log).strip()
        elif os.path.exists(log):
            with open(self.log, 'rt') as ff:
                line = next(ff).strip()
        else:
            raise ValueError('unknown file format')

        # parsing log file
        if line.startswith('#'):
            dd = self._bowtie_log()
        elif 'reads; of these' in line:
            dd = self._bowtie2_log()
        elif '|' in line:
            dd = self._star_log()
        else:
            raise ValueError('unknown file format')

        return dd



    def _tmp(self):
        """Create a temp file"""
        tmpfn = tempfile.NamedTemporaryFile(prefix='tmp',
                                            suffix='.json',
                                            delete=False)
        return tmpfn.name



    def saveas(self, _out=None):
        """Make a copy of statistics of mapping results"""
        log = self.log
        if _out is None:
            # _out = self._tmp()
            _out = os.path.splitext(log)[0] + '.json'

        dd = self.stat

        with open(_out, 'wt') as fo:
            json.dump(dd, fo, indent=4, sort_keys=True)

        return _out


## EOF
