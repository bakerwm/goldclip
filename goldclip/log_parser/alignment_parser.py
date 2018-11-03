#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import os
import io
import json
import glob
import tempfile
import pysam
import pandas as pd
from utils_parser import *


class Alignment_log(object):
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


    def __init__(self, log, unique_only=False):
        self.log = log
        self.unique_only = unique_only
        # stat
        if isinstance(log, Alignment_log):
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
                if not ':' in line or line.startswith('Warning'):
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
        # unique_only
        dd['unique'] = dd['map']
        dd['multiple'] = 0 #no multiple
        # if self.unique_only is True or 'multiple' in dd:
        #     pass
        # else:
        #     dd['multiple'] = 0
        dd['unmap'] = dd['total'] - dd['map']
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
        if self.unique_only is True:
            dd['map'] = dd['unique']
            # dd['multiple'] = 0
        else:
            # unique and multiple
            dd['map'] = dd['unique'] + dd['multiple']
        dd['unmap'] = dd['total'] - dd['map']
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
        if self.unique_only is True:
            dd['map'] = dd['unique']
            # dd['multiple'] = 0
        else:
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
        log_lines = []
        if isinstance(log, dict):
            return Alignment_log(log)
        elif isinstance(log, io.TextIOWrapper):
            for r in log:
                if r.startswith('Warning'): # skip warnings
                    continue
                log_lines.append(r.strip())
        elif os.path.exists(log):
            with open(self.log, 'rt') as ff:
                for r in ff:
                    if r.startswith('Warning'):
                        continue
                    log_lines.append(r.strip())
        else:
            raise ValueError('unknown file format')
        
        # parsing log file
        line = log_lines[0] # the first line
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



class Alignment_stat(object):
    """Parse mapping reads in directory
    1. for each rep bam, parse json files, 
    2. merged bam files, count bam lines
    """

    def __init__(self, path):
        self.path = path
        if isinstance(path, Alignment_stat):
            self.stat = path.stat
        elif isinstance(path, dict):
            self.stat = path
        elif os.path.exists(path):
            if not self._get_json_files() is False:
                self.stat = self.count_json_files()
            elif not self._get_bam_files() is False:
                self.stat = self.count_bam_files()
            else:
                raise ValueError('No bam and json files found: %s' % path)
        else:
            raise ValueError('unknown format')


    def _is_non_empty(self, fn):
        """Check if log file is empty"""
        if os.path.getsize(fn) > 0:
            return True
        else:
            return False


    def _get_json_files(self):
        path = self.path
        j_files = sorted(glob.glob(path + '/*.json'), key=len)
        j_files = [f for f in j_files if self._is_non_empty(f)] # not empty files
        if len(j_files) > 0:
            return j_files
        else:
            # raise ValueError('No json files detected in: %s' % path)
            return False


    # parse *.json files
    def count_json_files(self):
        path = self.path
        prefix = os.path.basename(path) # sample name
        j_files = self._get_json_files() # each group
        df = pd.DataFrame(columns=['name', 'group', 'count'])
        for j in j_files:
            dd = Json_file(j).stat # count
            # group            
            group = j.split('.')[-3] # group name, *map_genome.bowtie.json
            group = group.split('_')[1] # 
            # check spike-in
            if j_files.index(j) == 0 and group == 'genome' and len(j_files) > 1:
                group = 'spikein'
            num_map = dd['map']
            df = df.append(pd.DataFrame([[prefix, group, num_map]],
                           columns = ['name', 'group', 'count']),
                           ignore_index=True)
        # unmap reads
        dd = Json_file(j_files[-1]).stat
        unmap = dd['unmap']
        df = df.append(pd.DataFrame([[prefix, 'unmap', unmap]],
                       columns=['name', 'group', 'count']),
                       ignore_index=True)
        return df


    def _get_bam_files(self):
        path = self.path
        bam_files = sorted(glob.glob(path + '/*.bam'), key=len)
        bam_files = [f for f in bam_files if self._is_non_empty(f) 
                     and not os.path.islink(f)] # not empty files
        if len(bam_files) > 0:
            return bam_files
        else:
            raise ValueError('No .bam files found in: %s' % path)


    # count bam files
    def count_bam_files(self):
        path = self.path
        prefix = os.path.basename(path)
        bam_files = self._get_bam_files()
        df = pd.DataFrame(columns=['name', 'group', 'count'])
        for b in bam_files:
            b_cnt = pysam.AlignmentFile(b, 'rb').count()
            group = b.split('.')[-2] # group name*.map_genome.bam
            group = group.split('_')[1] # reference name
            if bam_files.index(b) == 0 and group == 'genome' and len(bam_files) > 1:
                group = 'spikein'
            df = df.append(pd.DataFrame([[prefix, group, b_cnt]],
                           columns=['name', 'group', 'count']),
                           ignore_index=True)
        # output
        return df


    def _tmp(self):
        """Create a temp file"""
        tmpfn = tempfile.NamedTemporaryFile(prefix='tmp',
                                            suffix='.csv',
                                            delete=False)
        return tmpfn.name


    def saveas(self, _out=None):
        """Make a copy of statistics of mapping results"""
        path = self.path
        if _out is None:
            prefix = os.path.basename(path)
            _out = os.path.join(os.path.dirname(path),
                                prefix + '.mapping_stat.csv')
        df = self.stat        

        default_kwargs = dict(sep='\t', header=False, index=False)
        if isinstance(_out, io.TextIOWrapper):
            print(df.to_string(index=False, header=False, justify='left'))
        else:
            df.to_csv(_out, **default_kwargs)
        
        return _out


## EOF
