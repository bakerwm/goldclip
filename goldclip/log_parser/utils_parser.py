#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import os
import io
import re
import json
import glob
import tempfile
import pysam
import pandas as pd


class Cutadapt_log(object):
    """Wrapper cutadapt log file"""

    def __init__(self, stat):
        self.stat = stat
        # stat
        if isinstance(stat, Cutadapt_log):
            self.stat = stat.stat
        elif isinstance(stat, dict):
            self.stat = stat
        elif isinstance(stat, io.TextIOWrapper):
            self.stat = self._log_parser()
        elif os.path.isfile(stat):
            self.stat = self._log_parser()
        else:
            raise ValueError('not supported file')


    def _log_parser(self):
        """Wrapper log file"""
        dd = {}
        with open(self.stat, 'rt') as ff:
            for line in ff:
                if line.startswith('This is cutadapt'):
                    sep = line.strip().split(' ')
                    dd['version'] = sep[3]
                    dd['python'] = sep[6]
                elif 'Command line parameters' in line:
                    dd['cmd'] = line.strip().split(':')[1]
                elif 'Total reads processed' in line:
                    value = line.strip().split(':')[1]
                    value = re.sub(',', '', value.strip())
                    dd['total'] = int(value)
                elif 'Reads written (passing filters)' in line:
                    value = line.strip().split(':')[1]
                    value = value.strip().split(' ')[0]
                    value = re.sub(',', '', value)
                    dd['clean'] = int(value)
                else:
                    continue
        pct = float(dd['clean']) / float(dd['total']) * 100
        dd['pct'] = '%.1f%%' % pct
        return dd



    def _tmp(self):
        """Create a temp file"""
        tmpfn = tempfile.NamedTemporaryFile(prefix='tmp',
                                            suffix='.json',
                                            delete=False)
        return tmpfn.name



    def saveas(self, _out=None):
        """Make a copy of statistics of mapping results"""
        if _out is None:
            _out = self._tmp()

        dd = self.stat

        with open(_out, 'wt') as fo:
            json.dump(dd, fo, indent=4)

        return _out



class Json_file(object):
    """Parsing Json and dict file"""

    def __init__(self, fn):
        self.fn = fn
        if isinstance(fn, Json_file):
            self.stat = fn.stat
        elif isinstance(fn, dict):
            self.stat = fn
        elif os.path.exists(fn):
            self.stat = self.json_reader()
        else:
            raise ValueError('unknown file format:')
          

    def _tmp(self):
        """Create a temp file"""
        tmpfn = tempfile.NamedTemporaryFile(prefix='tmp',
                                            suffix='.json',
                                            delete=False)
        return tmpfn.name
        

    def json_reader(self):
        """Load json file as dict"""
        fn = self.fn
        if os.path.isfile(fn) and os.path.getsize(fn) > 0:
            with open(fn, 'rt') as ff:
                return json.load(ff)



    def json_writer(self, to=None):
        """Write dict to file in json format"""
        fn = self.fn

        if to is None:
            to = self._tmp()

        if isinstance(fn, Json_file):
            fn = fn.fn
        elif isinstance(fn, dict):
            fn = fn
        with open(to, 'wt') as ff:
            json.dump(fn, ff, indent=4, sort_keys=True)



class Alignments(object):
    """Parse mapping reads in directory"""

    def __init__(self, path):
        self.path = path
        if isinstance(path, Alignments):
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



##--------------------##
## figure 1
def trim_wrapper(path, smp_name='demo'):
    """
    trimming and remove duplicates
    input: /path_out/input_reads/
    """
    json_files = sorted(glob.glob(path + '/*.cutadapt.json'))
    da = []
    for j in json_files:
        id = re.sub(r'.cutadapt.json', '', os.path.basename(j))
        nodup = os.path.join(os.path.dirname(j), id + '.reads.txt') # total reads, nodup
        d = json_reader(j)
        with open(nodup) as f:
            d['nodup'] = next(f).rstrip()
        tooshort = int(d['raw']) - int(d['clean'])
        dup = int(d['clean']) - int(d['nodup'])
        dn = pd.DataFrame({'group': ['raw', 'too_short', 'PCR_dup', 'no_dup'],
            id: [d['raw'], tooshort, dup, d['nodup']]})
        dn.set_index('group', inplace = True)
        da.append(dn)
    df = pd.concat(da, axis = 1)
    df = df.apply(pd.to_numeric)
    # add merge data
    df.insert(0, smp_name, df.sum(axis = 1))
    return df

# path_trim = os.path.join(path_out, 'input_reads')
# df = trim_wrapper(path_trim, smp_name)
# print(df)


## mapping pct
def map_wrapper(path, smp_name='demo'):
    """
    mapping to various genome
    input: /path_out/genome_mapping/
    """
    m_files = glob.glob(os.path.join(path, '*.mapping_stat.csv'))
    m_files = sorted(m_files, key=len)
    ma = []
    for m in m_files:
        # skip merge stat
        m_prefix = re.sub(r'.mapping_stat.csv', '', os.path.basename(m))
        if m_prefix == smp_name:
            continue
        dm = pd.read_csv(m, ',').filter(items=['group', 'read'])
        dm.set_index('group', inplace=True)
        dm2 = dm.rename(columns={'read': m_prefix})
        ma.append(dm2)
    df = pd.concat(ma, axis=1)
    # add merge data
    df.insert(0, smp_name, df.sum(axis=1))
    return df







