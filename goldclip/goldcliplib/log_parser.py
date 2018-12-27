#!/usr/bin/env python
# -*- encoding: utf-8 -*-


__author__ = 'Ming Wang <wangm08@hotmail.com>'
__copyright__ = '2018 by Ming Wang <wangm08@hotmail.com>'
__license__ = 'MIT'
__email__ = 'wangm08@hotmail.com'
__version__ = '0.0.1'


import os
# import sys
import re
import io
import glob
import json
import fnmatch
import tempfile
import shlex
import subprocess
import logging
import pandas as pd
# import pysam
# import pybedtools
# from operator import is_not
# from functools import partial
# import multiprocessing as mp
from goldclip.bin.bed_annotation import *
from goldclip.helper import *


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


class Cutadapt_log(object):
    """Wrapper cutadapt log file"""

    def __init__(self, log):
        self.log = log
        # stat
        if isinstance(log, Cutadapt_log):
            self.stat = stat.stat
        elif isinstance(log, dict):
            self.stat = stat
        elif isinstance(log, io.TextIOWrapper):
            self.stat = self._log_parser()
        elif os.path.isfile(log):
            self.stat = self._log_parser()
        else:
            raise ValueError('not supported file')


    def _log_parser(self):
        """Wrapper log file"""
        dd = {}
        with open(self.log, 'rt') as ff:
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
            _out = os.path.splitext(self.log)[0] + '.json'
            # _out = self._tmp()

        dd = self.stat

        with open(_out, 'wt') as fo:
            json.dump(dd, fo, indent=4)

        return _out


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
        elif os.path.isfile(log):
            self.stat = self._log_parser()
        else:
            raise ValueError('not supported file: %s' % log)


    def guess_aligner(self):
        """Guess the aligner of the log file:
        bowtie, bowtie2, STAR, ...
        """
        # read through log file
        log_lines = []
        with open(self.log, 'rt') as ff:
            for r in ff:
                if r.startswith('Warning'):
                    continue
                log_lines.append(r.strip())

        # parsing log file
        line = log_lines[0] # the first line
        if line.startswith('#'):
            log_parser = self._bowtie_parser
        elif 'reads; of these' in line:
            log_parser = self._bowtie2_parser
        elif '|' in line:
            log_parser = self._star_parser
        else:
            raise ValueError('unknown file format: %s' % self.log)
            pass
        return log_parser


    def _is_non_empty(self):
        """Check if log file is empty"""
        if os.path.getsize(self.log) > 0:
            return True
        else:
            return False


    def _bowtie_parser(self):
        """Wrapper bowtie log
        unique, multiple, unmap, map, total
        """
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
        dd['multiple'] = dd.get('multiple', 0) # default 0
        if self.unique_only:
            dd['map'] = dd['unique']
        else:
            dd['map'] = dd['unique'] + dd['multiple']
        dd['unmap'] = dd['total'] - dd['unique'] - dd['multiple']
        return dd


    def _bowtie2_parser(self):
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
        if self.unique_only:
            dd['map'] = dd['unique']
        else:
            dd['map'] = dd['unique'] + dd['multiple']
        dd['unmap'] = dd['total'] - dd['unique'] - dd['multiple']
        return dd


    def _star_parser(self):
        """Wrapper STAR *final.log"""
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
        else:
            dd['map'] = dd['unique'] + dd['multiple']
        dd['unmap'] = dd['total'] - dd['unique'] - dd['multiple']
        return dd


    def _log_parser(self):
        """Read log file as dict
        delimiter:
        bowtie:  ":"
        bowtie2:  ":"
        STAR:  "|"

        extra trimming
        1. trim "(10.00%)" 
        2. trim "blank" at both ends
        """
        log = self.log
        log_parser = self.guess_aligner()
        dd = log_parser()
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
            json.dump(dd, fo, indent=4, sort_keys=False)

        return _out


class Alignment_stat(object):
    """Parse mapping reads in directory
    1. for each rep bam, parse json files, 
    2. merge replicates
    """
    def __init__(self, path):
        self.path = path

        if isinstance(path, Alignment_stat):
            self.stat = path.stat
        elif isinstance(path, pd.DataFrame):
            self.stat = path
        elif os.path.isdir(path):
            json_files = self.json_files()
            bam_files = self.bam_files()
            if json_files is None: # no json files
                self.stat = self.merge_stat() 
            elif len(json_files) == 1:
                self.stat = self.single_stat()
            elif json_files:
                self.stat = self.rep_stat()
            elif bam_files:
                self.stat = self.merge_stat()
            else:
                raise ValueError('BAM or json files not found: %s' % path)
        else:
            raise ValueError('unknown format')


    def _is_non_empty(self, fn):
        """Check if log file is empty"""
        if os.path.getsize(fn) > 0:
            return True
        else:
            return False


    def findfiles(self, which, where='.'):
        """Returns list of filenames from `where` path matched by 'which'
        shell pattern. Matching is case-insensitive.
        # findfiles('*.ogg')
        """    
        # TODO: recursive param with walk() filtering
        rule = re.compile(fnmatch.translate(which), re.IGNORECASE)
        hits = [os.path.join(where, name) for name in os.listdir(where) if rule.match(name)]
        return hits


    def json_files(self):
        """Return all json files in path"""
        j = self.findfiles('*.json', self.path)
        j = [i for i in j if self._is_non_empty(i)] # exclude empty json
        if len(j) == 0:
            j = None
        return j


    def bam_files(self):
        """Return all BAM files in path
        skip symlink
        # directory of merged fastq, count bam file
        # to-do: summary replicates (except unmap)
        """
        bam_files = self.findfiles('*.bam', self.path)
        bam_files = [b for b in bam_files if self._is_non_empty(b) and not os.path.islink(b)]
        if len(bam_files) == 0:
            bam_files = None
        return bam_files


    def _json_log_reader(self, fn, to_dataframe=True):
        """Parse json file, save as DataFrame"""
        if fn is None:
            return None
        fn_name = os.path.basename(fn)
        group, aligner = fn.split('.')[-3:-1] # 
        group = re.sub('^map_', '', group)
        dd = Json_file(fn).stat # dict of count
        rpt = dd
        if to_dataframe:
            # total, unique, multiple, unmap
            df = pd.DataFrame(data=[list(dd.values())], columns=list(dd.keys()))
            df.pop('map') # drop "map" column
            df.index = [fn_name]
            rpt = df
        return rpt


    def single_stat(self):
        """Extract alignment log files in index_ext directory"""
        path = self.path.rstrip('/')
        json_files = self.json_files()
        json_files = sorted(json_files, key=len)
        prefix = os.path.basename(self.path)
        dd = self._json_log_reader(json_files[0], False)
        df = pd.DataFrame(dd, index=[prefix])
        return df


    def rep_stat(self):
        """Extract alignment log files in directory
        """
        path = self.path.rstrip('/')
        json_files = self.json_files()
        json_files = sorted(json_files, key=len)
        prefix = os.path.basename(self.path)
        # genome_rRNA, genome, sp_rRNA, sp, unmap, map, total
        if len(json_files) == 2 or len(json_files) == 4:
            pass
        else:
            raise ValueError('number of json files should be; 2|4, [%d]' % len(json_files))

        # genome
        # map rRNA, genome, spikein_rRNA, spikein
        json_genome_rRNA, json_genome = json_files[:2]
        json_sp_rRNA = json_sp = None
        if len(json_files) == 4:
            json_sp_rRNA, json_sp = json_files[2:]
        df1 = self._json_log_reader(json_genome_rRNA, False)
        df2 = self._json_log_reader(json_genome, False)
        df3 = self._json_log_reader(json_sp_rRNA, False)
        df4 = self._json_log_reader(json_sp, False)

        # genome
        n_total = df1['total']
        g_rRNA = df1['unique'] + df1['multiple']
        g_unique = df2['unique']
        g_multi = df2['multiple']
        n_unmap = df2['unmap']

        # spikein rRNA
        if isinstance(df3, dict):
            sp_rRNA = df3['unique'] + df3['multiple']
        else:
            sp_rRNA = 0

        # spikein
        if isinstance(df4, dict):
            sp_unique = df4['unique']
            sp_multi = df4['multiple']
            n_unmap = df4['unmap']
        else:
            sp_unique = sp_multi = 0

        # output
        data = {
            'genome_rRNA': g_rRNA,
            'genome_unique': g_unique,
            'genome_multiple': g_multi,
            'spikein_rRNA': sp_rRNA,
            'spikein_unique': sp_unique,
            'spikein_multi': sp_multi,
            'unmap': n_unmap,
            'total': n_total}
        df = pd.DataFrame(data, index=[prefix])

        # return
        return df


    def merge_stat(self):
        """Stat reads for merged sample
        combine all reads in each replicates
        no json files detedted in merged directory
        search *.csv files in up-level directory
        """
        merge_path_name = os.path.basename(self.path.rstrip('/'))
        parent_path = os.path.dirname(self.path.rstrip('/')) # 
        rep_path = [i for i in os.listdir(parent_path) if not i == merge_path_name]
        rep_csv_files = [os.path.join(parent_path, i + '.mapping_stat.csv') for i in rep_path]
        rep_csv_files = [f for f in rep_csv_files if os.path.isfile(f)]

        if len(rep_csv_files) > 0:
            frames = [pd.read_csv(i, '\t', index_col=0) for i in rep_csv_files]
            df = pd.concat(frames, axis=0)
            # merge
            df_merge = pd.DataFrame(data=[df.sum(axis=0)], columns=list(df.columns.values),
                index=[merge_path_name])
            return df_merge
        else:
            logging.error('%10s | not contain mapping files: %s' % ('failed', path))
            return None


    def saveas(self, _out=None):
        """Make a copy of statistics of mapping results"""
        path = self.path
        if _out is None:
            prefix = os.path.basename(path)
            _out = os.path.join(os.path.dirname(path),
                                prefix + '.mapping_stat.csv')
        df = self.stat
        df.to_csv(_out, sep='\t', header=True, index=True)
        
        return _out


##-------------------------------------------------------------------##
## functions for report

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

##--------------------##
## figure 2
# def anno_run(bed, genome, group, output, path_data=None):
#     """
#     annotate bed files
#     output : Queue (df)
#     """
#     df = bed_annotator(bed, genome, group, path_data)
#     output.put(df)


# def bed_anno(bed_files, genome, group, path_data=None):
#     """
#     intput: genome_mapping/ 
#     return the genome mapped bed files
#     """
#     #call peaks in parallel
#     output = mp.Queue()
#     processes = [mp.Process(target = anno_run, 
#         args = (b, genome, group, output)) for b in bed_files]
#     for p in processes:
#         p.start() #start process
#     for p in processes:
#         p.join() #exit completed process
#     results = [output.get() for p in processes]
#     df = pd.concat(results, axis = 1)
#     return df

##--------------------##
## figure 3
def bam_corr(fns, path_out, window=10000, multi_cores=8):
    """
    calculate the Pearson correlation between BAM files
    use window size: 500, 1k, 10k, 100k,
    optional: gene_level, transcript_level, ...
    output: cor, p_val
    #
    use deeptools to calculate bin-count
    """
    bam_ids = [os.path.splitext(os.path.basename(f))[0] for f in fns]
    para_bam = ' '.join(fns)
    # bam_ids = [os.path.splitext(os.path.basename(f))[0] for f in [bam1, bam2]]
    out_npz = os.path.join(path_out, 'results.npz')
    out_tab = os.path.join(path_out, 'results.tab')
    cor_png = os.path.join(path_out, 'cor_scatter.png')
    cor_tab = os.path.join(path_out, 'cor_matrix.tab')
    c1 = 'multiBamSummary bins --bamfiles {} -out {} --outRawCounts {} \
        -p {} --binSize {}'.format(para_bam, out_npz, out_tab, multi_cores,
        window)
    c2 = 'plotCorrelation -in {} --corMethod pearson --skipZeros \
        --removeOutliers -T {} --whatToPlot scatterplot -o {} \
        --outFileCorMatrix {}'.format(out_npz, 'Pearson_correlation', cor_png,
        cor_tab)
    cmd1 = shlex.split(c1)
    cmd2 = shlex.split(c2)
    p = subprocess.run(cmd1)
    p = subprocess.run(cmd2)
    # cal pearson value
    df = pd.read_csv(out_tab, '\t')
    cor_mat = df.drop(df.columns[:3], axis = 1).corr('pearson')
    return cor_mat

    # multiBamSummary bins 
    #     --bamfiles bam1 bam2
    #     --out results.npz
    #     --outRawCounts results.tab
    #     --binSize window
    #     --smartLabels
    #     --binSize window
    #     -p multi_cores
        

    # plotCorrelation 
    #     -in out_npz 
    #     --corMethod pearson 
    #     --skipZeros
    #     --plotTitle "Pearson Correlation"
    #     --whatToPlot scatterplot
    #     -o cor_scatter.png
    #     --outFileCorMatrix cor_matrix.tab


    # plotCorrelation 
    #     -in out_npz 
    #     --corMethod pearson 
    #     --skipZeros
    #     --plotTitle "Pearson Correlation"
    #     --whatToPlot heatmap
    #     --colorMap RdYlBu 
    #     --plotNumbers
    #     -o cor_heatmap.png
    #     --outFileCorMatrix cor_matrix.tab


    # plotPCA 
    #     -in results.npz 
    #     -o results_PCA.png
    #     --plotTitle 'PCA of read counts'


def bed_conservation(bed_in, bed_out, extend = 0, genome = 'hg19'):
    con_phylop = os.path.join(goldclip_home, 'data', genome, 'pub_data', \
                              'phyloP100way', genome + '.100way.phyloP100way.bw')
    # create temp file, BED-6
    tmp = tempfile.mkstemp(suffix = ".bed")[1]
    # fix BED6
    df_in = pd.read_csv(bed_in, '\t', header = None).iloc[:, 0:6]
    df_in[[4]] = df_in[[4]].astype(int) # convert score to int
    df_in = df_in.astype('str')
    df_in.loc[:, 3] = df_in.loc[:, 0] + df_in.loc[:, 1] + df_in.loc[:, 2] + df_in.loc[:, 3]
    df_in.to_csv(tmp, '\t', header = False, index = False)
    # cmd
    c = 'bigWigAverageOverBed {} {} {} -bedOut={}'.format(con_phylop, tmp, bed_out + '.tab', bed_out)
    if extend > 0:
        c += ' -sampleAroundCenter={}'.format(extend)
    cmd = shlex.split(c)
    p = subprocess.run(cmd)


# ## functions
# def json_reader(file):
#     """
#     load json file as dict
#     """
#     with open(file) as f:
#         return json.load(f)


# def json_writer(d, file):
#     """
#     write dict to file in json format
#     """
#     with open(file, 'w') as f:
#         json.dump(d, f)


# ## wrapper functions for cutadapt ##
# def cutadapt_log_parser(log):
#     """
#     Parsing the output of cutadpat
#     support multiple run output
#     save data in JSON format
#     """
#     logdict = {}
#     _cutadapt_ver = []
#     _python_ver = []
#     _cmd = []
#     _fname = []
#     _raw = []
#     _clean = []
#     outfile = os.path.splitext(log)[0] + ".json"
#     with open(log, 'r') as f:
#         for line in f.readlines():
#             if(len(re.findall('^This is cutadapt', line)) == 1):
#                 r1 = re.findall(r'(cutadapt|Python) (\d+\.\d+\.?\d+?)', line) # 1st line, tools
#                 for i in r1: logdict[i[0]] = i[1]
#             elif('Command line parameters' in line):
#                 _cmd.append(re.sub(r'Command line parameters: ', '', line).strip('\n'))
#                 _fname.append(os.path.basename(line.split()[-1])) # file name
#             elif('Total reads processed:' in line):
#                  _raw.append(re.sub(r',', '', line.split()[-1])) # first input
#             elif('Reads written (passing filters):' in line):
#                  _clean.append(re.sub(r',', '', line.split()[-2])) # output
#             else:
#                  continue
#     logdict['filename'] = _fname[0]
#     logdict['run_cutadapt_times'] = len(_raw)
#     logdict['command_lines'] = _cmd
#     logdict['raw'] = _raw[0]
#     logdict['clean'] = _clean[-1]
#     logdict['clean_pct'] = '{:.2f}%'.format(int(logdict['clean'])/int(logdict['raw'])*100)
#     with open(outfile, 'w') as fo:
#         json.dump(logdict, fo, indent = 4)
#     return outfile




# ## wrapper functions for bowtie mapping ##
# def bowtie_log_parser(path):
#     """
#     Parsing the log file of bowtie (to stderr)
#     fetch Input, mapped, unmapped reads
#     save in JSON format
#     return dict of all values
#     """
#     # !!! to-do
#     # convert to DataFrame, json, dict
#     # unmap = input - reported
#     logdict = {}
#     _input = []
#     _one_hit = []
#     _not_hit = []
#     _rpt = []
#     with open(path, 'r') as f:
#         for line in f.readlines():
#             if 'reads processed' in line:
#                 _num1 = re.sub(',', '', line.split()[-1])
#                 _input.append(int(_num1))
#             elif 'at least one reported' in line:
#                 _num2 = re.sub(',', '', line.split()[-2])
#                 _one_hit.append(int(_num2))
#             elif 'Reported' in line:
#                 _num4 = re.sub(',', '', line.split()[1])
#                 _rpt.append(int(_num4))
#             elif 'failed to align' in line:
#                 _num3 = re.sub(',', '', line.split()[-2])
#                 _not_hit.append(int(_num3))
#             else:
#                 continue
#     logdict['input_reads'] = _input[0] # first one
#     # logdict['mapped'] = sum(_one_hit) # sum all
#     logdict['mapped'] = sum(_rpt) # sum all
#     logdict['unmapped'] = int(_input[-1]) - int(_one_hit[-1]) # -m suppress
#     logdict['map_pct'] = '{:.2f}%'.\
#         format(int(logdict['mapped']) / int(logdict['input_reads'])*100)
#     json_out = os.path.splitext(path)[0] + '.json'
#     with open(json_out, 'w') as fo:
#         json.dump(logdict, fo, indent = 4)
#     return logdict


# def bowtie2_log_parser(path):
#     """
#     Parsing the log file of bowtie
#     fetch Input, unique, multiple, unmapped
#     save in JSON format
#     return dict of all values
#     """
#     logdict = {}
#     _input = []
#     _one_hit = []
#     _multi_hit = []
#     _not_hit = []
#     _rpt = []
#     with open(path, 'rt') as fi:
#         for line in fi.readlines():
#             line = line.strip()
#             _num = line.split(' ')[0]
#             _num = _num.strip('%')
#             _num = float(_num)
#             if line.endswith('reads; of these:'):
#                 _input.append(_num)
#             elif ') aligned 0 times' in line:
#                 _not_hit.append(_num)
#             elif ') aligned exactly 1 time' in line:
#                 _one_hit.append(_num)
#             elif ') aligned >1 times' in line:
#                 _multi_hit.append(_num)
#             else:
#                 continue
#     # save to dict
#     logdict['input_reads'] = int(_input[0]) # first one
#     logdict['mapped'] = int(sum(_one_hit + _multi_hit))
#     logdict['unique'] = int(sum(_one_hit))
#     logdict['multi'] = int(sum(_multi_hit))
#     logdict['unmapped'] = int(_not_hit[-1]) # -m suppress
#     logdict['map_pct'] = '{:.2f}%'.\
#         format(int(logdict['mapped']) / int(logdict['input_reads'])*100)
#     json_out = os.path.splitext(path)[0] + '.json'
#     with open(json_out, 'w') as fo:
#         json.dump(logdict, fo, indent=4)
#     return logdict


# def star_log_parser(path):
#     logdict = {}
#     with open(path) as f:
#         for line in f:
#             sep = line.strip().split('|')
#             if 'Number of input reads' in line:
#                 logdict['input_reads'] = int(sep[1].strip())
#             elif 'Uniquely mapped reads number' in line:
#                 logdict['unique'] = int(sep[1].strip())
#             elif 'Number of reads mapped to multiple loci' in line:
#                 logdict['multi'] = int(sep[1].strip())
#             else:
#                 pass
#     logdict['mapped'] = logdict['unique'] + logdict['multi']
#     logdict['unmapped'] = logdict['input_reads'] - logdict['mapped']
#     logdict['map_pct'] = '{:.2f}%'.format(logdict['mapped'] / logdict['input_reads'] * 100)
#     json_out = os.path.splitext(path)[0] + '.json'
#     with open(json_out, 'wt') as fo:
#         json.dump(logdict, fo, indent=4)
#     return logdict


# def rep_map_wrapper(path, save=True):
#     """
#     wrap all bowtie log files, [only for this script] namespace 
#     summarize mapping and RTStops 
#     header: name, group, read, RTStop
#     input: list of json files
#     output: pd.DataFrame
#     """
#     def _json_wrapper(fn): 
#         """
#         Only for bowtie map statistics
#         parsing only one json file
#         output: type, count
#         """
#         group = fn.split('.')[-3] # group name, *map_genome.bowtie.json
#         group = group.split('_')[1] # reference name
#         name = os.path.splitext(os.path.basename(fn))[0]
#         with open(fn, 'r') as f:
#             da = json.load(f)
#         df = [name, group, da['input_reads'], da['mapped'], 
#               da['unmapped']]
#         return df

#     # multiple json files
#     json_files = sorted(glob.glob(path + '/*.json'), key=len)
#     rep_prefix = os.path.basename(path)
#     df = pd.DataFrame(columns=['name', 'group', 'read'])
#     # 
#     for n in range(len(json_files)):
#         _, g, _, m1, _ = _json_wrapper(json_files[n])
#         if n == 0 and len(json_files) > 1 and g == 'genome':
#             g = 'spikein' # the first one - genome
#         df = df.append(pd.DataFrame([[rep_prefix, g, m1]], 
#             columns = ['name', 'group', 'read']), ignore_index=True)
#     _, gx, _, _, un = _json_wrapper(json_files[-1]) # unmap
#     df = df.append(pd.DataFrame([[rep_prefix, 'unmapped', un]], 
#                    columns=['name', 'group', 'read']), ignore_index=True)
#     save_csv = os.path.join(os.path.dirname(path), 
#                             rep_prefix + '.mapping_stat.csv')
#     if save:
#         df.to_csv(save_csv, ',', header=True, index=False)
#     return df



# def merge_map_wrapper(path, save=True):
#     """
#     count BAM files
#     Output: pd.DataFrame
#     """
#     bam_files = sorted(glob.glob(path + '/*.bam'), key=len) # bam files
#     merge_prefix = os.path.basename(path)
#     df = pd.DataFrame(columns=['name', 'group', 'read'])
#     bam_files = [f for f in bam_files if not os.path.islink(f)]
#     # iterate
#     for n in range(len(bam_files)):
#         b_cnt = pysam.AlignmentFile(bam_files[n], 'rb').count()
#         group = bam_files[n].split('.')[-2] # group name*.map_genome.bam
#         group = group.split('_')[1] # reference name
#         if n == 0 and len(bam_files) > 1 and group == 'genome':
#             group = 'spikein' # the first one - genome
#         dfx = pd.DataFrame([[merge_prefix, group, b_cnt]],
#                             columns=['name', 'group', 'read'])
#         df = df.append(dfx, ignore_index=True)
#     save_csv = os.path.join(os.path.dirname(path), 
#                             merge_prefix + '.mapping_stat.csv')
#     if save:
#         df.to_csv(save_csv, ',', header=True, index=False)
#     return df

