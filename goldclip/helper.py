#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
support functions for goldclip
"""

import os
import sys
import re
import datetime
import json
import glob
import argparse
import shlex
import subprocess
import pathlib
import warnings
import logging
import numpy as np
import pandas as pd
import pysam
import pybedtools
import binascii

from goldclip.configure import goldclip_home
from goldclip.helper import *

logging.basicConfig(format = '[%(asctime)s] %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S', 
                    level = logging.DEBUG)



##-------------------------------------------##
## formatter
def nested_dict_values(d):
    """
    get all values from nested dict
    """
    for v in d.values():
        if isinstance(v, dict):
            yield from nested_dict_values(v)
        else:
            yield v
            

def is_gz(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


def xopen(fn, mode='r', bgzip=False):
    """
    Read / Write regular and gzip file, also support stdin
    """
    assert isinstance(fn, str)
    if fn == '-':
        return sys.stdin if 'r' in mode else sys.stdout
    if fn.endswith('.gz') and mode.startswith('w') or is_gz(fn):
        return gzip.open(fn, mode)
    else:
        return open(fn, mode)


def get_time():
    """
    get current time in this format:
    2006-01-02 14:23:35
    """
    return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')


def is_path(path, create = True):
    """
    Check path, whether a directory or not
    if not, create it
    """
    assert isinstance(path, str)
    if os.path.exists(path):
        return True
    else:
        if create:
            try:
                os.makedirs(path)
                return True
            except IOError:
                logging.error('failed to create directories: %s' % path)
        else:
            return False


def seq_type(fn, top_n = 1000):
    """
    Check the top 1000 rows of fn
    identify @ for fastq, > for fasta, * unknown
    """
    assert isinstance(fn, str)
    tag = set()
    with xopen(fn, 'rt') as fi:
        for i, line in enumerate(fi):
            if i > top_n:
                break
            elif i % 4 == 0:
                b = line[0] # the first base
                if b.lower() in 'acgtn':
                    continue
                else:
                    tag.add(line[0])
            else:
                continue
    if tag ==  {'@'}:
        return 'fastq'
    elif tag ==  {'>'}:
        return 'fasta'
    else:
        return None


def is_fastq(fn):
    if seq_type(fn) == 'fastq':
        return True
    else:
        return False


def is_fasta(fn):
    if seq_type(fn) == 'fasta':
        return True
    else:
        return False


def file_row_counter(fn):
    """
    count the file rows
    count '\n' 
    from @glglgl on stackoverflow, modified
    https://stackoverflow.com/a/9631635/2530783
    """
    def blocks(files, size = 1024 * 1024):
        while True:
            b = files.read(size)
            if not b: break
            yield b
    freader = gzip.open if is_gz(fn) else open
    with freader(fn, 'rt', encoding="utf-8", errors='ignore') as fi:
        return sum(bl.count('\n') for bl in blocks(fi))


def str_common(strList, suffix = False):
    # extract longest prefix/suffix from list of strings
    # default: prefix
    # sort strings by len
    def iterStop(exp):
        if exp is False:
            raise StopIteration
        else:
            return True    

    def commonPrefix(s1, s2):
        # prefix
        return ''.join(list(val for i, val in enumerate(s1) 
                       if iterStop(s2[i] is val)))

    def fact(l):
        if len(l) ==  1:
            return l[0]
        else:
            la = l[0:2]
            lb = l[2:]
            s = commonPrefix(la[0], la[1])
            lb.insert(0, s)
            return fact(lb)

    ## empty or single item 
    if len(strList) ==  0:
        return ''
    elif len(strList) ==  1:
        return strList[0]
    else:
        ## save a copy of list
        L2 = sorted(strList, key = len)
        c = fact(L2)
    
    ## suffix, reverse strings
    if suffix is True:
        L2 = [i[::-1] for i in L2]
        c = fact(L2)
        c = c[::-1]

    return c # string 0-index


def file_prefix(fn, with_path = False):
    """
    extract the prefix of a file
    remove extensions
    .gz, .fq.gz
    """
    assert isinstance(fn, str)
    p1 = os.path.splitext(fn)[0]
    px = os.path.splitext(fn)[1]
    if px.endswith('gz') or px.endswith('.bz'):
        px = os.path.splitext(p1)[1] + px
        p1 = os.path.splitext(p1)[0]
    if not with_path:
        p1 = os.path.basename(p1)
    return [p1, px]


def rm_suffix1(fn):
    """
    simplify the name of bam files
    from: {name}.not_{}.not_{}.....map_{}
    to: {name}
    """
    p = os.path.splitext(fn)[0]
    px = os.path.splitext(fn)[1]
    if px.startswith('.not_') or px.startswith('.map_'):
        fn = p
    return filename_shorter(fn)


def filename_shorter(fn, with_path=False):
    """
    input: name1.not_spikein.not_mtrRNA.map_genome.bam 
           name2.not_spikein.not_mtrRNA.map_genome.bam
    output: name1.bam
            name2.bam
    """
    p1 = os.path.splitext(fn)[0]
    px = os.path.splitext(fn)[1]
    p2 = rm_suffix1(p1)
    if not with_path:
        p2 = os.path.basename(p2)
    return p2 + px


def bed_parser(fn, usecols = None):
    """
    read BED file as pandas DataFrame
    select specific columns, default all, (None)
    require at least 6 columns
    """
    if not pathlib.Path(fn).is_file() or os.path.getsize(fn) ==  0:
        df = pd.DataFrame(columns = ['chr', 'start', 'end', 'name', 'score', 
                                     'strand'])
        logging.warning('empty bed file: %s' % fn)
        return df
    else:
        df = pd.read_table(fn, '\t', usecols = usecols, header = None,
            dtype = {'0': np.str, '1': np.int64, '2': np.int64, '3': np.str, \
                '4': np.int64, '5': np.str})
        df = df.rename(index = str, columns = {0: 'chr', 1: 'start', 2: 'end', \
                3: 'name', 4: 'score', 5: 'strand'})
        return bed_fixer(df)



def bed_filter(fn, bed_exclude, bed_out, overlap = True, save = True):
    """
    remove records from fn that have overlap with bed_exclude, and 
    save to fn using pybedtools
    overlap, True: intersect, False: not intersect
    """
    assert pathlib.Path(fn).is_file()
    bed_out_path = os.path.dirname(bed_out)
    assert is_path(bed_out_path)
    a = pybedtools.BedTool(fn)
    b = pybedtools.BedTool(bed_exclude)
    if overlap is True:
        a_and_b = a.intersect(b, wa = True, u = True) # intersect with b
    elif overlap is False:
        a_and_b = a.intersect(b, wa = True, v = True) # exclude b
    else:
        logging.error('unknown overlap: %s' % overlap)
    if save is True:
        a_and_b.moveto(bed_out)
    else:
        return a_and_b # BedTool object


def bed_fixer(df):
    """
    filt BED records 
    1. start, end both are int
    2. start < end
    """
    dx = df[['start', 'end']].apply(pd.to_numeric)
    c = ((dx['start'] >=  0) & dx['end'] >=  0) & (dx['start'] < dx['end'])
    return df.loc[c, :]


##--------------------------------------------##
## virtualenv 
def is_venv():
    """
    determine if inside virtualenv
    """
    return (hasattr(sys, 'real_prefix') or 
        (hasattr(sys, 'base_prefix') and sys.base_prefix !=  sys.prefix))

# current, run
def _base_prefix():
    if hasattr(sys, 'real_prefix'):
        return sys.real_prefix
    elif hasattr(sys, 'base_prefix'):
        return sys.base_prefix
    else:
        return None


# enter env
def _venv_into(venv, out = False):
    venv_bin = os.path.join(venv, 'bin', 'activate_this.py')
    if not out:
        if sys.version_info[0:2] ==  (2, 7):
            execfile(venv_bin, dict(__file__ = venv_bin)) # python2        
        elif sys.version_info[0:1] >=  (3, ):
            exec(open(venv_bin).read(), {}, dict(__file__ = venv_bin)) #python3
        else:
            logging.error('unknown version of python: ' + sys.version)        
    else:
        pass
        #subprocess.run(['deactivate'])


def venv_checker(venv = '~/envs/py27', into_venv = True):
    """
    check virtualenv 
    if not, go into
    into_venv, into env or out env
    """
    venv_in = os.path.expanduser(venv)
    
    if into_venv: # go into env
        if is_venv() and venv_in ==  _base_prefix():
            return ('already in venv: ' + venv_in)
        elif os.path.exists(venv_in):
            _venv_into(venv_in)
        else:
            logging.error('virtualenv not exists - ' + venv)
    else: # exit env
        if is_venv() and venv_in ==  sys.prefix:
            _venv_into(venv_in, out = True)
        else:
            return ('not in venv: ' + venv_in)



##--------------------------------------------##
## config
def bam2bw(genome, pathout, bam, strandness = True, binsize = 1):
    """
    Convert bam to bigWig using deeptools
    effective length of genome:
    dm3=121,400,000,
    mm9=2,150,570,000,
    hg19=2,451,960,000
    Mappable sequence of a genome, see Table 1 in 
    url: https://www.nature.com/articles/nbt.1518.pdf
    """
    assert is_path(pathout, dir_create = True)
    effsize = {'dm3': 121400000,
               'mm9': 2150570000,
               'hg19': 2451960000,}
    gsize = effsize[genome]
    prefix = os.path.basename(os.path.splitext(bam)[0])
    if strandness:
        fwd_bw = os.path.join(pathout, prefix + '.fwd.bigWig')
        rev_bw = os.path.join(pathout, prefix + '.rev.bigWig')
        c2 = 'bamCoverage -b {} -o {} --binSize {} --filterRNAstrand forward --normalizeTo1x {}'.format(bam, fwd_bw, binsize, gsize)
        c3 = 'bamCoverage -b {} -o {} --binSize {} --filterRNAstrand reverse --normalizeTo1x {}'.format(bam, rev_bw, binsize, gsize)
        if os.path.exists(fwd_bw) and os.path.exists(rev_bw):
            logging.info('file exists, bigWig skipped ...')
        else:
            subprocess.run(shlex.split(c2), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            subprocess.run(shlex.split(c3), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            # os.system(c2)
            # os.system(c3)
    else:
        bw = os.path.join(pathout, prefix + '.bigWig')
        c1 = 'bamCoverage -b {} -o {} --binSize {} --normalizeTo1x {}'.format(bam, bw, binsize, gsize)
        if os.path.exists(bw):
            logging.info('file exists, bigWig skipped ...')
        else:
            subprocess.run(shlex.split(c1), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            #os.system(c1)


class genome_info():
    # genome data_dir
    def __init__(self, genome, genome_path = '/data/genome'):
        genome_path = os.path.join(genome_path, genome)
        self.genome_fa = os.path.join(genome_path, 'bigZips', genome + '.fa')
        self.genome_fa_size = os.path.join(genome_path, 'bigZips', genome + '.chrom.sizes')
        self.bowtie_index = os.path.join(genome_path, 'bowtie_index', 'genome')
        self.bowtie2_index = os.path.join(genome_path, 'bwotie2_index', 'genome')
        self.hisat2_index = os.path.join(genome_path, 'hisat2_index', 'genome')
        self.star_index = os.path.join(genome_path, 'STAR_index')
        self.phylop_100 = os.path.join(genome_path, 'phyloP100way', genome + '.100way.phyloP100way.bw')
        self.gene_bed = os.path.join(genome_path, 'annotation_and_repeats', genome + '.refseq.bed')
        self.rmsk_bed = os.path.join(genome_path, 'annotation_and_repeats', genome + '.rmsk.bed')

    def get_fa(self):
        assert os.path.exists(self.genome_fa)
        return self.genome_fa

    def get_fasize(self):
        assert os.path.exists(self.genome_fa_size)
        return self.genome_fa_size

    def get_bowtie_index(self):
        assert os.path.exists(self.bowtie_index + ".1.ebwt")
        return self.bowtie_index

    def get_bowtie2_index(self):
        assert os.path.exists(self.bowtie_index + ".1.bt2")
        return self.bowtie2_index

    def get_hisat2_index(self):
        assert os.path.exists(self.bowtie_index + ".1.ht2")
        return self.hisat2_index

    def get_star_index(self):
        assert os.path.exists(os.path.join(self.bowtie_index, "Genome"))
        return self.star_index

    def get_phylop100(self):
        assert os.path.exists(self.bowtie_index)
        return self.phylop100

    def get_gene(self):
        assert os.path.exists(self.gene_bed)
        return self.gene_bed

    def get_rmsk(self):
        assert os.path.exists(self.rmsk_bed)
        return self.rmsk_bed

    def get_anno(self):
        assert os.path.exists(self.anno)
        return self.anno


def bam_merge(bam_ins, bam_out):
    """
    merge multiple bam files
    input: list of bam files
    input: out.bam
    """
    # check input files
    bam_flag = []
    for b in bam_ins:
        if not os.path.exists(b) is True:
            bam_flag.append(b)
    if len(bam_flag) > 0:
        sys.exit('BAM files not exists:' + '\n'.join(bam_flag))
    # check output file
    if os.path.exists(bam_out) is True:
        pass
        # sys.exit('BAM exists:' + bam_out)
    else:
        # merge
        pysam.merge('-f', bam_out + '.unsorted.bam', *bam_ins) # overwrite output BAM
        pysam.sort('-o', bam_out, bam_out + '.unsorted.bam')
        pysam.index(bam_out)
        os.remove(bam_out + '.unsorted.bam')



##--------------------------------------------##
## virtual env
def aligner_index_validator(path, aligner = 'bowtie'):
    """validate the alignment index, bowtie, bowtie2, hisat2"""
    # bowtie index
    inspect = aligner + '-inspect'
    cmd = [inspect, '-s', path]
    n = subprocess.run(cmd, check = False, stdout = subprocess.PIPE).stdout
    if len(n) > 0:
        return True
    else:
        return False


def aligner_index_constructer(genome, name, aligner = 'bowtie'):
    """return the bowtie_index for specific type/name"""
    #bin_dir = os.path.join(os.path.dirname(__file__), '..', 'data')
    idx = os.path.join(goldclip_home, 'data', genome, aligner + '_index', name, name)
    if aligner_index_validator(idx, aligner):
        return idx
    else:
        return None


def aligner_index_picker(genome, spikein):
    """pick aligner index for specific groups"""
    groups = ['viral', 'repeatRNA', 'retroviral', 'MT_trRNA']
    # groups = ['MT_trRNA']
    g_idxes = [aligner_index_constructer(genome, i) for i in groups]
    spikein_idx = aligner_index_constructer(spikein, 'genome')
    genome_idx = aligner_index_constructer(genome, 'genome')
    if not genome_idx:
        sys.exit('genome index not detected - ' + genome_idx)
    if not spikein == genome and not spikein:
        sys.exit('spikein index not detected - ' + spikein_idx)
    if spikein == genome:
        spikein_idx = None
    idx_output = [spikein_idx] + g_idxes + [genome_idx]
    idx_output = [i for i in idx_output if i] # remove None values
    return idx_output

