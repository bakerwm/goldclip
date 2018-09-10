#!/usr/bin/env python
"""
Trimming reads
1. cut 3-adapter
2. trim low quality bases
3. remove N reads
4. trim N-bases from either ends of reads
5. limit read length

functions
"""

__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__date__ = "2018-03-21"
__version__ = "0.1"

import os
# import sys
# import re
import subprocess
import shlex
import logging
import goldclip
from goldclip.helper import *
from goldclip.goldcliplib.log_parser import *



def adapter_chopper(s, step=2, window=15):
    """
    chop the adapter by given length
    a series of adapters to trim
    """
    assert isinstance(s, str)
    assert isinstance(step, int)
    assert isinstance(window, int)
    p = []
    if len(s) <= window:
        p.append(s)
    else:
        for i in range(int(len(s) / step)):
            a = i * step
            b = a + window
            #b = b if b < len(s) else len(s)
            if b > len(s):
                continue
            p.append(s[a:b])
    return p



def dup_remover(fn, path_out, *, q=20, p=80):
    """
    Remove duplicate fastq sequences using fastx-collapser (fq to fa)
    convert fa to fastq
    optional:
    trim N-bases
    """
    pkg_dir, _ = os.path.split(goldclip.__file__)
    fa2fq = os.path.join(pkg_dir, 'bin', 'fasta_to_fastq.pl')
    path_out = os.path.dirname(fn) if path_out is None else path_out
    assert is_path(path_out)
    ## q, p, m
    fn_out_name = file_prefix(fn)[0] + '.nodup.fastq'
    fn_out_file = os.path.join(path_out, fn_out_name)
    if not os.path.exists(fn_out_file):
        freader = 'zcat' if is_gz(fn) else 'cat'
        c1 = '{} {}'.format(freader, fn)
        c2 = 'fastq_quality_filter -Q33 -q {} -p {}'.format(q, p)
        c3 = 'fastx_collapser -Q33'
        c4 = 'perl {} -'.format(fa2fq)
        cmd1 = shlex.split(c1)
        cmd2 = shlex.split(c2)
        cmd3 = shlex.split(c3)
        cmd4 = shlex.split(c4)
        p1 = subprocess.Popen(cmd1, stdout = subprocess.PIPE)
        p2 = subprocess.Popen(cmd2, stdin = p1.stdout, stdout = subprocess.PIPE)
        p3 = subprocess.Popen(cmd3, stdin = p2.stdout, stdout = subprocess.PIPE)
        with open(fn_out_file, 'wt') as fo:
            p4 = subprocess.Popen(cmd4, stdin = p3.stdout, stdout = fo)
            p5 = p4.communicate()
    return fn_out_file



def ends_trimmer(fn, path_out, cut_after_trim='0', len_min=15):
    """
    trim N-bases at either end of the read
    """
    assert os.path.exists(fn)
    assert isinstance(len_min, int)
    path_out = os.path.dirname(fn) if path_out is None else path_out
    assert is_path(path_out)
    fn_out_name = file_prefix(fn)[0] + '.cut.fastq'
    fn_out_file = os.path.join(path_out, fn_out_name)
    # trim either ends
    tmp = cutadapt_cut(cut_after_trim, False)
    if len(tmp) == 2:
        trim_5, trim_3 = tmp
    elif len(tmp) == 1:
        trim_5 = 0 if(int(tmp[0]) < 0) else tmp[0]
        trim_3 = tmp[0] if(int(tmp[0]) < 0) else 0
    else:
        trim_5 = trim_3 = 0
    # print('trim_5: %s' % trim_5)
    # print('trim_3: %s' % trim_3)
    with xopen(fn, 'rt') as fi, open(fn_out_file, 'wt') as fo:
        while True:
            try:
                fq_id, fq_seq, fq_plus, fq_qual = [next(fi).strip(), 
                                                   next(fi).strip(), 
                                                   next(fi).strip(), 
                                                   next(fi).strip(),]
                if len(fq_seq) < len_min + abs(trim_5) + abs(trim_3): 
                    continue # skip short reads
                fq_seq = fq_seq[trim_5:trim_3] if(trim_3 < 0) else fq_seq[trim_5:]
                fq_qual = fq_qual[trim_5:trim_3] if(trim_3 < 0) else fq_qual[trim_5:]
                fo.write('\n'.join([fq_id, fq_seq, fq_plus, fq_qual]) + '\n')
            except StopIteration:
                break
    return fn_out_file


def cutadapt_cut(s, cut_para=True):
    """
    recognize para: cut for cutadapt
    eg: cut=6, cut=-3, cut=6,-3
    """
    if ',' in s:
        n = s.split(',')
        if len(n) > 2:
            raise ValueError('illegal ad_cut: %s' % s)
        else:
            c1, c2 = (int(n[0]), int(n[1]))
            if c1 < 0 or c2 > 0:
                raise ValueError('illegal ad_cut: %s' % s)
        if cut_para:
            c_para = '--cut %s --cut %s' % (c1, c2)
        else:
            c_para = [c1, c2]
    else:
        if cut_para:
            c_para = '--cut %s' % int(s)
        else:
            c_para = [int(s), ]
    return c_para



def se_trimmer(fn, adapter3, path_out=None, adapter5=None, len_min=15, 
               double_trim=False, adapter_sliding=False, qual_min=20, 
               err_rate=0.1, multi_cores=1, rm_untrim=False, overlap=3,
               cut_before_trim=0, overwrite=False):
    """
    using cutadapt to trim SE read
    support one input only
    """
    assert os.path.exists(fn)
    assert isinstance(adapter3, str)
    path_out = os.path.dirname(read_in) if path_out is None else path_out
    assert is_path(path_out)
    fn_out_file = os.path.join(path_out, file_prefix(fn)[0] + '.clean.fastq')
    log_out_file = os.path.join(path_out, file_prefix(fn)[0] + '.cutadapt.log')
    # sliding adapter3
    ads = [adapter3, ]
    if adapter_sliding is True:
        ads.append(adapter_chopper(adapter3))
    para_ad = ' '.join(['-a {}'.format(i) for i in ads])
    # for adapter5
    if isinstance(adapter5, str):
        para_ad += ' -g %s' % adapter5    
    # for untrim
    if rm_untrim is True:
        fn_untrim_file = os.path.join(path_out, 
                                      file_prefix(fn)[0] + '.untrim.fastq')
        para_adx = ' --untrimmed-output=%s --cores=%s' % (fn_untrim_file, 1)
        # untrim, not support multiple threads
    else:
        para_adx = ' --cores=%s' % multi_cores
    # cut adapter *before* trimming
    if not cut_before_trim == 0:
        para_adx = '%s %s' % (para_adx, cutadapt_cut(cut_before_trim))

    ## command line
    if double_trim is True:
        c1 = 'cutadapt %s %s -m %s -q %s --overlap=%s --error-rate=%s \
              --times=1 --trim-n --max-n=0.1 %s' % (para_ad, para_adx, 
                                                    len_min, qual_min, 
                                                    overlap, err_rate, fn)
        c2 = 'cutadapt %s -m %s -q %s --overlap=%s --error-rate=%s --times=1 \
              --cores=%s --trim-n --max-n=0.1 -' % (para_ad, len_min, 
                                                    qual_min, overlap, 
                                                    err_rate, multi_cores)
        if not os.path.isfile(fn_out_file) or overwrite is True:
            with open(log_out_file, 'w') as fo1, open(log_out_file, 'a') as fo2, \
                 open(fn_out_file, 'wt') as fr:
                p1 = subprocess.Popen(shlex.split(c1), stdout=subprocess.PIPE, 
                                      stderr=fo1)
                p2 = subprocess.Popen(shlex.split(c2), stdin=p1.stdout, stdout=fr,
                                      stderr=fo2)
                px = p2.communicate()
                tmp1 = cutadapt_log_parser(log_out_file) # processing log
    else:
        c3 = 'cutadapt %s %s -m %s -q %s --overlap=%s --error-rate=%s \
              --times=1 --trim-n --max-n=0.1 %s' % (para_ad, para_adx, 
                                                    len_min, qual_min, 
                                                    overlap, err_rate, fn)
        if not os.path.isfile(fn_out_file) or overwrite is True:
            with open(log_out_file, 'w') as fo, open(fn_out_file, 'wt') as fr:
                p3 = subprocess.run(shlex.split(c3), stdout=fr, stderr=fo)
                tmp1 = cutadapt_log_parser(log_out_file) # processing log

    return fn_out_file



def trim(fns, adapter3, path_out=None, adapter5=None, len_min=15, 
         adapter_sliding=False,
         double_trim=False, qual_min=20, err_rate=0.1, overlap=3, 
         multi_cores=1, read12=1, rm_untrim=False, rm_dup=False, 
         cut_before_trim='0', cut_after_trim='0', overwrite=False):
    """
    processing reads
    """
    assert isinstance(fns, list)
    assert isinstance(adapter3, str)
    fn_out_files = []
    for fn in fns:
        assert os.path.exists(fn)
        logging.info('trimming reads: %s' % file_prefix(fn)[0])
        if rm_dup is True:
            if cut_after_trim == '0':
                fn_name = file_prefix(fn)[0] + '.clean.nodup.fastq'
            else:
                fn_name = file_prefix(fn)[0] + '.clean.nodup.cut.fastq'
        elif cut_after_trim == '0':
            fn_name = file_prefix(fn)[0] + '.clean.fastq'
        else:
            fn_name = file_prefix(fn)[0] + '.clean.cut.fastq'
        fn_out_check = os.path.join(path_out, fn_name)
        if os.path.exists(fn_out_check) and not overwrite:
            fn_out_files.append(fn_out_check)
            logging.info('file exists: %s' % fn_name)
            continue # next
        fn_out_file = se_trimmer(fn, adapter3, path_out, 
                                 len_min=len_min,
                                 double_trim=double_trim,
                                 adapter_sliding=adapter_sliding,
                                 qual_min=qual_min, 
                                 err_rate=err_rate, 
                                 multi_cores=multi_cores, 
                                 rm_untrim=rm_untrim,
                                 overlap=overlap,
                                 cut_before_trim=cut_before_trim,
                                 overwrite=overwrite)
        ## post-processing
        if rm_dup is True:
            if cut_after_trim == '0':
                f1 = dup_remover(fn_out_file, path_out, q=qual_min)
                os.remove(fn_out_file)
                fn_out_file = f1
            else:
                f2 = dup_remover(fn_out_file, path_out, q=qual_min)
                f3 = ends_trimmer(f2, path_out, 
                                  cut_after_trim=cut_after_trim, 
                                  len_min=len_min)
                os.remove(fn_out_file)                
                os.remove(f2)
                fn_out_file = f3
        else:
            if cut_after_trim == '0':
                pass
            else:
                f4 = ends_trimmer(fn_out_file, path_out, 
                                  cut_after_trim=cut_after_trim, 
                                  len_min=len_min)
                os.remove(fn_out_file)
                fn_out_file = f4
        ## post-stat
        fn_out_files.append(fn_out_file)
        fn_n = int(file_row_counter(fn_out_file) / 4)
        fn_n_file = os.path.join(path_out, file_prefix(fn)[0] + '.reads.txt')
        with open(fn_n_file, "wt") as f:
                f.write(str(fn_n) + '\n')
        ## report
        cutadapt_log = os.path.join(path_out, file_prefix(fn)[0] + '.cutadapt.json')
        with open(cutadapt_log, 'r') as f:
            dd = json.load(f)
        fn_raw = int(dd['raw'])
        logging.info('output: %d of %d (%.1f%%)' % (fn_n, fn_raw, fn_n / fn_raw * 100))
    return fn_out_files



## EOF
