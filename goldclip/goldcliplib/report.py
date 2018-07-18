#!/usr/bin/env python
"""
make a summary of the project

wrapper output of goldclip pipeline 
create stats and plots

## figure1.reads_mapping_stat.pdf
number of reads: raw, clean, no_dup, spikein, genome, unmap, ...

## figure2.reads_annotation_stat.pdf
number of reads: RNA categories, ...

## figure3.reads_correlation.pdf

## figure4.rtstop_correlation.pdf

## figure5.peak_number.pdf

## figure6.peak_length.pdf

## figure7.peak_annotation.pdf

## figure8.peak_motif.pdf

## figure9.peak_conservation.pdf

## figure10.hexmer_zscore.pdf

## figure11.peak_overlap.pdf

functions
"""

__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__date__ = "2018-03-21"
__version__ = "0.1"



import os
import sys
import re
import argparse
import json
import tempfile
import string
import random
import pybedtools
import pandas as pd
import numpy as np
import multiprocessing as mp
from goldclip.helper import *
from goldclip.bin.bed_annotation import *
from goldclip.bin.bed_fixer import Bed_parser
from goldclip.bin.bed_motif import *
from goldclip.goldcliplib.log_parser import *


def figure1_read_map(path, smp_name):
    """
    trim, PCR_dup
    mapped reads
    """
    logging.info('figure1 reads mapping')
    path_map = os.path.join(path, 'results', 'read_mapping')
    if not os.path.exists(path_map):
        os.makedirs(path_map)
    map_txt = os.path.join(path_map, 'read_mapping.txt')
    path_trim = os.path.join(path, 'input_reads')
    path_map = os.path.join(path, 'genome_mapping')
    df1 = trim_wrapper(path_trim, smp_name)
    df2 = map_wrapper(path_map, smp_name)
    # merge
    df = pd.concat([df1, df2], axis = 0)
    df.to_csv(map_txt, '\t', header = True, index = True)
    return df



def figure2_read_anno(path, genome, group = 'homer'):
    """
    genome mapped reads annotation
    """
    logging.info('figure2 reads annotation')
    path_anno = os.path.join(path, 'results', 'read_annotation')
    if not os.path.exists(path_anno):
        os.makedirs(path_anno)
    bed_files = locate_bam_files(path, bam_out = False)
    if len(bed_files) > 0:
        df = _bed_anno(bed_files, genome, group)
    anno_txt = os.path.join(path_anno, 'read_annotation.txt')
    df.to_csv(anno_txt, '\t', header = True, index = True)
    return(df)



def figure3_read_cor(path):
    logging.info('figure3 reads correlation')
    path_read_corr = os.path.join(path, 'results', 'read_correlation')
    if not os.path.exists(path_read_corr):
        os.makedirs(path_read_corr)
    bam_files = locate_bam_files(path, bam_out = True, rep_only = True)
    df = bam_corr(bam_files[0], bam_files[1], path_read_corr)
    return df


def figure4_rt_cor(path):
    """RTStop correlation between replicates"""
    logging.info('figure4 rtstop correlation')
    rt_dirs = next(os.walk(os.path.join(path, 'rtstops')))[1]
    rt_dirs = [i for i in rt_dirs if not re.search('_rep[0-9](_\w?\d?)?$', i)] #exclude rep
    rt_beds = [os.path.join(path, 'rtstops', d, d + '.RTStop.bed') for d in rt_dirs]
    if not len(rt_beds) == 1:
        sys.exit('only ONE rtstop file required: ' + ','.join(rt_beds))
    df = bed_parser(rt_beds[0])
    df2 = df.drop(df.columns[:6].tolist() + df.columns[-3:].tolist(), axis = 1)
    cor_ma = df2.apply(pd.to_numeric).corr(method = 'pearson')
    # save & report
    path_rt_corr = os.path.join(path, 'results', 'rtstop_correlation')
    if not os.path.exists(path_rt_corr):
        os.makedirs(path_rt_corr)
    cor_txt = os.path.join(path_rt_corr, 'cor_matrix.tab')
    cor_ma.to_csv(cor_txt, '\t', header = True, index = True)
    return cor_ma


def figure5_peak_num(path):
    """peak number"""
    logging.info('figure5 peak number')
    peak_files = glob.glob(os.path.join(path, 'peaks', '*', '*', "*.fixed.bed"))
    df = pd.DataFrame(columns = list('ABC'))
    for b in peak_files:
        t = b.split(r'/')
        tool, sample = t[-3:-1]
        cnt = file_row_counter(b)
        df2 = pd.DataFrame([[tool, sample, cnt]], columns = list('ABC'))
        df = df.append(df2, ignore_index = True)
    df.columns = ['tool', 'sample', 'count']
    # save & report
    path_peak_num = os.path.join(path, 'results', 'peak_number')
    if not os.path.exists(path_peak_num):
        os.makedirs(path_peak_num)
    peak_num_txt = os.path.join(path_peak_num, 'peak_number.txt')
    df.to_csv(peak_num_txt, '\t', header = True, index = True)
    return df


def figure6_peak_len(path):
    """peak length"""
    logging.info('figure6 peak length')
    peak_files = glob.glob(os.path.join(path, 'peaks', '*', '*', "*.fixed.bed"))
    # df = pd.DataFrame(columns = list('ABC'))
    a = []
    for b in peak_files:
        t = b.split(r'/')
        tool, sample = t[-3:-1]
        df = bed_parser(b)
        df2 = pd.DataFrame({'A': tool, 'B': sample,
            'C': df['end'] - df['start']})
        # df = df.append(df2, ignore_index = True)
        a.append(df2)
    df = pd.concat(a, axis = 0)
    df.columns = ['tool', 'sample', 'length']
    # save & report
    path_peak_len = os.path.join(path, 'results', 'peak_length')
    if not os.path.exists(path_peak_len):
        os.makedirs(path_peak_len)
    peak_len_txt = os.path.join(path_peak_len, 'peak_length.txt')
    df.to_csv(peak_len_txt, '\t', header = True, index = True)
    return df


def figure7_peak_anno(path, genome):
    """peak annotation"""
    logging.info('figure7 peak annotation')
    peak_files = glob.glob(os.path.join(path, 'peaks', '*', '*', "*.fixed.bed"))
    df = pd.DataFrame(columns = ['tool', 'sample', 'count'])
    for b in peak_files:
        t = b.split(r'/')
        tool, sample = t[-3:-1]
        dn = _bed_anno([b,], genome, group = 'homer')
        dn.columns = ['count']
        dn.insert(0, 'sample', sample)
        dn.insert(0, 'tool', tool)
        df = df.append(dn)
    # add group to column
    df.insert(0, 'group', df.index.tolist())
    # save & report
    path_peak_anno = os.path.join(path, 'results', 'peak_annotation')
    if not os.path.exists(path_peak_anno):
        os.makedirs(path_peak_anno)
    peak_anno_txt = os.path.join(path_peak_anno, 'peak_annotation.txt')
    df.to_csv(peak_anno_txt, '\t', header = True, index = True)
    return df


def figure8_motif_analysis(path, genome):
    """peak motif analysis, de novo analysis"""
    logging.info('figure8 motif analysis')
    peak_files = glob.glob(os.path.join(path, 'peaks', '*', '*', '*.fixed.bed'))
    path_peak_motif = os.path.join(path, 'results', 'peak_motif')
    if not os.path.exists(path_peak_motif):
        os.makedirs(path_peak_motif)
    for b in peak_files:
        t = b.split(r'/')
        tool, smp_id, smp_name = t[-3:]
        b_prefix = re.sub('.fixed.bed', '', smp_name)
        logging.info('motif - ' + tool + ' : ' + smp_id)
        b_path = os.path.join(path_peak_motif, tool, smp_id)
        bed2motif(b, genome, b_path, '4,5,6,7,8')



def figure9_peak_conservation(path, genome):
    logging.info('figure9 peak conservation')
    peak_files = glob.glob(os.path.join(path, 'peaks', '*', '*', "*.fixed.bed"))
    path_peak_conservation = os.path.join(path, 'results', 'peak_conservation')
    if not os.path.exists(path_peak_conservation):
        os.makedirs(path_peak_conservation)
    for b in peak_files:
        if file_row_counter(b) == 0:
            continue
        t = b.split(r'/')
        tool, smp_id, smp_name = t[-3:]
        b_con_name = re.sub(r'.bed', '', smp_name)
        path_b_con = os.path.join(path_peak_conservation, tool, smp_id)
        if not os.path.exists(path_b_con):
            os.makedirs(path_b_con)
        b_con = os.path.join(path_b_con, b_con_name)
        # bed_conservation(b, b_con + '.con.bed', 0, 'hg19')
        # bed_conservation(b, b_con + '.ext500.con.bed', 500, 'hg19')
        bed_conservation(b, b_con + '.ext1k.con.bed', 1000, genome)



def goldclip_report(path_out, smp_name, genome):
    """ make summary and report """
    ## figure
    df = figure1_read_map(path_out, smp_name)
    # df = figure2_read_anno(path_out, genome, 'homer')
    # df = figure3_read_cor(path_out)
    # df = figure4_rt_cor(path_out)
    # df = figure5_peak_num(path_out)
    # df = figure6_peak_len(path_out)
    # df = figure7_peak_anno(path_out, genome)
    # figure8_motif_analysis(path_out, genome)
    # figure9_peak_conservation(path_out, genome)
