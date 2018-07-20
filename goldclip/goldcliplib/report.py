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
import re
import random
import logging
import pandas as pd
import numpy as np
from goldclip.helper import *
from goldclip.bin.bed_fixer import Bed_parser
from goldclip.bin.bed_motif import *
from goldclip.bin.bed_annotation import *
from goldclip.goldcliplib.log_parser import *
# import sys
# import argparse
# import json
# import tempfile
# import string
# import pybedtools
# import multiprocessing as mp



class Goldclip_output:
    """
    record the output of Goldclip output
    """

    def __init__(self, path, smp_name, **kwargs):
        self.project_path = path
        self.project_name = smp_name


    def get_trim_stat(self):
        """
        return the reads processing of fastq files
        groups: raw, too_short, PCR_dup, no_dup
        """
        path_trim = os.path.join(self.project_path, 'input_reads')
        dfx = [] # list of pd.DataFrame
        with os.scandir(path_trim) as it:
            for entry in it:
                if not entry.name.endswith('.cutadapt.json'):
                    continue
                name = re.sub(r'.cutadapt.json', '', entry.name)
                fn = os.path.join(path_trim, entry.name)
                dd = json_reader(fn) # 
                # clean
                fn_nodup = os.path.join(path_trim, name + '.reads.txt')
                with open(fn_nodup) as fi:
                    dd['nodup'] = next(fi).strip()
                dd['too_short'] = int(dd['raw']) - int(dd['clean'])
                dd['dup'] = int(dd['clean']) - int(dd['nodup'])
                dx = pd.DataFrame({'group': ['raw', 'too_short', 'PCR_dup', 
                                             'no_dup'],
                                    name: [dd['raw'], dd['too_short'], dd['dup'],
                                         dd['nodup']]})
                dx.set_index('group', inplace=True)
                dfx.append(dx)
        df = pd.concat(dfx, axis=1)
        df = df.apply(pd.to_numeric)
        df.insert(0, self.project_name, df.sum(axis=1))
        return df


    def get_map_stat(self):
        """
        return the reads mapping 
        groups: spikein, MT_trRNA, genome
        """
        path_map = os.path.join(self.project_path, 'genome_mapping')
        dfx = []
        with os.scandir(path_map) as it:
            for entry in it:
                if not entry.name.endswith('.mapping_stat.csv'):
                    continue
                name = re.sub(r'.mapping_stat.csv', '', entry.name)
                fn = os.path.join(path_map, entry.name)
                # skip merged
                if name == self.project_name:
                    continue
                dx1 = pd.read_csv(fn, ',').filter(items=['group', 'read'])
                dx1.set_index('group', inplace=True)
                dx2 = dx1.rename(columns={'read': name})
                dfx.append(dx2)
        df = pd.concat(dfx, axis=1)
        df.insert(0, self.project_name, df.sum(axis=1))
        return df


    def get_bam_file(self, bam2bed=False, rep_only=False, merge_only=False):
        bam_files = []
        path_map = os.path.join(self.project_path, 'genome_mapping')
        with os.scandir(path_map) as it:
            for entry in it:
                if merge_only and not entry.name == self.project_name:
                    continue
                elif rep_only and entry.name == self.project_name:
                    continue
                else:
                    pass
                fx_ext = '.bed' if bam2bed else '.bam'
                fx = os.path.join(path_map, entry.name, entry.name + fx_ext)
                if os.path.exists(fx):
                    bam_files.append(fx)
        return bam_files


    def get_peak_file(self, rep_only=False, merge_only=False):
        peak_files = []
        # tools / samples
        path_peak = os.path.join(self.project_path, 'peaks')
        with os.scandir(path_peak) as tools:
            for tool in tools:
                path_tool = os.path.join(path_peak, tool.name)
                with os.scandir(path_tool) as smps:
                    for smp in smps:
                        if merge_only and not smp.name == self.project_name:
                            continue
                        elif rep_only and smp.name == self.project_name:
                            continue
                        else:
                            pass
                        fn = os.path.join(path_tool, smp.name, 
                                          smp.name + '.fixed.bed')
                        if os.path.exists(fn):
                            peak_files.append(fn)
        return peak_files


    def get_rtstop_file(self, rep_only=False, merge_only=False, rt_reads=False):
        rtstop_files = []
        path_rtstop = os.path.join(self.project_path, 'rtstops')
        with os.scandir(path_rtstop) as it:
            for entry in it:
                if merge_only and not entry.name == self.project_name:
                    continue
                elif rep_only and entry.name == self.project_name:
                    continue
                else:
                    pass
                fn_ext = '.RTRead.bed' if rt_reads else '.RTStop.bed'
                fn = os.path.join(path_rtstop, entry.name, entry.name + fn_ext)
                if os.path.exists(fn):
                    rtstop_files.append(fn)
        return rtstop_files






def figure1_read_map(path, smp_name):
    """
    trim, PCR_dup
    mapped reads
    """
    logging.info('figure1 reads mapping')
    figure1_path = os.path.join(path, 'results', 'read_mapping')
    figure1_txt = os.path.join(figure1_path, 'read_mapping.txt')
    assert os.path.exists(path)
    assert is_path(figure1_path)
    assert isinstance(smp_name, str)
    df1 = Goldclip_output(path, smp_name).get_trim_stat()
    df2 = Goldclip_output(path, smp_name).get_map_stat()
    df = pd.concat([df1, df2], axis=0, sort=False)
    df.to_csv(figure1_txt, '\t', header=True, index=True)
    return df


def figure2_read_anno(path, smp_name, genome, group='homer'):
    """
    genome mapped reads annotation
    # function
    df = bed_annotator(args.i.name, args.g, args.t, path_data)
    """
    logging.info('figure2 reads annotation')
    figure2_path = os.path.join(path, 'results', 'read_annotation')
    figure2_txt = os.path.join(figure2_path, 'read_annotation.txt')
    assert os.path.exists(path)
    assert is_path(figure2_path)
    assert isinstance(smp_name, str)
    bed_files = Goldclip_output(path, smp_name).get_bam_file(bam2bed=True)
    if len(bed_files) == 0:
        logging.error('failed, bed files not found: %s' % path)
        return None
    # annotate files
    dfx = [bed_annotator(bed, genome, group) for bed in bed_files]
    df = pd.concat(dfx, axis=1)
    df = df.reset_index()
    df = df.rename(columns={'index': 'type'})
    df.to_csv(figure2_txt, '\t', header=True, index=False)
    return(df)


def figure3_read_cor(path, smp_name, window=10000, multi_cores=8):
    logging.info('figure3 reads correlation')
    figure3_path = os.path.join(path, 'results', 'read_correlation')
    assert os.path.exists(path)
    assert is_path(figure3_path)
    assert isinstance(smp_name, str)
    bam_files = Goldclip_output(path, smp_name).get_bam_file(rep_only=True)
    if len(bam_files) < 2:
        logging.error('skipped, at least 2 bam files required: %s' % path)
        return None
    df = bam_corr(bam_files, figure3_path, window=window, multi_cores=multi_cores)
    return df


def figure4_rt_cor(path, smp_name):
    """RTStop correlation between replicates"""
    logging.info('figure4 rtstop correlation')
    figure4_path = os.path.join(path, 'results', 'rtstop_correlation')
    figure4_txt = os.path.join(figure4_path, 'cor_matrix.tab')
    assert os.path.exists(path)
    assert is_path(figure4_path)
    assert isinstance(smp_name, str)
    bed_merge = Goldclip_output(path, smp_name).get_rtstop_file(merge_only=True)
    if not os.path.exists(bed_merge[0]):
        logging.error('file not exists: %s' % bed_merge)
        return None
    df = Bed_parser(bed_merge[0]).bed
    # remove the first 6 columns, last 3 columns
    dfx = df.drop(df.columns[:6].tolist() + df.columns[-3:].tolist(), axis=1) #
    dfx_cor = dfx.apply(pd.to_numeric).corr(method='pearson')
    dfx_cor.to_csv(figure4_txt, '\t', header=True, index=False)
    return dfx_cor


def figure5_peak_num(path, smp_name):
    """peak number"""
    logging.info('figure5 peak number')
    figure5_path = os.path.join(path, 'results', 'peak_number')
    figure5_txt = os.path.join(figure5_path, 'peak_number.txt')
    assert os.path.exists(path)
    assert is_path(figure5_path)
    assert isinstance(smp_name, str)
    peak_files = Goldclip_output(path, smp_name).get_peak_file()
    df = pd.DataFrame(columns = ['tool', 'id', 'count'])
    for p in peak_files:
        sep = p.split('/') # only fo Linux path
        tool, sample = sep[-3:-1]
        num = file_row_counter(p)
        dfx = pd.DataFrame([[tool, sample, num]], 
                           columns=['tool', 'id', 'count'])
        df = df.append(dfx, ignore_index=True)
    df.to_csv(figure5_txt, '\t', header=True, index=False)
    return df


def figure6_peak_len(path, smp_name):
    """peak length"""
    logging.info('figure6 peak length')
    figure6_path = os.path.join(path, 'results', 'peak_length')
    figure6_txt = os.path.join(figure6_path, 'peak_length.txt')
    assert os.path.exists(path)
    assert is_path(figure6_path)
    assert isinstance(smp_name, str)
    peak_files = Goldclip_output(path, smp_name).get_peak_file()
    df = pd.DataFrame(columns = ['length', 'tool', 'id', 'count'])
    for p in peak_files:
        sep = p.split('/') # only for Linux path
        tool, sample = sep[-3:-1]
        dx = Bed_parser(p).bed
        dx2 = pd.DataFrame({'length': dx.end - dx.start,
                            'id': sample})
        s1 = dx2.groupby('length')['id'].nunique()
        dx3 = pd.DataFrame({'tool': tool, 
                            'id': sample,
                            'count': s1})
        dx3 = dx3.reset_index()
        dx3 = dx3.rename({'index': 'length'})
        df = df.append(dx3, ignore_index=True, sort=False)
    df.to_csv(figure6_txt, '\t', header=True, index=False)
    return df


def figure7_peak_anno(path, smp_name, genome, group='homer'):
    """peak annotation"""
    logging.info('figure7 peak annotation')
    figure7_path = os.path.join(path, 'results', 'peak_annotation')
    figure7_txt = os.path.join(figure7_path, 'peak_annotation.txt')
    assert os.path.exists(path)
    assert is_path(figure7_path)
    assert isinstance(smp_name, str)
    # peak_files = glob.glob(os.path.join(path, 'peaks', '*', '*', "*.fixed.bed"))
    peak_files = Goldclip_output(path, smp_name).get_peak_file()
    df = pd.DataFrame(columns = ['type', 'count', 'id', 'tool'])
    for p in peak_files:
        sep = p.split('/') # only for Linux path
        tool, sample = sep[-3:-1]
        dx = bed_annotator(p, genome, group) # index, id
        dx = dx.reset_index()
        dx.columns = ['type', 'count']
        dx['id'] = sample
        dx['tool'] = tool
        df = df.append(dx, ignore_index=True, sort=False)
    df.to_csv(figure7_txt, '\t', header=True, index=False)
    return df



def figure8_motif_analysis(path, smp_name, genome):
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



def figure9_peak_conservation(path, smp_name, genome):
    logging.info('figure9 peak conservation')
    figure9_path = os.path.join(path, 'results', 'peak_conservation')
    assert os.path.exists(path)
    assert is_path(figure9_path)
    assert isinstance(smp_name, str)
    peak_files = Goldclip_output(path, smp_name).get_peak_file()
    # peak_files = glob.glob(os.path.join(path, 'peaks', '*', '*', "*.fixed.bed"))
    # path_peak_conservation = os.path.join(path, 'results', 'peak_conservation')
    # if not os.path.exists(path_peak_conservation):
        # os.makedirs(path_peak_conservation)
    for peak in peak_files:
        if file_row_counter(peak) == 0:
            continue
        sep = peak.split('/') # Linux path
        tool, smp_id, smp_name = sep[-3:]
        peak_sub_path = os.path.join(figure9_path, tool, smp_id)
        assert is_path(peak_sub_path)
        peak_sub_prefix = re.sub('.bed$', '', smp_name)
        peak_sub_file = os.path.join(peak_sub_path, peak_sub_prefix)
        bed_conservation(peak, peak_sub_file + '.ext1k.con.bed', 1000, genome)


        # t = b.split(r'/')
        # tool, smp_id, smp_name = t[-3:]
        # b_con_name = re.sub(r'.bed', '', smp_name)
        # path_b_con = os.path.join(path_peak_conservation, tool, smp_id)
        # if not os.path.exists(path_b_con):
        #     os.makedirs(path_b_con)
        # b_con = os.path.join(path_b_con, b_con_name)
        # # bed_conservation(b, b_con + '.con.bed', 0, 'hg19')
        # # bed_conservation(b, b_con + '.ext500.con.bed', 500, 'hg19')
        # bed_conservation(b, b_con + '.ext1k.con.bed', 1000, genome)



def goldclip_report(path_out, smp_name, genome):
    """ make summary and report """
    ## figure
    # df = figure1_read_map(path_out, smp_name)
    # df = figure2_read_anno(path_out, smp_name, genome, 'homer')
    # df = figure3_read_cor(path_out, smp_name)
    # df = figure4_rt_cor(path_out, smp_name)
    # df = figure5_peak_num(path_out, smp_name)
    # df = figure6_peak_len(path_out, smp_name)
    # df = figure7_peak_anno(path_out, smp_name, genome)
    # figure8_motif_analysis(path_out, smp_name, genome)
    figure9_peak_conservation(path_out, smp_name, genome)
