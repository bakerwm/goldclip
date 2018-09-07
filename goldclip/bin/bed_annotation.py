#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
get annotation of bed
"""

__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__date__ = "2018-03-21"
__version__ = "0.1"


import os
import sys
import re
import argparse
import pathlib
import pybedtools
import pandas as pd
import numpy as np
import goldclip
from goldclip.bin.bed_fixer import Bed_parser

def get_args():
    ## parsing arguments
    parser = argparse.ArgumentParser(
        prog = 'bed_annotation',
        description = 'annotation for bed file',
        epilog = 'Example: ')
    parser.add_argument('-i', required = True, metavar='BED', 
        type = argparse.FileType('r'), help = 'BED file')
    parser.add_argument('-g', required = True, default = 'hg19', 
        metavar = 'GENOME', choices = ['dm3', 'hg19', 'hg38', 'mm10'],
        help = 'Reference genome : dm3, hg19, hg38, mm10')
    parser.add_argument('-t', required = False, default = 'homer', 
        choices = ['homer', 'basic'],
        metavar = 'Type', 
        help = 'Type of the annotation, basic|homer, default [homer]')
    parser.add_argument('--path_data', 
        help='The directory of genome files, default: \
        [$HOME/data/genome/]')
    args = parser.parse_args()
    return args


def anno_picker(genome, group='homer', path_data=None):
    """
    pick annotation files
    genome:
    group: basic, homer
    """
    #anno_dir = os.path.join(goldclip_home, 'data', genome, 'annotation')
    # bin_dir = os.path.join(os.path.dirname(__file__), '..', 'data')
    if path_data is None:
        path_data = os.path.join(pathlib.Path.home(), 'data', 'genome')
    anno_dir = os.path.join(path_data, genome, 'annotation_and_repeats')
    if group == 'basic':
        group_basic = [
            'genicPiRNA',
            'nonGenicPiRNA',
            'TE',
            'tRNA',
            'rRNA',
            'miRNA',
            'targetscan',
            'sncRNA',
            '3u',
            '5u',
            'exon',
            'intron', 
            'igr']
        anno = [os.path.join(anno_dir, 'basic', genome + '.'  + i + '.bed') for i in group_basic]
    if group == 'homer':
        group_homer = [
            'tts', 
            'rRNA',
            'pseudo', 
            'promoters', 
            'ncRNA', 
            'utr3',
            'utr5', 
            'coding', 
            'introns', 
            'intergenic']
        anno = [os.path.join(anno_dir, 'homer', 'ann.' + i + '.bed') for i in group_homer]
    # retrun only exists
    # anno = [f for f in anno if os.path.exists(f)]
    anno_fixed = []
    for f in anno:
        if os.path.exists(f):
            f_tmp = f + '.tmp'
            if not os.path.exists(f_tmp):
                os.replace(f, f_tmp)
                Bed_parser(f_tmp).bed_fixer().saveas(f)
            anno_fixed.append(f)
    if len(anno_fixed) == 0:
        raise ValueError('illegeal annotation files: %s' % anno_dir)
        # sys.exit('annotation bed files not found: %s' % genome)
    return anno_fixed


def bed_annotator(bed_in, genome, group='homer', path_data=None):
    """
    count reads in annos
    """
    if path_data is None:
        path_data = os.path.join(pathlib.Path.home(), 'data', 'genome')
    bed_prefix = os.path.splitext(os.path.basename(bed_in))[0]
    annos = anno_picker(genome, group, path_data=path_data)
    df = pd.DataFrame(columns=[bed_prefix])
    a = pybedtools.BedTool(bed_in)
    a_cnt = a.count() # all
    a_bed6 = True if a.to_dataframe().shape[1] == 6 else False # bed6 or bed3
    for n in annos:
        group = os.path.basename(n).split(r'.')[-2]
        b = pybedtools.BedTool(n)
        if a_bed6 is True:
            a_not_b = a.intersect(b, v=True, s=True)
        else:
            a_not_b = a.intersect(b, v=True)
        df.loc[group] = a.count() - a_not_b.count()
        a = a_not_b # convert to a
    # not in group
    df.loc['other'] = a_cnt - df.sum(axis=0)
    df['sample'] = bed_prefix
    return df


def main():
    args = get_args()
    path_data = args.path_data
    df = bed_annotator(args.i.name, args.g, args.t, path_data)
    print(df)



if __name__ == '__main__':
    main()



## EOF
