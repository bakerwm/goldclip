#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
HOMER motif analysis
input bed file
"""

__author__ = 'Ming Wang'
__email__ = "wangm08@hotmail.com"
__date__ = "2018-03-21"
__version__ = "0.1"

import os
import sys
import glob
import tempfile
import subprocess
import shlex
import pybedtools
import random_peak # custom package
from goldclip.helper import *

# import logging
# logging.basicConfig(format = '[%(asctime)s] %(message)s', 
#                     datefmt = '%Y-%m-%d %H:%M:%S', 
#                     level = logging.DEBUG)


def get_args():
    """parsing arguments"""
    parser = argparse.ArgumentParser(
        prog = 'bed2motif.homer',
        description = '',
        epilog = 'Example:')
    parser.add_argument('--bed', required = True, metavar = 'BED', 
        type = argparse.FileType('r'),
        help = 'peak files in BED format')
    parser.add_argument('--genome', required = True, metavar = 'GENOME',
        help = 'The reference genome for the BED file.')
    parser.add_argument('--out', required = True, metavar = 'output',
        help = 'The directory to save results.')
    parser.add_argument('--motif_width', default = '4,5,6,7', 
        help = 'The width of motifs, default: 4,5,6,7')
    parser.add_argument('--num', default = 1, 
        help = 'The number of random bed to generate, default: 1')
    parser.add_argument('--tool', default='homer',
        help='Which tool to use for motif analysis, default: HOMER')
    args = parser.parse_args()
    return args



def bed_to_fa(bed, genome):
    """
    convert bed to fasta using bedtools
    """
    assert isinstance(bed, str)
    assert os.path.exists(bed)
    g_fa = Genome(genome).get_fa()
    fa_out = tempfile.mkstemp(suffix = '.fa')[1]
    p = pybedtools.BedTool(bed).sequence(fi=g_fa, fo=fa_out, s=True, name=True)
    return fa_out


def motif_width_formater(w):
    """expect input numbers, 4,5,6"""
    a = [i for i in w.split(r',') if i > 3 and i < 10] # range 3:10
    a = list(set(a)) # unique
    if len(a) > 0:
        b = ','.join(sorted(a))
        return b
    else:
        return None


def get_random_peak(bed, genome, path_out=None, num=1):
    """generate random peaks using random_peak"""
    assert isinstance(bed, str)
    assert os.path.exists(bed)
    gene_bed = Genome(genome).gene_bed()
    gene_rmsk = Genome(genome).gene_bed(rmsk=True)
    if path_out is None:
        path_out = tempfile.TemporaryDirectory()
    assert is_path(path_out)
    c = 'random_peak -i {} -g {} -m {} -o {} -f {} -n {}'.format(bed, 
        gene_bed, gene_rmsk, path_out, 'random', num)
    p = subprocess.run(shlex.split(c), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    random_bed = os.path.join(path_out, 'random1')
    return random_bed


def get_random_peak2(bed, genome, path_out = None, num = 1):
    """get random peaks using bedtools shuffle"""
    gene_bed = Genome(genome).gene_bed()
    random_bed = os.path.join(path_out, 'random.bed')
    b = pybedtools.BedTool(bed).shuffle(genome=genome, incl=gene_bed, 
        chrom=True, seed=1, noOverlapping=True).saveas(random_bed)
    return random_bed


def homer_success(path):
    """check homer.html file"""
    result_html = os.path.join(path, 'homerResults.html')
    result_motif = glob.glob(os.path.join(path, 'homerResults', 'motif*.motif'))
    if os.path.exists(result_html) and len(result_motif) > 0:
        return True
    else:
        return False


def fa_to_motif(fa_bed, fa_bg, path_out, motif_width='4,5,6,7', remove_fa=False):
    """
    de nove motif analysis using Homer
    """
    homer_para = '-p 40 -rna -S 20 -noconvert -nogo'
    c = 'findMotifs.pl {} fasta {} -fasta {} -len {} {}'.format(fa_bed, path_out,
        fa_bg, motif_width, homer_para)
    p = subprocess.run(shlex.split(c), stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)
    if remove_fa:
        [os.remove(f) for f in [fa_bed, fa_bg]]
    return homer_success(path_out)


def bed2motif(bed_in, genome, path_out, motif_width='4,5,6,7'):
    """ Homer motif analysis """
    # randome peak directory
    path_random = tempfile.TemporaryDirectory().name
    assert is_path(path_random)
    bed_bg = get_random_peak(bed_in, genome, path_random, 1)
    fa_bg = bed_to_fa(bed_bg, genome)
    fa_in = bed_to_fa(bed_in, genome)
    if file_row_counter(fa_in) > 0 and file_row_counter(fa_bg):
        s = fa_to_motif(fa_in, fa_bg, path_out, motif_width, remove_fa=False)
        if not s:
            logging.warning('error: Homer failed, exiting...')
    else:
        logging.warning('failed, empty files: %s %s' % (fa_in, fa_bg))


def main():
    args = get_args()
    bed2motif(args.bed.name, args.genome, args.out, args.motif_width)


if __name__ == '__main__':
    main()


## EOF