#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
parsing arguments for sub-commands
"""

__author__ = 'Ming Wang <wangm08@hotmail.com>'
__copyright__ = '2018 by Ming Wang <wangm08@hotmail.com>'
__license__ = 'MIT'
__email__ = 'wangm08@hotmail.com'
__version__ = '0.0.1'


import argparse


def add_demx_args(parser):
    """
    processing GoldCLIP illumina datasets
    only one of the PE reads
    ## type1: goldclip_version_1
    read1: {NNN} - {bc} - {NN} - <insert>

    ## type2: goldclip_version_2
    read1: {N10} - <insert> - A{barcode}
    read2: {barcode}A - <insert> - {N10}
    
    ## type3: eCLIP
    read1: {barcode} - <insert> - {N10}
    read2: {N10} - <insert> - {bracode}
    """
    parser.add_argument('--fq1', required=True, metavar='READ1', 
        type=argparse.FileType('r'),
        help='Read file contains barcode at the 5-prime end, in FASTQ format')
    parser.add_argument('--fq2', default=None, 
        metavar='READ2', type=argparse.FileType('r'),
        help='Another read of PE, does not contain barcodes (optional)')
    parser.add_argument('--bc_file', required=True, metavar='BARCODE', 
        type=argparse.FileType('r'),
        help='file of barcodes, only demultiplex P7 or barcode: \
        2-column <barcode> <sample_name> \
        if both P7 and barcode, require 3-column:\
        <p7-index> <barcode> <sample_name>')
    parser.add_argument('--out', required=True, metavar='OUTPUT', 
        help='The directory to save results')
    parser.add_argument('--p7_and_bc', action='store_true',
        help='if specified, demx P7 and barcode at the same time, \
        ignore --p7')
    parser.add_argument('--p7', action='store_true',
        help = 'if specified, demx P7 index, in fastq comment-field')
    parser.add_argument('--n_left', default=3, metavar='N-LEFT', type=int,
        help='Number of randomers at the left of read. default: 3')
    parser.add_argument('--n_right', default=2, metavar='N-RIGHT', type=int,
        help='Number of randomers at the right of read. default: 2')
    parser.add_argument('--n_mismatch', default=0, metavar = 'Mismatches', 
        type = int,
        help = 'Number of mismatches allowed in barcode/index, default: 0')
    parser.add_argument('--cut', action='store_true',
        help = 'if specified, cut the barcode from read')
    parser.add_argument('--bioawk', action='store_true',
        help='use bioawk to demulplex, only support one read')

    return parser


def add_trim_args(parser):
    """
    processing SE read 
    - remove 3' adapter(s) (default: TruSeq RNA-Seq)
    - remove 5' adatper
    - trim low-quality bases on both 5 and 3 end
    - trim N reads
    - cut N-bases at either end of read
    """
    parser.add_argument('-i', nargs='+', required=True, metavar='file', 
        type=argparse.FileType('r'),
        help='reads in FASTQ format, support (*.gz), 1-4 files.')
    parser.add_argument('-a', 
        default=None, metavar='adapter', type=str,
        help='Adapter sequence in the 3 prime end of read, \
        default [AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG].')
    parser.add_argument('-o', default=None, metavar='output', 
        help='The directory to save results.')
    parser.add_argument('--read12', type=int, default=1, metavar='read12',
        help='which one of PE reads, 1=read1, 2=read2, default: 1')
    parser.add_argument('-m', default=15, metavar='len_min', 
        type=int, help='Minimum length of reads after trimming, defualt [15]')
    parser.add_argument('-p', default=80, metavar='percent', 
        type=int,
        help='minimum percent of bases that must have -q quality, default [80]')
    parser.add_argument('-q', default=20, metavar='quality', 
        type=int,
        help='The cutoff of base quality, default [20]')    
    parser.add_argument('-e', default=0.1, metavar='err_rate', 
        type=float,
        help='Maximum allowed error rate, default [0.1]')
    parser.add_argument('-O', default=1, metavar='overlap',
        help='Required N bases overlap between reads and adapter, default [1]')
    parser.add_argument('--rm_untrim', action='store_false',
        help='if specified, discard reads without adapter')
    parser.add_argument('--rm_dup', action='store_true',
        help='if specified, remove duplicated reads' )
    parser.add_argument('--cut_before_trim', default='0', metavar='cut1', 
        help='cut bases before trimming adapter, Number of bases to cut from each \
        read, plus on 5-prime end, minus on 3-prime end, could be \
        single, or double numbers, eg: 3 or -4 or 3,-4, default [0]')
    parser.add_argument('--cut_after_trim', default='0', metavar='cut2', 
        help='cut bases after trimming adapter, Number of bases to cut from each \
        read, plus on 5-prime end, minus on 3-prime end, , could be \
        single, or double numbers, eg: 3 or -4 or 3,-4, default [0]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    parser.add_argument('--threads', default=1, 
        metavar='threads', type=int,
        help='Number of threads to launch, default [1]')

    return parser


def add_map_args(parser):
    """
    Mapping SE read or one of PE reads to reference genome
    using bowtie, STAR, ... (universal)
    """
    parser.add_argument('-i', nargs='+', required=True, metavar='INPUT', 
        type=argparse.FileType('r'),
        help='CLIP reads in FASTQ format, (not *.gz), 1-4 files.')
    parser.add_argument('-n', required=True, metavar='NAME',
        help='Name of the experiment')
    parser.add_argument('-o', default=None, 
        metavar='OUTPUT',  help='The directory to save results.')
    parser.add_argument('-g', default='hg19', 
        metavar='GENOME', choices=['dm3', 'hg19', 'hg38', 'mm10'],
        help='Reference genome : dm3, hg19, GRCh38, mm10, GRCm38')
    parser.add_argument('-k', default='hg19', 
        metavar='Spike-in', choices=['dm3', 'hg19', 'hg38', 'mm10'],
        help='Spike-in genome : dm3, hg19, hg38, mm10')
    parser.add_argument('--threads', default=1, 
        metavar='THREADS', type=int, 
        help='Number of threads to launch, default [1].')
    parser.add_argument('--aligner', default='bowtie', 
        choices=['bowtie', 'bowtie2', 'star'],
        help='Choose which aligner to use. default: bowtie')
    return parser



def add_peak_args(parser):
    """
    call peaks using clipper, pyicoclip
    """
    parser.add_argument('-tool', required=True, metavar='TOOL', 
        choices=['clipper', 'pyicoclip'], 
        help='Peak-caller, clipper|pyicoclip')
    parser.add_argument('-i', nargs='+', required=True, metavar='BAM', 
        type=argparse.FileType('r'),
        help='mapped BAM files, sorted, 1-4 files.')
    parser.add_argument('-n', required=True, metavar='NAME',
        help='Name of the experiment')
    parser.add_argument('-g', default='GRCh38', 
        metavar='GENOME', choices=['dm3', 'hg19', 'GRCh38', 'mm10', 'GRCm38'],
        help='Reference genome : hg19, GRCh38, mm10, GRCm38')
    parser.add_argument('-o', default=None, 
        metavar='OUTPUT', help='The directory to save results.')
    parser.add_argument('--threads', default=1, 
        metavar='THREAD', type=int,
        help='Number of threads to launch, default [1].')

    return parser



def add_rtstop_args(parser):
    """
    call RT-Stops
    """
    pass



def add_report_args(parser):
    """
    create report, advance analysis based on previous results
    """
    pass


