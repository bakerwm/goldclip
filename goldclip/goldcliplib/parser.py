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
    parser.add_argument('--read1', required=True, metavar='READ1', 
        type=argparse.FileType('r'), dest='fq1',
        help='Read1 of PE reads, FASTQ format, support *.gz')
    parser.add_argument('--read2', default=None, metavar='READ2', 
        type=argparse.FileType('r'), dest='fq2',
        help='Read2 of PE reads, FASTQ format, support *.gz, (optional)')
    parser.add_argument('--bc-file', required=True, metavar='BARCODE', 
        dest='bc_file', type=argparse.FileType('r'),
        help='file of barcodes, only demultiplex P7 or barcode: \
        2-column <barcode> <sample_name> \
        if both P7 and barcode, require 3-column:\
        <p7-index> <barcode> <sample_name>')
    parser.add_argument('--out', required=True, metavar='OUTPUT', 
        help='The directory to save results')
    parser.add_argument('--bc-in-read12', default=1, metavar='READ12',
        dest='bc_in_read12', choices=[1, 2], type=int,
        help='Barcode located in read1 or read2, [1, 2], default: 1')
    parser.add_argument('--p7-and-bc', action='store_true',
        dest='p7_and_bc',
        help='if specified, demx P7 and barcode at the same time, \
        ignore --p7, --bc')
    parser.add_argument('--bc', action='store_true', dest='bc_only',
        help='if specified, demx barcode only.')
    parser.add_argument('--p7', action='store_true', dest='p7_only',
        help = 'if specified, demx P7 index, in fastq comment-field')
    parser.add_argument('--n-left', default=3, metavar='N-LEFT', type=int,
        dest='n_left',
        help='Number of randomers at the left of read. default: 3')
    parser.add_argument('--n-right', default=2, metavar='N-RIGHT', type=int,
        dest='n_right',
        help='Number of randomers at the right of read. default: 2')
    parser.add_argument('--n-mismatch', default=0, metavar = 'Mismatches', 
        dest='n_mismatch', type = int,
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
    parser.add_argument('--rm-untrim', action='store_false', dest='rm_untrim',
        help='if specified, discard reads without adapter')
    parser.add_argument('--rm-dup', action='store_true', dest='rm_dup',
        help='if specified, remove duplicated reads' )
    parser.add_argument('--cut-before-trim', default='0', metavar='cut1', 
        dest='cut_before_trim',
        help='cut bases before trimming adapter, Number of bases to cut from each \
        read, plus on 5-prime end, minus on 3-prime end, could be \
        single, or double numbers, eg: 3 or -4 or 3,-4, default [0]')
    parser.add_argument('--cut-after-trim', default='0', metavar='cut2', 
        dest='cut_after_trim',
        help='cut bases after trimming adapter, Number of bases to cut from each \
        read, plus on 5-prime end, minus on 3-prime end, , could be \
        single, or double numbers, eg: 3 or -4 or 3,-4, default [0]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    parser.add_argument('--threads', default=1, 
        metavar='threads', type=int,
        help='Number of threads to launch, default [1]')
    return parser



def add_align_args(parser):
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
        metavar='OUTPUT',  help='The directory to save results, default, \
        the same as input fastq files.')
    parser.add_argument('-g', required=True, default='hg19', 
        metavar='GENOME', choices=['dm3', 'hg19', 'hg38', 'mm10'],
        help='Reference genome : dm3, hg19, hg39, mm10, default: hg19')
    parser.add_argument('-k', default='hg19', 
        metavar='Spike-in', choices=['dm3', 'hg19', 'hg38', 'mm10'],
        help='Spike-in genome : dm3, hg19, hg38, mm10, default: None')
    parser.add_argument('--unique-only', action='store_true', 
        dest='unique_only',
        help='if specified, keep unique mapped reads only')
    parser.add_argument('--align-to-rRNA', action='store_true',
        dest='align_to_rRNA', 
        help='if specified, remove rRNA before mapping to reference genome')
    parser.add_argument('--threads', default=1, 
        metavar='THREADS', type=int, 
        help='Number of threads to launch, default: 1.')
    parser.add_argument('--aligner', default='bowtie', 
        choices=['bowtie', 'bowtie2', 'star'],
        help='Choose which aligner to use. default: bowtie')
    parser.add_argument('--path_data', 
        help='The directory of genome files, default: \
        [$HOME/data/genome/]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    return parser



def add_peak_args(parser):
    """
    call peaks using clipper, pyicoclip
    """
    parser.add_argument('--tool', required=True, metavar='TOOL', 
        choices=['clipper', 'pyicoclip'], 
        help='Peak-caller, clipper|pyicoclip')
    parser.add_argument('-i', nargs='+', required=True, metavar='BAM', 
        type=argparse.FileType('r'),
        help='BAM files, sorted, 1-4 files.')
    parser.add_argument('-g', default='hg19', 
        metavar='GENOME', choices=['dm3', 'hg19', 'hg38', 'mm10'],
        help='Reference genome, support: dm3, hg19, hg38, mm10, default: hg19')
    parser.add_argument('-o', default=None, 
        metavar='OUTPUT', help='The directory to save results.')
    parser.add_argument('--threads', default=1, 
        metavar='THREAD', type=int,
        help='Number of threads to launch, default 1.')
    return parser



def add_rtstop_args(parser):
    """
    call RT-Stops
    """
    parser.add_argument('-i', nargs = '+', required = True, metavar = 'BED', 
        type = argparse.FileType('r'),
        help = 'BED files to call RTStops, 1-4 files.')
    parser.add_argument('-n', required = True, metavar = 'name',
        help = 'Name of the experiment')
    parser.add_argument('-o', required = False, default = None, 
        metavar = 'output',  help = 'The directory to save results.')
    parser.add_argument('-t', required = False, default = 1, 
        choices = list(range(1, 4)), metavar = 'threshold', type = int, 
        help = 'The threshold to filt RTStops, default [1].')
    parser.add_argument('-c', required = False, default = 0, choices = [0, 1], 
        metavar = 'intersect', type = int,
        help = 'how to merge 0=union, 1=intersect, default [0]')
    parser.add_argument('-f', required = False, action = "store_true",
        help = 'Overwrite the output files if exist')
    return parser




def add_report_args(parser):
    """
    create report, advance analysis based on previous results
    """
    pass


def add_run_args(parser):
    """
    run goldclip program from fastq to peaks
    """
    parser.add_argument('-i', nargs='+', required=True, metavar='file', 
        type=argparse.FileType('r'),
        help='reads in FASTQ format, support (*.gz), 1-4 files.')
    parser.add_argument('-o', default=None, metavar='output', 
        help='The directory to save results.')
    parser.add_argument('-g', required=True, default='hg19', 
        metavar='GENOME', choices=['dm3', 'hg19', 'hg38', 'mm10'],
        help='Reference genome : dm3, hg19, hg39, mm10, default: hg19')
    parser.add_argument('-n', required=True, metavar='NAME',
        help='Name of the experiment')
    parser.add_argument('-k', default='hg19', 
        metavar='Spike-in', choices=['dm3', 'hg19', 'hg38', 'mm10'],
        help='Spike-in genome : dm3, hg19, hg38, mm10, default: None')
    parser.add_argument('--trimmed', action='store_true',
        help='if input fastq files were clean reads, specify this option')

    parser.add_argument('-a', default=None, metavar='adapter',
        help='Adapter sequence in the 3 prime end of read, \
        default [AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG].')
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
    parser.add_argument('--rm-untrim', action='store_false', dest='rm_untrim',
        help='if specified, discard reads without adapter')
    parser.add_argument('--rm-dup', action='store_true', dest='rm_dup',
        help='if specified, remove duplicated reads' )
    parser.add_argument('--cut-before-trim', default='0', metavar='cut1', 
        dest='cut_before_trim',
        help='cut bases before trimming adapter, Number of bases to cut from each \
        read, plus on 5-prime end, minus on 3-prime end, could be \
        single, or double numbers, eg: 3 or -4 or 3,-4, default [0]')
    parser.add_argument('--cut-after-trim', default='0', metavar='cut2', 
        dest='cut_after_trim',
        help='cut bases after trimming adapter, Number of bases to cut from each \
        read, plus on 5-prime end, minus on 3-prime end, , could be \
        single, or double numbers, eg: 3 or -4 or 3,-4, default [0]')
    parser.add_argument('--aligner', default='bowtie', 
        choices=['bowtie', 'bowtie2', 'star'],
        help='Choose which aligner to use. default: bowtie')
    parser.add_argument('--unique-only', action='store_true', 
        dest='unique_only',
        help='if specified, keep unique mapped reads only')
    parser.add_argument('--align-to-rRNA', action='store_true',
        dest='align_to_rRNA', 
        help='if specified, remove rRNA before mapping to reference genome')
    parser.add_argument('-t', required = False, default = 1, 
        choices = list(range(1, 4)), metavar = 'threshold', type = int, 
        help = 'The threshold to filt RTStops, default [1].')
    parser.add_argument('-c', required = False, default = 0, choices = [0, 1], 
        metavar = 'intersect', type = int,
        help = 'how to merge 0=union, 1=intersect, default [0]')
    parser.add_argument('--path_data', 
        help='The directory of genome files, default: \
        [$HOME/data/genome/]')
    parser.add_argument('--threads', default=1, 
        metavar='threads', type=int,
        help='Number of threads to launch, default [1]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    return parser


def add_report_args(parser):
    parser.add_argument('--path', required=True, metavar='goldclip_output',
        help = 'The directory of goldclip output')
    parser.add_argument('--name', required=True, metavar='project name',
        help = 'The name of the project')
    parser.add_argument('--genome', required=True, metavar='GENOME',
        help = 'The reference genome of the project')
    return parser



