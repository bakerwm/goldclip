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


def add_trim_args(parser):
    """
    processing SE read 
    - remove 3' adapter(s) (default: TruSeq RNA-Seq)
    - trim low-quality bases on both 5 and 3 end
    - trim N reads
    - cut N-bases at either end of read
    """
    parser.add_argument('-i', nargs='+', required=True, metavar='file', 
        type=argparse.FileType('r'),
        help='reads in FASTQ format, support (*.gz), 1-4 files.')
    parser.add_argument('-a', default='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC', 
        metavar='adapter', type=str,
        help='3-Adapter, default: [AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC].')
    parser.add_argument('-o', default=None, metavar='output', 
        help='The directory to save results.')
    parser.add_argument('-g', default='', metavar='adapter-5', type=str,
        help='5-Adapter, default: None')
    parser.add_argument('-m', default=15, metavar='len_min', 
        type=int, help='Minimum length of reads after trimming, defualt [15]')
    parser.add_argument('--read12', type=int, default=1, metavar='read12',
        help='which one of PE reads, 1=read1, 2=read2, default: 1')
    
    ## global arguments    
    parser.add_argument('-q', default=20, metavar='quality', type=int,
        help='The cutoff of base quality, default [20]')    
    parser.add_argument('-e', default=0.1, metavar='err_rate', type=float,
        help='Maximum allowed error rate, default [0.1]')
    parser.add_argument('-O', default=3, metavar='overlap', type=int,
        help='Required N bases overlap between reads and adapter, default [3]')
    parser.add_argument('-p', default=80, metavar='percent', type=int,
        help='minimum percent of bases that must have -q quality, default [80]')
    parser.add_argument('--rm-untrim', action='store_true', dest='rm_untrim',
        help='if specified, discard reads without adapter')
    parser.add_argument('--threads', default=1, metavar='threads', type=int,
        help='Number of threads to launch, default [1]')
    parser.add_argument('--overwrite', action='store_true',
        help='if spcified, overwrite exists file')
    parser.add_argument('--keep-name', action='store_true', dest='keep_name',
        help='if specified, do not change file names')
    
    ## extra arguments
    parser.add_argument('--adapter-sliding', dest='adapter_sliding', 
        action='store_true',
        help='Trim reads by sliding windows on adapter')
    parser.add_argument('--trim-times', dest='trim_times', type=int,
        default=1, help='Trim adapter from reads by N times, default:1')
    parser.add_argument('--double-trim', action='store_true', 
        dest='double_trim', help='if specified, trim adapters twice')
    parser.add_argument('--rm-dup', action='store_true', dest='rm_dup',
        help='if specified, remove duplicated reads' )
    parser.add_argument('--cut-before-trim', default='0', metavar='cut1', 
        dest='cut_before_trim',
        help='cut bases before trimming adapter, Number of bases to cut \
              from each read, plus on 5-prime end, minus on 3-prime end, \
              could be single, or double numbers, eg: 3 or -4 or 3,-4, \
              default [0]')
    parser.add_argument('--cut-after-trim', default='0', metavar='cut2', 
        dest='cut_after_trim',
        help='cut bases after trimming adapter, Number of bases to cut \
              from each read, plus on 5-prime end, minus on 3-prime end, \
              could be single, or double numbers, eg: 3 or -4 or 3,-4, \
              default [0]')
    parser.add_argument('--trim-to-length', default=0, metavar='max-length',
        dest='trim_to_length', type=int,
        help='trim reads from right, save the specific length of reads. \
              default: [0], 0=the full length')

    ## PE arguments
    parser.add_argument('--fq2', nargs='+', default=None, 
        help='The read2 of pair-end reads')
    parser.add_argument('-A', default='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
        help='The 3 adapter of read2, default: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT')
    parser.add_argument('-G', default=None,
        help='The 5 adapter of read1, default: None')
    return parser


def add_align_args(parser):
    """
    Mapping SE read or one of PE reads to reference genome
    using bowtie, STAR, ... (universal)
    """
    parser.add_argument('-i', nargs='+', required=True, metavar='INPUT', 
        type=argparse.FileType('r'),
        help='CLIP reads in FASTQ format, (not *.gz), 1-4 files.')
    parser.add_argument('-o', default=None, 
        metavar='OUTPUT',  help='The directory to save results, default, \
        current working directory.')
    parser.add_argument('-n', required=True, metavar='NAME',
        help='Name of the experiment')
    parser.add_argument('-g', required=True, default='hg19', 
        metavar='GENOME', choices=['dm3', 'hg19', 'hg38', 'mm10', 'mm9'],
        help='Reference genome : dm3, hg19, hg39, mm10, default: hg19')
    parser.add_argument('-k', default=None, 
        metavar='Spike-in', choices=[None, 'dm3', 'hg19', 'hg38', 'mm10'],
        help='Spike-in genome : dm3, hg19, hg38, mm10, default: None')
    parser.add_argument('-x', nargs='+', metavar='align_index',
        help='Provide alignment index(es) for alignment, support multiple\
        indexes. if specified, ignore -g, -k')
    parser.add_argument('--threads', default=8, 
        metavar='THREADS', type=int, 
        help='Number of threads to launch, default: 8.')
    parser.add_argument('--unique-only', action='store_true',
        dest='unique_only',
        help='if specified, keep unique mapped reads only')
    parser.add_argument('--n-map', dest='n_map', type=int, default=0,
        help='Report up to N alignments per read. use -k for bowtie and \
        bowtie2 (default 1), --outFilterMultimapNmax for STAR \
        (default 20).')
    parser.add_argument('--aligner', default='bowtie', 
        choices=['bowtie', 'bowtie2', 'STAR'],
        help='Choose which aligner to use. default: bowtie')
    parser.add_argument('--align-to-rRNA', dest='align_to_rRNA',
        action='store_true',
        help='if specified, align to rRNA before genome')
    parser.add_argument('--repeat-masked-genome', dest='repeat_masked_genome',
        action='store_true',
        help='map to repeat masked reference genome, data from EnsEMBL')
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
    parser.add_argument('--peak-caller', required=True, 
        dest='peak_caller',
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
    parser.add_argument('--project-path', required=True, 
        dest='project_path', metavar='project_output',
        help = 'The directory of goldclip output')
    parser.add_argument('--project-name', required=True, 
        dest='project_name', metavar='project name',
        help = 'The name of the project')
    parser.add_argument('-g', required=True, metavar='genome',
        help = 'The reference genome of the project')
    ## to-do
    ## automaticly extract the information from the directory
    return parser


def add_run_args(parser):
    """
    run goldclip program from fastq to peaks
    """
    parser.add_argument('-i', nargs='+', required=True, 
        metavar='input fastq files',
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

    ## alignment
    parser.add_argument('--aligner', default='bowtie', 
        choices=['bowtie', 'bowtie2', 'STAR'],
        help='Choose which aligner to use. default: bowtie')
    parser.add_argument('--align-to-rRNA', dest='align_to_rRNA',
        action='store_true',
        help='if specified, align to rRNA before genome')
    parser.add_argument('-x', nargs='+', metavar='align_index',
        help='Provide alignment index(es) for alignment, support multiple\
        indexes. if specified, ignore -g, -k')
    parser.add_argument('--unique-only', action='store_true',
        dest='unique_only',
        help='if specified, keep unique mapped reads only')
    parser.add_argument('--n-map', dest='n_map', type=int, default=0,
        help='Report up to N alignments per read. use -k for bowtie and \
        bowtie2 (default 1), --outFilterMultimapNmax for STAR \
        (default 20).')
    parser.add_argument('--repeat-masked-genome', dest='repeat_masked_genome',
        action='store_true',
        help='map to repeat masked reference genome, data from EnsEMBL')
    parser.add_argument('--path_data', 
        help='The directory of genome files, default: \
        [$HOME/data/genome/]')

    ## call peaks
    parser.add_argument('--peak-caller', required=True, dest='peak_caller',
        choices=['clipper', 'pyicoclip'], 
        help='Peak-caller, clipper|pyicoclip')

    ## call rtstop
    parser.add_argument('-t', required = False, default = 1, 
        choices = list(range(1, 4)), metavar = 'threshold', type = int, 
        help = 'The threshold to filt RTStops, default [1].')
    parser.add_argument('-c', required = False, default = 0, choices = [0, 1], 
        metavar = 'intersect', type = int,
        help = 'how to merge 0=union, 1=intersect, default [0]')
    parser.add_argument('--overwrite', required = False, action = "store_true",
        help = 'Overwrite the output files if exist')

    return parser



