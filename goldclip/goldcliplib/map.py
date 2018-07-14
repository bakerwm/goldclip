#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Mapping reads to reference genome
1. to spike-in genome
2. to repbase
3. to rRNA, tRNA, MT
4. to genome
"""
# TO-DO
# count RT in rep1/2


__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__date__ = "2018-03-21"
__version__ = "0.1"


import os
import sys
import re
import shlex
import subprocess
import pathlib
import logging
import pysam
import pybedtools

from goldclip.helper import *

logging.basicConfig(format = '[%(asctime)s] %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S', 
                    level = logging.DEBUG)


def bowtie_se(fn, idx, path_out, para='-v 2 -k 2', multi_cores=1, overwrite=False):
    """
    Mapping SE reads to idx using Bowtie
    """
    assert isinstance(fn, str)
    assert os.path.exists(fn)
    assert idx_checker(idx, 'bowtie')
    path_out = os.path.dirname(fn) if path_out is None else path_out
    assert is_path(path_out)
    ## parameters
    fn_type = seq_type(fn)
    if not para is None:
        para_bowtie += ' ' + para
    if seq_type(fn) == 'fasta':
        para_bowtie += ' -f'
    elif seq_type(fn) == 'fastq':
        para_bowtie += ' -q'
    else:
        raise ValueError('unknown type of reads')
    ## prefix
    fn_prefix = file_prefix(fn)[0]
    fn_prefix = re.sub('\.clean|\.nodup|\.cut\-?\d+?', '', fn_prefix)
    idx_name = os.path.basename(idx)
    fn_unmap_file = os.path.join(path_out, '%s.not_%s.%s' % (fn_prefix, dix_name, fn_type))
    fn_map_prefix = os.path.join(path_out, fn_prefix)
    fn_map_bam = fn_map_prefix + '.bam'
    fn_map_bed = fn_map_prefix + '.bed'
    fn_map_log = fn_map_prefix + '.bowtie.log'
    if os.path.exists(fn_map_bam) and os.path.exists(fn_unmap_file) and overwrite is False:
        logging.info('file exists: %s' % fn_prefix)
    else:
        c1 = 'bowtie %s -p %s --mm --best --sam --no-unal --un %s %s %s' % (para_bowtie,
            multi_cores, fn_unmap_file, idx, fn)
        c2 = 'samtools view -bhS -F 0x4 -@ %s -' % multi_cores
        c3 = 'samtools sort -@ %s -o %s -' % (multi_cores, fn_map_bam)
        with open(fn_map_log, 'wt') as fo:
            p1 = subprocess.Popen(shlex.split(c1), stdout=subprocess.PIPE, stderr=fo)
            p2 = subprocess.Popen(shlex.split(c2), stdin=p1.stdout, stdout=subprocess.PIPE)
            p3 = subprocess.Popen(shlex.split(c3), stdin=p2.stdout)
            p4 = p3.communicate()
        pysam.index(fn_map_bam)
        pybedtools.BedTool(fn_map_bam).bam_to_bed().saveas(fn_map_bed)
    ## statistics
    d = bowtie_log_parser(fn_map_log)

    return [fn_map_bam, fn_unmap_file]



def bowtie2_se(fn, idx, path_out, para='-v 2 -k 2', multi_cores=1, overwrite=False):
    """
    Mapping SE reads to idx using Bowtie
    """
    assert isinstance(fn, str)
    assert os.path.exists(fn)
    assert idx_checker(idx, 'bowtie2')
    path_out = os.path.dirname(fn) if path_out is None else path_out
    assert is_path(path_out)
    ## parameters
    fn_type = seq_type(fn)
    if not para is None:
        para_bowtie2 += ' ' + para
    if seq_type(fn) == 'fasta':
        para_bowtie2 += ' -f'
    elif seq_type(fn) == 'fastq':
        para_bowtie2 += ' -q'
    else:
        raise ValueError('unknown type of reads')
    ## prefix
    fn_prefix = file_prefix(fn)[0]
    fn_prefix = re.sub('\.clean|\.nodup|\.cut\-?\d+?', '', fn_prefix)
    idx_name = os.path.basename(idx)
    fn_unmap_file = os.path.join(path_out, '%s.not_%s.%s' % (fn_prefix, dix_name, fn_type))
    fn_map_prefix = os.path.join(path_out, fn_prefix)
    fn_map_bam = fn_map_prefix + '.bam'
    fn_map_bed = fn_map_prefix + '.bed'
    fn_map_log = fn_map_prefix + '.bowtie2.log'
    if os.path.exists(fn_map_bam) and os.path.exists(fn_unmap_file) and overwrite is False:
        logging.info('file exists: %s' % fn_prefix)
    else:
        c1 = 'bowtie2 %s -p %s --mm --best --sam --no-unal --un %s -x %s %s' % (para_bowtie2,
              multi_cores, fn_unmap_file, idx, fn)
        c2 = 'samtools view -bhS -F 0x4 -@ %s -' % multi_cores
        c3 = 'samtools sort -@ %s -o %s -' % (multi_cores, fn_map_bam)
        with open(fn_map_log, 'wt') as fo:
            p1 = subprocess.Popen(shlex.split(c1), stdout=subprocess.PIPE, stderr=fo)
            p2 = subprocess.Popen(shlex.split(c2), stdin=p1.stdout, stdout=subprocess.PIPE)
            p3 = subprocess.Popen(shlex.split(c3), stdin=p2.stdout)
            p4 = p3.communicate()
        pysam.index(fn_map_bam)
        pybedtools.BedTool(fn_map_bam).bam_to_bed().saveas(fn_map_bed)
    ## statistics
    d = bowtie2_log_parser(fn_map_log)

    return [fn_map_bam, fn_unmap_file]



def star_se(fn, idx, path_out, para=1, multi_cores=1):
    """
    mapping single read to one index using STAR
    Input: fastq|a
    Output: bam (sorted), unmapped reads
    """
    assert isinstance(fn, str)
    assert os.path.exists(fn)
    assert idx_checker(idx, 'star')
    path_out = os.path.dirname(fn) if path_out is None else path_out
    assert is_path(path_out)
    ## prefix
    fn_prefix = file_prefix(fn)[0]
    fn_prefix = re.sub('\.clean|\.nodup|\.cut\-?\d+?', '', fn_prefix)
    idx_name = os.path.basename(idx)
    fn_unmap_file = os.path.join(path_out, '%s.not_%s.%s' % (fn_prefix, dix_name, fn_type))
    fn_map_prefix = os.path.join(path_out, fn_prefix)
    fn_map_bam = fn_map_prefix + '.bam'
    fn_map_bed = fn_map_prefix + '.bed'
    fn_map_log = fn_map_prefix + '.bowtie2.log'
    ## skip exist files
    if os.path.exists(fn_map_bam) and os.path.exists(fn_unmap_file) and overwrite is False:
        logging.info('file exists: %s' % fn_prefix)
    else:
        c1 = 'STAR --runMode alignReads --genomeDir %s --readFilesIn %s \
              --readFilesCommand %s --outFileNamePrefix %s \
              --runThreadN %s ' % (idx, fn, freader, fn_map_prefix, multi_cores)
        c1 += ' --limitOutSAMoneReadBytes 1000000 \
             --genomeLoad LoadAndKeep \
             --limitBAMsortRAM 10000000000 \
             --outSAMtype BAM SortedByCoordinate \
             --outReadsUnmapped Fastx \
             --outFilterMismatchNoverLmax 0.05 \
             --seedSearchStartLmax 20'
        p1 = subprocess.run(shlex.split(c1))
        # rename exists file
        os.rename(fn_map_prefix + 'Aligned.sortedByCoord.out.bam', fn_map_bam)
        os.rename(fn_map_prefix + 'Unmapped.out.mate1', fn_unmap_file)
        os.rename(fn_map_prefix + 'Log.final.out', fn_map_log)
        pysam.index(fn_map_bam)
        d = star_log_parser(log_out)

    return [fn_map_bam, fn_unmap_file]



def bowtie_se_batch(fn, idxes, path_out, para='-v 2 -k 1', multi_cores=1):
    """
    mapping fastq to multiple indexes
    """
    assert isinstance(fn, str)
    assert os.path.exists(fn)
    assert isinstance(idxes, list)
    path_out = os.path.dirname(fn) if path_out is None else path_out
    assert is_path(path_out)
    # iterate index
    fn_bam_files = []
    fn_input = fn
    for idx in idxes:
        para_bowtie = '-v 2 -m 1' if idx is idxes[-1] else '-v 2 -k 1'
        fn_bam_idx, fn_unmap_idx = bowtie_se(fn_input, idx, path_out, 
                                             para_bowtie,
                                             multi_cores=multi_cores)
        fn_input = fn_unmap_idx
        fn_bam_files.append(fn_bam_idx)
    return fn_bam_files



def map_se(fns, path_out, genome, spikein, smp_name, genome_index, multi_cores=1):
    """
    mapping reads to multiple indexes, one-by-one
    call RT stops
    merge RT stops 
    filter merged_rt by snoRNA, miRNA, ...
    Input: read1, read2, ... [1 to 4]
    Input: idx1, idx2, ...
    Output: RTStops
    """
    assert isinstance(smp_name, str)
    assert isinstance(fns, list)
    idxes = aligner_index_picker(genome, spikein)
    fn_bam_files = []
    logging.info('mapping reads')
    ## mapping 
    for fn in fns:
        fn_prefix = file_prefix(fn)
        fn_prefix = re.sub('\.clean|\.nodup|\.cut\-?\d+?', '', read_prefix)
        path_out_fn = os.path.join(path_out, fn_prefix)
        b = bowtie_se_batch(fn, idxes, path_out_fn, multi_cores=multi_cores)
        fn_bam_files.append(b)
        se_map_wrapper(b)

    # merge bam files
    path_out_merge = os.path.join(path_out, smp_name)
    assert is_path(path_out_merge)
    merge_bam_files = []
    if len(fn_bam_files) > 1:
        for i in range(len(fn_bam_files[0])):
            se_bam_files = [b[i] for  b in fn_bam_files]
            merge_bam_name = smp_name + str_common(se_bam_files, suffix=True)
            merge_bam_file = os.path.join(path_out_merge, merge_bam_name)
            merge_bed_file = re.sub('.bam$', '.bed', merge_bam_file)
            if os.path.exists(merge_bam_file) and overwrite is False:
                logging.info('file exists: %s' % merge_bam_name)
            else:
                tmp = bam_merge(se_bam_files, merge_bam_file)
                pybedtools.BedTool(merge_bam_file).bam_to_bed().saveas(merge_bed_file)
            merge_bam_files.append(merge_bam_file)
        merge_map_wrapper(path_out_merge)
        fn_bam_files.append(merge_bam_files)

    ## get genome mapping files (the last one)
    genome_bam_files = [f[-1] for f in fn_bam_files]

    ## rename genome bam, to a shorter name
    genome_ids = [filename_shorter(f) for f in genome_bam_files]
    gbam_files = []
    gbed_files = []
    for i in range(len(genome_bam_files)):
        bam_from = genome_bam_files[i]
        bam_to = os.path.join(os.path.dirname(bam_from), genome_ids[i])
        if not os.path.exists(bam_to):
            os.symlink(os.path.basename(bam_from), bam_to)
        if not os.path.exists(bam_to + '.bai'):
            if not path.exists(bam_from + '.bai'):
                pysam.index(bam_from)
            os.symlink(os.path.basename(bam_from) + '.bai', bam_to + '.bai')
        gbam_files.append(bam_to)
        # rename .bed
        bed_from = re.sub('.bam$', '.bed', bam_from)
        bed_to = re.sub('.bam$', '.bed', bam_to)
        if os.path.exists(bed_from) and not os.path.exists(bed_to):
            os.symlink(os.path.basename(bed_from), bed_to)
        gbed_files.append(bed_to)

    return [gbam_files, gbed_files]



def run():
    fqs = [f.name for f in self.kwargs['i']]
    smp_name = self.kwargs['n']
    out_path = self.kwargs['o']
    genome = self.kwargs['g']
    spikein = self.kwargs['k']
    genome_index = self.kwargs['x']
    multi_cores = self.kwargs['threads']
    aligner = self.kwargs['aligner']
    tmp = map_se(fqs, smp_name, out_path, genome, spikein, genome_index,
                  multi_cores)

## EOF
