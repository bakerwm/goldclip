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
# import pathlib
import logging
import pysam
import pybedtools
from goldclip.helper import *
from goldclip.goldcliplib.log_parser import *

logging.basicConfig(format = '[%(asctime)s] %(message)s', 
                    datefmt = '%Y-%m-%d %H:%M:%S', 
                    level = logging.DEBUG)


# def bowtie_se(fn, idx, path_out, para=1, multi_cores=1, overwrite=False):
def bowtie_se(fn, idx, path_out, multi_cores=1, unique_only=False, overwrite=False):
    """
    Mapping SE reads to idx using Bowtie
    """
    assert isinstance(fn, str)
    assert os.path.exists(fn)
    assert is_idx(idx, 'bowtie')
    path_out = os.path.dirname(fn) if path_out is None else path_out
    assert is_path(path_out)
    ## parameters
    if unique_only is True:
        para_bowtie = '-v 2 -m 1' # report only unique map reads
    else:
        para_bowtie = '-v 2 -k 1' # report only one hit for each read
    fn_type = seq_type(fn)
    if fn_type == 'fasta':
        para_bowtie += ' -f'
    elif fn_type == 'fastq':
        para_bowtie += ' -q'
    else:
        raise ValueError('unknown type of reads')
    ## prefix
    fn_prefix = file_prefix(fn)[0]
    fn_prefix = re.sub('\.clean|\.nodup|\.cut', '', fn_prefix)
    idx_name = os.path.basename(idx)
    fn_unmap_file = os.path.join(path_out, '%s.not_%s.%s' % (fn_prefix, idx_name, fn_type))
    fn_map_prefix = os.path.join(path_out, fn_prefix)
    fn_map_bam = fn_map_prefix + '.map_%s.bam' % idx_name
    fn_map_bed = fn_map_prefix + '.map_%s.bed' % idx_name
    fn_map_log = fn_map_prefix + '.map_%s.bowtie.log' % idx_name
    if os.path.exists(fn_map_bam) and os.path.exists(fn_unmap_file) and overwrite is False:
        logging.info('file exists: %s' % fn_map_bam)
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



def bowtie2_se(fn, idx, path_out, multi_cores=1, unique_only=False, overwrite=False):
    """
    Mapping SE reads to idx using Bowtie2
    """
    assert isinstance(fn, str)
    assert os.path.exists(fn)
    assert is_idx(idx, 'bowtie2')
    path_out = os.path.dirname(fn) if path_out is None else path_out
    assert is_path(path_out)
    ## parameters
    if unique_only is True:
        para_unique = '-q 10' # samtools MAPQ 10
    else:
        para_unique = ''
    fn_type = seq_type(fn)
    if seq_type(fn) == 'fasta':
        para_bowtie2 = '-f'
    elif seq_type(fn) == 'fastq':
        para_bowtie2 = '-q'
    else:
        raise ValueError('unknown type of reads')
    ## prefix
    fn_prefix = file_prefix(fn)[0]
    fn_prefix = re.sub('\.clean|\.nodup|\.cut', '', fn_prefix)
    idx_name = os.path.basename(idx)
    fn_unmap_file = os.path.join(path_out, '%s.not_%s.%s' % (fn_prefix, idx_name, fn_type))
    fn_map_prefix = os.path.join(path_out, fn_prefix)
    fn_map_bam = fn_map_prefix + '.map_%s.bam' % idx_name
    fn_map_bed = fn_map_prefix + '.map_%s.bed' % idx_name
    fn_map_log = fn_map_prefix + '.map_%s.bowtie2.log' % idx_name
    if os.path.exists(fn_map_bam) and os.path.exists(fn_unmap_file) and overwrite is False:
        logging.info('file exists: %s' % fn_map_bam)
    else:
        c1 = 'bowtie2 %s -p %s --mm --no-unal --un %s -x %s %s' % (para_bowtie2,
              multi_cores, fn_unmap_file, idx, fn)
        c2 = 'samtools view -bhS -F 0x4 -@ %s %s -' % (multi_cores, para_unique)
        c3 = 'samtools sort -@ %s -o %s -' % (multi_cores, fn_map_bam)
        with open(fn_map_log, 'wt') as fo:
            p1 = subprocess.Popen(shlex.split(c1), stdout=subprocess.PIPE, stderr=fo)
            p2 = subprocess.Popen(shlex.split(c2), stdin=p1.stdout, stdout=subprocess.PIPE)
            p3 = subprocess.Popen(shlex.split(c3), stdin=p2.stdout)
            p4 = p3.communicate()
        pysam.index(fn_map_bam)
        # pybedtools.BedTool(fn_map_bam).bam_to_bed().saveas(fn_map_bed)
    ## statistics
    d = bowtie2_log_parser(fn_map_log)

    return [fn_map_bam, fn_unmap_file]



def star_se(fn, idx, path_out, multi_cores=1, unique_only=False, overwrite=False):
    """
    mapping single read to one index using STAR
    Input: fastq|a
    Output: bam (sorted), unmapped reads
    #
    unique mapped:
    --outFilterMismatchNoverLmax 0.07 \
    --outFilterMultimapNmax 1 \

    --outFilterMismatchNoverLMax 0.05 \
    --seedSearchStartLmax 20

    """
    assert isinstance(fn, str)
    assert os.path.exists(fn)
    assert is_idx(idx, 'star')
    path_out = os.path.dirname(fn) if path_out is None else path_out
    assert is_path(path_out)
    ## parameters
    if unique_only is True:
        para_star = '--outFilterMismatchNoverLmax 0.07 --outFilterMultimapNmax 1'
    else:
        para_star = '--outFilterMismatchNoverLmax 0.05 --seedSearchStartLmax 20'
    ## prefix
    fn_prefix = file_prefix(fn)[0]
    fn_prefix = re.sub('\.clean|\.nodup|\.cut', '', fn_prefix)
    idx_name = os.path.basename(idx)
    freader = 'zcat' if is_gz(fn) else '-'
    fn_type = seq_type(fn)
    fn_unmap_file = os.path.join(path_out, '%s.not_%s.%s' % (fn_prefix, idx_name, fn_type))
    fn_map_prefix = os.path.join(path_out, fn_prefix)
    fn_map_bam = fn_map_prefix + '.map_%s.bam' % idx_name
    fn_map_bed = fn_map_prefix + '.map_%s.bed' % idx_name
    fn_map_log = fn_map_prefix + '.map_%s.star.log' % idx_name
    ## skip exist files
    if os.path.exists(fn_map_bam) and overwrite is False:
        logging.info('file exists: %s' % fn_map_bam)
    else:
        c1 = 'STAR --runMode alignReads \
              --genomeDir %s \
              --readFilesIn %s \
              --readFilesCommand %s \
              --outFileNamePrefix %s \
              --runThreadN %s \
              --limitOutSAMoneReadBytes 1000000 \
              --genomeLoad LoadAndRemove \
              --limitBAMsortRAM 10000000000 \
              --outSAMtype BAM SortedByCoordinate \
              --outReadsUnmapped Fastx\
              %s' % (idx, fn, freader, fn_map_prefix, multi_cores, para_star)
        p1 = subprocess.run(shlex.split(c1))
        # rename exists file
        os.rename(fn_map_prefix + 'Aligned.sortedByCoord.out.bam', fn_map_bam)
        os.rename(fn_map_prefix + 'Unmapped.out.mate1', fn_unmap_file)
        os.rename(fn_map_prefix + 'Log.final.out', fn_map_log)
        pysam.index(fn_map_bam)
        d = star_log_parser(fn_map_log)

    return [fn_map_bam, fn_unmap_file]



def align_se_batch(fn, idxes, path_out, aligner='STAR', multi_cores=1, 
                   unique_only=False, overwrite=False):
    """
    mapping fastq to multiple indexes
    """
    assert isinstance(fn, str)
    assert os.path.exists(fn)
    assert isinstance(idxes, list)
    path_out = os.path.dirname(fn) if path_out is None else path_out
    assert is_path(path_out)
    # aligner
    if aligner.lower() == 'star':
        aligner_se = star_se
    elif aligner.lower() == 'bowtie2':
        aligner_se = bowtie2_se
    elif aligner.lower() == 'bowtie':
        aligner_se = bowtie_se
    else:
        raise ValueError('unknown aligner: %s' % aligner)
    # iterate index
    fn_bam_files = []
    fn_input = fn
    for idx in idxes:
        # logging.error('idx: %s' % idx)
        fn_bam_idx, fn_unmap_idx = aligner_se(fn_input, idx, path_out,
                                              multi_cores=multi_cores,
                                              unique_only=unique_only,
                                              overwrite=overwrite)
        fn_input = fn_unmap_idx
        fn_bam_files.append(fn_bam_idx)
    return fn_bam_files



def align(fns, smp_name, path_out, genome, spikein=None, multi_cores=1, 
          aligner='bowtie', path_data=None, unique_only=False, 
          align_to_rRNA=False, overwrite=False):
    """
    mapping multiple reads to multiple indexes, one-by-one
    merge
    """
    assert isinstance(fns, list)
    assert isinstance(genome, str)
    assert isinstance(smp_name, str)
    # get indexes
    sp = idx_picker(spikein, path_data=path_data, aligner=aligner) #
    if align_to_rRNA is True:
        sg = idx_grouper(genome, path_data=path_data, aligner=aligner) #
        idxes = sg if spikein == genome else sg.append(sp)
    else:
        sg = idx_picker(genome, path_data=path_data, aligner=aligner) #
        idxes = [sg, ] if spikein == genome else [sg, sp]
    idxes = list(filter(None.__ne__, idxes)) # idxes
    # remove duplicate records, keep orders
    idxes = list(dict.fromkeys(idxes))
    if len(idxes) == 0:
        raise ValueError('genome index not exists: ' + path_data)
    # mapping se reads
    fn_bam_files = []
    for fn in fns:
        logging.info('mapping file: %s' % fn)
        fn_prefix = file_prefix(fn)[0]
        fn_prefix = re.sub('\.clean|\.nodup|\.cut', '', fn_prefix)
        path_out_fn = os.path.join(path_out, fn_prefix)
        b = align_se_batch(fn, idxes, path_out_fn, multi_cores=multi_cores,
                           aligner=aligner, unique_only=unique_only, 
                           overwrite=overwrite)
        fn_bam_files.append(b) # bam files
        rep_map_wrapper(path_out_fn)

    # merge bam files
    path_out_merge = os.path.join(path_out, smp_name)
    merge_bam_files = []
    if len(fn_bam_files) > 1:
        assert is_path(path_out_merge)
        for i in range(len(fn_bam_files[0])): # merge each sub-index
            se_bam_files = [b[i] for  b in fn_bam_files]
            merge_suffix = str_common(se_bam_files, suffix=True)
            merge_suffix = re.sub('^_[12]|_R[12]', '', merge_suffix)
            merge_bam_name = smp_name + merge_suffix
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

    # get genome mapping files (the last one)
    genome_bam_files = [f[-1] for f in fn_bam_files]

    # rename genome bam, to a shorter name
    # remove "not_index." suffix
    gbam_files = []
    gbed_files = []
    for i in range(len(genome_bam_files)):
        bam_from = genome_bam_files[i]
        bam_to = os.path.join(os.path.dirname(bam_from), 
                              filename_shorter(bam_from))
        if not os.path.exists(bam_to):
            os.symlink(os.path.basename(bam_from), bam_to)
        if not os.path.exists(bam_to + '.bai'):
            if not os.path.exists(bam_from + '.bai'):
                pysam.index(bam_from)
            os.symlink(os.path.basename(bam_from) + '.bai', 
                       bam_to + '.bai')
        gbam_files.append(bam_to)
        # rename .bed
        bed_from = re.sub('.bam$', '.bed', bam_from)
        bed_to = re.sub('.bam$', '.bed', bam_to)
        if os.path.exists(bed_from) and not os.path.exists(bed_to):
            os.symlink(os.path.basename(bed_from), bed_to)
        gbed_files.append(bed_to)

    # return [gbam_files, gbed_files]
    return gbam_files


## EOF
