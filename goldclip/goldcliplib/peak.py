#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
call peaks using CLIPper or Pyicoclip
"""


import os
import sys
import re
import shlex
import logging
import multiprocessing as mp
import numpy as np
import pysam
import subprocess
import goldclip
from goldclip.helper import *



def adjust_score(bed):
    """
    fix clipper output 
    template: 
    https://github.com/YeoLab/gscripts/blob/master/gscripts/clipseq/fix_scores.py
    """
    peak_center = str((int(bed[6]) + int(bed[7])) / 2)
    qValue = bed.score
    pValue = str(-1) #to fix
    signalValue = str(-1)
    if float(bed.score) != 0.0:
        bed.score = str(min(int(-10 * np.log10(float(bed.score))), 1000))
    else:
        bed.score = '0'
    bed[6] = signalValue
    bed[7] = pValue
    bed.append(qValue)
    bed.append(peak_center)
    return bed
       


def run_clipper(genome, bam_in, path_out, output):
    """
    call peaks using clipper, only support hg19 now
    """
    pkg_dir, _ = os.path.split(goldclip.__file__)
    bam_to_peak_clipper = os.path.join(pkg_dir, 'bin', 'bam2peaks.clipper.sh')
    if not os.path.exists(bam_to_peak_clipper):
        sys.exit(bam_to_peak_clipper + ' - script not found.')
    # Check output peak file
    peak_prefix = os.path.basename(os.path.splitext(bam_in)[0])
    path_peak = os.path.join(path_out, peak_prefix)
    peak_bed = os.path.join(path_peak, peak_prefix + '.bed')
    peak_fix = os.path.join(path_peak, peak_prefix + '.fixed.bed')
    c1 = 'bash {} {} {} {}'.format(bam_to_peak_clipper, genome, bam_in, path_out)
    os.system(c1)    
    if os.path.exists(peak_bed):
        pybedtools.BedTool(peak_bed).each(adjust_score).saveas(peak_fix)
    if os.path.exists(peak_fix):
        output.put(peak_fix)
    return peak_fix



def run_pyicoclip(genome, bam_in, path_out, output):
    """
    call peaks using pyicpclip, support hg19, dm3 now
    output : Queue
    """
    pkg_dir, _ = os.path.split(goldclip.__file__)
    bam_to_peak_pyicoclip = os.path.join(pkg_dir, 'bin', 'bam2peaks.pyicoclip.sh')
    if not os.path.exists(bam_to_peak_pyicoclip):
        sys.exit(bam_to_peak_clipper + ' - script not found.')
    # Check output peak file
    peak_prefix = os.path.basename(os.path.splitext(bam_in)[0])
    path_peak = os.path.join(path_out, peak_prefix)
    peak_fix = os.path.join(path_peak, peak_prefix + '.fixed.bed')
    c2 = 'bash {} {} {} {}'.format(bam_to_peak_pyicoclip, genome, bam_in, path_out)
    os.system(c2)
    if os.path.exists(peak_fix):
        output.put(peak_fix)
    return peak_fix



def call_peak(genome, bam_ins, path_out, tool='clipper'):
    """
    call peaks using clipper, pyicoclip
    multiple BAM files supported
    """
    logging.info('calling peaks - %s' % tool)
    path_out = path_out if path_out else os.getcwd()
    assert is_path(path_out)
    if tool == 'clipper':
        if not genome in ['hg19', 'mm9']:
            # raise ValueError('CLIPper not support: %s' % genome)
            logging.error('CLIPper not support: %s' % genome)
            return None
    if tool == 'pyicoclip':
        if not genome in ['dm3', 'mm9', 'mm10', 'hg19', 'hg38']:
            # raise ValueError('Pyicoclip not support: %s' % genome)
            logging.error('Pyicoclip not support: %s' % genome)
            return None
    output = mp.Queue() # run in parallel
    if tool == 'clipper':
        processes = [mp.Process(target = run_clipper, 
                                args = (genome, b, path_out, output))
                    for b in bam_ins]
    elif tool == 'pyicoclip':
        processes = [mp.Process(target = run_pyicoclip, 
                                args = (genome, b, path_out, output)) 
                    for b in bam_ins]
    else:
        sys.exit(tool + ' : unknown tool, [clipper|pyicoclip]')        
    # Run processes
    for p in processes:
        p.start()
    # Exit completed process
    for p in processes:
        p.join()
    results = [output.get() for p in processes]
    return results


## EOF
