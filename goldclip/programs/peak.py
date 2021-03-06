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
import pysam
import subprocess

# from goldclip.configure import goldclip_home
from goldclip.helper import *


class Peak:

    def __init__(self, *args, **kwargs):
        self.kwargs = kwargs


    def adjust_score(read):
        """
        fix clipper output 
        template: 
        https://github.com/YeoLab/gscripts/blob/master/gscripts/clipseq/fix_scores.py
        """
        import numpy as np
        peak_center = str((int(read[6]) + int(read[7])) / 2)
        qValue = read.score
        pValue = str(-1) #to fix
        signalValue = str(-1)
        if float(read.score) != 0.0:
            read.score = str(min(int(-10 * np.log10(float(read.score))), 1000))
        else:
            read.score = '0'
        read[6] = signalValue
        read[7] = pValue
        read.append(qValue)
        read.append(peak_center)
        return read
           


    def run_clipper(genome, bam_in, path_out, output):
        """
        call peaks using clipper, only support hg19 now
        """
        pkg_dir, _ = os.path.split(__file__)
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
        pkg_dir, _ = os.path.split(__file__)
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



    def call_peak(genome, bam_ins, path_out, tool = 'clipper'):
        """
        call peaks using clipper, pyicoclip
        multiple BAM files supported
        """
        logging.info('calling peaks')
        path_out = path_out if path_out else os.getcwd()
        assert is_path(path_out)
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


    def run():
        tool = self.kwargs['tool']
        genome = self.kwargs['g']
        bam_files = [f.name for f in self.kwargs['i']]
        path_out = self.kwargs['o']
        peak_file = call_peak(genome, bam_files, path_out, tool)
        return peak_file


## EOF
