#!/usr/bin/env python
"""
Run goldclip pipeline
1. trimming reads
2. mapping reads
3. call peaks
4. call RTStops
5. report

"""

__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__date__ = "2018-12-25"
__version__ = "0.1"


import os
import sys
from goldclip.helper import BAM
from goldclip.goldcliplib.trim import Trimmer
from goldclip.goldcliplib.alignment import Alignment
from goldclip.goldcliplib.peak import call_peak
from goldclip.goldcliplib.rtstop import *
from goldclip.goldcliplib.report import *


def run_goldclip(fqs, path_out, genome, smp_name, lib_type=1, **kwargs):
    """
    run goldclip pipeline for all
    """
    assert is_path(path_out)
    assert isinstance(genome, str)
    assert isinstance(smp_name, str)

    args = Argument(lib_type).all()
    args = {**args, **kwargs} # update parameters

    # trim
    path_trim = os.path.join(path_out, 'input_reads')
    if is_trimmed:
        clean_fq_files = fq_files
    else:
        clean_fq_files = trim(fq_files, path_trim, **args)

    # # map
    # path_map = os.path.join(path_out, 'genome_mapping')

    # map_bam_files = Alignment(clean_fq_files,
    #         path_map,
    #         smp_name=smp_name,
    #         genome=genome, **args).run()

    # # bigWig
    # path_bw = os.path.join(path_out, 'bigWig')
    # for bam in map_bam_files:
    #     logging.info('making bigWig: %s' % os.path.basename(bam))
    #     bam2bw(bam, genome, path_bw, strandness=True, binsize=10)

    # # peak
    # path_peak1 = os.path.join(path_out, 'peaks', 'clipper')
    # peak_clipper_files = call_peak(genome, map_bam_files, path_peak1, 
    #                                tool='clipper')
    # path_peak2 = os.path.join(path_out, 'peaks', 'pyicoclip')
    # peak_pyicoclip_files = call_peak(genome, map_bam_files, path_peak2, 
    #                                  tool='pyicoclip')

    # # rtstop
    # path_rtstop = os.path.join(path_out, 'rtstops')
    # rtstop_files = call_rtstop(map_bed_files, path_rtstop,  smp_name,
    #                            threshold=threshold, intersect=intersect,
    #                            overwrite=overwrite)

    # # report
    # tmp = Goldclip_output(path_out, smp_name, genome).get_all_figures()
    # goldclip_report(path_out, smp_name, genome)



