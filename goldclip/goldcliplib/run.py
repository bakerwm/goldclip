#!/usr/bin/env python
"""
Run goldclip pipeline
1. trimming reads
2. mapping reads
3. call peaks
4. call RTStops
5. report

functions
"""

__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__date__ = "2018-03-21"
__version__ = "0.1"


import os
import sys
from goldclip.goldcliplib.trim import *
from goldclip.goldcliplib.map import *
from goldclip.goldcliplib.peak import *
from goldclip.goldcliplib.rtstop import *
from goldclip.goldcliplib.report import *

def run_goldclip(fq_files, path_out, genome, smp_name, spikein=None,
                 is_trimmed=False, ad3=None, read12=1, len_min=15, qual_pct=80,
                 qual_min=20, err_rate=0.1, overlap=1, rm_untrim=False, 
                 rm_dup=False, cut_before_trim=0, cut_after_trim=0, 
                 aligner='bowtie', threshold=1, intersect=1, path_data=None, 
                 threads=1, overwrite=False):
    """
    run goldclip pipeline for all
    """
    assert is_path(path_out)
    assert isinstance(genome, str)
    assert isinstance(smp_name, str)

    # trim
    path_trim = os.path.join(path_out, 'input_reads')
    if is_trimmed:
        clean_fq_files = fq_files
    else:
        clean_fq_files = trim(fq_files, path_trim, adapter3=ad3, 
                              len_min=len_min, qual_min=qual_min, 
                              err_rate=err_rate, multi_cores=threads,
                              rm_untrim=rm_untrim, overlap=overlap, 
                              read12=read12, overwrite=overwrite, 
                              rm_dup=rm_dup, cut_before_trim=cut_before_trim,
                              cut_after_trim=cut_after_trim)

    # # map
    # path_map = os.path.join(path_out, 'genome_mapping')
    # map_bam_files, map_bed_files = map(clean_fq_files, smp_name, path_map, 
    #                                    genome, spikein, multi_cores=threads, 
    #                                    aligner=aligner, path_data=path_data, 
    #                                    overwrite=overwrite)

    # # # bigWig
    # # path_bw = os.path.join(path_out, 'bigWig')
    # # for bam in map_bam_files:
    # #     logging.info('making bigWig: %s' % os.path.basename(bam))
    # #     bam2bw(bam, genome, path_bw, strandness=True, binsize=10)

    # # # peak
    # # path_peak1 = os.path.join(path_out, 'peaks', 'clipper')
    # # peak_clipper_files = call_peak(genome, map_bam_files, path_peak1, 
    # #                                tool='clipper')
    # path_peak2 = os.path.join(path_out, 'peaks', 'pyicoclip')
    # peak_pyicoclip_files = call_peak(genome, map_bam_files, path_peak2, 
    #                                  tool='pyicoclip')

    # # rtstop
    # path_rtstop = os.path.join(path_out, 'rtstops')
    # rtstop_files = call_rtstop(map_bed_files, path_rtstop,  smp_name,
    #                            threshold=threshold, intersect=intersect,
    #                            overwrite=overwrite)

    # report
    tmp = Goldclip_output(path_out, smp_name, genome).get_all_figures()
    # goldclip_report(path_out, smp_name, genome)



