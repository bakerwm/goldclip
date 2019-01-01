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
from goldclip.goldcliplib.default_arguments import Argument


# def run_goldclip(fqs, path_out, genome, smp_name, library_type=1, **kwargs):

# class Goldclip_all_in_one(object):
#     """Run goldclip pipeline in one-step
#     """

#     def __init__(self, fqs, path_out, genome, smp_name, library_type=1, **kwargs):
#         """Fetch all arguments for goldclip in one step"""
#         self.fqs = fqs
#         self.path_out = path_out,
#         self.genome = genome
#         self.smp_name = smp_name
#         self.library_type = library_type
#         self.kwargs = kwargs

#     ## run
#     def run(self):


#         ##Trim adapters
#         ##Trim and Collapse
#         args = Argument(self.library_type).all()
#         args = {**args, **self.kwargs} # update arguments
#         print('' %self.pat_out)
#         path_trim = os.path.join(self.path_out, '01.trimming')
#         if args['is_trimmed']:
#             # make links
#             trim_fq_files = self.fqs
#         else:
#             trim_fq_files = Trimmer(self.fqs, path_out=path_trim, **args).run()

#         ##Map reads
#         logging.info('02.Alignment')
#         path_map = os.path.join(self.path_out, '02.genome_mapping')
#         bam_files, _ = Alignment(trim_fq_files,
#             path_map,
#             self.smp_name,
#             self.genome,
#             **args).run()

#         ##Call peaks
#         logging.info('03.Call Peaks')
#         path_peak1 = os.path.join(self.path_out, '03.call_peaks', 'clipper')
#         path_peak2 = os.path.join(self.path_out, '03.call_peaks', 'pyicoclip')
#         peak_files1 = call_peak(self.genome, bam_files, path_peak1, peak_caller='clipper')
#         peak_files1 = call_peak(self.genome, bam_files, path_peak1, peak_caller='pyicoclip')

#         ##Call rtstops
#         logging.info('04.Call RT-Stops')
#         path_rtstop = os.path.join(self.path_out, '04.call_rtstops')
#         bed_files = []
#         # convert bam to bed
#         for bam in bam_files:
#             bed = os.path.splitext(bam)[0] + '.bed'
#             if not os.path.exists(bed):
#                 logging.info('convert BAM to BED: %s' % bam)
#                 bed = BAM(bam).to_bed()
#             bed_files.append(bed)
#         rtstop_files = call_rtstop(
#             bed_files=bed_files, 
#             path_out=path_rtstop,
#             smp_name=self.smp_name, 
#             threshold=args['threshold'], 
#             intersect=args['intersect'],
#             overwrite=args['overwrite'])

#         ## report
#         logging.info('05.Generate report')
#         Goldclip_output(
#             project_path=self.path_out,
#             project_name=self.smp_name,
#             genome=self.genome,
#             threads=self.threads).get_all_figures()

#     # # map
#     # path_map = os.path.join(path_out, 'genome_mapping')

#     # map_bam_files = Alignment(clean_fq_files,
#     #         path_map,
#     #         smp_name=smp_name,
#     #         genome=genome, **args).run()

#     # # bigWig
#     # path_bw = os.path.join(path_out, 'bigWig')
#     # for bam in map_bam_files:
#     #     logging.info('making bigWig: %s' % os.path.basename(bam))
#     #     bam2bw(bam, genome, path_bw, strandness=True, binsize=10)

#     # # peak
#     # path_peak1 = os.path.join(path_out, 'peaks', 'clipper')
#     # peak_clipper_files = call_peak(genome, map_bam_files, path_peak1, 
#     #                                tool='clipper')
#     # path_peak2 = os.path.join(path_out, 'peaks', 'pyicoclip')
#     # peak_pyicoclip_files = call_peak(genome, map_bam_files, path_peak2, 
#     #                                  tool='pyicoclip')

#     # # rtstop
#     # path_rtstop = os.path.join(path_out, 'rtstops')
#     # rtstop_files = call_rtstop(map_bed_files, path_rtstop,  smp_name,
#     #                            threshold=threshold, intersect=intersect,
#     #                            overwrite=overwrite)

#     # # report
#     # tmp = Goldclip_output(path_out, smp_name, genome).get_all_figures()
#     # goldclip_report(path_out, smp_name, genome)



