#!/usr/bin/env python
# -*- encoding: utf-8 -*-


__author__ = 'Ming Wang <wangm08@hotmail.com>'
__copyright__ = '2018 by Ming Wang <wangm08@hotmail.com>'
__license__ = 'MIT'
__email__ = 'wangm08@hotmail.com'
__version__ = '0.0.2'


import sys
# from goldclip.goldcliplib.demx import *
# from goldclip.goldcliplib.aligner import *
from goldclip.helper import BAM
from goldclip.goldcliplib.alignment import Alignment
from goldclip.goldcliplib.trim import Trimmer
from goldclip.goldcliplib.peak import *
from goldclip.goldcliplib.rtstop import *
from goldclip.goldcliplib.run import *
from goldclip.goldcliplib.report import *


# class Demx:

#     """
#     processing GoldCLIP illumina datasets
#     only one of the PE reads
#     ## type1: goldclip_version_1
#     read1: {NNN} - {bc} - {NN} - <insert>

#     ## type2: goldclip_version_2
#     read1: {N10} - <insert> - A{barcode}
#     read2: {barcode}A - <insert> - {N10}

#     ## type3: eCLIP
#     read1: {barcode} - <insert> - {N10}
#     read2: {N10} - <insert> - {bracode}

#     """

#     def __init__(self, *args, **kwargs):
#         self.kwargs = kwargs


#     def run(self):
#         """run demx"""
#         r1 = self.kwargs['fq1']
#         r2 = self.kwargs['fq2']
#         barcode = self.kwargs['bc_file']
#         bc_in_read12 = self.kwargs['bc_in_read12']
#         path_out = self.kwargs['out']
#         n_left = self.kwargs['n_left']
#         n_right = self.kwargs['n_right']
#         is_bioawk = self.kwargs['bioawk']
#         bc = self.kwargs['bc_only']
#         p7 = self.kwargs['p7_only']
#         p7_and_bc = self.kwargs['p7_and_bc']
#         mm = self.kwargs['n_mismatch']
#         cut = self.kwargs['cut']
#         # demx p7, then barcode
#         read1 = r1.name
#         barcode_file = barcode.name
#         assert is_path(path_out)
#         if p7_and_bc: # demx both p7 and barcode
#             if r2:
#                 logging.info('demx P7 and barcode, PE reads')
#                 read2 = r2.name
#                 tmp = p7_bc_demx_pe(read1, read2, barcode_file, path_out,
#                                     n_left, n_right, 
#                                     bc_in_read12=bc_in_read12,cut=cut, mm=mm)
#             else:
#                 logging.info('demx P7 and barcode, SE reads')
#                 tmp = p7_bc_demx_se(read1, barcode_file, path_out, n_left, n_right,
#                                     cut=cut, mm=mm)
#         elif p7: # require demx p7, in fastq-comment-field
#             if r2:
#                 logging.info('demx P7, PE reads')
#                 read2 = r2.name
#                 tmp = p7_demx_pe(read1, read2, barcode_file, path_out, 
#                                  bc_in_read12=bc_in_read12, mm=mm)
#             else:
#                 logging.info('demx P7, SE reads')
#                 tmp = p7_demx_se(read1, barcode_file, path_out, mm)
#         else: # only barcode
#             if r2:
#                 logging.info('demx barcode, PE reads')
#                 read2 = r2.name
#                 tmp = bc_demx_pe(read1, read2, barcode_file, path_out, n_left, 
#                                  n_right, bc_in_read12=bc_in_read12, cut=cut, 
#                                  mm=mm)
#             else:
#                 if is_bioawk:
#                     logging.info('demx barcode, SE reads - bioawk')
#                     tmp = demx_se_bioawk(read1, barcode_file, path_out, n_left,
#                                          n_right)
#                 else:
#                     logging.info('demx barcode, SE reads')
#                     tmp = bc_demx_se(read1, barcode_file, path_out, n_left,
#                                      n_right, cut=cut, mm=mm)
#         logging.info('demx finish!')


class Trim:

    def __init__(self, *args, **kwargs):
        self.kwargs = kwargs

    def run(self):
        args = self.kwargs
        fq1_files = [f.name for f in args['i']]

        ## SE mode
        if args['fq2'] is None: 
            for fq1 in fq1_files:
                tmp = Trimmer(fq1, 
                    adapter3=args['a'], 
                    path_out=args['o'], 
                    len_min=args['m'],
                    adapter5=args['g'], 
                    read12=args['read12'], 
                    qual_min=args['q'], 
                    error_rate=args['e'], 
                    overlap=args['O'],
                    rm_untrim=args['rm_untrim'], 
                    threads=args['threads'], 
                    overwrite=args['overwrite'], 
                    keep_name=args['keep_name'],
                    adapter_sliding=args['adapter_sliding'], 
                    trim_times=args['trim_times'],
                    double_trim=args['double_trim'],
                    rm_dup=args['rm_dup'],
                    cut_before_trim=args['cut_before_trim'],
                    cut_after_trim=args['cut_after_trim'],
                    trim_to_length=args['trim_to_length']).run()
        ## PE mode
        else:
            fq2_files = [f.name for f in args['fq2']]
            for fq1, fq2 in zip(fq1_files, fq2_files):
                tmp = Trimmer(fq1, 
                    adapter3=args['a'], 
                    path_out=args['o'], 
                    len_min=args['m'],
                    adapter5=args['g'], 
                    read12=args['read12'], 
                    fq2=fq2, 
                    AD3=args['A'], 
                    AD5=args['G'],
                    qual_min=args['q'], 
                    error_rate=args['e'], 
                    overlap=args['O'],
                    rm_untrim=args['rm_untrim'], 
                    threads=args['threads'], 
                    overwrite=args['overwrite'], 
                    keep_name=args['keep_name'],
                    adapter_sliding=args['adapter_sliding'], 
                    trim_times=args['trim_times'],
                    double_trim=args['double_trim'],
                    rm_dup=args['rm_dup'],
                    cut_before_trim=args['cut_before_trim'],
                    cut_after_trim=args['cut_after_trim'],
                    trim_to_length=args['trim_to_length']).run()
        logging.info('trimming finish!')


class Align:
    """
    Mapping SE reads to reference genome
    specify: fq, path_out, index, parameters, tools
    """

    def __init__(self, **kwargs):
        self.kwargs = kwargs


    def run(self):
        args = self.kwargs
        fqs = [f.name for f in self.kwargs['i']]
        tmp = Alignment(fqs,
            args['o'],
            smp_name=args['n'],
            genome=args['g'],
            spikein=args['k'], 
            index_ext=args['x'],
            threads=args['threads'], 
            unique_only=args['unique_only'],
            n_map=args['n_map'],
            aligner=args['aligner'],
            align_to_rRNA=args['align_to_rRNA'],
            repeat_masked_genome=args['repeat_masked_genome'],
            path_data=args['path_data'],
            overwrite=args['overwrite']).run()
        logging.info('mapping finish!')
        return tmp[0]


class Peak:
    """
    call peaks using CLIPper, pyicoclip
    """
    def __init__(self, *args, **kwargs):
        self.kwargs = kwargs

    def run(self):
        logging.info('peak-calling start')
        bam_files = [f.name for f in self.kwargs['i']]
        genome = self.kwargs['g']
        path_out = self.kwargs['o']
        tool = self.kwargs['peak_caller']
        peak_files = call_peak(genome, bam_files, path_out, tool)
        logging.info('peak-calling finish')
        return peak_files


class Rtstop:
    """
    call RT-stops from BAM files
    """
    def __init__(self, *args, **kwargs):
        self.kwargs = kwargs

    def run(self):
        args = self.kwargs

        logging.info('RTStop-calling start')
        bed_files = []

        # convert bam to bed
        for i in args['i']:
            ifile = i.name
            if ifile.endswith('bam'):
                logging.info('convert BAM to BED: %s' % ifile)
                ibed = BAM(ifile).to_bed()
                bed_files.append(ibed)
            elif ifile.endswith('bed'):
                bed_files.append(ifile)
            else:
                continue

        tmp = call_rtstop(bed_files=bed_files, 
            path_out=args['o'],
            smp_name=args['n'],
            threshold=args['t'],
            intersect=args['c'],
            overwrite=args['f'])
        logging.info('RTstop-calling finish')

        return tmp


class Run_all:
    """
    call RT-stops from BAM files
    """
    def __init__(self, *args, **kwargs):
        self.kwargs = kwargs


    def run(self):
        logging.info('GoldCLIP start')

        tmp = run_goldclip(
            fq_files = [f.name for f in self.kwargs['i']],
            path_out = self.kwargs['o'],
            genome = self.kwargs['g'],
            smp_name = self.kwargs['n'],
            spikein = self.kwargs['k'],
            is_trimmed = self.kwargs['trimmed'],
            ad3 = self.kwargs['a'],
            read12 = self.kwargs['read12'],
            len_min = self.kwargs['m'],
            qual_pct = self.kwargs['p'],
            qual_min = self.kwargs['q'],
            err_rate = self.kwargs['e'],
            overlap = self.kwargs['O'],
            rm_untrim = self.kwargs['rm_untrim'],
            rm_dup = self.kwargs['rm_dup'],
            cut_before_trim = self.kwargs['cut_before_trim'],
            cut_after_trim = self.kwargs['cut_after_trim'],
            aligner = self.kwargs['aligner'],
            threshold = self.kwargs['t'],
            intersect = self.kwargs['c'],
            path_data = self.kwargs['path_data'],
            threads = self.kwargs['threads'],
            overwrite = self.kwargs['overwrite'])

        logging.info('GoldCLIP finish')

        return tmp


class Report:
    """
    create report of goldclip
    """
    def __init__(self, *args, **kwargs):
        self.kwargs = kwargs

    def run(self):
        path_out = self.kwargs['path']
        smp_name = self.kwargs['name']
        genome = self.kwargs['genome']
        tmp = goldclip_report(path_out, smp_name, genome)
        return tmp

