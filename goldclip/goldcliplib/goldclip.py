#!/usr/bin/env python
# -*- encoding: utf-8 -*-


__author__ = 'Ming Wang <wangm08@hotmail.com>'
__copyright__ = '2018 by Ming Wang <wangm08@hotmail.com>'
__license__ = 'MIT'
__email__ = 'wangm08@hotmail.com'
__version__ = '0.0.2'


import sys
from goldclip.helper import BAM
from goldclip.goldcliplib.trim import Trimmer
from goldclip.goldcliplib.alignment import Alignment
from goldclip.goldcliplib.peak import *
from goldclip.goldcliplib.rtstop import *
from goldclip.goldcliplib.run import *
from goldclip.goldcliplib.report import *
from goldclip.goldcliplib.default_arguments import Argument


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
    specify: fq, path_out, index, parameters, 
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
        peak_caller = self.kwargs['peak_caller']
        peak_files = call_peak(genome, bam_files, path_out, peak_caller)
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


class Report:
    """
    create report of goldclip

    args : project_path, the directory of goldclip output
    args : project_name, the smp_name of the project, -n in Alignment
    args : g, the reference genome of the project
    """
    def __init__(self, **kwargs):
        self.kwargs = kwargs

    def run(self):
        args = self.kwargs
        Goldclip_output(
            project_path=args['project_path'],
            project_name=args['project_name'],
            genome=args['g'],
            threads=8).get_all_figures()


class Run_all:
    """
    call RT-stops from BAM files
    """
    def __init__(self,  **kwargs):
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

