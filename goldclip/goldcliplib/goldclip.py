#!/usr/bin/env python
# -*- encoding: utf-8 -*-


__author__ = 'Ming Wang <wangm08@hotmail.com>'
__copyright__ = '2018 by Ming Wang <wangm08@hotmail.com>'
__license__ = 'MIT'
__email__ = 'wangm08@hotmail.com'
__version__ = '0.0.2'


import os
import sys
from goldclip.helper import BAM, logging
from goldclip.goldcliplib.arguments import args_init
from goldclip.goldcliplib.trim import Trimmer
from goldclip.goldcliplib.alignment import Alignment
from goldclip.goldcliplib.peak import call_peak
from goldclip.goldcliplib.rtstop import call_rtstop
from goldclip.goldcliplib.report import Goldclip_report


class Trim:
    def __init__(self, **kwargs):
        self.kwargs = args_init(kwargs, trim=True)

    def run(self):
        logging.info('trimming start')
        args = self.kwargs
        fq1_files = args.pop('fq1', None) # remove 'fq1' from args

        ## SE mode
        if args['fq2'] is None: 
            for fq1 in fq1_files:
                tmp = Trimmer(fq1=fq1, **args).run()
        ## PE mode
        ## !!!! goldclip only works on SE reads !!!! ##
        else:
            fq2_files = args.pop('fq2', None) # remove 'fq2' from args
            for fq1, fq2 in zip(fq1_files, fq2_files):
                tmp = Trimmer(fq1=fq1, fq2=fq2, **args).run()
        logging.info('trimming finish!')


class Align:
    """
    Mapping SE reads to reference genome
    specify: fq, path_out, index, parameters, 
    """
    def __init__(self, **kwargs):
        self.kwargs = args_init(kwargs, align=True)


    def run(self):
        tmp, _ = Alignment(**self.kwargs).run() # genome_bam, extra_bam
        logging.info('mapping finish!')
        return tmp[0]


class Peak:
    """Call peaks using CLIPper, pyicoclip
    """
    def __init__(self, **kwargs):
        self.kwargs = args_init(kwargs, call_peak=True)

    def run(self):
        args = self.kwargs

        logging.info('peak-calling start')
        peak_files = call_peak(
            genome=args['genome'], 
            bam_ins=args['bam_files'],
            path_out=args['path_out'], 
            peak_caller=args['peak_caller'])
        logging.info('peak-calling finish')
        return peak_files


class Rtstop:
    """
    call RT-stops from BAM files
    """
    def __init__(self, **kwargs):
        self.kwargs= args_init(kwargs, call_peak=True)

    def run(self):
        args = self.kwargs

        logging.info('RTStop-calling start')
        bed_files = []
        # convert bam to bed
        for b in args['bed_files']:
            if b.endswith('.bam'):
                logging.info('convert BAM to BED: %s' % b)
                bed = BAM(b).to_bed()
                bed_files.append(bed)
            elif b.endswith('.bed'):
                bed_files.append(b)
            else:
                continue

        tmp = call_rtstop(
            bed_files=bed_files, 
            path_out=args['path_out'],
            smp_name=args['smp_name'],
            threshold=args['threshold'],
            intersect=args['intersect'],
            overwrite=args['overlap'])
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
        self.kwargs = args_init(kwargs, demx=True, trim=True, align=True, call_peak=True)

    def run(self):
        args = self.kwargs
        Goldclip_report(**args).get_all_figures()


class Goldclip_all_in_one(object):
    """Run goldclip pipeline in one-step
    01.trimming
    02.genome_mapping
    03.call_peaks
    04.call_rtstops
    05.report

    """
    def __init__(self, **kwargs):
        """Fetch all arguments for goldclip in one step"""
        self.kwargs = args_init(kwargs, demx=True, trim=True, align=True, call_peak=True)

    ## run
    def run(self):
        args = self.kwargs

        ## Trim adapters
        logging.info('01.Trimming')
        args_trim = args.copy()
        path_trim = os.path.join(args_trim['path_out'], '01.trimming')
        fq1_files = args_trim.pop('fq1', None) # remove 'fq1'
        args_trim['path_out'] = path_trim


        if args['trimmed']:
            # make links
            trim_fq_files = fq1_files
        else:
            trim_fq_files = []
            for fq1 in fq1_files:
                tmp = Trimmer(fq1=fq1, **args_trim).run()
                trim_fq_files.append(tmp)

        ## Map reads
        logging.info('02.Alignment')
        path_map = os.path.join(args['path_out'], '02.genome_mapping')
        ## update arguments
        args_map = args.copy()
        args_map['fqs'] = trim_fq_files
        args_map['path_out'] = path_map
        bam_files, _ = Alignment(**args_map).run()

        ## Call peaks
        logging.info('03.Call Peaks')
        ## clipper
        path_peak1 = os.path.join(args['path_out'], '03.call_peaks', 'clipper')
        peak_files1 = call_peak(args['genome'], bam_files, path_peak1, peak_caller='clipper')
        
        ## pyicoclip
        path_peak2 = os.path.join(args['path_out'], '03.call_peaks', 'pyicoclip')
        peak_files1 = call_peak(args['genome'], bam_files, path_peak2, peak_caller='pyicoclip')

        ## Call rtstops
        logging.info('04.Call RT-Stops')
        path_rtstop = os.path.join(args['path_out'], '04.call_rtstops')
        bed_files = []
        # convert bam to bed
        for bam in bam_files:
            bed = os.path.splitext(bam)[0] + '.bed'
            if not os.path.exists(bed):
                logging.info('convert BAM to BED: %s' % bam)
                bed = BAM(bam).to_bed()
            bed_files.append(bed)
        rtstop_files = call_rtstop(
            bed_files=bed_files, 
            path_out=path_rtstop,
            smp_name=args['smp_name'], 
            threshold=args['threshold'], 
            intersect=args['intersect'],
            overwrite=args['overwrite'])

        ## Report
        logging.info('05.Generate report')
        Goldclip_report(
            project_path=args['path_out'],
            project_name=args['smp_name'],
            genome=args['genome'],
            threads=args['threads']).get_all_figures()


