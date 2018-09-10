#!/usr/bin/env python
# -*- encoding: utf-8 -*-


__author__ = 'Ming Wang <wangm08@hotmail.com>'
__copyright__ = '2018 by Ming Wang <wangm08@hotmail.com>'
__license__ = 'MIT'
__email__ = 'wangm08@hotmail.com'
__version__ = '0.0.1'



from goldclip.goldcliplib.demx import *
from goldclip.goldcliplib.trim import *
from goldclip.goldcliplib.aligner import *
from goldclip.goldcliplib.peak import *
from goldclip.goldcliplib.rtstop import *
from goldclip.goldcliplib.run import *
from goldclip.goldcliplib.report import *


class Demx:

    """
    processing GoldCLIP illumina datasets
    only one of the PE reads
    ## type1: goldclip_version_1
    read1: {NNN} - {bc} - {NN} - <insert>

    ## type2: goldclip_version_2
    read1: {N10} - <insert> - A{barcode}
    read2: {barcode}A - <insert> - {N10}

    ## type3: eCLIP
    read1: {barcode} - <insert> - {N10}
    read2: {N10} - <insert> - {bracode}

    """

    def __init__(self, *args, **kwargs):
        self.kwargs = kwargs


    def run(self):
        """run demx"""
        r1 = self.kwargs['fq1']
        r2 = self.kwargs['fq2']
        barcode = self.kwargs['bc_file']
        bc_in_read12 = self.kwargs['bc_in_read12']
        path_out = self.kwargs['out']
        n_left = self.kwargs['n_left']
        n_right = self.kwargs['n_right']
        is_bioawk = self.kwargs['bioawk']
        bc = self.kwargs['bc_only']
        p7 = self.kwargs['p7_only']
        p7_and_bc = self.kwargs['p7_and_bc']
        mm = self.kwargs['n_mismatch']
        cut = self.kwargs['cut']
        # demx p7, then barcode
        read1 = r1.name
        barcode_file = barcode.name
        assert is_path(path_out)
        if p7_and_bc: # demx both p7 and barcode
            if r2:
                logging.info('demx P7 and barcode, PE reads')
                read2 = r2.name
                tmp = p7_bc_demx_pe(read1, read2, barcode_file, path_out,
                                    n_left, n_right, 
                                    bc_in_read12=bc_in_read12,cut=cut, mm=mm)
            else:
                logging.info('demx P7 and barcode, SE reads')
                tmp = p7_bc_demx_se(read1, barcode_file, path_out, n_left, n_right,
                                    cut=cut, mm=mm)
        elif p7: # require demx p7, in fastq-comment-field
            if r2:
                logging.info('demx P7, PE reads')
                read2 = r2.name
                tmp = p7_demx_pe(read1, read2, barcode_file, path_out, 
                                 bc_in_read12=bc_in_read12, mm=mm)
            else:
                logging.info('demx P7, SE reads')
                tmp = p7_demx_se(read1, barcode_file, path_out, mm)
        else: # only barcode
            if r2:
                logging.info('demx barcode, PE reads')
                read2 = r2.name
                tmp = bc_demx_pe(read1, read2, barcode_file, path_out, n_left, 
                                 n_right, bc_in_read12=bc_in_read12, cut=cut, 
                                 mm=mm)
            else:
                if is_bioawk:
                    logging.info('demx barcode, SE reads - bioawk')
                    tmp = demx_se_bioawk(read1, barcode_file, path_out, n_left,
                                         n_right)
                else:
                    logging.info('demx barcode, SE reads')
                    tmp = bc_demx_se(read1, barcode_file, path_out, n_left,
                                     n_right, cut=cut, mm=mm)
        logging.info('demx finish!')



class Trim:

    def __init__(self, *args, **kwargs):
        self.kwargs = kwargs

    def run(self):
        logging.info('trimming files')
        fq_files = [f.name for f in self.kwargs['i']]
        ad3 = self.kwargs['a']
        path_out = self.kwargs['o']
        len_min = self.kwargs['m']
        qual_pct = self.kwargs['p']
        qual_min = self.kwargs['q']
        overlap = self.kwargs['O']
        err_rate = self.kwargs['e']
        threads = self.kwargs['threads']
        read12 = self.kwargs['read12']
        double_trim = self.kwargs['double_trim']
        rm_untrim = self.kwargs['rm_untrim'],
        rm_dup = self.kwargs['rm_dup']
        cut_before_trim = self.kwargs['cut_before_trim']
        cut_after_trim = self.kwargs['cut_after_trim']
        overwrite = self.kwargs['overwrite']
        rm_untrim = rm_untrim[0] # is tuple, not bool? !!!! why?
        tmp = trim(fq_files, adapter3=ad3, path_out=path_out, len_min=len_min,
                   double_trim=double_trim, qual_min=qual_min, 
                   err_rate=err_rate, overlap=overlap, multi_cores=threads, 
                   read12=read12, rm_untrim=rm_untrim, rm_dup=rm_dup, 
                   cut_before_trim=cut_before_trim, 
                   cut_after_trim=cut_after_trim,
                   overwrite=overwrite,)
        logging.info('trimming finish!')



class Align:
    """
    Mapping SE reads to reference genome
    specify: fq, path_out, index, parameters, tools
    """

    def __init__(self, *args, **kwargs):
        self.kwargs = kwargs


    def run(self):
        logging.info('mapping files')
        fqs = [f.name for f in self.kwargs['i']]
        smp_name = self.kwargs['n']
        path_out = self.kwargs['o']
        genome = self.kwargs['g']
        spikein = self.kwargs['k']
        unique_only = self.kwargs['unique_only']
        align_to_rRNA = self.kwargs['align_to_rRNA']
        multi_cores = self.kwargs['threads']
        aligner = self.kwargs['aligner']
        path_data = self.kwargs['path_data']
        overwrite = self.kwargs['overwrite']
        tmp = align(fqs, smp_name, path_out, genome, spikein, 
                    unique_only=unique_only, 
                    align_to_rRNA=align_to_rRNA,
                    multi_cores=multi_cores,
                    aligner=aligner, 
                    path_data=path_data, 
                    overwrite=overwrite)
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
        tool = self.kwargs['tool']
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
        logging.info('RTStop-calling start')
        bed_files = [f.name for f in self.kwargs['i']]
        path_out = self.kwargs['o']
        smp_name = self.kwargs['n']
        threshold = self.kwargs['t']
        intersect = self.kwargs['c']
        overwrite = self.kwargs['f']
        tmp = call_rtstop(bed_files, path_out, smp_name, threshold,
                          intersect, overwrite)
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
        fq_files = [f.name for f in self.kwargs['i']]
        path_out = self.kwargs['o']
        genome = self.kwargs['g']
        smp_name = self.kwargs['n']
        spikein = self.kwargs['k']
        is_trimmed = self.kwargs['trimmed']
        ad3 = self.kwargs['a']
        read12 = self.kwargs['read12']
        len_min = self.kwargs['m']
        qual_pct = self.kwargs['p']
        qual_min = self.kwargs['q']
        err_rate = self.kwargs['e']
        overlap = self.kwargs['O']
        rm_untrim = self.kwargs['rm_untrim'],
        rm_dup = self.kwargs['rm_dup']
        cut_before_trim = self.kwargs['cut_before_trim']
        cut_after_trim = self.kwargs['cut_after_trim']
        aligner = self.kwargs['aligner']
        unique_only = self.kwargs['unique_only']
        align_to_rRNA = self.kwargs['align_to_rRNA']
        threshold = self.kwargs['t']
        intersect = self.kwargs['c']
        path_data = self.kwargs['path_data']
        threads = self.kwargs['threads']
        overwrite = self.kwargs['overwrite']
        tmp = run_goldclip(fq_files, path_out, genome, smp_name, spikein,
                           is_trimmed, ad3, read12, len_min, qual_pct,
                           qual_min, err_rate, overlap, rm_untrim, rm_dup,
                           cut_before_trim, cut_after_trim, aligner, threshold,
                           intersect, path_data, threads, overwrite)
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

