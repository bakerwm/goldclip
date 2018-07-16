#!/usr/bin/env python
# -*- encoding: utf-8 -*-


__author__ = 'Ming Wang <wangm08@hotmail.com>'
__copyright__ = '2018 by Ming Wang <wangm08@hotmail.com>'
__license__ = 'MIT'
__email__ = 'wangm08@hotmail.com'
__version__ = '0.0.1'



from goldclip.goldcliplib.demx import *
from goldclip.goldcliplib.trim import *
from goldclip.goldcliplib.map import *
from goldclip.goldcliplib.peak import *



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
        path_out = self.kwargs['out']
        n_left = self.kwargs['n_left']
        n_right = self.kwargs['n_right']
        is_bioawk = self.kwargs['bioawk']
        p7 = self.kwargs['p7']
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
                                    n_left, n_right, cut=cut, mm=mm)
            else:
                logging.info('demx P7 and barcode, SE reads')
                tmp = p7_bc_demx_se(read1, barcode_file, path_out, n_left, n_right,
                                    cut=cut, mm=mm)
        elif p7: # require demx p7, in fastq-comment-field
            if r2:
                logging.info('demx P7, PE reads')
                read2 = r2.name
                tmp = p7_demx_pe(read1, read2, barcode_file, path_out, mm)
            else:
                logging.info('demx P7, SE reads')
                tmp = p7_demx_se(read1, barcode_file, path_out, mm)
        else: # only barcode
            if r2:
                logging.info('demx barcode, PE reads')
                read2 = r2.name
                tmp = bc_demx_pe(read1, read2, barcode_file, path_out, n_left, 
                                 n_right, cut=cut, mm=mm)
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
        rm_untrim = self.kwargs['rm_untrim'],
        rm_dup = self.kwargs['rm_dup']
        cut_before_trim = self.kwargs['cut_before_trim']
        cut_after_trim = self.kwargs['cut_after_trim']
        overwrite = self.kwargs['overwrite']
        tmp = trim(fq_files, path_out, adapter3=ad3, len_min=len_min, 
                   qual_min=qual_min, 
                   err_rate=err_rate, multi_cores=threads,
                   rm_untrim=rm_untrim,
                   overlap=overlap, read12=read12, overwrite=overwrite,
                   rm_dup=rm_dup, cut_before_trim=cut_before_trim,
                   cut_after_trim=cut_after_trim)
        logging.info('trimming finish!')



class Map:
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
        multi_cores = self.kwargs['threads']
        aligner = self.kwargs['aligner']
        path_data = self.kwargs['path_data']
        overwrite = self.kwargs['overwrite']
        tmp = map(fqs, smp_name, path_out, genome, spikein, 
                     multi_cores=multi_cores, aligner=aligner,
                     path_data=path_data, overwrite=overwrite)
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




