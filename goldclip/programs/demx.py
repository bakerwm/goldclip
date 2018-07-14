#!/usr/bin/env python
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

__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__date__ = "2018-04-01"
__version__ = "0.1"


from goldclip.goldcliplib.demx import *


class Demx:

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

## EOF
