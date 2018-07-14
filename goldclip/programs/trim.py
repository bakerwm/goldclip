#!/usr/bin/env python
"""
Trimming reads
1. cut 3-adapter
2. trim low quality bases
3. remove N reads
4. trim N-bases from either ends of reads
5. limit read length
"""

__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__date__ = "2018-03-21"
__version__ = "0.1"


from goldclip.goldcliplib.trim import *


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

## EOF
