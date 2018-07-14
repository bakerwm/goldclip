#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Mapping reads to reference genome
1. to spike-in genome
2. to repbase
3. to rRNA, tRNA, MT
4. to genome
"""
# TO-DO
# count RT in rep1/2


__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__date__ = "2018-03-21"
__version__ = "0.1"


from goldclip.goldcliplib.map import *


class Map:
    """
    Mapping SE reads to reference genome
    specify: fq, path_out, index, parameters, tools
    """

    def __init__(self, *args, **kwargs):
        self.kwargs = kwargs


    def run():
        fqs = [f.name for f in self.kwargs['i']]
        out_path = self.kwargs['o']
        genome = self.kwargs['g']
        genome_index = self.kwargs['x']
        smp_name = self.kwargs['n']
        multi_cores = self.kwargs['threads']
        aligner = self.kwargs['aligner']

        tmp = map_se(fqs, smp_name, out_path, genome, spikein, genome_index,
                      multi_cores)



## EOF
