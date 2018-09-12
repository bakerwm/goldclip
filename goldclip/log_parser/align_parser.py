#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import os
import re
import multiprocessing as mp
from goldclip.helper import *
from goldclip.bin.bed_annotation import *



class Bowtie_log(object):
	"""
	processing log of bowtie
	report: total reads, unique mapped reads, multiple mapped reads
	"""

    def __init__(self):
        assert isinstance(genome, str)
        self.kwargs = kwargs
       

    def count(self, group='total'):
    	"""Return the number of reads for each categories
    	group: total, unmapped, unique, multiple
    	"""
    	assert isinstance(group, str)
    	if group in self.log:
    		_num1 = re.sub(',', '', line.split()[-1])

    		

    def .is_file(self):
    	"""Check the log file is exists, not empty
    	"""
    	pass




        



## wrapper functions for bowtie mapping ##
def bowtie_log_parser(path):
    """
    Parsing the log file of bowtie (to stderr)
    fetch Input, mapped, unmapped reads
    save in JSON format
    return dict of all values
    """
    # !!! to-do
    # convert to DataFrame, json, dict
    # unmap = input - reported
    logdict = {}
    _input = []
    _one_hit = []
    _not_hit = []
    _rpt = []
    with open(path, 'r') as f:
        for line in f.readlines():
            if 'reads processed' in line:
                _num1 = re.sub(',', '', line.split()[-1])
                _input.append(int(_num1))
            elif 'at least one reported' in line:
                _num2 = re.sub(',', '', line.split()[-2])
                _one_hit.append(int(_num2))
            elif 'Reported' in line:
                _num4 = re.sub(',', '', line.split()[1])
                _rpt.append(int(_num4))
            elif 'failed to align' in line:
                _num3 = re.sub(',', '', line.split()[-2])
                _not_hit.append(int(_num3))
            else:
                continue
    logdict['input_reads'] = _input[0] # first one
    # logdict['mapped'] = sum(_one_hit) # sum all
    logdict['mapped'] = sum(_rpt) # sum all
    logdict['unmapped'] = int(_input[-1]) - int(_one_hit[-1]) # -m suppress
    logdict['map_pct'] = '{:.2f}%'.\
        format(int(logdict['mapped']) / int(logdict['input_reads'])*100)
    json_out = os.path.splitext(path)[0] + '.json'
    with open(json_out, 'w') as fo:
        json.dump(logdict, fo, indent = 4)
    return logdict


def bowtie2_log_parser(path):
    """
    Parsing the log file of bowtie
    fetch Input, unique, multiple, unmapped
    save in JSON format
    return dict of all values
    """
    logdict = {}
    _input = []
    _one_hit = []
    _multi_hit = []
    _not_hit = []
    _rpt = []
    with open(path, 'rt') as fi:
        for line in fi.readlines():
            line = line.strip()
            _num = line.split(' ')[0]
            _num = _num.strip('%')
            _num = float(_num)
            if line.endswith('reads; of these:'):
                _input.append(_num)
            elif ') aligned 0 times' in line:
                _not_hit.append(_num)
            elif ') aligned exactly 1 time' in line:
                _one_hit.append(_num)
            elif ') aligned >1 times' in line:
                _multi_hit.append(_num)
            else:
                continue
    # save to dict
    logdict['input_reads'] = int(_input[0]) # first one
    logdict['mapped'] = int(sum(_one_hit + _multi_hit))
    logdict['unique'] = int(sum(_one_hit))
    logdict['multi'] = int(sum(_multi_hit))
    logdict['unmapped'] = int(_not_hit[-1]) # -m suppress
    logdict['map_pct'] = '{:.2f}%'.\
        format(int(logdict['mapped']) / int(logdict['input_reads'])*100)
    json_out = os.path.splitext(path)[0] + '.json'
    with open(json_out, 'w') as fo:
        json.dump(logdict, fo, indent=4)
    return logdict


def star_log_parser(path):
    logdict = {}
    with open(path) as f:
        for line in f:
            sep = line.strip().split('|')
            if 'Number of input reads' in line:
                logdict['input_reads'] = int(sep[1].strip())
            elif 'Uniquely mapped reads number' in line:
                logdict['unique'] = int(sep[1].strip())
            elif 'Number of reads mapped to multiple loci' in line:
                logdict['multi'] = int(sep[1].strip())
            else:
                pass
    logdict['mapped'] = logdict['unique'] + logdict['multi']
    logdict['unmapped'] = logdict['input_reads'] - logdict['mapped']
    logdict['map_pct'] = '{:.2f}%'.format(logdict['mapped'] / logdict['input_reads'] * 100)
    json_out = os.path.splitext(path)[0] + '.json'
    with open(json_out, 'wt') as fo:
        json.dump(logdict, fo, indent=4)
    return logdict