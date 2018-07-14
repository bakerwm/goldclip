#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import os
import re
from goldclip.helper import *

__author__ = 'Ming Wang <wangm08@hotmail.com>'
__copyright__ = '2018 by Ming Wang <wangm08@hotmail.com>'
__license__ = 'MIT'
__email__ = 'wangm08@hotmail.com'
__version__ = '0.0.1'



## wrapper functions for cutadapt ##
def cutadapt_log_parser(log):
    """
    Parsing the output of cutadpat
    support multiple run output
    save data in JSON format
    """
    logdict = {}
    _cutadapt_ver = []
    _python_ver = []
    _cmd = []
    _fname = []
    _raw = []
    _clean = []
    outfile = os.path.splitext(log)[0] + ".json"
    with open(log, 'r') as f:
        for line in f.readlines():
            if(len(re.findall('^This is cutadapt', line)) == 1):
                r1 = re.findall(r'(cutadapt|Python) (\d+\.\d+\.?\d+?)', line) # 1st line, tools
                for i in r1: logdict[i[0]] = i[1]
            elif('Command line parameters' in line):
                _cmd.append(re.sub(r'Command line parameters: ', '', line).strip('\n'))
                _fname.append(os.path.basename(line.split()[-1])) # file name
            elif('Total reads processed:' in line):
                 _raw.append(re.sub(r',', '', line.split()[-1])) # first input
            elif('Reads written (passing filters):' in line):
                 _clean.append(re.sub(r',', '', line.split()[-2])) # output
            else:
                 continue
    logdict['filename'] = _fname[0]
    logdict['run_cutadapt_times'] = len(_raw)
    logdict['command_lines'] = _cmd
    logdict['raw'] = _raw[0]
    logdict['clean'] = _clean[-1]
    logdict['clean_pct'] = '{:.2f}%'.format(int(logdict['clean'])/int(logdict['raw'])*100)
    with open(outfile, 'w') as fo:
        json.dump(logdict, fo, indent = 4)
    return outfile


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



def rep_map_wrapper(path, save = True):
    """
    wrap all bowtie log files, [only for this script] namespace 
    summarize mapping and RTStops 
    header: name, group, read, RTStop
    input: list of json files
    output: pd.DataFrame
    """
    def _json_wrapper(file):
        """
        parsing only one json file
        output: type, count
        """
        try:
            group = file.split('.')[-2]
            group = re.sub('map_', '', group)
            name = os.path.basename(file).split('.')[0]
            with open(file, 'r') as f:
                da = json.load(f)
            return [name, group, da['input_reads'], da['mapped'], 
                    da['unmapped']]
        except IOError:
            print('json file faild: ' + file)
    ## sort files by len
    json_files = sorted(glob.glob(path + '/*.json'), key = len)
    rep_prefix = os.path.basename(path)
    df = pd.DataFrame(columns = ['name', 'group', 'read'])
    #for ff in json_files:
    for n in range(len(json_files)):
        _, g, _, c, _ = _json_wrapper(json_files[n])
        # the first one - genome
        g = 'spikein' if g == 'genome' and n == 0 else g
        df = df.append(pd.DataFrame([[rep_prefix, g, c]], 
            columns = ['name', 'group', 'read']), ignore_index = True)
    n1, g1, _, _, c1 = _json_wrapper(json_files[-1]) # unmap
    df = df.append(pd.DataFrame([[n1, 'unmapped', c1]], 
            columns = ['name', 'group', 'read']), ignore_index = True)
    save_csv = os.path.join(os.path.dirname(path), rep_prefix + '.mapping_stat.csv')
    if save:
        try:
            df.to_csv(save_csv, ',', header = True, index = False)
        except IOError:
            print('[failed] saving data to file: ' + save_csv)
    return df



def merge_map_wrapper(path, save = True):
    """
    count BAM files
    Output: pd.DataFrame
    """
    b_files = sorted(glob.glob(path + '/*.bam'), key = len) # bam files
    b_prefix = os.path.basename(path)
    df = pd.DataFrame(columns = ['name', 'group', 'read'])
    #for b in b_files:
    for n in range(len(b_files)):
        b_cnt = pysam.AlignmentFile(b_files[n], 'rb').count()
        group = os.path.basename(b_files[n]).split('.')[-2] #
        group = re.sub('^map_', '', group)
        group = 'spikein' if n == 0 and group == 'genome' else group
        df = df.append(pd.DataFrame([[b_prefix, group, b_cnt]],
                       columns = ['name', 'group', 'read']),
                       ignore_index = True)
    save_csv = os.path.join(os.path.dirname(path), b_prefix + '.mapping_stat.csv')
    if save:
        try:
            df.to_csv(save_csv, ',', header = True, index = False)
        except IOError:
            print('[failed] saving data to file: ' + save_csv)
    return df


