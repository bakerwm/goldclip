#!/usr/bin/env python
"""
Functions for demx
to-do: barcode demx only
"""

__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__date__ = "2018-04-01"
__version__ = "0.1"


import os
import sys
import json
import gzip
import subprocess
import shlex
import shutil
from collections import defaultdict
from goldclip.helper import *


def str_mismatch(a, b):
    """
    calculate the mismatches between two strings
    query: a
    subject: b
    if b == null, return 0 mismatch
    """
    assert isinstance(a, str)
    assert isinstance(b, str)
    m = sum(list(map(lambda x, y : 0 if x == y else 1, a, b)))
    if b.lower() == 'null':
        m = 1000
    return m


def bc_parser(fn):
    """
    parse the barcode and sample name
    file format:
    <barcode> <name>
    """
    bc_dict = {}
    with open(fn, 'rt') as f:
        for n in f:
            tabs = n.strip().split('\t') # <barcode> <name>
            bc_dup = False if not tabs[0] in bc_dict else True
            bc_dict[tabs[0]] = tabs[1]
    if bc_dup:
        raise ValueError('duplicate [barcodes] detected, exiting...')
    _name = sorted(list(bc_dict.values()))
    _name_uniq = sorted(list(set(_name)))
    if not _name == _name_uniq:
        raise ValueError('duplicate [names] detected, exiting...')
    # barcode in same length
    len_uniq = set(list(map(len, list(bc_dict.keys()))))
    if len(len_uniq) > 1:
        raise ValueError('barcode not in same length')
    return bc_dict


def bc_validater(s, bc_dict, mm=0):
    """
    check whether barcode exists in list
    allow no more than {mm} mismatche(s)
    return name of the barcode
    """
    n = [x for x in list(bc_dict.keys()) if 
          str_mismatch(s, x) <= mm]
    if len(n) == 1:
        return n[0] # bc (from list)
    elif len(n) > 1:
        return 'multi'
    else:
        return None # undemx


##-----------------------------------------##
## demx P7 index   
def p7_split(seq_unit, p7_dict, p7_len, mm=0):
    """
    extract P7 index from comment of fastq
    input: seq_unit : [name, seq, +, qual]
    name: @ST-E00310:586:HJH7JCCXY:7:1101:27275:1555 1:N:0:NCAACAAT
    """
    _name, _seq, _flag, _qual = seq_unit
    s = _name.split(':')[-1]
    p7_query = s[:p7_len]
    m = bc_validater(p7_query, p7_dict, mm=mm)
    return m


##-----------------------------------------##
## demx barcode
def bc_split(seq_unit, bc_dict, bc_length, n_left=3, n_right=2, 
             cut=False, mm=0):
    """
    extract barcode and randomer
    seq_unit : [name, seq, +, qual]
    append the randomer to the name
    """
    _name, _seq, _flag, _qual = seq_unit
    _x = (n_left + bc_length)
    _y = (n_left + bc_length + n_right)
    bc_query = _seq[n_left:_x]
    bc_random = _seq[0:n_left] + _seq[_x:_y]
    # check bc_dict
    m = bc_validater(bc_query, bc_dict, mm=mm) # barcode in list
    _name_list = _name.split(" ")
    _name_list.insert(1, bc_random)
    _name = " ".join(_name_list)
    if cut and m:
        _seq = _seq[_y:]
        _qual = _qual[_y:]
    return [m, [_name, _seq, _flag, _qual]]


##-----------------------------------------##
## demx P7 and barcode at the same time
def get_keys(d, levels=1):
    """
    return all keys of nested dict
    support N levels
    """
    # for levels
    if isinstance(d, dict) and levels == 1:
        return list(d.keys())
    elif levels == 2:
        k1 = list(d.keys())
        kn = []
        for i in k1:
            if isinstance(d[i], dict):
                kn += list(d[i].keys())
            else:
                kn.append(i)
        levels -= 1
        return kn
    elif levels > 2:
        pass
    else:
        pass


def p7_bc_parser(fn):
    """
    parse the barcode and sample name
    file format:
    <p7> <barcode> <name>
    """
    d = defaultdict(dict)
    with open(fn, 'rt') as f:
        for n in f:
            tabs = n.strip().split('\t') # <p7> <barcode> <name>
            if len(tabs) < 3:
                continue
            if tabs[0] in d and tabs[1] in d[tabs[0]]:
                raise ValueError('Duplicate barcodes detected: %s' % n)
            d[tabs[0]][tabs[1]] = tabs[2]
    _name = sorted(list(nested_dict_values(d)))
    _name_uniq = sorted(list(set(_name)))
    if not _name == _name_uniq:
        raise ValueError('Duplicate [names] detected')
    p7_len = set(list(map(len, list(d.keys()))))
    if len(p7_len) > 1:
        raise ValueError('P7 index not in same length')
    bc_list = get_keys(d, 2)
    bc_len = set([len(i) for i in bc_list if not i.lower() == 'null'])
    if len(bc_len) > 1:
        raise ValueError('barcode not in same length')
    return d



def p7_bc_validater(p7, bc, d, mm=0):
    """
    check whether p7 exists, barcode exists
    allow no more than {mm} mismatche(s)
    return the name
    """
    p7_list = list(d.keys()) # for p7 list
    bc_list = get_keys(d, 2) # for barcode list
    a = [x for x in p7_list if str_mismatch(p7, x) <= mm]
    b1 = [x for x in bc_list if str_mismatch(bc, x) <= mm]
    b2 = [x for x in bc_list if str_mismatch(bc, x) > 1000] # check null, 1kb
    b = b2 if len(b2) == 1 else b1
    if len(a) == 1 and len(b) == 1:
        return d[a[0]][b[0]] # name
    elif len(a) > 1 or len(b) > 1:
        return 'multi'
    else:
        return None # undemx



def p7_bc_split(seq_unit, p7_bc_dict, n_left=3, n_right=2, cut=False, mm=0):
    """
    extract P7 index from comment of fastq
    input: seq_unit : [name, seq, +, qual]
    name: @ST-E00310:586:HJH7JCCXY:7:1101:27275:1555 1:N:0:NCAACAAT
    """
    _name, _seq, _flag, _qual = seq_unit
    p7_list = list(p7_bc_dict.keys())
    p7_len = int(sum(list(map(len, p7_list))) / len(p7_list))
    bc_list = get_keys(p7_bc_dict, 2)
    bc_list.remove('undemx') #remove last 2 elements, undemx, multi
    bc_list.remove('multi')
    bc_len = int(sum(list(map(len, bc_list))) / len(bc_list))
    # check p7 index
    s = _name.split(':')[-1]
    p7_query = s[:p7_len]
    # check barcode
    _x = (n_left + bc_len)
    _y = (n_left + bc_len + n_right)
    bc_query = _seq[n_left:_x]
    bc_random = _seq[0:n_left] + _seq[_x:_y]
    # validate barcodes
    n = p7_bc_validater(p7_query, bc_query, p7_bc_dict, mm=mm)
    # output sequence
    _name_list = _name.split(" ")
    _name_list.insert(1, bc_random)
    _name = " ".join(_name_list)
    if cut and n:
        _seq = _seq[_y:]
        _qual = _qual[_y:]
    return [n, [_name, _seq, _flag, _qual]]


##------------------------------------------##
## demx P7 index
def p7_demx_se(fn, p7_file, path_out, cut=False, mm=0):
    """
    demultiplex reads by P7 index
    """
    p7_dict = bc_parser(p7_file)
    p7_len = int(sum(len(x) for x in list(p7_dict.keys())) / len(p7_dict))
    p7_dict['undemx'] = 'undemx' # undemx file
    p7_dict['multi'] = 'multi' # barcode match multi hits
    # writer
    p7_count = {}
    p7_writer = {}
    for idx in p7_dict:
        p7_count[idx] = {}
        p7_count[idx]['count'] = 0
        p7_count[idx]['name'] = p7_dict[idx]
        p7_writer[idx] = gzip.open(os.path.join(path_out, 
                                   p7_dict[idx] + '.fq.gz'), 'wt')
    fq_reader = gzip.open if is_gz(fn) else open
    with fq_reader(fn, 'rt') as fi:
        while True:
            try:
                seq_unit = [next(fi).strip(), 
                            next(fi).strip(), 
                            next(fi).strip(), 
                            next(fi).strip(),]
                tag = p7_split(seq_unit, p7_dict, p7_len, mm)
                tag = tag if tag else 'undemx' # None to 'undemx'
                p7_writer[tag].write('\n'.join(seq_unit) + '\n')
                p7_count[tag]['count'] += 1
            except StopIteration:
                break
    for idx in p7_writer: 
        p7_writer[idx].close() # close writers
    # save report
    report_file = os.path.join(path_out, "report_demx.json")
    with open(report_file, "w") as fo:
        json.dump(p7_count, fo, indent = 4)
    return p7_count


def p7_demx_pe(fn1, fn2, p7_file, path_out, bc_in_read12=1, mm=0):
    """
    demultiplex PE reads, index in NAME
    """
    p7_dict = bc_parser(p7_file)
    p7_len = int(sum(len(x) for x in list(p7_dict.keys())) / len(p7_dict))
    p7_dict['undemx'] = 'undemx' # undemx file
    p7_dict['multi'] = 'multi' # barcode match multi hits
    r1_suffix = '_1.fq.gz'
    r2_suffix = '_2.fq.gz'
    if bc_in_read12 == 2:
        (fn1, fn2) = (fn2, fn1) # switch read 1/2
        (r1_suffix, r2_suffix) = (r2_suffix, r1_suffix) # switch read 1/2

    ###################
    ## create writer ##
    ###################
    p7_count = {}
    p7_writer = {}
    for idx in p7_dict:
        p7_count[idx] = {}
        p7_count[idx]['count'] = 0
        p7_count[idx]['name'] = p7_dict[idx]
        p7_writer[idx] = [gzip.open(os.path.join(path_out, 
                                    p7_dict[idx] + r1_suffix), 'wt'),
                          gzip.open(os.path.join(path_out, 
                                    p7_dict[idx] + r2_suffix), 'wt'),]
    fq_reader1 = gzip.open if is_gz(fn1) else open
    fq_reader2 = gzip.open if is_gz(fn2) else open
    with fq_reader1(fn1, 'rt') as f1, fq_reader2(fn2, 'rt') as f2:
        while True:
            try:
                seq_unit1 = [next(f1).strip(), 
                             next(f1).strip(), 
                             next(f1).strip(), 
                             next(f1).strip(),]
                seq_unit2 = [next(f2).strip(), 
                             next(f2).strip(), 
                             next(f2).strip(), 
                             next(f2).strip(),]
                if not seq_unit1[0].split(" ")[0] == seq_unit2[0].split(" ")[0]:
                    raise NameError(seq_unit1[0] + "\n" + seq_unit2[0])
                tag = p7_split(seq_unit1, p7_dict, p7_len, mm)
                tag = tag if tag else 'undemx' # None to 'undemx'
                p7_writer[tag][0].write('\n'.join(seq_unit1) + '\n')
                p7_writer[tag][1].write('\n'.join(seq_unit2) + '\n')
                p7_count[tag]['count'] += 1
            except StopIteration:
                break
    for idx in p7_writer: 
        p7_writer[idx][0].close() # close writers
        p7_writer[idx][1].close() # close writers
    # save report
    report_file = os.path.join(path_out, "report_demx.json")
    with open(report_file, "w") as fo:
        json.dump(p7_count, fo, indent=4)
    return p7_count



##-----------------------------------------##
## demx barcode
def bc_demx_se(fn, bc_file, path_out, n_left=3, n_right=2, 
               cut=True, mm=0):
    """
    demultiplex SE reads, 
    inline-barcode at the beginning of the read
    """
    bc_dict = bc_parser(bc_file)
    bc_len = int(sum(len(x) for x in list(bc_dict.keys())) / len(bc_dict))
    bc_dict['undemx'] = 'undemx' # add undemx
    bc_dict['multi'] = 'multi' # barcode match multi hits
    # writer
    bc_count = {}
    bc_writer = {}
    for bc in bc_dict:
        bc_count[bc] = {}
        bc_count[bc]['count'] = 0
        bc_count[bc]['name'] = bc_dict[bc]
        bc_writer[bc] = gzip.open(os.path.join(path_out, 
                                               bc_dict[bc] + '.fq.gz'), 'wt')
    fq_reader = gzip.open if is_gz(fn) else open
    with fq_reader(fn, 'rt') as fi:
        while True:
            try:
                seq_unit = [next(fi).strip(), 
                            next(fi).strip(), 
                            next(fi).strip(), 
                            next(fi).strip(),]
                tag, seq_unit_new = bc_split(seq_unit, 
                                             bc_dict, 
                                             bc_len, 
                                             n_left=n_left, 
                                             n_right=n_right, 
                                             cut=cut, 
                                             mm=mm)
                tag = tag if tag else 'undemx' # None
                bc_writer[tag].write('\n'.join(seq_unit_new) + '\n')
                bc_count[tag]['count'] += 1
            except StopIteration:
                break
    for bc in bc_writer:
        bc_writer[bc].close() # close writers
    # save report
    report_file = os.path.join(path_out, "report_demx.json")
    with open(report_file, "w") as fo:
        json.dump(bc_count, fo, indent = 4)
    return bc_count


def bc_demx_pe(fn1, fn2, bc_file, path_out, n_left=3, n_right=2, 
               bc_in_read12=1, cut=True, mm=0):
    """
    demultiplex PE reads
    barcode located in the 5-prime end of fn1
    the second file: fn2 does not contain barcode
    """
    bc_dict = bc_parser(bc_file)
    bc_list = list(bc_dict.keys())
    bc_len = int(sum(map(len, bc_list)) / len(bc_dict))
    bc_dict['undemx'] = 'undemx' # add undemx
    bc_dict['multi'] = 'multi' # barcode match multi hits
    r1_suffix = '_1.fq.gz'
    r2_suffix = '_2.fq.gz'
    if bc_in_read12 == 2:
        (fn1, fn2) = (fn2, fn1) # switch read 1/2
        (r1_suffix, r2_suffix) = (r2_suffix, r1_suffix) # switch read 1/2

    ###################
    ## create writer ##
    ###################
    bc_count = {}
    bc_writer = {}
    for bc in bc_dict:
        bc_count[bc] = {}
        bc_count[bc]['count'] = 0
        bc_count[bc]['name'] = bc_dict[bc]
        bc_writer[bc] = [gzip.open(os.path.join(path_out, 
                                   bc_dict[bc] + r1_suffix), 'wt'),
                         gzip.open(os.path.join(path_out, 
                                   bc_dict[bc] + r2_suffix), 'wt'),]
    ###################
    ## iterate reads ##
    ###################
    fq_reader1 = gzip.open if is_gz(fn1) else open
    fq_reader2 = gzip.open if is_gz(fn2) else open
    with fq_reader1(fn1, 'rt') as f1, fq_reader2(fn2, 'rt') as f2:
        while True:
            try:
                seq_unit1 = [next(f1).strip(), 
                             next(f1).strip(), 
                             next(f1).strip(), 
                             next(f1).strip(),]
                seq_unit2 = [next(f2).strip(), 
                             next(f2).strip(), 
                             next(f2).strip(), 
                             next(f2).strip(),]
                if not seq_unit1[0].split(" ")[0] == seq_unit2[0].split(" ")[0]:
                    raise NameError(seq_unit1[0] + "\n" + seq_unit2[0])
                tag, seq_unit1 = bc_split(seq_unit1, bc_dict, bc_len, 
                                          n_left=n_left, n_right=n_right,
                                          cut=cut, mm=mm)
                tag = tag if tag else 'undemx' # None
                bc_writer[tag][0].write('\n'.join(seq_unit1) + '\n')
                bc_writer[tag][1].write('\n'.join(seq_unit2) + '\n')
                bc_count[tag]['count'] += 1
            except StopIteration:
                break
    for bc in bc_writer: 
        bc_writer[bc][0].close()
        bc_writer[bc][1].close()
    # save report
    report_file = os.path.join(path_out, "report_demx.json")
    with open(report_file, "w") as fo:
        json.dump(bc_count, fo, indent = 4)
    return bc_count


##-----------------------------------------##
## demx P7 and barcode at the same time
def p7_bc_demx_se(fn, p7_bc_file, path_out, n_left=3, n_right=2, 
                  cut=False, mm=0):
    """
    demultiplex reads by P7 index and barcode
    """
    p7_bc_dict = p7_bc_parser(p7_bc_file)
    # bc_list = get_keys(p7_bc_dict, 2)
    # bc_len = int(sum(list(map(len, bc_list))) / len(bc_list))
    p7_bc_dict['undemx'] = 'undemx' # undemx file
    p7_bc_dict['multi'] = 'multi' # barcode match multi hits
    # writer
    p7_bc_count = {}
    p7_bc_writer = {}
    for p7 in p7_bc_dict: #
        # assert isinstance(p7_bc_dict[p7], dict)
        if isinstance(p7_bc_dict[p7], dict):
            for bc in p7_bc_dict[p7]: #
                f_name = p7_bc_dict[p7][bc]
                f_key = p7 + bc
                p7_bc_count[f_name] = {}
                p7_bc_count[f_name]['count'] = 0
                p7_bc_count[f_name]['name'] = f_name
                p7_bc_writer[f_name] = gzip.open(os.path.join(path_out, 
                                                 f_name + '.fq.gz'), 'wt')
        elif isinstance(p7_bc_dict[p7], str):
            f_name = p7_bc_dict[p7]
            f_key = p7
            p7_bc_count[f_name] = {}
            p7_bc_count[f_name]['count'] = 0
            p7_bc_count[f_name]['name'] = f_name
            p7_bc_writer[f_name] = gzip.open(os.path.join(path_out, 
                                             f_name + '.fq.gz'), 'wt')
        else:
            continue
    fq_reader = gzip.open if is_gz(fn) else open
    with fq_reader(fn, 'rt') as fi:
        while True:
            try:
                seq_unit = [next(fi).strip(), 
                            next(fi).strip(), 
                            next(fi).strip(), 
                            next(fi).strip(),]
                tag, seq_unit_new = p7_bc_split(seq_unit, 
                                                p7_bc_dict, 
                                                n_left=n_left,
                                                n_right=n_right,
                                                cut=cut,
                                                mm=mm) 
                tag = tag if tag else 'undemx' # None to 'undemx'
                p7_bc_writer[tag].write('\n'.join(seq_unit_new) + '\n')
                p7_bc_count[tag]['count'] += 1
            except StopIteration:
                break
    for idx in p7_bc_writer: 
        p7_bc_writer[idx].close() # close writers
    # save report
    report_file = os.path.join(path_out, "report_demx.json")
    with open(report_file, "w") as fo:
        json.dump(p7_bc_count, fo, indent = 4)
    return p7_bc_count



def p7_bc_demx_pe(fn1, fn2, p7_bc_file, path_out, n_left=3, n_right=2, 
                  bc_in_read12=1, cut=False, mm=0):
    """
    demultiplex PE reads, by p7 index and barcode
    p7 in comment of fastq
    barcode in fn1
    """
    p7_bc_dict = p7_bc_parser(p7_bc_file)
    p7_bc_dict['undemx'] = 'undemx' # undemx file
    p7_bc_dict['multi'] = 'multi' # barcode match multi hits
    r1_suffix = '_1.fq.gz'
    r2_suffix = '_2.fq.gz'
    if bc_in_read12 == 2:
        (fn1, fn2) = (fn2, fn1) # switch read 1/2
        (r1_suffix, r2_suffix) = (r2_suffix, r1_suffix) # switch read 1/2

    ###################
    ## create writer ##
    ###################
    p7_bc_count = {}
    p7_bc_writer = {}
    for p7 in p7_bc_dict: #
        if isinstance(p7_bc_dict[p7], dict):
            # one p7 -> multiple barcode
            for bc in p7_bc_dict[p7]: #
                f_name = p7_bc_dict[p7][bc]
                f_key = p7 + bc
                p7_bc_count[f_name] = {}
                p7_bc_count[f_name]['count'] = 0
                p7_bc_count[f_name]['name'] = f_name
                p7_bc_writer[f_name] = [gzip.open(os.path.join(path_out, 
                                        f_name + r1_suffix), 'wt'),
                                        gzip.open(os.path.join(path_out, 
                                        f_name + r2_suffix), 'wt'),]
        elif isinstance(p7_bc_dict[p7], str):
            # one p7 -> one barcode
            f_name = p7_bc_dict[p7]
            f_key = p7
            p7_bc_count[f_name] = {}
            p7_bc_count[f_name]['count'] = 0
            p7_bc_count[f_name]['name'] = f_name
            p7_bc_writer[f_name] = [gzip.open(os.path.join(path_out, 
                                    f_name + r1_suffix), 'wt'),
                                    gzip.open(os.path.join(path_out, 
                                    f_name + r2_suffix), 'wt'),]
        else:
            continue
    fq_reader1 = gzip.open if is_gz(fn1) else open # contain barcode in first file
    fq_reader2 = gzip.open if is_gz(fn2) else open
    with fq_reader1(fn1, 'rt') as f1, fq_reader2(fn2, 'rt') as f2:
        while True:
            try:
                seq_unit1 = [next(f1).strip(), 
                             next(f1).strip(), 
                             next(f1).strip(), 
                             next(f1).strip(),]
                seq_unit2 = [next(f2).strip(), 
                             next(f2).strip(), 
                             next(f2).strip(), 
                             next(f2).strip(),]
                if not seq_unit1[0].split(" ")[0] == seq_unit2[0].split(" ")[0]:
                    raise NameError(seq_unit1[0] + "\n" + seq_unit2[0])
                tag, seq_unit1_new = p7_bc_split(seq_unit1, 
                                                 p7_bc_dict,
                                                 n_left=n_left,
                                                 n_right=n_right,
                                                 cut=cut,
                                                 mm=mm) 
                tag = tag if tag else 'undemx' # None to 'undemx'
                p7_bc_writer[tag][0].write('\n'.join(seq_unit1_new) + '\n')
                p7_bc_writer[tag][1].write('\n'.join(seq_unit2) + '\n')
                p7_bc_count[tag]['count'] += 1
            except StopIteration:
                break
    for idx in p7_bc_writer: 
        p7_bc_writer[idx][0].close() # close writers
        p7_bc_writer[idx][1].close() # close writers
    # save report
    report_file = os.path.join(path_out, "report_demx.json")
    with open(report_file, "w") as fo:
        json.dump(p7_bc_count, fo, indent = 4)
    return p7_bc_count


##------------------------------------------##
## demx by bioawk
def demx_se_bioawk(fn, bc_file, outdir, n_left, n_right):
    """
    demultiplexing SE reads using bioawk
    ~ 7x faster than demx_se()
    """
    if not shutil.which("bioawk"):
        sys.exist("bioawk not found in your PATH")
    if not os.path.exists(outdir) is True:
        os.makedirs(outdir)
    fn = os.path.abspath(fn)
    cwd = os.getcwd()
    tmpdir = os.path.join(outdir, "undemx")
    if not os.path.exists(tmpdir) is True:
        os.makedirs(tmpdir)
    # # counter
    bc_dict = bc_parser(bc_file)
    bc_list = list(bc_dict.keys())
    bc_tmp = list(set([len(i) for i in bc_list]))
    bc_length = bc_tmp[0] # choose only one
    bc_dict["undemx"] = "undemx" # add undemx
    # if len(bc_tmp) > 1:
    #     sys.exit("--barcode, barcodes are not in same length")
    # # run
    # bioawk -c fastx -v lft=3 -v rgt=2 -v n=4 '{sum+=1; s=lft+1; e=lft+n+1; d=lft+n+rgt+1; g=len($seq); \
    # bc=substr($seq,s,n); r1=substr($seq,1,lft); r2=substr($seq,e,rgt); \
    # id="@"$name" "r1""r2" "$comment; sq=substr($seq,d,g); qa=substr($qual,d,g); \
    # bcfile=bc".fq"; print id"\n"sq"\n+\n"qa > bcfile}END{print sum > report_count.txt}' in.fq
    # the command line
    c1 = "bioawk -c fastx -v lft={} -v rgt={} -v n={} ".format(n_left, n_right, bc_length)
    c2 = "'{sum+=1; s=lft+1; e=lft+n+1; d=lft+n+rgt+1; g=length($seq); bc=substr($seq,s,n); r1=substr($seq,1,lft); r2=substr($seq,e,rgt); "
    c3 = "id=\"@\"$name\" \"r1\"\"r2\" \"$comment; sq=substr($seq,d,g); qa=substr($qual,d,g); "
    c4 = "bcfile=bc\".fq\"; print id\"\\n\"sq\"\\n+\\n\"qa > bcfile}END{print sum > \"report_count.txt\"}' "
    cmd1 = shlex.split(c1 + c2 + c3 + c4 + fn)
    os.chdir(tmpdir)
    p1 = subprocess.Popen(cmd1)
    p1.communicate()
    os.chdir(cwd)
    # rname files, count reads
    bc_count = {}
    _tmp_demx = 0
    for bc in bc_dict:
        bc_count[bc] = {}
        bc_count[bc]["count"] = 0
        bc_count[bc]["name"] = bc_dict[bc]
        bc_from = os.path.join(outdir, "undemx", bc + ".fq")
        bc_to = os.path.join(outdir, bc_dict[bc] + ".fq")
        if os.path.exists(bc_from):
            _tmp = file_row_counter(bc_from) / 4 # fq
            _tmp_demx += _tmp
            bc_count[bc]["count"] = _tmp
            shutil.move(bc_from, bc_to)
            cmd2 = shlex.split("gzip -f {}".format(bc_to))
            subprocess.run(cmd2)
            # with open(bc_to, 'rb') as read, gzip.open(bc_to + ".gz", 'wb') as write:
            #     shutil.copyfileobj(read, write)
    # total
    with open(os.path.join(outdir, "undemx", "report_count.txt"), "r") as f:
        _tmp_total = f.read().strip()
        _tmp_total = int(_tmp_total)
    # undemx
    bc_count["undemx"] = {}
    bc_count["undemx"]["count"] = _tmp_total - _tmp_demx
    bc_count["undemx"]["name"] = "undemx"
    # save report
    report_file = os.path.join(outdir, "report_demx.json")
    with open(report_file, "w") as fo:
        json.dump(bc_count, fo, indent = 4)
    return bc_count

## EOF
