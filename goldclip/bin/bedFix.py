#!/usr/bin/env python
# fix BED format
#
# Date: 2018-04-16 switch to python code
# Date: 2017-08-04

import os, sys, re
import pandas as pd
import numpy as np

# start >= 0
# end >= 1
# start < end
# len(id) < 40
# score is int
# strand + - *
# extra add

# def bed_fix(bed):
#     """fix bed format"""
#     bed[1] = str(abs(int(bed[1])))
#     bed[1] = str(bed[1]) if int(bed[1]) > 0 else 0
#     bed[2] = str(abs(int(bed[2])))
#     bed[2] = str(bed[2]) if int(bed[2]) > 1 else 1
#     bed[1], bed[2] = (bed[2], bed[1]) if bed[1] > bed[2] else bed[1:3]
#     bed[3] = bed[3][0:40] # length at most
#     bed[5] = bed[5] if re.fullmatch("^[+-]$", bed[5]) else "*"
#     return bed


def coord_fixer(n):
    """
    fix chromStart, chromEnd
    """
    m = re.sub(',|\s+', '', n) # 1,103,274 -> 1103274
    m = 0 if int(m) < 0 else int(m)
    return int(m) # convert minus to plus


def coord_valid(s, e):
    """
    should be m < n
    """
    s = coord_fixer(s) # 0-index
    e = coord_fixer(e) # 1-index
    e = 1 if e == 0 else e # 
    s, e = [s, e] if (int(s) < int(e)) else [e, s]
    return [str(s), str(e)]


def bed_fix(bed):
    """
    fix bed format errors, simple bed3, bed6 format
    1. chromStart < chromEnd, int
    2. name : length <= 40
    3. score : int, 0-1000
    4. strand : +, -, *
    input: list
    """
    if isinstance(bed, list) and len(bed) >= 6:
        chrom, chromStart, chromEnd, name, score, strand = bed[:6]
        chromStart, chromEnd = coord_valid(chromStart, chromEnd)
        if len(name) > 40:
            name = name[0:40] # maximum width: 40
        if re.fullmatch('-?\d+(\.0+)?', score):
            pass
        else:
            bed.append(score)
            score = '100'
        strand = strand if re.fullmatch('[+-]', strand) else '*'
        return [chrom, chromStart, chromEnd, name, score, strand] + bed[6:]
    elif isinstance(bed, list) and len(bed) >= 3:
        chrom, chromStart, chromEnd = bed[:3]
        chromStart, chromEnd = coord_valid(chromStart, chromEnd)
        return [chrom, chromStart, chromEnd] + bed[3:]
    else:
        return None

## parameters
if len(sys.argv) < 2:
    sys.exit("Usage: bedFix.py <in.bed>")
bed_in = sys.argv[1]

with open(bed_in, "r") as f:
    for line in f:
        line = line.rstrip()
        b = bed_fix(line.split('\t'))
        if b:
            print('\t'.join(b))

## version-2017-01-18
##
## 1. fix the coordinates of BED file
## 2. chop the id_field of BED (limit: 40 characters)
## save original file to *.bak
#
#[[ $# -lt 1 ]] && echo "Usage: $0 <in.bed>" && exit 0
#in=$1
#bak="${in}.bak"
#cp -rf ${in} ${bak}
#cat ${bak} | awk '{FS=OFS="\t"}{$4=substr($4,0,20); if($2>$3){t=$2;$2=$3;$3=t} print $0}' > ${in}


## EOF
