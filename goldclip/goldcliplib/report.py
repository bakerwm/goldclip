#!/usr/bin/env python
"""
make a summary of the project

wrapper output of goldclip pipeline 
create stats and plots

## figure1.reads_mapping_stat.pdf
number of reads: raw, clean, no_dup, spikein, genome, unmap, ...

## figure2.reads_annotation_stat.pdf
number of reads: RNA categories, ...

## figure3.reads_correlation.pdf

## figure4.rtstop_correlation.pdf

## figure5.peak_number.pdf

## figure6.peak_length.pdf

## figure7.peak_annotation.pdf

## figure8.peak_motif.pdf

## figure9.peak_conservation.pdf

## figure10.hexmer_zscore.pdf

## figure11.peak_overlap.pdf

functions
"""

__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__date__ = "2018-03-21"
__version__ = "0.1"



import os
import sys
import re
import argparse
import json
import tempfile
import string
import random
import pybedtools
import pandas as pd
import numpy as np
import multiprocessing as mp
from goldclip.helper import *


from .helper import *
from .bed_annotation import bed_annotator
from .bed2motif_homer import bed2motif
from .settings import goldclip_home



if path_data is None:
        path_data = os.path.join(pathlib.Path.home(), 'data', 'genome')