#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""
The main purpose of this package is to process GoldCLIP datasets, 
including, demultiplexing, trimming, mapping, peak-calling, 
and also further reporting.

This package was designed to be easy maintained and human-readable.

"""

__author__ = 'Ming Wang <wangm08@hotmail.com>'
__copyright__ = '2018 by Ming Wang <wangm08@hotmail.com>'
__license__ = 'MIT'
__email__ = 'wangm08@hotmail.com'
__version__ = '0.0.1'


import sys
import argparse
from goldclip.goldcliplib.parser import *
from goldclip.goldcliplib.goldclip import Trim, Align, Peak, Rtstop, Report, Goldclip_all_in_one


parser = argparse.ArgumentParser(prog='goldclip', description='Goldclip commands.')
subparsers = parser.add_subparsers(help='sub-commands.', dest='subcommand')
subparsers.required = True

## options for sub-programs
subparser_demx = subparsers.add_parser('all-in-one', help='running goldclip pipeline from fastq to peaks.')
subparser_demx = add_all_in_one_args(subparser_demx)
subparser_trim = subparsers.add_parser('trim', help='trimming reads.')
subparser_trim = add_trim_args(subparser_trim)
subparser_align = subparsers.add_parser('align', help='alignment reads to reference genome.')
subparser_align = add_align_args(subparser_align)
subparser_peak = subparsers.add_parser('call-peak', help='calling peaks.')
subparser_peak = add_peak_args(subparser_peak)
subparser_rtstop = subparsers.add_parser('call-rtstop', help='calling RT-Stops.')
subparser_rtstop = add_rtstop_args(subparser_rtstop)
subparser_report = subparsers.add_parser('report', help='making a report for goldclip project.')
subparser_report = add_report_args(subparser_report)


def main(flag=None):
    # define subprocess, classes
    sub_process = {
        'all-in-one': Goldclip_all_in_one,
        'trim': Trim,
        'align': Align,
        'call-peak': Peak,
        'call-rtstop': Rtstop,
        'report': Report
    }
    if flag is None:
        flag = sys.argv[1:]
    args = vars(parser.parse_args(flag))
    #
    subcommand = args['subcommand'].lower()
    if subcommand in sub_process:
        tmp = sub_process[subcommand](**args).run()
    else:
        raise Exception('unknown subcommand: %s' % subcommand)

if __name__ == '__main__':
    main()

## EOF
