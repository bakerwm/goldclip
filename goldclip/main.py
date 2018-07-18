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
from goldclip.goldcliplib.goldclip import *


parser = argparse.ArgumentParser(prog='goldclip', description='Sub-commands of goldclip.')
subparsers = parser.add_subparsers(help='sub-commands.', dest='subcommand')
subparsers.required = True

## options for sub-programs
subparser_demx = subparsers.add_parser('run', help='running goldclip from fastq to peaks.')
subparser_demx = add_run_args(subparser_demx)
subparser_demx = subparsers.add_parser('demx', help='demultiplexing illumina reads.')
subparser_demx = add_demx_args(subparser_demx)
subparser_trim = subparsers.add_parser('trim', help='trimming reads.')
subparser_trim = add_trim_args(subparser_trim)
subparser_map = subparsers.add_parser('map', help='mapping reads to reference.')
subparser_map = add_map_args(subparser_map)
subparser_peak = subparsers.add_parser('peak', help='calling peaks.')
subparser_peak = add_peak_args(subparser_peak)
subparser_rtstop = subparsers.add_parser('rtstop', help='calling RT-Stops.')
subparser_rtstop = add_rtstop_args(subparser_rtstop)
subparser_report = subparsers.add_parser('report', help='making a report for the project.')
subparser_report = add_report_args(subparser_report)


class goldclip:

    def __init__(self, *args, **kwargs):
        self.kwargs = kwargs

    def start_process(self):
        if 'demx' in self.kwargs['subcommand'].lower():
            d = Demx(**self.kwargs).run()
        elif 'trim' in self.kwargs['subcommand'].lower():
            t = Trim(**self.kwargs).run()
        elif 'map' in self.kwargs['subcommand'].lower():
            m = Map(**self.kwargs).run()
        elif 'peak' in self.kwargs['subcommand'].lower():
            p = Peak(**self.kwargs).run()
        elif 'rtstop'in self.kwargs['subcommand'].lower():
            r = Rtstop(**self.kwargs).run()
        elif 'report' in self.kwargs['subcommand'].lower():
            r = Report(**self.kwargs).run()
        elif 'run' in self.kwargs['subcommand'].lower():
            a = Run_all(**self.kwargs).run()
        else:
            raise ValueError('Unknown subcommand: %s' % self.kwargs['subcommand'])


def main(flags=None):
    if flags == None:
        flags = sys.argv[1:]
    flags = vars(parser.parse_args(flags))
    goldclip(**flags).start_process()


if __name__ == '__main__':
    main()

## EOF
