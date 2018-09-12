#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import os
import io
import re
import json
import tempfile


class Cutadapt_log(object):
    """Wrapper cutadapt log file"""

    def __init__(self, stat):
        self.stat = stat
        # stat
        if isinstance(stat, Cutadapt_log):
            self.stat = stat.stat
        elif isinstance(stat, dict):
            self.stat = stat
        elif isinstance(stat, io.TextIOWrapper):
            self.stat = self._log_parser()
        elif os.path.isfile(stat):
            self.stat = self._log_parser()
        else:
            raise ValueError('not supported file')


    def _log_parser(self):
        """Wrapper log file"""
        dd = {}
        with open(self.stat, 'rt') as ff:
            for line in ff:
                if line.startswith('This is cutadapt'):
                    sep = line.strip().split(' ')
                    dd['version'] = sep[3]
                    dd['python'] = sep[6]
                elif 'Command line parameters' in line:
                    dd['cmd'] = line.strip().split(':')[1]
                elif 'Total reads processed' in line:
                    value = line.strip().split(':')[1]
                    value = re.sub(',', '', value.strip())
                    dd['total'] = int(value)
                elif 'Reads written (passing filters)' in line:
                    value = line.strip().split(':')[1]
                    value = value.strip().split(' ')[0]
                    value = re.sub(',', '', value)
                    dd['clean'] = int(value)
                else:
                    continue
        pct = float(dd['clean']) / float(dd['total']) * 100
        dd['pct'] = '%.1f%%' % pct
        return dd



    def _tmp(self):
        """Create a temp file"""
        tmpfn = tempfile.NamedTemporaryFile(prefix='tmp',
                                            suffix='.json',
                                            delete=False)
        return tmpfn.name



    def saveas(self, _out=None):
        """Make a copy of statistics of mapping results"""
        if _out is None:
            _out = self._tmp()

        dd = self.stat

        with open(_out, 'wt') as fo:
            json.dump(dd, fo, indent=4)

        return _out




class Json_file(object):
    """Parsing Json and dict file"""

    def __init__(self, fn):
        if isinstance(fn, Json_file):
            self.fn = fn.fn
        elif isinstance(fn, dict):
            self.fn = fn
        elif os.path.exists(fn):
            self.fn = json_reader()
        else:
            raise ValueError('unknown file format:')
          

    def _tmp(self):
        """Create a temp file"""
        tmpfn = tempfile.NamedTemporaryFile(prefix='tmp',
                                            suffix='.json',
                                            delete=False)
        return tmpfn.name
        

    def json_reader(self):
        """Load json file as dict"""
        fn = self.fn
        if os.path.isfile(fn) and os.path.getsize(fn) > 0:
            with open(fn, 'rt') as ff:
                return json.load(ff)



    def json_writer(self, to=None):
        """Write dict to file in json format"""
        fn = self.fn

        if to is None:
            to = self._tmp()

        if isinstance(fn, Json_file):
            fn = fn.fn
        elif isinstance(fn, dict):
            fn = fn
            with open(to, 'wt') as ff:
                json.dump(fn, to)















