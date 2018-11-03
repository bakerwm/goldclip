#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Call peaks using multiple tools:
1. CLIPper
2. Pyicoclip
3. Pyiraha


# cmd1
clipper -b <bam> -s <genome> -o <out.bed> <para> 

# cmd2
pyicoclip --stranded -p-value 0.001 --region <region.bed> -f bed <in.bed> <out.bed> 

# cmd3
pyraha

# cmd4
CLAM

"""

import os
import sys
import pathlib
import logging
import shutil
import shlex
import subprocess
import pybedtools
import multiprocessing as mp
from helper import *



class Peaktools(object):
    """Call peaks using CLIPper (python 2.7)
    input: <BAM>, <genome> 
    output: <peak.bed>, <fixed.bed>
    """

    def __init__(self, bam, genome, path_out, overwrite=False):
        self.bam = bam
        self.genome = genome
        self.path_out = path_out
        self.overwrite = overwrite



    # current, run
    def _get_base(self):
        """Get base prefix of current path"""
        if hasattr(sys, 'real_prefix'):
            return sys.real_prefix
        elif hasattr(sys, 'base_prefix'):
            return sys.base_prefix
        else:
            return None



    def _within_venv(self):
        """
        determine if inside virtualenv
        """
        return (hasattr(sys, 'real_prefix') or 
            (hasattr(sys, 'base_prefix') and 
            sys.base_prefix !=  sys.prefix))



    def _to_venv(self, venv, out=False):
        """Enter a venv"""
        venv_bin = os.path.join(venv, 'bin', 'activate_this.py')
        assert os.path.exists(venv_bin)
        if sys.version_info[0:2] ==  (2, 7):
            execfile(venv_bin, dict(__file__ = venv_bin)) # python2
            assert self._within_venv()
        elif sys.version_info[0:1] >=  (3, 5):
            exec(open(venv_bin).read(), {}, dict(__file__ = venv_bin)) #python3
            assert self._within_venv()
        else:
            logging.error('unknown version of python: %s' % sys.version)



    def _check_venv(self, venv, enter_venv=True):
        """Check venv, if not , switch to it"""
        venv = os.path.expanduser(venv)

        if self._within_venv() and venv == self._get_base():
            logging.info('already in venv: %s' % venv)
            return True
        elif os.path.exists(venv):
            if enter_venv is True:
                logging.info('switch to venv: %s' % venv)
                self._to_venv(venv)
            return True
        else:
            logging.error('virtualenv not exists: %s' % venv)
            return False



    def adjust_score(self, bed):
        """
        fix clipper output 
        template: 
        https://github.com/YeoLab/gscripts/blob/master/gscripts/clipseq/fix_scores.py
        """
        peak_center = str((int(bed[6]) + int(bed[7])) / 2)
        qValue = bed.score
        pValue = str(-1) #to fix
        signalValue = str(-1)
        if float(bed.score) != 0.0:
            bed.score = str(min(int(-10 * np.log10(float(bed.score))), 1000))
        else:
            bed.score = '0'
        bed[6] = signalValue
        bed[7] = pValue
        bed.append(qValue)
        bed.append(peak_center)
        return bed



    def run_clipper(self):
        """Run CLIPper"""
        clipper_env = os.path.join(str(pathlib.Path.home()), 'envs', 'clipper_env')
        self._check_venv(clipper_env)
        flag = which('clipper')
        if flag is None:
            raise ValueError('clipper command not found')

        bam = self.bam
        genome = self.genome
        path_out = self.path_out
        overwrite = self.overwrite
        assert os.path.exists(bam)
        assert is_path(path_out)

        prefix = file_prefix(bam)[0]
        peak_name = os.path.join(path_out, prefix)
        peak_bed = peak_name + '.bed'
        peak_fixed = peak_name + '.fixed.bed'
        peak_log = peak_name + '.clipper.log'
        clipper_para = '--bonferroni --superlocal --threshold-method binomial \
                        --save-pickle'
        c1 = 'clipper %s -b %s -s %s -o %s' % (clipper_para, bam, genome, peak_bed)
        # overwrite
        if os.path.exists(peak_fixed) and overwrite is False:
            logging.info('peak file exists: %s' % peak_fixed)
        else:
            with open(peak_log, 'wt') as fo:
                p1 = subprocess.run(shlex.split(c1), stdout=fo)
            pybedtools.BedTool(peak_bed).each(self.adjust_score).saveas(peak_fixed)
        return peak_fixed



    def _tmp(self):
        """Create a temp file"""
        tmpfn = tempfile.NamedTemporaryFile(prefix='tmp',
                                            suffix='.bed',
                                            delete=False)
        return tmpfn.name



    def _get_homer_anno(self, anno_dir=None):
        """Get annotation files in BED format"""
        genome = self.genome
        if anno_dir is None:
            home = str(pathlib.Path.home())
        anno_dir = os.path.join(home, 'data', 'genome', genome, 'annotation', 
            'homer')
        groups = ['tts', 'rRNA', 'pseudo', 'promoters', 'ncRNA',
            'utr3', 'utr5', 'coding', 'introns', 'intergenic']
        anno_files = [os.path.join(anno_dir, a + '.ann.bed') for a in groups]
        anno_files = [f for f in anno_files if os.path.exists(f)]
        if len(anno_files) < 8:
            logging.error('Annotation files not found: %s' % genome)
            raise ValueError('Annotation files not found')
        return anno_files



    def bed2peak_pyicoclip(self, bed, bed_anno, path_out):
        """Run single annotation for bed
        return before_*bedpk
        """
        pyicoclip_env = os.path.join(pathlib.Path.home(), 'envs', 'pyicoclip_env')
        self._check_venv(pyicoclip_env)
        flag = which('pyicoclip')
        if flag is None:
            raise ValueError('pyicoclip  not found')

        assert os.path.exists(bed)
        assert os.path.exists(bed_anno)
        assert is_path(path_out)
        prefix = file_prefix(bed)[0]
        anno_prefix = os.path.basename(bed_anno).split('.')[0]
        peak_name = os.path.join(path_out, prefix + '.' + anno_prefix)
        peak_bed = peak_name + '.bed'
        peak_log = peak_name + '.log'
        print(peak_bed)
        with open(peak_log, 'wt') as fo:
            c1 = 'pyicoclip --stranded --p-value 0.001 --region %s \
                  -f bed %s %s' % (bed, bed_anno, peak_bed)
            p1 = subprocess.run(shlex.split(c1), stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
        os.remove('pyicoclip.log') # log
        return peak_bed



    def run_pyicoclip(self):
        """Run Pyicoclip
        input: bed, anno
        output: bed
        """
        bam = self.bam
        genome = self.genome
        path_out = self.path_out
        overwrite = self.overwrite
        anno_files = self._get_homer_anno()

        prefix = file_prefix(bam)[0]
        peak_name = os.path.join(path_out, prefix)
        peak_unique = peak_name + '.unique.bed'
        peak_rmdup_log = peak_name + '.rmdup.log'
        peak_fixed = peak_name + '.fixed.bed'

        if os.path.exists(peak_fixed) and overwrite is False:
            logging.info('peak exists: %s' % peak_fixed)
            return True

        # convert BAM to bed ?
        if bam.lower().endswith('.bed'):
            bed_in = bam
        elif bam.lower().endswith('.bam'):
            bed_in = os.path.split(bam)[0] + '.bed'
            if os.path.exists(bed_in) and overwrite is False:
                pass
            else:
                pybedtools.BedTool(bam).bam_to_bed().saveas(bed_in)

        # call peaks for each region
        path_tmp = os.path.join(path_out, 'temp_peaks')
        assert is_path(path_tmp)
       
        # run in parallel
        pool = mp.Pool(processes=10)
        peak_sub_files = [pool.apply(self.bed2peak_pyicoclip, 
                          args=(bed_in, anno_file, path_tmp)) 
                          for anno_file in anno_files]
        print(peak_sub_files)
        peak_dup = self._tmp()
        with open(peak_dup, 'wb') as wfd:
            for f in peak_sub_files:
                with open(f, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd, 1024*1024*10)
                    # 10MB per writing chunk to aovid reading large file into memory

        # remove duplicates
        c1 = 'pyicos remduplicates %s %s' % (peak_dup, peak_unique)
        with open(peak_rmdup_log, 'wt') as fo:
            subprocess.run(shlex.split(c1), stderr=fo, stdout=fo)

        # fix bed
        Bed_parser(peak_unique).bed_fixer().saveas(peak_fixed)
        return peak_fixed






