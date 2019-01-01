#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Mapping fastq to reference genome
1. rRNA, spikein, optional
2. genome
"""

import os
import sys
import re
import io
import glob
import json
import fnmatch
import tempfile
import shlex
import subprocess
import logging
import pandas as pd
import pysam
import pybedtools
from operator import is_not
from functools import partial

from goldclip.helper import *
from goldclip.goldcliplib.log_parser import *


class Alignment(object):
    """Run alignment for SE reads using bowtie/bowtie2/STAR"""

    def __init__(self, fqs, path_out, smp_name, genome, **kwargs):
        """parse arguments
        required arguments: fqs, path_out, smp_name, genome, genome, spikein,
        index_ext, threads, unique_only, n_map, aligner, align_to_rRNA,
        genome_path, overwrite
        """
        self.fqs = fqs
        self.path_out = path_out
        self.smp_name = smp_name
        self.genome = genome
        self.kwargs = kwargs
        self.args = self._args_init()


    def _args_init(self):
        """Inititate the arguments, assign the default values to arg
        """
        args = self.kwargs
        args['fqs'] = self.fqs
        args['path_out'] = self.path_out
        args['smp_name'] = self.smp_name
        args['genome'] = self.genome
        args['spikein'] = args.get('spikein', None)
        args['index_ext'] = args.get('index_ext', None)
        args['threads'] = args.get('threads', 1)
        args['unique_only'] = args.get('unique_only', False)
        args['n_map'] = args.get('n_map', 0)
        args['aligner'] = args.get('aligner', 'bowtie')
        args['align_to_rRNA'] = args.get('align_to_rRNA', True)
        args['repeat_masked_genome'] = args.get('repeat_masked_genome', False)
        args['merge_rep'] = args.get('merge_rep', True)
        args['genome_path'] = args.get('genome_path', None)
        args['overwrite'] = args.get('overwrite', False)
        # check
        if args['spikein'] == self.genome:
            args['spikein'] = None #
        return args


    def _path_init(self, fq, index, reference=None, align_path=None):
        """Create folders for the alignment, 
        Alignment, genome versions
        1.genome_rRNA
        2.genome
        3.spikein_rRNA
        4.spikein

        return files:
        prefix, bam, bed, log, unmap
        """
        args = self.args

        if not reference:
            reference = args['genome'] # default is reference genome
        if not align_path:
            align_path = args['path_out']

        fq_prefix = file_prefix(fq)[0]
        fq_prefix = re.sub('\.clean|\.nodup|\.cut', '', fq_prefix)
        fq_type = seq_type(fq)

        map_prefix = os.path.join(align_path, '%s.map_%s' % (fq_prefix, reference))
        unmap_prefix = os.path.join(align_path, '%s.not_%s' % (fq_prefix, reference))
        map_bam = map_prefix + '.bam'
        # map_bed = map_prefix + '.bed'
        map_log = map_prefix + '.%s.log' % args['aligner']
        unmap_fq = unmap_prefix + '.%s' % fq_type

        return [fq_prefix, map_bam, map_log, unmap_fq, reference]

    
    def _index_builder(self, rRNA=False, genome=None):
        """Return the genome index
        """
        args = self.args
        aligner = args['aligner']

        if not genome:
            genome = self.genome
        # check aligner
        if aligner == 'bowtie':
            index = Genome(genome, repeat_masked_genome=args['repeat_masked_genome']).bowtie_index(rRNA=rRNA)
        elif aligner == 'bowtie2':
            index = Genome(genome, repeat_masked_genome=args['repeat_masked_genome']).bowtie2_index(rRNA=rRNA)
        elif aligner == 'STAR':
            index = Genome(genome, repeat_masked_genome=args['repeat_masked_genome']).star_index(rRNA=rRNA)
        else:
            logging.error('unknown aligner: %s' % aligner)
            index = None # unknonwn aligner
        return index


    def _index_list(self):
        """List the align index (es) for the job
        rRNA, reference, genome
        """
        args = self.args

        # aligner index
        idx = {
            'genome_rRNA': self._index_builder(rRNA=True),
            'genome': self._index_builder(rRNA=False),
            'sp_rRNA': self._index_builder(rRNA=True, genome=args['spikein']),
            'sp': self._index_builder(rRNA=False, genome=args['spikein'])}

        # determine
        if not args['align_to_rRNA']:
            idx['genome_rRNA'] = None

        # save in dict
        return idx # dictionary


    def wrap_log(self, log):
        """Wrapper alignment log file, save as json"""
        args = self.args
        j_file = Alignment_log(log, args['unique_only']).saveas() # save as json


    def bowtie_se(self, fq, index, reference=None, unique_map=False, align_path=None):
        """Run bowtie

        arguments:
        reference: genome, genome_rRNA, spikein, spikein_rRNA
        """
        args = self.args

        bowtie_exe = which('bowtie')
        if not align_path:
            align_path = args['path_out']

        # output directory
        prefix, map_bam, map_log, unmap, reference = self._path_init(fq, index, 
            reference, align_path)

        # determine parameters
        n_map = args['n_map']
        if n_map < 1:
            n_map = 1 # default
        if unique_map:
            para_unique = '-m 1'
        else:
            para_unique = '-v 2 -k %s' % n_map # default: 1

        if seq_type(fq) == 'fasta':
            para_fq = '-f'
        else:
            para_fq = '-q'

        # file exists
        if os.path.exists(map_bam) and args['overwrite'] is False:
            logging.info('bam file exists: %s' % map_bam)
        else:
            c1 = '%s %s %s -p %s --mm --best --sam --no-unal --un %s %s \
                %s' % (bowtie_exe, para_fq, para_unique, args['threads'], 
                    unmap, index, fq)
            c2 = 'samtools view -bhS -F 0x4 -@ %s -' % args['threads']
            c3 = 'samtools sort -@ %s -o %s -' % (args['threads'], map_bam)
            with open(map_log, 'wt') as ff:
                p1 = subprocess.Popen(shlex.split(c1), stdout=subprocess.PIPE,
                                      stderr=ff)
                p2 = subprocess.Popen(shlex.split(c2), stdin=p1.stdout,
                                      stdout=subprocess.PIPE)
                p3 = subprocess.Popen(shlex.split(c3), stdin=p2.stdout)
                px = p3.communicate()

        # process log file
        self.wrap_log(map_log)
        return [map_bam, unmap]


    def bowtie2_se(self, fq, index, reference=None, unique_map=False, align_path=None):
        """Run bowtie2

        arguments:
        reference: genome, genome_rRNA, spikein, spikein_rRNA
        """
        args = self.args

        bowtie2_exe = which('bowtie2')
        if not align_path:
            align_path = args['path_out']

        # output directory
        prefix, map_bam, map_log, unmap, reference = self._path_init(fq, index, 
            reference, align_path)
        
        # determine parameters
        if unique_map:
            para_unique = '-q 10'
        else:
            para_unique = '-q 0'
        
        # multi map
        n_map = args['n_map']
        if n_map == 0:
            # n_map = 1 # default 1, report 1 hit for each read
            # default: #look for multiple alignments, report best, with MAPQ
            para_fq = ''
        else:
            para_fq = '-k %s' % n_map

        # fq type
        if seq_type(fq) == 'fasta':
            para_fq = para_fq + ' -f'
        else:
            para_fq = para_fq + ' -q'

        # file exists
        if os.path.exists(map_bam) and args['overwrite'] is False:
            logging.info('bam file exists: %s' % map_bam)
        else:
            c1 = '%s %s -p %s --very-sensitive-local --mm --no-unal --un %s -x %s -U %s' % (bowtie2_exe, 
                para_fq, args['threads'], unmap, index, fq)
            # c1 = '%s --very-sensitive-local --mm --no-unal -p %s --un %s -x %s -U %s' % (bowtie2_exe, 
            #     args['threads'], unmap, index, fq)
            c2 = 'samtools view -bhS -F 0x4 -@ %s %s -' % (args['threads'], para_unique)
            c3 = 'samtools sort -@ %s -o %s -' % (args['threads'], map_bam)
            with open(map_log, 'wt') as ff:
                p1 = subprocess.Popen(shlex.split(c1), stdout=subprocess.PIPE,
                                      stderr=ff)
                p2 = subprocess.Popen(shlex.split(c2), stdin=p1.stdout,
                                      stdout=subprocess.PIPE)
                p3 = subprocess.Popen(shlex.split(c3), stdin=p2.stdout)
                px = p3.communicate()

        # process log file
        self.wrap_log(map_log)
        return [map_bam, unmap]    


    def star_se(self, fq, index, reference, unique_map=False, align_path=None):
        """Run STAR, default kwargs
        
        args['unique_only'] is TRUE, unique_map=True:
        """
        args = self.args
        star_exe = which('STAR')
        if not align_path:
            align_path = args['path_out']

        # output directory
        prefix, map_bam, map_log, unmap, reference = self._path_init(fq, index, 
            reference, align_path)
        
        # determine parameters
        n_map = args['n_map']
        if n_map > 1:
            n_map = n_map # n_map default: 0
        else:
            n_map = 10 # STAR default: 10
        para_unique = '--outFilterMultimapNmax %s' % n_map

        fr = 'zcat' if is_gz(fq) else '-'
        # file exists
        map_prefix = os.path.join(align_path, prefix)
        if os.path.exists(map_bam) and args['overwrite'] is False:
            logging.info('bam file exists: %s' % map_bam)
        else:
            c1 = 'STAR --runMode alignReads \
              --genomeDir %s \
              --readFilesIn %s \
              --readFilesCommand %s \
              --outFileNamePrefix %s \
              --runThreadN %s \
              --limitOutSAMoneReadBytes 1000000 \
              --genomeLoad NoSharedMemory  \
              --limitBAMsortRAM 10000000000 \
              --outSAMtype BAM SortedByCoordinate \
              --outFilterMismatchNoverLmax 0.07 \
              --seedSearchStartLmax 20 \
              --outReadsUnmapped Fastx %s %s' % (index, fq, fr, map_prefix, 
                args['threads'], unmap, para_unique)
            p1 = subprocess.run(shlex.split(c1))
            
            # filter unique mapped reads
            if unique_map: # only unique mapped reads, -q 10
                pysam.view('-bhS', '-q', '10', '-@', str(args['threads']),
                    '-o', map_bam, map_prefix + 'Aligned.sortedByCoord.out.bam',
                    catch_stdout=False)
            else:
                os.rename(map_prefix + 'Aligned.sortedByCoord.out.bam', map_bam)
            os.rename(map_prefix + 'Unmapped.out.mate1', unmap)
            os.rename(map_prefix + 'Log.final.out', map_log)

        # process log file
        self.wrap_log(map_log)
        return [map_bam, unmap]


    def align_se_batch(self, fq, align_path=None):
        """Align reads to multiple indexes in specific order,
        return align.stat, map_bam, unmap_reads
        determine the index order
        return bam_files
        """
        args = self.args

        if not align_path:
            align_path = args['path_out']

        # define aligner
        aligner_dict = {
            'bowtie': self.bowtie_se,
            'bowtie2': self.bowtie2_se,
            'STAR': self.star_se}

        aligner_exe = aligner_dict.get(args['aligner'], None) # determine aligner
        if not aligner_exe:
            raise ValueError('unknown aligner: %s' % args['aligner'])

        # get all index in order
        index_dict = self._index_list() # genome_rRNA, genome, sp_rRNA, sp
        bam_files = []
        fq_input = fq

        # 1. genome_rRNA (rRNA: both unique, multiple)
        idx1 = index_dict['genome_rRNA']
        if idx1 is None:
            raise ValueError('genome_rRNA index not found: %s' % args['genome'])
        reference = self.genome + '_rRNA'
        bam_idx1, unmap_idx1 = aligner_exe(fq=fq_input, index=idx1, 
            reference=reference, unique_map=False, 
            align_path=align_path)
        fq_input = unmap_idx1

        # 2. genome
        idx2 = index_dict['genome']
        if idx2 is None:
            raise ValueError('genome index not found: %s' % args['genome'])

        reference = self.genome
        bam_idx2, unmap_idx2 = aligner_exe(fq=fq_input, index=idx2, 
            reference=reference, unique_map=args['unique_only'], 
            align_path=align_path)
        fq_input = unmap_idx2

        if args['spikein']: # add spikein
            # 3. sp_rRNA  (rRNA: both unique, multiple)
            idx3 = index_dict['sp_rRNA']
            reference = args['spikein'] + '_rRNA'
            bam_idx3, unmap_idx3 = aligner_exe(fq=fq_input, index=idx3, 
                reference=reference, unique_map=False,
                align_path=align_path)
            fq_input = unmap_idx3

            # 4. sp (optional)
            idx4 = index_dict['sp']
            reference = args['spikein']
            bam_idx4, unmap_idx4 = aligner_exe(fq=fq_input, index=idx4, 
                reference=reference, unique_map=args['unique_only'],
                align_path=align_path)
            fq_input = unmap_idx4

            bam_files = [bam_idx1, bam_idx2, bam_idx3, bam_idx4]
        else: # no spikein
            bam_idx3 = bam_idx4 = None

        bam_files = [bam_idx1, bam_idx2, bam_idx3, bam_idx4]
        bam_files = list(filter(partial(is_not, None), bam_files)) # remove None

        return bam_files


    def align_extra(self, fq, align_path=None):
        """Align reads to extra index
        such as GFP, white, firefly, transposon
        return bam_files
        """
        args = self.args

        if not align_path:
            align_path = args['path_out']

        # define aligner
        aligner_dict = {
            'bowtie': self.bowtie_se,
            'bowtie2': self.bowtie2_se,
            'STAR': self.star_se}

        aligner_exe = aligner_dict.get(args['aligner'], None) # determine aligner
        if not aligner_exe:
            raise ValueError('unknown aligner: %s' % args['aligner'])

        # get all index in order
        # index_ext = args['index_ext']
        bam_ext_list = []

        for ext in args['index_ext']:
            if index_validator(ext, args['aligner']):
                reference = os.path.basename(ext)
                bam_ext, unmap_ext = aligner_exe(fq=fq, 
                    index=ext, reference=reference, unique_map=True, 
                    align_path=align_path)
            else:
                bam_ext = None
            bam_ext_list.append(bam_ext)
        return bam_ext_list


    def run_extra(self):
        """Run the alignment for specific fastq file onto extra index
        1. run alignment for each replicate
        2. merge replicates
        3. run log parser, in json format
        4. organize the log files, saved in one report, including the following groups:
        """
        args = self.args

        bam_out = []

        # run alignment for replicates
        for fq in args['fqs']:
            logging.info('alignment: %s' % fq)
            fq_prefix = file_prefix(fq)[0]
            fq_prefix = re.sub('\.clean|\.nodup|\.cut', '', fq_prefix)
            fq_path = os.path.join(args['path_out'], 'extra_mapping', fq_prefix)
            assert is_path(fq_path)
            logging.info('align to index_ext: %s' % fq_prefix)
            bam_files = self.align_extra(fq, fq_path)
            bam_out.append(bam_files) #
            # stat alignment
            Alignment_stat(fq_path).saveas()


        # merge bam files
        if args['merge_rep']: 
            merged_path = os.path.join(args['path_out'], 'extra_mapping', args['smp_name'])
            merged_files = []
            if len(bam_out) > 1: # for multiple bam files
                assert is_path(merged_path)
                for i in range(len(bam_out[0])):
                    rep_bam_files = [b[i] for b in bam_out]
                    merged_suffix = str_common(rep_bam_files, suffix=True)
                    merged_suffix = re.sub('^_[12]|_R[12]', '', merged_suffix)
                    merged_bam_name = args['smp_name'] + merged_suffix
                    merged_bam_file =  os.path.join(merged_path, merged_bam_name)
                    if os.path.exists(merged_bam_file) and args['overwrite'] is False:
                        logging.info('file exists: %s' % merged_bam_file)
                    else:
                        tmp = bam_merge(rep_bam_files, merged_bam_file)
                    merged_files.append(merged_bam_file)
                Alignment_stat(merged_path).saveas()
                bam_out.append(merged_files)

        return bam_out


    def run(self):
        """Run the alignment for specific fastq file
        1. run alignment for each replicate
        2. merge replicates
        3. run log parser, in json format
        4. organize the log files, saved in one report, including the following groups:
        genome_rRNA, genome_unique, genome_multi, sp_rRNA, sp_unique, sp_multi, unmap
        """
        args = self.args

        bam_out = []

        # run alignment for replicates
        for fq in args['fqs']:
            logging.info('alignment: %s' % fq)
            fq_prefix = file_prefix(fq)[0]
            fq_prefix = re.sub('\.clean|\.nodup|\.cut', '', fq_prefix)
            fq_path = os.path.join(args['path_out'], fq_prefix)
            assert is_path(fq_path)
            bam_files = self.align_se_batch(fq, fq_path)
            bam_out.append(bam_files) # 
            # stat alignment
            Alignment_stat(fq_path).saveas()

        # merge bam files
        if args['merge_rep']: 
            merged_path = os.path.join(args['path_out'], args['smp_name'])
            merged_files = []
            if len(bam_out) > 1: # for multiple bam files
                assert is_path(merged_path)
                for i in range(len(bam_out[0])):
                    rep_bam_files = [b[i] for b in bam_out]
                    merged_suffix = str_common(rep_bam_files, suffix=True)
                    merged_suffix = re.sub('^_[12]|_R[12]', '', merged_suffix)
                    merged_bam_name = args['smp_name'] + merged_suffix
                    merged_bam_file =  os.path.join(merged_path, merged_bam_name)
                    if os.path.exists(merged_bam_file) and args['overwrite'] is False:
                        logging.info('file exists: %s' % merged_bam_file)
                    else:
                        tmp = bam_merge(rep_bam_files, merged_bam_file)
                    merged_files.append(merged_bam_file)
                Alignment_stat(merged_path).saveas()
                bam_out.append(merged_files)

        # make short names for genome bam files
        genome_bam_files = []
        for b in bam_out: # nested array
            bam_from = b[1]
            bam_to = os.path.join(os.path.dirname(bam_from),
                filename_shorter(bam_from))
            if not os.path.exists(bam_to):
                os.symlink(os.path.basename(bam_from), bam_to)
            if not os.path.exists(bam_to + '.bai'):
                if not os.path.exists(bam_from + '.bai'):
                    if os.path.getsize(bam_from) < 1000: # !!!! empty bam files
                        continue
                    else:
                        pysam.index(bam_from) # empty bam
                os.symlink(os.path.basename(bam_from) + '.bai', bam_to + '.bai')
            genome_bam_files.append(bam_to)


        # run extra index mapping
        ext_bam_files = None
        if not args['index_ext'] is None:
            ext_bam_files = self.run_extra()
        
        return [genome_bam_files, ext_bam_files]

