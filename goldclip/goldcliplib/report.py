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
__date__ = "2018-12-25"
__version__ = "0.1"



import os
import re
import random
import logging
import pandas as pd
# import numpy as np
from goldclip.helper import *
from goldclip.bin.bed_fixer import Bed_parser
from goldclip.bin.bed_motif import bed2motif
from goldclip.bin.bed_annotation import bed_annotator
from goldclip.goldcliplib.log_parser import *
from goldclip.goldcliplib.log_parser import Json_file



class Goldclip_report(object):
    """Record the output of Goldclip output

    directory structure of goldclip output

    01.trimming
    02.genome_mapping
    03.call_peaks
    04.call_rtstops
    05.report

    """

    def __init__(self, project_path, project_name, genome, group='homer',
                 window=10000, threads=8, **kwargs):
        self.project_path = project_path
        self.project_name = project_name
        self.genome = genome
        # self.group = group
        self.window = window # bam correlation
        self.threads = threads
        ## directory structure
        path_keys = [
            'trim_path', 
            'align_path', 
            'peak_path', 
            'rtstop_path', 
            'report_path']
        path_values = [
            '01.trimming',
            '02.genome_mapping',
            '03.call_peaks',
            '04.call_rtstops',
            '05.report']
        path_values = [os.path.join(self.project_path, i) for i in path_values]
        self.project_subpath = dict(zip(path_keys, path_values))


    def get_trim_stat(self):
        """Return the reads processing of fastq files
        groups: raw, too_short, PCR_dup, no_dup
        """
        path_trim = self.project_subpath['trim_path']
        dfx = [] # list of pd.DataFrame
        if os.path.exists(path_trim):
            with os.scandir(path_trim) as it:
                for entry in it:
                    if not entry.name.endswith('.cutadapt.json'):
                        continue
                    name = re.sub(r'.cutadapt.json', '', entry.name)
                    fn = os.path.join(path_trim, entry.name)
                    dd = Json_file(fn).json_reader()
                    # clean
                    fn_nodup = os.path.join(path_trim, name + '.clean_reads.txt')
                    with open(fn_nodup) as fi:
                        dd['nodup'] = next(fi).strip()
                    dd['too_short'] = int(dd['total']) - int(dd['clean'])
                    dd['dup'] = int(dd['clean']) - int(dd['nodup'])
                    dx = pd.DataFrame({'group': ['raw_count', 
                                                 'too_short', 
                                                 'PCR_dup', 
                                                 'no_dup'],
                                        name: [dd['total'], 
                                               dd['too_short'], 
                                               dd['dup'], 
                                               dd['nodup']]})
                    dx.set_index('group', inplace=True)
                    dfx.append(dx)
        df = pd.DataFrame()
        if len(dfx) > 0:
            df = pd.concat(dfx, axis=1)
            df = df.apply(pd.to_numeric)
            df.insert(0, self.project_name, df.sum(axis=1))
        return df


    def get_map_stat(self):
        """Return the reads mapping 
        categories
        """
        path_map = self.project_subpath['align_path']
        dfx = []
        with os.scandir(path_map) as it:
            for entry in it:
                if not entry.name.endswith('.mapping_stat.csv'):
                    continue
                name = re.sub(r'.mapping_stat.csv', '', entry.name)
                fn = os.path.join(path_map, entry.name)
                # skip merged
                if name == self.project_name:
                    continue
                dx1 = pd.read_csv(fn, '\t')
                dx1 = dx1.rename(index=str, columns={'Unnamed: 0': 'name'})
                dx2 = pd.melt(dx1, id_vars=['name'], var_name='group', value_name='count')
                dx3 = dx2[['group', 'count']].set_index('group').rename(index=str, columns={'count': name})
                dfx.append(dx3)
        df = pd.concat(dfx, axis=1)
        df.insert(0, self.project_name, df.sum(axis=1))
        return df


    def get_bam_file(self, bam2bed=False, rep_only=False, merge_only=False):
        path_map = self.project_subpath['align_path']
        bam_files = []
        with os.scandir(path_map) as it:
            for entry in it:
                if merge_only and not entry.name == self.project_name:
                    continue
                elif rep_only and entry.name == self.project_name:
                    continue
                else:
                    pass
                bam = os.path.join(path_map, entry.name, entry.name + '.bam')
                bed = os.path.join(path_map, entry.name, entry.name + '.bed')
                if not os.path.exists(bam):
                    continue
                if bam2bed:
                    if not os.path.exists(bed):
                        BAM(bam).to_bed()
                    bam_files.append(bed)
                else:
                    bam_files.append(bam)
        return bam_files


    def get_peak_file(self, rep_only=False, merge_only=False):
        path_peak = self.project_subpath['peak_path']
        peak_files = []
        with os.scandir(path_peak) as tools:
            for tool in tools:
                path_tool = os.path.join(path_peak, tool.name)
                with os.scandir(path_tool) as smps:
                    for smp in smps:
                        if merge_only and not smp.name == self.project_name:
                            continue
                        elif rep_only and smp.name == self.project_name:
                            continue
                        else:
                            pass
                        fn = os.path.join(path_tool, smp.name, 
                                          smp.name + '.fixed.bed')
                        if os.path.exists(fn):
                            peak_files.append(fn)
        return peak_files


    def get_rtstop_file(self, rep_only=False, merge_only=False, rt_reads=False):
        path_rtstop = self.project_subpath['rtstop_path']
        rtstop_files = []
        with os.scandir(path_rtstop) as it:
            for entry in it:
                if merge_only and not entry.name == self.project_name:
                    continue
                elif rep_only and entry.name == self.project_name:
                    continue
                else:
                    pass
                fn_ext = '.RTRead.bed' if rt_reads else '.RTStop.bed'
                fn = os.path.join(path_rtstop, entry.name, entry.name + fn_ext)
                if os.path.exists(fn):
                    rtstop_files.append(fn)
        return rtstop_files


    ##---------------------------------##
    ## statistics for figures
    ##---------------------------------##
    def figure1_trim_map(self):
        """Trim, PCR_dup, mapped reads
        """
        logging.info('figure1 reads mapping')
        path_report = self.project_subpath['report_path']
        project_name = self.project_name
        figure1_path = os.path.join(path_report, 'read_mapping')
        figure1_txt = os.path.join(figure1_path, 'read_mapping.txt')
        assert is_path(figure1_path)
        df1 = self.get_trim_stat()
        df2 = self.get_map_stat()
        df = pd.concat([df1, df2], axis=0, sort=False)
        df.to_csv(figure1_txt, '\t', header=True, index=True)
        return df


    def figure2_read_anno(self, genome, group='homer'):
        """Return categories of mapped reads
        # function
        df = bed_annotator(args.i.name, args.g, args.t, path_data)
        """
        logging.info('figure2 reads annotation')
        path_report = self.project_subpath['report_path']
        project_name = self.project_name        
        figure2_path = os.path.join(path_report, 'read_annotation')
        figure2_txt = os.path.join(figure2_path, 'read_annotation.txt')
        assert is_path(figure2_path)
        bed_files = self.get_bam_file(bam2bed=True)
        if len(bed_files) == 0:
            logging.error('failed, bed files not found: %s' % path)
            return None
        # annotate files
        dfx = []
        for bed in bed_files:
            df_anno = bed_annotator(bed, genome, group).drop(columns=['sample'])
            dfx.append(df_anno)
        df = pd.concat(dfx, axis=1).reset_index().rename(columns={'index': 'type'})
        df.to_csv(figure2_txt, '\t', header=True, index=False)
        return(df)


    def figure3_read_cor(self):
        """Return the correlation between replicates, window with fixed width
        using deeptools
        """
        logging.info('figure3 reads correlation')
        path_report = self.project_subpath['report_path']
        project_name = self.project_name        
        figure3_path = os.path.join(path_report, 'read_correlation')
        assert is_path(figure3_path)
        bam_files = self.get_bam_file(rep_only=True)
        if len(bam_files) < 2:
            logging.error('skipped, at least 2 bam files required: %s' % path_report)
            return None
        df = bam_corr(bam_files, figure3_path, window=self.window, threads=self.threads)
        return df


    def figure4_rtstop_cor(self):
        """Return the correlation between replicates, using rtstops counts
        using pandas
        """
        logging.info('figure4 rtstop correlation')
        path_report = self.project_subpath['report_path']
        project_name = self.project_name
        figure4_path = os.path.join(path_report, 'rtstop_correlation')
        figure4_txt = os.path.join(figure4_path, 'cor_matrix.tab')
        assert is_path(figure4_path)
        bed_merge = self.get_rtstop_file(merge_only=True)
        if not os.path.exists(bed_merge[0]):
            logging.error('file not exists: %s' % bed_merge)
            return None
        df = Bed_parser(bed_merge[0]).bed
        # remove the first 6 columns, last 3 columns
        # column7, column8, ...
        dfx = df.drop(df.columns[:6].tolist() + df.columns[-3:].tolist(), axis=1) #
        dfx_cor = dfx.apply(pd.to_numeric).corr(method='pearson')
        dfx_cor.to_csv(figure4_txt, '\t', header=True, index=True)
        return dfx_cor


    def figure5_peak_count(self):
        """
        Number of peaks for each sample
        """
        logging.info('figure5 peak number')
        path_report = self.project_subpath['report_path']
        project_name = self.project_name
        figure5_path = os.path.join(path_report, 'peak_number')
        figure5_txt = os.path.join(figure5_path, 'peak_number.txt')
        assert is_path(figure5_path)
        peak_files = self.get_peak_file()
        df = pd.DataFrame(columns = ['tool', 'id', 'count'])
        for p in peak_files:
            sep = p.split('/') # only fo Linux path
            tool, sample = sep[-3:-1]
            num = file_row_counter(p)
            dfx = pd.DataFrame([[tool, sample, num]], 
                               columns=['tool', 'id', 'count'])
            df = df.append(dfx, ignore_index=True)
        df.to_csv(figure5_txt, '\t', header=True, index=False)
        return df


    def figure6_peak_length(self):
        """peak length"""
        logging.info('figure6 peak length')
        path_report = self.project_subpath['report_path']
        project_name = self.project_name
        figure6_path = os.path.join(path_report, 'peak_length')
        figure6_txt = os.path.join(figure6_path, 'peak_length.txt')
        assert is_path(figure6_path)
        peak_files = self.get_peak_file()
        df = pd.DataFrame(columns = ['length', 'tool', 'id', 'count'])
        for p in peak_files:
            tmp1, sample = os.path.split(p) # basename
            tmp2 = os.path.split(tmp1)[0] # tmp
            tool = os.path.split(tmp2)[1] # pyicoclip/clipper
            dx = Bed_parser(p).bed
            dx['length'] = dx['end'] - dx['start']
            dx2 = dx['length'].value_counts().to_frame().reset_index()
            dx2['id'] = sample
            dx2['tool'] = tool
            dx2 = dx2.rename(index=str, columns={'index': 'length', 'length': 'count'})
            df = df.append(dx2, sort=False)
        df.to_csv(figure6_txt, '\t', header=True, index=False)
        return df


    def figure7_peak_anno(self, genome, group='homer'):
        """peak annotation"""
        logging.info('figure7 peak annotation')
        path_report = self.project_subpath['report_path']
        project_name = self.project_name
        figure7_path = os.path.join(path_report, 'peak_annotation')
        figure7_txt = os.path.join(figure7_path, 'peak_annotation.txt')
        assert is_path(figure7_path)
        peak_files = self.get_peak_file()
        df = pd.DataFrame(columns = ['type', 'count', 'id', 'tool'])
        for p in peak_files:
            tmp1, sample = os.path.split(p) # basename
            tmp2 = os.path.split(tmp1)[0] # tmp
            tool = os.path.split(tmp2)[1] # pyicoclip/clipper
            dx = bed_annotator(p, genome, group) # index, id
            dx = dx.reset_index()
            dx.columns = ['type', 'count', 'id']
            dx['tool'] = tool
            df = df.append(dx, ignore_index=True, sort=False)
        df.to_csv(figure7_txt, '\t', header=True, index=False)
        return df


    def figure8_peak_motif(self, genome):
        """peak motif analysis, de novo analysis"""
        logging.info('figure8 motif analysis')
        path_report = self.project_subpath['report_path']
        project_name = self.project_name
        figure8_path = os.path.join(path_report, 'peak_motif')
        assert is_path(figure8_path)
        peak_files = self.get_peak_file()
        for p in peak_files:
            tmp1, project_name = os.path.split(p) # basename
            tmp2, smp_id = os.path.split(tmp1) # tmp
            tool = os.path.split(tmp2)[1] # pyicoclip/clipper
            project_name = re.sub('.fixed.bed', '', project_name)
            logging.info('find motifs: %s - %s' % (tool, smp_id))
            peak_path = os.path.join(figure8_path, tool, smp_id)
            assert is_path(peak_path)
            bed2motif(p, genome, peak_path, '4,5,6,7,8')


    def figure9_peak_conservation(self, genome):
        """
        conservation of peaks
        phyloP1000
        """
        path = self.project_path
        project_name = self.project_name
        logging.info('figure9 peak conservation')
        figure9_path = os.path.join(path, 'results', 'peak_conservation')
        assert os.path.exists(path)
        assert is_path(figure9_path)
        assert isinstance(project_name, str)
        peak_files = self.get_peak_file()
        for peak in peak_files:
            if file_row_counter(peak) == 0:
                continue
            sep = peak.split('/') # Linux path
            tool, smp_id, project_name = sep[-3:]
            peak_sub_path = os.path.join(figure9_path, tool, smp_id)
            assert is_path(peak_sub_path)
            peak_sub_prefix = re.sub('.bed$', '', project_name)
            peak_sub_file = os.path.join(peak_sub_path, peak_sub_prefix)
            bed_conservation(peak, peak_sub_file + '.ext1k.con.bed', 1000, genome)


    def get_all_figures(self):
        """
        report all figures data
        """
        path = self.project_path
        project_name = self.project_name
        genome = self.genome
        window = self.window
        threads = self.threads
        df = self.figure1_trim_map()
        df = self.figure2_read_anno(genome=genome)
        df = self.figure3_read_cor()
        df = self.figure4_rtstop_cor()
        df = self.figure5_peak_count()
        df = self.figure6_peak_length()
        df = self.figure7_peak_anno(genome=genome)
        df = self.figure8_peak_motif(genome=genome)
        # df = self.figure9_peak_conservation(genome=genome)
