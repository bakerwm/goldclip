#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
call RT-Stops
"""

__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__date__ = "2018-03-21"
__version__ = "0.1"


import os
import sys
import pandas as pd
import goldclip
from goldclip.helper import *



def _df_merge(dfs, how='inner'):
    """
    merge a list of DataFame, only for RTStops
    header: ['chr', 'RTStop', 'strand', 'count']
    """
    if len(dfs) == 1:
        return dfs[0]
    else:
        df = pd.merge(dfs[0], dfs[1], how=how,
            on = ['chr', 'start', 'end', 'name', 'score', 'strand'])
        return _df_merge([df] + dfs[2:], how=how)



def _rt_extracter(df):
    """
    Extract RTStops, 1-nt upstream of reads
    positive_strand: start - 1
    negative_strand: end + 1
    Input: DataFrame (BED)
    Output: DataFrame (RTStop)
    """
    if len(df) ==  0:
        df_rt = df # empty DataFrame
    else:
        df_bed6 = df[df.columns[:6]]
        df_ext = df[df.columns[6:]]
        df_bed6['RTStop'] = df_bed6.apply(lambda row: row.start - 1 
                            if(row.strand == '+') else row.end + 1, axis = 1) #
        df_bed6 = df_bed6[['chr', 'RTStop', 'strand', 'name', 'score']]
        # columns in ext not in bed6
        m = [i for i in df_ext.columns if i not in df_bed6.columns]
        df_ext2 = df_ext.loc[:, m]
        # merge
        df_rt = pd.concat([df_bed6, df_ext2], axis=1)
    return df_rt



def _rt_counter(df, threshold=1):
    """
    convert RTReads to RTStops
    group by 'chr', 'start', 'strand'
    Input: DataFrame (RTReads)
    Output: DataFrame (grouped)
    """
    if len(df) == 0:
        return df # empty file
    else:
        df_cnt = df.groupby(['chr', 'RTStop', 'strand']).size()
        df_cnt = df_cnt[df_cnt >= threshold]
        return df_cnt.reset_index()



def _rt_to_bed(df):
    """
    convert RTRead / RTStop to BED format
    Input: ['chr', 'RTStop', 'strand', ...]
    Output: ['chr', 'start', 'end', 'name', 'score', 'strand']
    """
    df_rt = df[df.columns[:3]]
    df_ext = df[df.columns[3:]]
    ##
    name = df_ext['name'] if 'name' in df_ext.columns else 'RTStop'
    score = df_ext['score'] if 'score' in df_ext.columns else 255
    df_rt.insert(1, 'start', df_rt['RTStop'] - 1)
    df_rt.insert(3, 'name', name)
    df_rt.insert(4, 'score', score)
    # columns dup in ext
    m = [i for i in df_ext.columns if i not in df_rt.columns]
    df_ext2 = df_ext[m]
    df_bed = pd.concat([df_rt, df_ext2], axis=1)
    # rename 'RTStop' to 'end'
    df_bed2 = df_bed.rename(index=str, columns={'RTStop': 'end'})
    return df_bed2



def rt_merge(rt_stops, intersect = 1):
    """
    merge and filt RTStops 
    filter by threshold for each replicate
    Input: RTStop in DataFrame (bed6+)
    """
    if len(rt_stops) == 0:
        sys.exit('no bed file detected')
    elif len(rt_stops) == 1:
        m = rt_stops[0].rename(index=int, columns={6: 'sum'})
    else:
        if intersect ==  1:
            m = _df_merge(rt_stops, how='inner') # intersect
        else:
            m = _df_merge(rt_stops, how='outer') # union
            m = m.fillna(0) # replace NaN by 0
        # rename header, rep1 to repN
        m.columns = list(m.columns[:6]) + \
            ['rep' + str(i + 1) for i in range(len(m.columns) - 6)]
    
        # sum replicates
        md = m[m.columns[:6]] # id
        mc = m[m.columns[6:]].astype('int') # count
        rep_sum = mc.sum(axis=1).astype('int')
        rep_mean = mc.mean(axis=1).map('{:.4f}'.format)
        rep_std = mc.std(axis=1).map('{:.4f}'.format)
        mc = mc.assign(sum=rep_sum)
        mc = mc.assign(mean=rep_mean)
        mc = mc.assign(std=rep_std)
        m = pd.concat([md, mc], axis=1) # rtstops
    return m



def _bed_to_rtstop(bed, path_out, threshold):
    """
    convert BED to rtstops, filter by threshold
    """
    rep_rt_read = '' # pd.DataFrame, rtstop format
    rep_rt_stop = '' # pd.DataFrame, rtstop format
    bed_prefix = os.path.basename(os.path.splitext(bed)[0]) #
    path_sub = os.path.join(path_out, bed_prefix)
    assert is_path(path_sub)

    # load BED
    bed_df = bed_parser(bed) # convert to DataFrame
        
    # call RTRead
    rep_rt_read = _rt_extracter(bed_df) # rt reads [full version]
    rep_rt_read_bed = _rt_to_bed(rep_rt_read) # to BED

    # call RTStop
    rt_stop = os.path.join(path_sub, bed_prefix + '.RTStop.bed') # rt stops
    rep_rt_stop = _rt_counter(rep_rt_read, threshold) # rt stops
    rep_rt_stop_bed = _rt_to_bed(rep_rt_stop) # to BED
    rep_rt_stop_bed.to_csv(rt_stop, '\t', header=False, index=False)

    # recover RTStop to RTRead
    # 6th column = RTStop count
    rt_read = os.path.join(path_sub, bed_prefix + '.RTRead.bed') # rt reads
    df = rep_rt_stop_bed.copy()
    df = df.loc[df.index.repeat(df[0])].reset_index() #repeat rows
    rep_rt_read_stop = df.drop(['index'], axis=1)
    rep_rt_read_stop.to_csv(rt_read, '\t', header=False, index=False)
    del df # remove tmp name
        
    return [rep_rt_read_stop, rep_rt_stop_bed]



def call_rtstop(bed_files, path_out, smp_name, threshold, intersect, overwrite=False):
    """
    merge rt reads and rt stops
    """
    logging.info('calling RTStops')
    if len(bed_files) == 0:
        sys.exit('bed_files not exists')
    # loading rep read and stops
    rep_rt_reads = []
    rep_rt_stops = []
    report_rt_stops = []
    for bed in bed_files:
        if file_row_counter(bed) == 0:
            continue # skip empty files
        bed_prefix = os.path.basename(os.path.splitext(bed)[0]) #
        path_sub = os.path.join(path_out, bed_prefix)
        rt_read = os.path.join(path_sub, bed_prefix + '.RTRead.bed') # rt reads
        rt_stop = os.path.join(path_sub, bed_prefix + '.RTStop.bed') # rt stops
        if not os.path.exists(rt_read) or not os.path.exists(rt_stop) or overwrite:
            rep_rt_read, rep_rt_stop = _bed_to_rtstop(bed, path_out, threshold)
        else:    
            rep_rt_read = bed_parser(rt_read)
            rep_rt_stop = bed_parser(rt_stop)
        rep_rt_reads.append(rep_rt_read)
        rep_rt_stops.append(rep_rt_stop)
        report_rt_stops.append(rt_stop)

    # merge
    path_merge = os.path.join(path_out, smp_name)
    merge_rt_read = os.path.join(path_merge, smp_name + '.RTRead.bed')
    merge_rt_stop = os.path.join(path_merge, smp_name + '.RTStop.bed')
    if not os.path.exists(merge_rt_stop) or overwrite:
        if len(rep_rt_reads) == 0:
            return None
        assert is_path(path_merge)

        # merge RTRead
        merge_rt_read_bed = pd.concat(rep_rt_reads, axis=0, ignore_index=True)

        # merge RTstop
        merge_rt_stop_bed = rt_merge(rep_rt_stops, intersect)
        merge_rt_stop_bed.to_csv(merge_rt_stop, '\t', header=False, index=False)

        # recover RTStop to RTRead
        # 'sum' column = RTStop count
        df = merge_rt_stop_bed.copy()
        if len(rep_rt_stops) == 1:  # only one input file
            df['sum'] = df[0]
        df = df.loc[df.index.repeat(df['sum'])].reset_index() #repeat rows
        merge_rt_read_stop = df.drop(['index'], axis = 1)
        merge_rt_read_stop.to_csv(merge_rt_read, '\t', header=False, index=False)
        del df # remove tmp name 

    ## to report
    report_rt_stops.append(merge_rt_stop)
    return report_rt_stops


## EOF