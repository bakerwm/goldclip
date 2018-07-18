#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import os
import re
from goldclip.helper import *
from goldclip.bin.bed_annotation import *

__author__ = 'Ming Wang <wangm08@hotmail.com>'
__copyright__ = '2018 by Ming Wang <wangm08@hotmail.com>'
__license__ = 'MIT'
__email__ = 'wangm08@hotmail.com'
__version__ = '0.0.1'


## functions
def json_reader(file):
    """
    load json file as dict
    """
    with open(file) as f:
        return json.load(f)


def json_writer(d, file):
    """
    write dict to file in json format
    """
    with open(file, 'w') as f:
        json.dump(d, f)


## wrapper functions for cutadapt ##
def cutadapt_log_parser(log):
    """
    Parsing the output of cutadpat
    support multiple run output
    save data in JSON format
    """
    logdict = {}
    _cutadapt_ver = []
    _python_ver = []
    _cmd = []
    _fname = []
    _raw = []
    _clean = []
    outfile = os.path.splitext(log)[0] + ".json"
    with open(log, 'r') as f:
        for line in f.readlines():
            if(len(re.findall('^This is cutadapt', line)) == 1):
                r1 = re.findall(r'(cutadapt|Python) (\d+\.\d+\.?\d+?)', line) # 1st line, tools
                for i in r1: logdict[i[0]] = i[1]
            elif('Command line parameters' in line):
                _cmd.append(re.sub(r'Command line parameters: ', '', line).strip('\n'))
                _fname.append(os.path.basename(line.split()[-1])) # file name
            elif('Total reads processed:' in line):
                 _raw.append(re.sub(r',', '', line.split()[-1])) # first input
            elif('Reads written (passing filters):' in line):
                 _clean.append(re.sub(r',', '', line.split()[-2])) # output
            else:
                 continue
    logdict['filename'] = _fname[0]
    logdict['run_cutadapt_times'] = len(_raw)
    logdict['command_lines'] = _cmd
    logdict['raw'] = _raw[0]
    logdict['clean'] = _clean[-1]
    logdict['clean_pct'] = '{:.2f}%'.format(int(logdict['clean'])/int(logdict['raw'])*100)
    with open(outfile, 'w') as fo:
        json.dump(logdict, fo, indent = 4)
    return outfile


## wrapper functions for bowtie mapping ##

def bowtie_log_parser(path):
    """
    Parsing the log file of bowtie (to stderr)
    fetch Input, mapped, unmapped reads
    save in JSON format
    return dict of all values
    """
    # !!! to-do
    # convert to DataFrame, json, dict
    # unmap = input - reported
    logdict = {}
    _input = []
    _one_hit = []
    _not_hit = []
    _rpt = []
    with open(path, 'r') as f:
        for line in f.readlines():
            if 'reads processed' in line:
                _num1 = re.sub(',', '', line.split()[-1])
                _input.append(int(_num1))
            elif 'at least one reported' in line:
                _num2 = re.sub(',', '', line.split()[-2])
                _one_hit.append(int(_num2))
            elif 'Reported' in line:
                _num4 = re.sub(',', '', line.split()[1])
                _rpt.append(int(_num4))
            elif 'failed to align' in line:
                _num3 = re.sub(',', '', line.split()[-2])
                _not_hit.append(int(_num3))
            else:
                continue
    logdict['input_reads'] = _input[0] # first one
    # logdict['mapped'] = sum(_one_hit) # sum all
    logdict['mapped'] = sum(_rpt) # sum all
    logdict['unmapped'] = int(_input[-1]) - int(_one_hit[-1]) # -m suppress
    logdict['map_pct'] = '{:.2f}%'.\
        format(int(logdict['mapped']) / int(logdict['input_reads'])*100)
    json_out = os.path.splitext(path)[0] + '.json'
    with open(json_out, 'w') as fo:
        json.dump(logdict, fo, indent = 4)
    return logdict



def rep_map_wrapper(path, save = True):
    """
    wrap all bowtie log files, [only for this script] namespace 
    summarize mapping and RTStops 
    header: name, group, read, RTStop
    input: list of json files
    output: pd.DataFrame
    """
    def _json_wrapper(file):
        """
        parsing only one json file
        output: type, count
        """
        try:
            group = file.split('.')[-2]
            group = re.sub('map_', '', group)
            name = os.path.basename(file).split('.')[0]
            with open(file, 'r') as f:
                da = json.load(f)
            return [name, group, da['input_reads'], da['mapped'], 
                    da['unmapped']]
        except IOError:
            print('json file faild: ' + file)
    ## sort files by len
    json_files = sorted(glob.glob(path + '/*.json'), key = len)
    rep_prefix = os.path.basename(path)
    df = pd.DataFrame(columns = ['name', 'group', 'read'])
    #for ff in json_files:
    for n in range(len(json_files)):
        _, g, _, c, _ = _json_wrapper(json_files[n])
        # the first one - genome
        g = 'spikein' if g == 'genome' and n == 0 else g
        df = df.append(pd.DataFrame([[rep_prefix, g, c]], 
            columns = ['name', 'group', 'read']), ignore_index = True)
    n1, g1, _, _, c1 = _json_wrapper(json_files[-1]) # unmap
    df = df.append(pd.DataFrame([[n1, 'unmapped', c1]], 
            columns = ['name', 'group', 'read']), ignore_index = True)
    save_csv = os.path.join(os.path.dirname(path), rep_prefix + '.mapping_stat.csv')
    if save:
        try:
            df.to_csv(save_csv, ',', header = True, index = False)
        except IOError:
            print('[failed] saving data to file: ' + save_csv)
    return df



def merge_map_wrapper(path, save = True):
    """
    count BAM files
    Output: pd.DataFrame
    """
    b_files = sorted(glob.glob(path + '/*.bam'), key = len) # bam files
    b_prefix = os.path.basename(path)
    df = pd.DataFrame(columns = ['name', 'group', 'read'])
    #for b in b_files:
    for n in range(len(b_files)):
        b_cnt = pysam.AlignmentFile(b_files[n], 'rb').count()
        group = os.path.basename(b_files[n]).split('.')[-2] #
        group = re.sub('^map_', '', group)
        group = 'spikein' if n == 0 and group == 'genome' else group
        df = df.append(pd.DataFrame([[b_prefix, group, b_cnt]],
                       columns = ['name', 'group', 'read']),
                       ignore_index = True)
    save_csv = os.path.join(os.path.dirname(path), b_prefix + '.mapping_stat.csv')
    if save:
        try:
            df.to_csv(save_csv, ',', header = True, index = False)
        except IOError:
            print('[failed] saving data to file: ' + save_csv)
    return df


##-------------------------------------------------------------------##
## functions for report

##--------------------##
## figure 1
def trim_wrapper(path, smp_name='demo'):
    """
    trimming and remove duplicates
    input: /path_out/input_reads/
    """
    j_files = sorted(glob.glob(path + '/*.cutadapt.json'))
    da = []
    for j in j_files:
        id = re.sub(r'.cutadapt.json', '', os.path.basename(j))
        d = json_reader(j)
        nodup = os.path.join(os.path.dirname(j), id + '.clean.nodup.txt')
        with open(nodup) as f:
            d['nodup'] = next(f).rstrip()
        tooshort = int(d['raw']) - int(d['clean'])
        dup = int(d['clean']) - int(d['nodup'])
        dn = pd.DataFrame({'group': ['raw', 'too_short', 'PCR_dup', 'no_dup'],
            id: [d['raw'], tooshort, dup, d['nodup']]})
        dn.set_index('group', inplace = True)
        da.append(dn)
    df = pd.concat(da, axis = 1)
    df = df.apply(pd.to_numeric)
    # add merge data
    df.insert(0, smp_name, df.sum(axis = 1))
    return df

# path_trim = os.path.join(path_out, 'input_reads')
# df = trim_wrapper(path_trim, smp_name)
# print(df)


## mapping pct
def map_wrapper(path, smp_name = 'demo'):
    """
    mapping to various genome
    input: /path_out/genome_mapping/
    """
    m_files = sorted(glob.glob(os.path.join(path, '*.mapping_stat.csv')))
    ma = []
    for m in m_files:
        id = re.sub(r'.mapping_stat.csv', '', os.path.basename(m))
        if id == smp_name:
            continue
        dm = pd.read_csv(m, ',').filter(items = ['group', 'read'])
        dm.set_index('group', inplace = True)
        #dm2 = dm.transpose().rename(index = {'read': id})
        dm2 = dm.rename(columns = {'read': id})
        ma.append(dm2)
    df = pd.concat(ma, axis = 1)
    # add merge data
    df.insert(0, smp_name, df.sum(axis = 1))
    return df

##--------------------##
## figure 2
def anno_run(bed, genome, group, output):
    """
    annotate bed files
    output : Queue (df)
    """
    df = bed_annotator(bed, genome, group)
    output.put(df)


def _bed_anno(bed_files, genome, group):
    """
    intput: genome_mapping/ 
    return the genome mapped bed files
    """
    #call peaks in parallel
    output = mp.Queue()
    processes = [mp.Process(target = anno_run, 
        args = (b, genome, group, output)) for b in bed_files]
    for p in processes:
        p.start() #start process
    for p in processes:
        p.join() #exit completed process
    results = [output.get() for p in processes]
    df = pd.concat(results, axis = 1)
    return df

##--------------------##
## figure 3
def bam_corr(bam1, bam2, path_out, window = 10000, multi_cores = 16):
    """
    calculate the Pearson correlation between BAM files
    use window size: 500, 1k, 10k, 100k,
    optional: gene_level, transcript_level, ...
    output: cor, p_val
    #
    use deeptools to calculate bin-count
    """
    bam_ids = [os.path.splitext(os.path.basename(f))[0] for f in [bam1, bam2]]
    out_npz = os.path.join(path_out, 'results.npz')
    out_tab = os.path.join(path_out, 'results.tab')
    cor_png = os.path.join(path_out, 'cor_scatter.png')
    cor_tab = os.path.join(path_out, 'cor_matrix.tab')
    c1 = 'multiBamSummary bins --bamfiles {} {} -out {} --outRawCounts {} \
        -p {} --binSize {}'.format(bam1, bam2, out_npz, out_tab, multi_cores,
        window)
    c2 = 'plotCorrelation -in {} --corMethod pearson --skipZeros \
        --removeOutliers -T {} --whatToPlot scatterplot -o {} \
        --outFileCorMatrix {}'.format(out_npz, 'Pearson_correlation', cor_png,
        cor_tab)
    cmd1 = shlex.split(c1)
    cmd2 = shlex.split(c2)
    p = subprocess.run(cmd1)
    p = subprocess.run(cmd2)
    # cal pearson value
    df = pd.read_csv(out_tab, '\t')
    cor_mat = df.drop(df.columns[:3], axis = 1).corr('pearson')
    return cor_mat

    # multiBamSummary bins 
    #     --bamfiles bam1 bam2
    #     --out results.npz
    #     --outRawCounts results.tab
    #     --binSize window
    #     --smartLabels
    #     --binSize window
    #     -p multi_cores
        

    # plotCorrelation 
    #     -in out_npz 
    #     --corMethod pearson 
    #     --skipZeros
    #     --plotTitle "Pearson Correlation"
    #     --whatToPlot scatterplot
    #     -o cor_scatter.png
    #     --outFileCorMatrix cor_matrix.tab


    # plotCorrelation 
    #     -in out_npz 
    #     --corMethod pearson 
    #     --skipZeros
    #     --plotTitle "Pearson Correlation"
    #     --whatToPlot heatmap
    #     --colorMap RdYlBu 
    #     --plotNumbers
    #     -o cor_heatmap.png
    #     --outFileCorMatrix cor_matrix.tab


    # plotPCA 
    #     -in results.npz 
    #     -o results_PCA.png
    #     --plotTitle 'PCA of read counts'



def bed_conservation(bed_in, bed_out, extend = 0, genome = 'hg19'):
    con_phylop = os.path.join(goldclip_home, 'data', genome, 'pub_data', \
                              'phyloP100way', genome + '.100way.phyloP100way.bw')
    # create temp file, BED-6
    tmp = tempfile.mkstemp(suffix = ".bed")[1]
    # fix BED6
    df_in = pd.read_csv(bed_in, '\t', header = None).iloc[:, 0:6]
    df_in[[4]] = df_in[[4]].astype(int) # convert score to int
    df_in = df_in.astype('str')
    df_in.loc[:, 3] = df_in.loc[:, 0] + df_in.loc[:, 1] + df_in.loc[:, 2] + df_in.loc[:, 3]
    df_in.to_csv(tmp, '\t', header = False, index = False)
    # cmd
    c = 'bigWigAverageOverBed {} {} {} -bedOut={}'.format(con_phylop, tmp, bed_out + '.tab', bed_out)
    if extend > 0:
        c += ' -sampleAroundCenter={}'.format(extend)
    cmd = shlex.split(c)
    p = subprocess.run(cmd)






