#!/usr/bin/env python3

"""
pre-defined arguments for specific pipeline

1. GoldCLIP - NSR

2. GoldCLIP - eCLIP type

3. GoldCLIP - iCLIP type


"""



class Argument(object):
    """Pre-defined arguments for goldclip analysis
    triming
    alignment
    peak-calling
    rtstop-calling
    report
    """

    def __init__(self, mode=1):
        """Get the default values for arguments
        all 5 modules
        return dict
        """
        self.mode = mode

    def trim(self):
        ## illumina TruSeq
        args_trim_basic = {
            'adapter3' : 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
            'adapter5' : None,
            'len_min' : 15,
            'read12' : 1,
            'qual_min' : 20,
            'err_rate' : 0.1,
            'AD3' : 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
            'AD5' : None,
            'overlap' : 3,
            'threads': 8,
            'overwrite' : False,
            'rm_untrim' : False,
            'keep_name' : True,
            'adapter_sliding' : False,
            'rm_dup' : False,
            'cut_before_trim' : 0,
            'trim_to_length' : 0,
        }

        ##
        cut_args = {
            1 : { 'cut_after_trim' : '7,-7' }, # NSR
            2 : { 'cut_after_trim' : '10,-7' }, # read1: N{10}-----{7-nt}
            3 : { 'cut_after_trim' : '9'} # read1: NNN{bc-4-nt}NN
        }

        args_trim_cut = cut_args[self.mode]

        args_trim = {**args_trim_basic, **args_trim_cut}

        return args_trim


    def align(self):
        # required
        # -i, -o, -g, -n, 
        args_align = {
            'spikein' : None,
            'index_ext' : None,
            'threads': 8,
            'unique_only' : False,
            'n_map' : 0,
            'aligner' : 'STAR',
            'align_to_rRNA' : True,
            'repeat_masked_genome' : False,
            'merge_rep' : True,
            'overwrite' : False
        }
        return args_align


    def peak(self):
        args_peak = {
            'peak_caller' : 'pyicoclip',
            'threads' : 8,
            'overwrite' : False,
        }
        return args_peak


    def rtstop(self):
        args_rtstop = {
            'threshold' : 1, # threshold
            'intersect' : 0, # intersect, 0, 1
            'overwrite' : False
        }
        return args_rtstop


    def all(self):
        args_all = {**self.trim(), **self.align(), **self.peak(), **self.rtstop()}
        return args_all
