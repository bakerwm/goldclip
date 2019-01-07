
"""
Prepare arguments for goldclip pipeline
"""

def args_init(args=None, demx=False, trim=False, align=False, call_peak=False):
        """Inititate the arguments, assign the default values to arg
        """
        if isinstance(args, dict):
            pass
        elif args is None:
            args = {} # init dictionary
        else:
            raise Exception('unknown argument: args=%s' % args)

        args['fq1'] = args.get('fq1', None)
        args['fq2'] = args.get('fq2', None)
        args['path_out'] = args.get('path_out', None)

        ## optional
        args['genome_path'] = args.get('genome_path', None)
        args['overwrite'] = args.get('overwrite', False)
        args['threads'] = args.get('threads', 8)

        ## demx
        if demx:
            args['demx_type'] = args.get('demx_type', 'p7') # p7, barcode, both
            args['n_mismatch'] = args.get('n_mismatch', 0)
            args['bc_n_left'] = args.get('bc_n_left', 3)
            args['bc_n_right'] = args.get('bc_n_right', 2)
            args['bc_in_read'] = args.get('bc_in_read', 1)
            args['cut'] = args.get('cut', False)

        ## trimming
        if trim:
            args['adapter3'] = args.get('adapter3', 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC')
            args['len_min']  = args.get('len_min', 15)
            args['adapter5'] = args.get('adapter5', '')
            args['qual_min'] = args.get('qual_min', 20)
            args['error_rate'] = args.get('error_rate', 0.1)
            args['overlap'] = args.get('overlap', 3)
            args['percent'] = args.get('percent', 80)
            args['rm_untrim'] = args.get('rm_untrim', False)
            args['keep_name'] = args.get('keep_name', True)
            args['adapter_sliding'] = args.get('adapter_sliding', False)
            args['trim_times'] = args.get('trim_times', 1)
            args['double_trim'] = args.get('double_trim', False)
            args['rm_dup'] = args.get('rm_dup', True)
            args['cut_before_trim'] = args.get('cut_before_trim', '0')
            args['cut_after_trim'] = args.get('cut_after_trim', '7,-7') # NSR
            args['trim_to_length'] = args.get('trim_to_length', 0)

            ## PE trimming options
            args['AD3'] = args.get('AD3', 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT')
            args['AD5'] = args.get('AD5', None)

        ## alignment
        if align:
            args['spikein'] = args.get('spikein', None)
            args['index_ext'] = args.get('index_ext', None)
            args['unique_only'] = args.get('unique_only', True) # unique map
            args['aligner'] = args.get('aligner', 'bowtie') # bowtie alignment
            args['align-to-rRNA'] = args.get('align-to-rRNA', True)
            args['n_map'] = args.get('n_map', 0)
            args['align_to_rRNA'] = args.get('align_to_rRNA', True)
            args['repeat_masked_genome'] = args.get('repeat_masked_genome', False)
            args['merge_rep'] = args.get('merge_rep', True)

        ## peak-calling
        if call_peak:
            args['peak_caller'] = args.get('peak_caller', 'pyicoclip')

            ## rtstop-calling
            args['threshold'] = args.get('threshold', 1)
            args['intersect'] = args.get('intersect', 0)
            args['threads'] = args.get('threads', 8)

        return args
