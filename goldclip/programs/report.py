
import os
import sys
sys.path.append('..') # 



class Report:

    def __init__(self, *args, **kwargs):
        self.kwargs = kwargs

    def run(self):
        f = self.kwargs['i']
        g = self.kwargs['g']
        print('input_file: %s' % f)
        print('gtf_file: %s' % g)

