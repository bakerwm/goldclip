
__author__ = "Ming Wang"
__email__ = "wangm08@hotmail.com"
__date__ = "2018-03-21"
__version__ = "0.1"


import os
import sys
from goldclip.helper import *

class Run_all:

    def __init__(self, *args, **kwargs):
        self.kwargs = kwargs


    def run(self):
        f = self.kwargs['i']
        g = self.kwargs['g']
        print('input_file: %s' % f)
        print('gtf_file: %s' % g)