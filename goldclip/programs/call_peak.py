
import os
import sys
sys.path.append('..') # 



class diff:

    def __init__(self, *args, **kwargs):
        self.kwargs = kwargs

    def run(self):
        a = self.kwargs['a']
        b = self.kwargs['b']
        print('matrix_a: %s' % a)
        print('matrix_b: %s' % b)

