"""
settings.py
global variables
"""

import os
import sys

global goldclip_home
goldclip_home = os.environ.get("GOLDCLIP_PATH")
if not goldclip_home:
    sys.exit('Error - GOLDCLIP_PATH variable not defined, please set it to the \
              goldclip source directory')
