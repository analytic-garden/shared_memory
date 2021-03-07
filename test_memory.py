#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_memory.py

@author: Bill Thompson
@license: GPL 3
@copyright: 2020_03_05
"""
import os
import psutil
import gc
from Bio import AlignIO
from Bio import SeqIO
import numpy as np
import pandas as pd

def print_memory_use():
    """
    Print percent of memory used.
    From https://stackoverflow.com/questions/276052/how-to-get-current-cpu-and-ram-usage-in-python

    Returns
    -------
    None.

    """    
    pid = os.getpid()
    py = psutil.Process(pid)
    memoryUse = py.memory_info()[0]/2.**30  # memory use in GB...I think
    print('memory use:', memoryUse)
    
def read_alignment(filename):
    # align = AlignIO.read(filename, 'fasta')
    align = list(SeqIO.parse(filename, 'fasta'))
    # n = np.ones((78193, 32280), dtype = int)
    # p = pd.DataFrame(n)
    
def main():
    align_file = 'msa_usa.fasta'

    print('Initial:')
    print_memory_use()
    for r in range(5):
        print('Round:', r + 1)
        read_alignment(align_file)
        gc.collect()
        print_memory_use()
        
if __name__ == "__main__":
    main()
