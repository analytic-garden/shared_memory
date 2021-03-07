#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_multi.py

@author: Bill Thompson
@license: GPL 3
@copyright: %(date)
"""
import time
import os
import psutil
from collections import Counter
from Bio import AlignIO
import multiprocessing as mp

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
    
def col_variation(align, col):
    return Counter(align[:, col])

def gisaid_get_columns_variation(align, num_processes = 10):
    """
    Get counts of each nucleotide type in each column of MSA.

    Parameters
    ----------
    align : Bio.Align,MultipleSeqAlignmnet object
        Alignment returned by Bio.AlignIO

    Returns
    -------
    column_variation : a list
        A list counts for each column of the MSA.

    """
    pool = mp.Pool(processes = num_processes)
    column_variation = pool.starmap(col_variation, [(align, col)
                                                    for col in range(align.get_alignment_length())])

    return column_variation

def main():
    align_file = 'msa_usa.fasta'
    
    t1 =time.time()
    align = AlignIO.read(align_file, 'fasta')
    print('Load time:', time.time() - t1)
    print_memory_use()
    
    num_processors = mp.cpu_count()
    
    t1 =time.time()
    column_variation = gisaid_get_columns_variation(align, num_processors - 1)
    print('Calculation time:', time.time() - t1)
    print_memory_use()
    
if __name__ == "__main__":
    main()
