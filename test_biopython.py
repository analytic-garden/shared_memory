#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_biopython.py

@author: Bill Thompson
@license: GPL 3
@copyright: 2020_03_05
"""
import time
import os
import psutil
from collections import Counter
from Bio import AlignIO

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

def gisaid_get_columns_variation(align):
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
    column_variation = list()
    for col in range(align.get_alignment_length()):
        c = Counter(align[:, col])
        column_variation.append(c)

    return column_variation
           
def main():
    align_file = 'msa_usa.fasta'
    
    t1 =time.time()
    align = AlignIO.read(align_file, 'fasta')
    print('Load time:', time.time() - t1)
    print_memory_use()
    
    t1 =time.time()
    column_variation = gisaid_get_columns_variation(align)
    print('Calculation time:', time.time() - t1)
    print_memory_use()

if __name__ == "__main__":
    main()
