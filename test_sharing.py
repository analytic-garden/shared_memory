#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_sharing.py

@author: Bill Thompson
@license: GPL 3
@copyright: %(date)
"""
import time
import os
import psutil
import subprocess
import numpy as np
from collections import Counter
import multiprocessing as mp
from multiprocessing import shared_memory
from Bio import SeqIO

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
    
def read_alignment_to_numpy(filename):
    """
    Read a MSA in fasta format into a numpy array.
    
    Sequences are stored in a (seq_count, align_length) numpy array
    of characters. The array is stored in column major order to speed
    up column access.
    
    Parameters
    ----------
    filename : str
        The name of an MSA file.

    Returns
    -------
    None.
    
    Requires
    -------
    All sequences must have the same length.

    """         
    # First we need to figure out how much memory we will need
    # Make a pass through the file and count number of sequences and 
    # width of each record.
    
    def get_msa_size():    
        res = subprocess.check_output(['/usr/bin/grep', '>', filename]).split(b'\n')
        res.remove(b'')
        rec = next(SeqIO.parse(filename, 'fasta'))
        seq_count = len(res)
        align_length = len(rec.seq)
        
        return (seq_count, align_length)

    seq_count, align_length = get_msa_size()
    
    # Reserve some memory for sequence data.                .
    shm = shared_memory.SharedMemory(create = True, 
                                      size = seq_count * align_length)
    align_array = np.ndarray((seq_count, align_length), 
                              dtype = '|S1', 
                              order = 'F',
                              buffer = shm.buf)
    
    # load the sequences
    seq_count = 0
    curr_line = ''
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line[0] == '>':                    
                if len(curr_line) > 0:
                    align_array[seq_count, :] = list(curr_line)
                    seq_count += 1
                curr_line = ''
            elif len(line) > 0:
                curr_line += line

    # Add the last sequence.                
    if len(curr_line) > 0:
        align_array[seq_count, :] = list(curr_line)
       
    return (align_array, shm)

def get_col_variation(shm_name, shape, col):
    """
    Get count of each nucleotide in each column.
    Get's the shared memory object as a numpy array and counts
    the nucleotides in the column.

    Parameters
    ----------
    shm_name : str
        The id of the shared memory block.
    shape : a tuple of ints.
        (number of sequences, alignmnet width)
    col : int
        The column number.

    Returns
    -------
    TYPE
        A count of each nucleotide in the column.
        
    Requires
    --------
    The np.ndarray must have the same shape and format as when
    created by read_alignment_to_numpy.

    """
    shm = shared_memory.SharedMemory(name=shm_name)
    align = np.ndarray(shape, 
                       dtype = '|S1', 
                       order = 'F',
                       buffer = shm.buf)
    
    return Counter(align[:, col])

def gisaid_get_columns_variation(shm_name, shape, num_processes = 10):
    """
    Get counts of each nucleotide type in each column of MSA.

    Parameters
    ----------
    shm_name : str
        The assigned name of the shared memory block.
    shape : a tuple of ints
        (number of sequences, alignment width)
    num_processes : int
        The number of subprocesses to start.

    Returns
    -------
    column_variation : a list
        A list counts for each column of the MSA.

    """
    pool = mp.Pool(processes = num_processes)
    column_variation = pool.starmap(get_col_variation, [(shm_name, shape, col)
                                                        for col in range(shape[1])])
        
    return column_variation

def main():
    align_file = 'msa_usa.fasta'
    
    t1 =time.time()
    align_array, shm = read_alignment_to_numpy(align_file)
    print('Load time:', time.time() - t1)
    print_memory_use()
    
    num_processors = mp.cpu_count()
    
    t1 =time.time()
    column_variation = gisaid_get_columns_variation(shm.name, 
                                                    align_array.shape,
                                                    num_processors - 1)
    print('Calculation time:', time.time() - t1)
    print_memory_use()
    
    shm.close()
    shm.unlink()

if __name__ == "__main__":
    main()
