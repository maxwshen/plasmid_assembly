# debruijn.py
# 

import sys, string, datetime, random, copy, os
import numpy as np
from Bio import SeqIO
from collections import defaultdict

def main():
  reads_fn = sys.argv[1]
  _k = sys.argv[2]

  build_debruijn_graph(reads_fn, _k)

  return

def build_debruijn_graph(reads_fn, _k):
  for record in SeqIO.parse(reads_fn, 'fastq'):
    curr_read = record.seq

  return

# Initiates program and records total time
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start