# debruijn.py
# 

import sys, string, datetime, random, copy, os
import numpy as np
# from Bio import SeqIO
from collections import defaultdict

def main():
  reads_fn = sys.argv[1]
  _k = int(sys.argv[2])

  kmer_freq(reads_fn, _k)

  return

def kmer_freq(reads_fn, _k):
  kmers = dict()
  with open(reads_fn) as f:
    for i, line in enumerate(f):
      if i % 4 == 1:
        curr_read = line.strip()
        for i in range(len(curr_read) - _k + 1):
          kmer = curr_read[i : i + _k]
          if 'N' not in kmer:
            if kmer not in kmers:
              kmers[kmer] = 1
            else:
              kmers[kmer] +=1

  for k in kmers:
    print kmers[k]
  return

class Debruijn_Graph:
  def __init__(self, reads_fn, _k):
    # Constructs the de bruijn graph
    self.nodes = []
    self.edges = []

    # for record in SeqIO.parse(reads_fn, 'fastq'):
    #   curr_read = record.seq
    #   for i in range(len(curr_read) - _k + 1):
    #     kmer = curr_read[i : i + _k]
    #     pref = kmer[:-1]
    #     suff = kmer[1:]

class Node:
  # Nodes for the de bruijn graph
  def __init__(self):
    self.inc = []
    self.out = []
    self.k1mer = ''

class Edge:
  def __init__(self):
    self.inc = []
    self.out = []
    self.kmer = ''


# Initiates program and records total time
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start