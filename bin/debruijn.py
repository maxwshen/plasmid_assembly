# debruijn.py
# 

import sys, string, datetime, random, copy, os
import numpy as np
# from Bio import SeqIO
from collections import defaultdict

def main():
  reads_fn = sys.argv[1]
  _k = int(sys.argv[2])

  db = Debruijn_Graph(reads_fn, _k)

  # kmer_freq(reads_fn, _k)

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
    self.nodes = dict()   # Key = (k-1)-mer, Val = edge obj
    self.edges = dict()   # Key = kmer, Val = edge obj

    with open(reads_fn) as f:
      for i, line in enumerate(f):
        if i % 4 == 1:
          print i
          curr_read = line.strip()
          for i in range(len(curr_read) - _k + 1):
            kmer = curr_read[i : i + _k]
            if 'N' not in kmer:
              self.add_to_graph(kmer)

  def add_to_graph(self, kmer):
    # Inserts node and edge objects into self.nodes/edges, and updates internal class vars 
    pref = kmer[:-1]
    suff = kmer[1:]
    if pref not in self.nodes:
      self.nodes[pref] = Node(pref)
    if suff not in self.nodes:
      self.nodes[suff] = Node(suff)
    if kmer not in self.edges:
      self.edges[kmer] = Edge(kmer)
    
    pn = self.nodes[pref]
    sn = self.nodes[suff]
    e = self.edges[kmer]

    e.cov += 1
    pn.add_outedge(e)
    sn.add_inedge(e)
    e.add_nodes(pn, sn)
    return

  def allcov(self):
    return [s.cov for s in self.edges]

class Node:
  # Nodes for the de bruijn graph
  def __init__(self, k1mer):
    self.inc = []
    self.out = []
    self.k1mer = k1mer

  def add_inedge(self, e):
    if e not in self.inc:
      self.inc.append(e)

  def add_outedge(self, e):
    if e not in self.out:
      self.out.append(e)

class Edge:
  def __init__(self, kmer):
    self.inc = []
    self.out = []
    self.kmer = kmer
    self.cov = 0

  def add_nodes(self, pn, sn):
    if pn not in self.inc:
      self.inc.append(pn)
    if sn not in self.out:
      self.out.append(sn)

# Initiates program and records total time
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start