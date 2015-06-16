# find_coverage.py
#
# Finds the coverage of a genome given a set of reads.
# Used to test the assumption that NGS provides highly uniform coverage 

import sys, string, datetime, random, copy, os
import numpy as np
from collections import defaultdict

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

def main():
  reads_fn = sys.argv[1]
  # genome_fn = sys.argv[2]

  r = 0
  r = SeqIO.convert(reads_fn, 'fastq', reads_fn.split('.')[0] + '.fasta', 'fasta')

  print r
  # coverage(reads_fn, genome_fn)
  return

def coverage(reads_fn, genome_fn):
  stdout, stderr = NcbiblastnCommandline(cmd = 'blastn', out = 'out.xml', outfmt = 5, query = reads_fn, db = genome_fn)()
  return

# Initiates program and records total time
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start