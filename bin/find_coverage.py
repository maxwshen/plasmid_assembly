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
  # reads_fn = sys.argv[1]
  # genome_fn = sys.argv[2]
  out_fn = sys.argv[1]


  # r = 0
  # r = SeqIO.convert(reads_fn, 'fastq', reads_fn.split('.')[0] + '.fasta', 'fasta')

  # coverage(reads_fn, genome_fn)
  interpretXML(out_fn)

  return

def coverage(reads_fn, genome_fn):
  stdout, stderr = NcbiblastnCommandline(cmd = 'blastn', out = 'out.xml', outfmt = 5, query = reads_fn, db = genome_fn)()
  return

def interpretXML(out_fn):
  g = [0 for i in range(6000000)]

  with open(out_fn) as f:
    get = False
    for i, line in enumerate(f):
      if get:
        get = False
        r2 = int(line.split('>')[1].split('<')[0])
        if r1 < r2:
          for i in range(r1, r2 + 1):
            g[i] += 1
        else:
          for i in range(r2, r1 + 1):
            g[i] += 1
      key = line.split('<')[1].split('>')[0]
      if key == 'Hsp_hit-from':
        get = True
        r1 = int(line.split('>')[1].split('<')[0])

  print '\n'.join([str(s) for s in g])


# Initiates program and records total time
if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start