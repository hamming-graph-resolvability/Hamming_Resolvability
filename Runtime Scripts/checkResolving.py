#Python code to determine whether or not a set is resolving on a given hamming graph

from hammingGrobner import check_resolving_grobner
from itertools import product, repeat, islice#, izip
import numpy as np
import multiprocessing as mp
import argparse
import time
#from memory_profiler import profile

##################
### READ/WRITE ###
##################
#
def readTSV(f):
  L = []
  with open(f, 'r') as f:
    f.readline() #header
    for line in f:
      l = line.strip().split('\t')
      R = l[0].split(',')
      alphabet = l[1].split(',')
      k = int(l[2])
      L.append((R, alphabet, k))
  return L

def writeTSV(results, o):
  with open(o, 'w') as o:
    o.write('\t'.join(['R', 'alphabet', 'k', 'resolving?'])+'\n') #header
    for (R, alphabet, k, res) in results:
      o.write('\t'.join([','.join(map(str, R)), ','.join(map(str, alphabet)), str(k), str(res[0]), ','.join(map(str, res[1]))])+'\n')

######################
### PARALLEL BRUTE ###
######################
#
def parallelBrute(R, alphabet, k, dictSize=100000, chunkSize=100000, procs=1):
  start = time.clock()
  tags = {}
  kmers = product(alphabet, repeat=k)
  dictJobs = zip(kmers, repeat(R))
  dictChunk = list(islice(dictJobs, 0, dictSize))
  d = 0
  while len(dictChunk)>0:
    print('Brute Force: dict '+str(d))
    d += len(dictChunk)
    pool = mp.Pool(processes=procs)
    results = pool.map_async(genTag, dictChunk)
    results = results.get()
    pool.close()
    pool.join()
    
    curDict = {}
    for (tag, kmer) in results:
#      print(''.join(kmer), tag)
      if tag in curDict: return (False, (''.join(curDict[tag]), ''.join(kmer)))
      curDict[tag] = kmer

    checkKmers = product(alphabet, repeat=k)
    checkJobs = zip(checkKmers, repeat(R))
    chunk = list(islice(checkJobs, d, chunkSize+d))
    i = 0
    while len(chunk)>0:
      print('   Chunk: dict '+str(d)+' chunk '+str(i))
      i += len(chunk)
      pool = mp.Pool(processes=procs)
      results = pool.map_async(genTag, chunk)
      results = results.get()
      pool.close()
      pool.join()

      for (tag, kmer) in results:
#        print('   ', ''.join(kmer), tag)
        if tag in curDict: return (False, (''.join(curDict[tag]), ''.join(kmer)))

      chunk = list(islice(checkJobs, 0, chunkSize))

    dictChunk = list(islice(dictJobs, 0, dictSize))

  elapsed = time.clock() - start
  print('elapsed time', elapsed)
  return (True, [])

def genTag(args):
  (kmer, R) = args
  tag = ';'.join(map(str, [hammingDist(kmer, r) for r in R]))
  return (tag, kmer)

def hammingDist(a, b):
  return sum(1 for (x,y) in zip(a,b) if x!=y)

########################
### PARALLEL GROBNER ###
########################
#
#@profile
def grobnerPar(R, alphabet, k, procs=1):
  print('******************************************')
  print('Parallel Grobner:', 'R', R, 'alphabet', alphabet, 'k', k)
  start = time.clock()
  res = check_resolving_grobner(R, k, len(alphabet), alphabet=alphabet, procs=procs)
  elapsed = time.clock() - start
  print('elapsed time', elapsed)
  return (res, [])

############
### MAIN ###
############
#
if __name__=='__main__':
  parser = argparse.ArgumentParser(description='Check Resolvability Timing')
  parser.add_argument('--i', type=str, default='check_res_sets.tsv', required=False,
                      help='file from which to read potential resolving sets')
  parser.add_argument('--o', type=str, default='check_res_results.tsv', required=False,
                      help='file where results will be printed')
  parser.add_argument('--method', type=int, default=1, required=False,
                      help='method to use for testing resolvability (1) parallel Grobner basis (2) brute force')
  parser.add_argument('--dictSize', type=int, default=100000, required=False,
                      help='the size of each dictionary to check if using brute force')
  parser.add_argument('--chunkSize', type=int, default=100000, required=False,
                      help='the size of chunks to check against each dictionary if using brute force')
  parser.add_argument('--procs', type=int, default=1, required=False,
                      help='number of processes to use')
  parser.add_argument('--test', action='store_true',
                      help='flag if set test against information in res_set_data_3.tsv')
  args = parser.parse_args()

  if args.test:
    print('TESTING')
    with open('res_set_data_3.tsv', 'r') as f:
      f.readline() #header
      i = 0
      for line in f:
        i += 1
        l = line.strip().split('\t')
        k = int(l[0])
        a = int(l[1])
        alphabet = list(map(str, range(a)))
        R = l[2].split(',')
        isResolving = (l[3]=='True')
        if a*k<=25:
          (res, R) = parallelBrute(R, alphabet, k, dictSize=args.dictSize, chunkSize=args.chunkSize, procs=args.procs) if args.method==2 else grobnerPar(R, alphabet, k, procs=args.procs)
          if res!=isResolving:
            print('FAIL', res, isResolving)
            print(k, a, R)
            exit(0)
          else: print('SUCCESS', k, a, R, isResolving, 'line', i)
    exit(0)

  print('READ TSV')
  RSets = readTSV(args.i)

  results = []
  for (R, alphabet, k) in RSets:
    if args.method==1: results.append((R, alphabet, k, grobnerPar(R, alphabet, k, procs=args.procs)))
    elif args.method==2: results.append((R, alphabet, k, parallelBrute(R, alphabet, k, dictSize=args.dictSize, chunkSize=args.chunkSize, procs=args.procs)))

  writeTSV(results, args.o)

  print('DONE')
  
  
