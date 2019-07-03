###############################################################################################
# Python code to run Grobner basis resolving set timing tests in parallel
# Author: Lucas Laird, Richard Carter Tillquist
# Last Updated: June 2019
###############################################################################################

from hammingGrobner import check_resolving_grobner, brute_force_resolve
import multiprocessing as mp
from itertools import product
import numpy as np
import matplotlib.pyplot as plt
import argparse
import glob
import pickle
import time
import os

def writeDict(d, outFile):
  with open(outFile, 'wb') as o:
    pickle.dump(d, o, protocol=pickle.HIGHEST_PROTOCOL)

def readDict(inFile):
  d = {}
  with open(inFile, 'rb') as f:
    d = pickle.load(f)
  return d

def readTSV(inFile):
  L = []
  with open(inFile, 'r') as f:
    f.readline() #header
    for line in f:
      l = line.strip().split('\t')
      k = int(l[0])
      a = int(l[1])
      R = l[2].split(',')
      isResolving = l[3]=='True'
      L.append((k, a, R, isResolving))
  return L

def readFiles(prefix=''):
  data = {}
  if len(glob.glob(prefix+'brute_force_times_all.dict'))>0: data = mergeDicts(data, readDict(prefix+'brute_force_times_all.dict'))
  if len(glob.glob(prefix+'grobner_times_all.dict'))>0: data = mergeDicts(data, readDict(prefix+'grobner_times_all.dict'))
  if len(glob.glob(prefix+'parallel_grobner_times_all.dict'))>0: data = mergeDicts(data, readDict(prefix+'parallel_grobner_times_all.dict'))
  if len(data)>0: return data
  for f in [f for f in glob.glob(prefix+'brute_force_*.dict')+glob.glob(prefix+'grobner_*.dict')+glob.glob(prefix+'parallel_grobner_*.dict') if 'all' not in f]:
    data = mergeDicts(data, readDict(f))
  return data

def runAll(examples, prefix='', intermediate='', repeats=1, procs=1):
  data = readFiles(prefix=prefix)
  data = {1:data.get(1, {}), 2:data.get(2, {}), 3:data.get(3, {})}
  jobs = [(a, k, R, isResolving, funcNum, intermediate) for a in examples for k in examples[a] for (R, isResolving) in examples[a][k] for funcNum in [1,2] for _ in range(repeats-len([1 for (r,_,_) in data[funcNum].get(a, {}).get(k, {}) if r==R]))]
  pool = mp.Pool(processes=procs)
  results = pool.map_async(runJob, jobs)
  results = results.get()
  pool.close()
  pool.join()

  jobs = [(a, k, R, isResolving) for a in examples for k in examples[a] for (R, isResolving) in examples[a][k] for _ in range(repeats-len([1 for (r,_,_) in data[3].get(a, {}).get(k, {}) if r==R]))]
  results += [runGrobnerPar(a, k, R, isResolving, intermediate=prefix+'intermediate_parallel_grobner_'+str(a)+'.tsv', procs=procs) for (a, k, R, isResolving) in jobs]

  for (funcNum, a, k, R, isResolving, elapsed) in results:
    if a not in data[funcNum]: data[funcNum][a] = {}
    if k not in data[funcNum][a]: data[funcNum][a][k] = []
    data[funcNum][a][k].append((R, isResolving, elapsed))

  writeDict({1:data[1]}, prefix+'brute_force_times_all.dict')
  writeDict({2:data[2]}, prefix+'grobner_times_all.dict')
  writeDict({3:data[3]}, prefix+'parallel_grobner_times_all.dict')

  if intermediate and len(glob.glob(prefix+'intermediate_*.tsv'))>0:
    for f in glob.glob(prefix+'intermediate_*.tsv'): os.remove(f)

def runSize(a, funcNum, examples, prefix='', intermediate='', repeats=1, procs=1):
  data = readFiles(prefix=prefix) #funcNum -> a -> k -> (R, isResolving, time)
  data = {funcNum:{a:data.get(funcNum, {a:{}}).get(a, {})}}
  jobs = [(a, k, R, isResolving, funcNum, intermediate) for k in examples[a] for (R, isResolving) in examples[a][k] for _ in range(repeats-len([1 for (r,_,_) in data[funcNum][a].get(k, {}) if r==R]))]

  results = []
  if funcNum==3: results = [runGrobnerPar(a, k, R, isResolving, intermediate=intermediate, procs=procs) for (a, k, R, isResolving, _, _) in jobs]
  else:
    pool = mp.Pool(processes=procs)
    results = pool.map_async(runJob, jobs)
    results = results.get()
    pool.close()
    pool.join()

  for (funcNum, a, k, R, isResolving, elapsed) in results:
    if k not in data[funcNum][a]: data[funcNum][a][k] = []
    data[funcNum][a][k].append((R, isResolving, elapsed))

  outFile = ('brute_force' if funcNum==1 else ('grobner' if funcNum==2 else 'parallel_grobner'))+'_times_'+str(a)+'.dict'
  writeDict({funcNum:data[funcNum]}, prefix+outFile)

  if intermediate and len(glob.glob(prefix+'intermediate_'+('brute_force' if funcNum==1 else ('grobner' if funcNum==2 else 'parallel_grobner'))+'_'+str(a)+'.tsv'))>0:
    os.remove(prefix+'intermediate_'+('brute_force' if funcNum==1 else ('grobner' if funcNum==2 else 'parallel_grobner'))+'_'+str(a)+'.tsv')

def runJob(arg): #((a, k, R, isResolving, funcNum)):
  (a, k, R, isResolving, funcNum, intermediate) = arg
  print('run', 'a', a, 'k', k, 'func', funcNum)
  res = -1
  numRuns = 0
  start = time.clock()
  while (time.clock()-start)<2.: #run for at least 2 seconds
    res = brute_force_resolve(R, k, a) if funcNum==1 else check_resolving_grobner(R, k, a)
    numRuns += 1
  elapsed = (time.clock() - start) / float(numRuns)
  if res!=isResolving:
    print('res != isResolving')
    elapsed = 'FAIL'
  if intermediate:
    with open(intermediate, 'a+') as o:
      o.write('\t'.join(map(str, [funcNum, a, k, ','.join(R), isResolving, elapsed]))+'\n')
  return (funcNum, a, k, R, isResolving, elapsed)

def runGrobnerPar(a, k, R, isResolving, intermediate='', procs=1):
  print('run', 'a', a, 'k', k, 'parallel grobner procs', procs)
  res = -1
  numRuns = 0
  start = time.clock()
  while (time.clock()-start)<2.:
    res = check_resolving_grobner(R, k, a, procs=procs)
    numRuns += 1
  elapsed = (time.clock() - start) / float(numRuns)
  if res!=isResolving:
    print('res != isResolving')
    elapsed = 'FAIL'
  if intermediate:
    with open(intermediate, 'a+') as o:
      o.write('\t'.join(map(str, [3, a, k, ','.join(R), isResolving, elapsed]))+'\n')
  return (3, a, k, R, isResolving, elapsed)

def combineFiles(prefix=''):
  data = {}
  for f in [f for f in glob.glob(prefix+'brute_force_*.dict')+glob.glob(prefix+'grobner_*.dict')+glob.glob(prefix+'parallel_grobner_*.dict') if '_all' not in f]:
    b = readDict(f)
    data = mergeDicts(data, b)

  #only when the *_all.dict was not updated at the end... how to check this? delete intermediate files when run completes
  for f in glob.glob(prefix+'intermediate_*'):
    with open(f, 'r') as g:
      for line in g:
        l = line.strip().split('\t')
        funcNum = int(l[0])
        a = int(l[1])
        k = int(l[2])
        R = l[3].split(',')
        isResolving = (l[4]=='True')
        elapsed = float(l[5])
        if funcNum not in data: data[funcNum] = {}
        if a not in data[funcNum]: data[funcNum][a] = {}
        if k not in data[funcNum][a]: data[funcNum][a][k] = []
        data[funcNum][a][k].append((R, isResolving, elapsed))
    os.remove(f)

  if 1 in data: writeDict({1:data[1]}, prefix+'brute_force_times_all.dict')
  if 2 in data: writeDict({2:data[2]}, prefix+'grobner_times_all.dict')
  if 3 in data: writeDict({3:data[3]}, prefix+'parallel_grobner_times_all.dict')

def mergeDicts(x,y):
  for funcNum in y:
    if funcNum not in x: x[funcNum] = {}
    for a in y[funcNum]:
      if a not in x[funcNum]: x[funcNum][a] = {}
      for k in y[funcNum][a]:
        if k not in x[funcNum][a]: x[funcNum][a][k] = []
        x[funcNum][a][k] += y[funcNum][a][k]
  return x

#show how many repeats have been collected
def dictTests(examples, repeats, maxAK=-1, prefix=''):
  data = readFiles(prefix=prefix)

  funcNames = {1:'Brute Force', 2:'Grobner', 3:'Parallel Grobner'}
  for funcNum in sorted(data.keys()):
    print(str(funcNum)+': '+funcNames[funcNum])
    for a in sorted(data[funcNum].keys()):
      for k in sorted(data[funcNum][a].keys()):
        if a*k<=maxAK:
          (res, nonRes) = (0, 0)
          for (R, isResolving, t) in data[funcNum][a][k]:
            if isResolving: res += 1
            else: nonRes += 1
          exRes = len(list(filter(lambda x: x[1]==True, examples[a][k])))*repeats
          exNonRes = len(list(filter(lambda x: x[1]==False, examples[a][k])))*repeats
        print('   a: '+str(a)+', k: '+str(k)+', size: '+str(a*k)+', resolving: '+str(res)+'/'+str(exRes)+', non-resolving: '+str(nonRes)+'/'+str(exNonRes)+', total: '+str(res+nonRes)+'/'+str(exRes+exNonRes)+('*****' if (res+nonRes!=exRes+exNonRes) else ''))

############
### MAIN ###
############
#
if __name__=='__main__':
  parser = argparse.ArgumentParser(description='Check Resolvability Timing')
  parser.add_argument('--a', type=int, default=2, required=False,
                      help='alphabet size (assumed to be <= 10)')
  parser.add_argument('--funcNum', type=int, default=1, required=False,
                      help='method to test. (1) for brute force (2) for grobner basis approach (3) for parallel Grobner')
  parser.add_argument('--data', type=str, default='res_set_data_3.tsv', required=False,
                      help='a file containing test instances')
  parser.add_argument('--procs', type=int, default=1, required=False,
                      help='number of processes to use')
  parser.add_argument('--repeats', type=int, default=1, required=False,
                      help='number of times to run each example')
  parser.add_argument('--all', action='store_true',
                      help='flag if set runs all combinations in a single job')
  parser.add_argument('--combine', action='store_true',
                      help='flag if set combine files')
  parser.add_argument('--test', action='store_true',
                      help='flag if set tet number of repeats recorded in full dicts')
  parser.add_argument('--maxAK', type=int, default=25, required=False,
                      help='the maximum a*k values to examine from the given data set')
  parser.add_argument('--intermediate', action='store_true',
                      help='flag if set save runs of each process. note that some writes on linux on atomic (not guaranteed here)')
  parser.add_argument('--prefix', type=str, default='', required=False,
                      help='a prefix to include on saved files')
  args = parser.parse_args()

  if args.a>10:
    print('a value too large, must be <=10. a:'+str(args.a))
    exit(0)

  print('READ TSV FILES')
  files = [args.data]
  intermediate = ''
  if args.intermediate: intermediate = args.prefix+'intermediate_'+('brute_force' if args.funcNum==1 else ('grobner' if args.funcNum==2 else 'parallel_grobner'))+'_'+str(args.a)+'.tsv'
  examples = {}
  for f in files:
    L = readTSV(f)
    for (k, a, R, isResolving) in L:
      if a*k<=args.maxAK:
        if a not in examples: examples[a] = {}
        if k not in examples[a]: examples[a][k] = []
        examples[a][k].append((R, isResolving))

  print('START RUNS')
  if args.combine: combineFiles(prefix=args.prefix)
  elif args.test: dictTests(examples, args.repeats, maxAK=args.maxAK, prefix=args.prefix)
  elif args.all: runAll(examples, prefix=args.prefix, intermediate=intermediate, repeats=args.repeats, procs=args.procs)
  else: runSize(args.a, args.funcNum, examples, prefix=args.prefix, intermediate=intermediate, repeats=args.repeats, procs=args.procs)

  print('DONE')
  
  
