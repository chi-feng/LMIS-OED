import multiprocessing as mp
from multiprocessing.pool import Pool
import subprocess
import random
import numpy as np
import os, sys, time, datetime

executable = './build/Driver'
if not os.path.isfile(executable[2:]):
  print "could not find executable"
  exit()

'''
if sys.argv[1] == '2':
  dim, sigeps = 2, 0.1

if sys.argv[1] == '4':
  dim, sigeps = 4, 0.4
'''

dim, sigeps = 4, 0.4

npts = 51
designs = np.linspace(0, 1, npts)

N = 1000
M1 = 100
M2 = 100

#algorithms = ['prior','is','mis','post']

algorithms = ['prior','mis']
trials = {'prior':1000,'mis':1}

options = {"-dim":dim,
  "-sigeps":sigeps,
  "-poi":1,
  "-design":0,
  "-seed":0,
  "-dumpfile":"none",
  "-outfile":"none",
  "-N":N, "-M1":M1, "-M2":M2,
  "-useIS":0,
  "-useMIS":0,
  "-useReverseLikelihood":1,
  "-useMinSampleDistance":0,
  "-sort":1,
  "-useExactPosterior":0,
  "-maxComponents":min(N,M1),
  "-biasingDistributionType":"MVT",
  "-dof":2.5,
  "-nugget":1e-3}

if 'joint' in sys.argv:
  options['-poi'] = dim;

runs = []

for algo in algorithms:
  options['-useExactPosterior'] = 0
  options['-biasingDistributionType'] = 'MVT'
  options['-N'] = N;
  options['-M1'] = M1;
  options['-M2'] = M2;

  if algo is 'prior':
    options['-useIS'], options['-useMIS'] = 0, 0
    options['-N'] = M1;
    options['-M1'] = N;
    options['-M2'] = N;

  if algo is 'is':
    options['-useIS'], options['-useMIS'] = 1, 0
  if algo is 'mis':
    options['-useIS'], options['-useMIS'] = 0, 1
  if algo is 'post':
    options['-useIS'], options['-useMIS'] = 1, 0
    options['-useExactPosterior'] = 1
    options['-biasingDistributionType'] = 'MVN'

  for i_d, d in enumerate(designs):
    options['-design'] = d
    isjoint = '_joint' if 'joint' in sys.argv else ''
    filename = 'out/sweep/sweep%dd_%s%s_%03d_of_%03d' % (dim, isjoint, algo, i_d, npts)
    arguments = [executable]
    for k, v in options.iteritems():
      arguments.append('%s' % k)
      arguments.append('%s' % v)
    if os.path.isfile(filename):
      print "Warning: file already exists: %s" % filename
    for trial in xrange(trials[algo]):
      runs.append({"id": len(runs), "filename":filename, "arguments":arguments})

random.shuffle(runs)

def execute(run):
  output = subprocess.check_output(run['arguments'])
  f = open(run['filename'], 'a')
  f.write(output)
  f.close()

pool = Pool(processes=mp.cpu_count())
result = pool.map_async(execute, runs, chunksize=1)
start_time = time.time()
while not result.ready():
  completed = len(runs) - result._number_left
  if completed > 0:
    eta = (time.time() - start_time) / completed * result._number_left / 60
    print '%6d left, ETA: %4.1f minutes' % (result._number_left, eta)
  sys.stdout.flush()
  time.sleep(1)

pool.close()
#nice parallel -j 6 --"command1""command2"