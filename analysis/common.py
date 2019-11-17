# -*- coding: utf-8 -*-
# encoding: utf-8
import multiprocessing
from multiprocessing.pool import Pool
import subprocess
import numpy as np
import random
import os, sys, time, traceback
import sqlite3
import pprint

pp = pprint.PrettyPrinter(indent=2)

executable = './build/Driver';
if not os.path.isfile(executable[2:]):
  print('could not find %s' % executable)
  # exit()

database = 'default.sqlite';
# database = 'mossbauer_sweep3d.sqlite'

if 'database' in sys.argv:
  database = sys.argv[sys.argv.index('database') + 1]

mossbauer_defaults = {
  '-experiment': 'mossbauer',
  "-dim":3,
  "-sigeps": 0.1,
  "-poi":1,
  "-design": '-1.0,0.0,1.0',
  "-seed":0,
  "-dumpfile":"none",
  "-outfile":"none",
  "-N":100,
  "-M1":100,
  "-M2":100,
  "-useIS":0,
  "-useMIS":0,
  "-useMarginal":0,
  "-useReverseLikelihood":1,
  "-useMinSampleDistance":0,
  "-sort":1,
  "-useExactPosterior":0,
  "-maxComponents":100,
  "-biasingDistributionType":"MVT",
  "-dof":'2.5',
  "-nugget":'0.001',
  "-simulate":0,
  "-tag":''}

linear_defaults = {
  '-experiment': 'linear',
  "-dim":4,
  "-sigeps": '0.4',
  "-poi":1,
  "-design": '0.8',
  "-index": 0,
  "-seed":0,
  "-dumpfile":"none",
  "-outfile":"none",
  "-N":0,
  "-M1":0,
  "-M2":0,
  "-useIS":0,
  "-useMIS":0,
  "-useMarginal":0,
  "-useReverseLikelihood":1,
  "-useMinSampleDistance":0,
  "-sort":1,
  "-useExactPosterior":0,
  "-maxComponents":0,
  "-biasingDistributionType":"MVT",
  "-dof":'2.5',
  "-nugget":'0.001',
  "-simulate":0,
  "-tag":''}

columns = [
  {'name':'experiment','type':'text'},
  {'name':'dim','type':'int'},
  {'name':'sigeps','type':'real'},
  {'name':'poi','type':'int'},
  {'name':'index','type':'int'},
  {'name':'design','type':'text'},
  {'name':'N','type':'int'},
  {'name':'M1','type':'int'},
  {'name':'M2','type':'int'},
  {'name':'useMIS','type':'int'},
  {'name':'biasingDistributionType','type':'text'},
  {'name':'dof','type':'real'},
  {'name':'nugget','type':'real'},
  {'name':'useMarginal','type':'int'},
  {'name':'tag','type':'text'},
]

conn = sqlite3.connect(database)
cursor = conn.cursor()
colspec = ['id integer PRIMARY KEY']
colspec += ['`%s` %s' % (col['name'], col['type']) for col in columns];
colspec += ['eig','timestamp real','runtime real','command']
sql = 'CREATE TABLE IF NOT EXISTS runs (%s);' % (','.join(colspec))
cursor.execute(sql);
conn.commit()
conn.close()

def getMossbauerNM(algo, evaluations):
  if algo == 'mis' or algo == 'mis_marg':
    M_values = map(int, np.floor(np.sqrt(evaluations / 10) / 2 + 0.5))
    N_values = map(int, np.floor(np.sqrt(evaluations / 10) * 10 + 0.5))
    return (N_values, M_values)
  elif algo == 'prior':
    N_values = map(int, np.floor(np.sqrt(evaluations / 20) + 0.5))
    M_values = map(int, 10 * np.sqrt(evaluations / 20) + 0.5);
    return (N_values, M_values)
  else:
    print('unknown algo %s' % algo)
    exit()

def get_linear_NM(algo, evaluations):
  ratio = 100;
  def get_samples(budget, ratio):
    N = 0.25 * (-ratio + np.sqrt(ratio) * np.sqrt(ratio + 8 * budget))
    M = (-ratio + np.sqrt(ratio) * np.sqrt(ratio + 8 * budget)) / (4 * ratio)
    if (budget < 10**3):
      return int(N+0.5), int(M+0.5)
    else:
      return int(N), int(M)
  if algo == 'mis' or algo == 'mis_marg' or algo == 'exact':
    N = [get_samples(evals, ratio)[0] for evals in evaluations]
    M = [get_samples(evals, ratio)[1] for evals in evaluations]
    return (N,M)
  elif algo == 'prior':
    N = [get_samples(evals, 1.0/ratio)[0] for evals in evaluations]
    M = [get_samples(evals, 1.0/ratio)[1] for evals in evaluations]
    return (N,M)


def format_options(options):
  topsep = '┌' + ''.join(['─']*27) + '┬' + ''.join(['─']*14) + '┐' + '\n'
  bottomsep = '└' + ''.join(['─']*27) + '┴' + ''.join(['─']*14) + '┘'
  entries = [ ]
  for k, v in options.iteritems():
    entries.append('│ %25s │ %12s │\n' % (k, v));
  return topsep + ''.join(entries) + bottomsep;

def get_eig(cursor, options):
  where = ["`%s`='%s'" % (key,value) for key, value in options.iteritems()]
  sql = 'SELECT * FROM runs WHERE %s' % ' AND '.join(where);
  cursor.execute(sql);
  values = []
  runtimes = []
  for row in cursor:
    values.append(float(row['eig']))
    runtimes.append(float(row['runtime']))
  print(format_options(options))
  print('got %d entries. average runtime %0.1f' % (len(values), np.mean(np.array(runtimes))))
  if len(values) > 0:
    print('%0.4f +- %0.4f [%0.4f, %0.4f]\n' % (np.mean(values), np.std(values), np.min(values), np.max(values)))
  return np.array(values)

def read_binary(filename):
  f = open(filename, 'rb')
  rows = struct.unpack('i', f.read(4))[0]
  cols = struct.unpack('i', f.read(4))[0]
  m = np.empty([rows, cols])
  for i in xrange(rows):
    for j in xrange(cols):
      m[i, j] = struct.unpack('d', f.read(8))[0]
  return m

def read_mixture_indices(filename):
  f = open(filename, 'r')
  num_indices = []
  for line in f:
    num_indices.append(len(line.split()) -1)
  return np.array(num_indices)


def execute(options):
  try:
    arguments = [executable]
    simulate = False
    for k, v in options.iteritems():
      if k == '-simulate':
        simulate = v
        continue
      arguments.append('%s' % k)
      arguments.append('%s' % v)
    start = time.time();
    # print format_options(options)
    output = subprocess.check_output(arguments)
    runtime = time.time() - start;
    for line in output.split('\n'):
      if line[0:3] == 'EIG':
        eig = float(line[6:]);
        if simulate:
          continue
        # record to database
        conn = sqlite3.connect(database)
        cursor = conn.cursor()
        names = [col['name'] for col in columns]
        values = []
        for col in columns:
          key = '-%s' % col['name']
          values.append(options[key]);
        names += ['eig','timestamp','runtime','command']
        values += [eig,start,runtime,' '.join(list(sys.argv))]
        values = tuple(values);
        placeholder = ','.join(['?'] * len(values));
        sql = 'INSERT INTO runs (%s) VALUES (%s)' % (','.join(['`%s`' % name for name in names]), placeholder);
        cursor.execute(sql, values)
        conn.commit()
        conn.close()
  except Exception as e:
    print(traceback.format_exception(*sys.exc_info()))
    raise # reraises the exception
  return 0

def progressbar(fraction):
  barwidth = 70
  pos = int(barwidth * fraction)
  bar = '╞'
  for i in range(barwidth):
    if i < pos:
      bar += '═'
    elif i == pos:
      bar += '═'
    else:
      bar += ' '
  bar += '╡ '
  return bar

def execute_all(runs):
  print('starting runs')
  nproc = multiprocessing.cpu_count()
  pool = Pool(processes=nproc)
  result = pool.map_async(execute, runs, chunksize=1)
  start_time = time.time()
  while not result.ready():
    completed = len(runs) - result._number_left
    if completed > 0:
      rate = 60.0 * completed / (time.time() - start_time)
      percent = 100.0 * completed / len(runs)
      print('%04d/%04d %4.1f%% %s %4d/min  ETA: %3.1f\n' % (completed, len(runs), percent, progressbar(percent / 100), rate, 1.0 / rate * result._number_left))
    sys.stdout.flush()
    time.sleep(1)
  print('\ncompleted', len(runs) - result._number_left)
  pool.close()

def eig_marginal(dim, d, s=0.4):
  # print 'dim = %d, design = %f, sigma=%f' % (dim, d, s)
  if dim == 2:
    return -0.5 * np.log((s**2 * (26. + 25. * (-2. + d) * d + s**2)) / ((1. + 25. * (-1. + d) * d)**2 + (27. + 50. * (-1. + d) * d) * s**2 + s**4))
  if dim == 4:
    return -0.5 * np.log((s**2 * (52. + 5. * (-14. + 5. * d) * d + s**2)) / ((3. + 5. * (-7. + 5. * d) * d)**2 + 5. * (11. + 2. * d * (-7. + 5. * d)) * s**2 + s**4))
  if dim == 8:
    return -0.5 * np.log((4.0*(3204.0+125.0*d*(-22.0+5.0*d)))/(44141.0+125.0*d*(-11.0+5.0*d)*(358.0+125.0*d*(-11.0+5.0*d))))

def mean(array):
  return np.mean(array), np.std(array) / np.sqrt(array.size)

def var(array):
  return np.var(array), np.var(array) * np.sqrt(2.0 / (array.size - 1))

def mse(array, exact):
  bias, bias_err = mean(array)
  bias = bias - exact
  vari, vari_err = var(array)
  alt = np.sum((array - exact)**2)/len(array);
  return vari + bias**2, np.sqrt((vari_err/vari)**2 + (2 * bias_err/bias)**2) * (vari + bias**2)
