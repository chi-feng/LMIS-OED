import numpy as np
from numpy import *
import os, sys
from common import *
from plotting import *

if sys.argv[1] == '2':
  dim, sigma = 2, 0.1

if sys.argv[1] == '4':
  dim, sigma = 4, 0.4

design = 0.8
algorithms =  ['prior', 'is', 'mis']
grid = array([50, 77, 118, 181, 277, 425, 652, 1000])
xmin, xmax = 40, 1100
exact = eig_marginal(dim, design, sigma)

for algo in algorithms:
  trends[algo] = {}
  for trend in ['mean','mean_err','bias','bias_err','var','var_err','mse','mse_err']:
    trends[algo][trend] = zeros((len(grid), len(grid)))
  for i, N in enumerate(grid):
    for j, M in enumerate(grid):
      filename = 'out/mesh/mesh%dd_%s_N%04d_M%04d' % (dim, algo, N, M)
      data = get_eig(filename)
      trends[algo]['mean'][i,j], trends[algo]['mean_err'][i,j] = mean(data)
      trends[algo]['var'][i,j], trends[algo]['var_err'][i,j]   = var(data)
      trends[algo]['mse'][i,j], trends[algo]['mse_err'][i,j]   = mse(data, exact)
  trends[algo]['bias'] = trends[algo]['mean'] - exact
  trends[algo]['bias_err'] = trends[algo]['mean_err']

def contourplot(algo, trend, log10min, log10max, overlay=False):
  print 'plotting', (algo, trend)
  N, M = meshgrid(grid, grid)
  values = trends[algo][trend].T;
  fig = plt.figure()
  fig.set_size_inches(5, 4)
  ax = plt.gca()
  ax.set_xscale('log')
  ax.set_yscale('log')
  powers = range(int(ceil(log10min)), int(ceil(log10max)))
  levels = np.logspace(log10min, log10max, 21)
  ticks = [10**power for power in powers]
  ticklabels = ['$10^{%d}$' % power for power in powers]
  cs = plt.contourf(N, M, values, levels=levels, norm=matplotlib.colors.LogNorm(), cmap=matplotlib.cm.PuBu_r)
  cb = plt.colorbar(ticks=ticks)
  cb.ax.set_yticklabels(ticklabels)
  if overlay:
    plt.plot([50, 1e4/100],[1e4/100, 50], 'k-', alpha=0.5)
    plt.plot([50, 1e5/100],[1e5/100, 50], 'k-', alpha=0.5)
    plt.text(75, 75, '$10^4$ evals', size=9, rotation=-45, ha="center", va="center")
    plt.text(240, 240, '$10^5$ evals', size=9, rotation=-45, ha="center", va="center")
  plt.xlabel(r'$N$ (outer)')
  plt.ylabel(r'$M_1=M_2$ (inner)')
  save_fig('plot/mesh%dd_%s_%s.pdf' % (dim, algo, trend))

options = []
for algo in algorithms:
  options.append((algo, 'bias', -3, 1, True))
  options.append((algo, 'var', -3, 1, True))
  options.append((algo, 'mse', -4, 1.5, True))

def dowork(option):
  contourplot(*option)

import multiprocessing as mp
from multiprocessing.pool import Pool
pool = Pool(processes=mp.cpu_count())
result = pool.map(dowork, options)
