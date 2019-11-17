import numpy as np
from numpy import *
import os, sys
from common import *
from plotting import *

pdN = {}
pdM = {}

dim, sigma = 8, 0.4 

N_values = [50, 100, 200, 400, 800, 1600]
M_values = [50, 100, 200, 400, 800, 1600]

N_fixed = 800
M_fixed = 800

xmin, xmax = 40, 2000

design = 0.8

trials = 1000

algos =  ['prior', 'mis', 'exact']

trends = {
  'prior':{'name':'Prior'},
  # 'is':{'name':'LIS'},
  'mis':{'name':'LMIS'},
  'exact':{'name':'Exact'}
}

colors = {'prior':'#e41a1c', 'is':'#377eb8', 'mis':'#4daf4a', 'exact':'#984ea3'} # from colorbrewer

exact = eig_marginal(dim, design, sigma)

for i, algo in enumerate(algos):

  trends[algo]['N'] = {'mean': np.zeros(len(N_values)), 'mean_err': np.zeros(len(N_values)),
                       'var': np.zeros(len(N_values)),  'var_err': np.zeros(len(N_values))}
  for j, N in enumerate(N_values):
    data = get_eig('out/convergence%dd_%s_N%04d_M%04d' % (dim, algo, N, M_fixed))
    trends[algo]['N']['mean'][j], trends[algo]['N']['mean_err'][j] = mean(data)
    trends[algo]['N']['var'][j], trends[algo]['N']['var_err'][j] = var(data)

  trends[algo]['M'] = {'mean': np.zeros(len(M_values)), 'mean_err': np.zeros(len(M_values)),
                       'var': np.zeros(len(M_values)),  'var_err': np.zeros(len(M_values))}
  for j, M in enumerate(M_values):
    data = get_eig('out/convergence%dd_%s_N%04d_M%04d' % (dim, algo, N_fixed, M))
    trends[algo]['M']['mean'][j], trends[algo]['M']['mean_err'][j] = mean(data)
    trends[algo]['M']['var'][j], trends[algo]['M']['var_err'][j] = var(data)

  pdN[(dim,algo,'values')] = N_values;
  pdN[(dim,algo,'bias')] = np.abs(trends[algo]['N']['mean'] - exact)
  pdN[(dim,algo,'bias_err')] = 2 * trends[algo]['N']['mean_err']
  pdN[(dim,algo,'var')] = np.abs(trends[algo]['N']['var'])
  pdN[(dim,algo,'var_err')] = 2 * trends[algo]['N']['var_err']
  pdM[(dim,algo,'values')] = M_values;
  pdM[(dim,algo,'bias')] = np.abs(trends[algo]['M']['mean'] - exact)
  pdM[(dim,algo,'bias_err')] = 2 * trends[algo]['M']['mean_err']
  pdM[(dim,algo,'var')] = np.abs(trends[algo]['M']['var'])
  pdM[(dim,algo,'var_err')] = 2 * trends[algo]['M']['var_err']

styles4d = {}
for algo in algos:
  styles4d[algo] = plt.Line2D((0,1),(0,0), color=colors[algo], linestyle='-', linewidth=1)
styles8d = {}
for algo in algos:
  styles8d[algo] = plt.Line2D((0,1),(0,0), color=colors[algo], linestyle='--', linewidth=1)

def make_legend():
  labels = []
  label_styles = []
  for algo in algos:
    labels.append(r'%s 4D' % trends[algo]['name'])
    labels.append(r'%s 8D' % trends[algo]['name'])
    label_styles.append(styles4d[algo])
    label_styles.append(styles8d[algo])
  plt.legend(label_styles, labels,
    bbox_to_anchor=(0, 1.05, 1, 1.02), loc=3, ncol=3, mode="expand", borderaxespad=0., numpoints=1, prop={'size':8})

for c, index in enumerate(['N', 'M']):
  x_values = N_values if index is 'N' else M_values
  label = 'Outer' if index is 'N' else 'Inner'
  fig = plt.figure(c)
  fig.set_size_inches(4, 5)
  plt.subplot(211)
  ax = plt.gca()
  for i, algo in enumerate(algos):
    if algo == 'exact':
      continue
    bias = trends[algo][index]['mean'] - exact
    print algo, trends[algo][index]['mean'], np.mean(trends[algo][index]['mean'])
    plt.errorbar(x_values, np.abs(bias), 2 * (np.array(trends[algo][index]['mean_err'])), lw=1, fmt='--', color=colors[algo], elinewidth=0.4, capthick=0.4, capsize=2, barsabove=True)
    for j, value in enumerate(bias):
      if value < 0:
        plt.scatter((x_values[j]),(np.abs(value)), s=80, facecolors='none', edgecolors='r')
  plt.subplot(212)
  ax = plt.gca()
  ax.set_xscale('log')
  ax.set_yscale('log', nonposy='clip')
  for i, algo in enumerate(algos):
    variance = trends[algo][index]['var']
    plt.errorbar(x_values, variance, 2 * (np.array(trends[algo][index]['var_err'])), lw=1, fmt='--', color=colors[algo], elinewidth=0.4, capthick=0.4, capsize=2, barsabove=True)

dim, sigma = 4, 0.4 
design = 0.8

exact = eig_marginal(dim, design, sigma)

for i, algo in enumerate(algos):

  trends[algo]['N'] = {'mean': np.zeros(len(N_values)), 'mean_err': np.zeros(len(N_values)),
                       'var': np.zeros(len(N_values)),  'var_err': np.zeros(len(N_values))}
  for j, N in enumerate(N_values):
    data = get_eig('out/convergence%dd_%s_N%04d_M%04d' % (dim, algo, N, M_fixed))
    trends[algo]['N']['mean'][j], trends[algo]['N']['mean_err'][j] = mean(data)
    trends[algo]['N']['var'][j], trends[algo]['N']['var_err'][j] = var(data)

  trends[algo]['M'] = {'mean': np.zeros(len(M_values)), 'mean_err': np.zeros(len(M_values)),
                       'var': np.zeros(len(M_values)),  'var_err': np.zeros(len(M_values))}
  for j, M in enumerate(M_values):
    data = get_eig('out/convergence%dd_%s_N%04d_M%04d' % (dim, algo, N_fixed, M))
    trends[algo]['M']['mean'][j], trends[algo]['M']['mean_err'][j] = mean(data)
    trends[algo]['M']['var'][j], trends[algo]['M']['var_err'][j] = var(data)

  pdN[(dim,algo,'values')] = N_values;
  pdN[(dim,algo,'bias')] = np.abs(trends[algo]['N']['mean'] - exact)
  pdN[(dim,algo,'bias_err')] = 2 * trends[algo]['N']['mean_err']
  pdN[(dim,algo,'var')] = np.abs(trends[algo]['N']['var'])
  pdN[(dim,algo,'var_err')] = 2 * trends[algo]['N']['var_err']
  pdM[(dim,algo,'values')] = M_values;
  pdM[(dim,algo,'bias')] = np.abs(trends[algo]['M']['mean'] - exact)
  pdM[(dim,algo,'bias_err')] = 2 * trends[algo]['M']['mean_err']
  pdM[(dim,algo,'var')] = np.abs(trends[algo]['M']['var'])
  pdM[(dim,algo,'var_err')] = 2 * trends[algo]['M']['var_err']

for c, index in enumerate(['N', 'M']):

  x_values = N_values if index is 'N' else M_values
  label = 'Outer' if index is 'N' else 'Inner'
  fig = plt.figure(c)
  plt.subplot(211)
  ax = plt.gca()
  ax.set_xscale('log')
  ax.set_yscale('log', nonposy='clip')
  ymin, ymax = 1e8, 0
  for i, algo in enumerate(algos):
    if algo == 'exact':
      continue
    bias = trends[algo][index]['mean'] - exact
    print algo, trends[algo][index]['mean'], np.mean(trends[algo][index]['mean'])
    plt.errorbar(x_values, np.abs(bias), 2 * (np.array(trends[algo][index]['mean_err'])), lw=1, fmt='-', color=colors[algo], elinewidth=0.4, capthick=0.4, capsize=2, barsabove=True)
    for j, value in enumerate(bias):
      if value < 0:
        plt.scatter((x_values[j]),(np.abs(value)), s=80, facecolors='none', edgecolors='r')
  if index == 'N':
    ymin = 1e-4
    ymax = 1e0
  else:
    ymin = 10**(-3.5)
    ymax = 1e1
  ymin, ymax = 1e-4, 1e1
  print [ymin,ymax]    
  plot_loglines(log10(xmin), log10(xmax), log10(ymin), log10(ymax))
  plt.xlim([xmin, xmax])
  plt.ylim([ymin, ymax])
  plt.ylabel(r'Estimator bias')
  # plt.xlabel(r'%s Monte Carlo samples $%s$' % (label, index))
  make_legend()
  plt.subplot(212)
  ax = plt.gca()
  ax.set_xscale('log')
  ax.set_yscale('log', nonposy='clip')
  for i, algo in enumerate(algos):
    variance = trends[algo][index]['var']
    plt.errorbar(x_values, variance, 2 * (np.array(trends[algo][index]['var_err'])), lw=1, fmt='-', color=colors[algo], elinewidth=0.4, capthick=0.4, capsize=2, barsabove=True)
  print [ymin,ymax]    
  plot_loglines(log10(xmin), log10(xmax), log10(ymin), log10(ymax))
  plt.xlim([xmin, xmax])
  plt.ylim([ymin, ymax])
  plt.ylabel(r'Estimator variance')
  plt.xlabel(r'%s Monte Carlo samples $%s$' % (label, index))
  filename = 'plot/convergence%s_sweep_overlay.pdf' % (index)
  # save_fig(filename)
  # os.system('gnome-open %s' % filename)


for dim in [4,8]:
  columns = []
  matN = np.zeros((len(pdN[(dim,'prior','values')]),5*len(algos)))
  matM = np.zeros((len(pdM[(dim,'prior','values')]),5*len(algos)))
  for i, algo in enumerate(algos):
    columns = columns + [
      '%s%ddvalues' % (algo,dim),
      '%s%ddbias' % (algo,dim),
      '%s%ddbiaserr' % (algo,dim),
      '%s%ddvar' % (algo,dim),
      '%s%ddvarerr' % (algo,dim),
    ]
    matN[:,0+i*5] = pdN[(dim,algo,'values')]
    if algo != 'exact':
      matN[:,1+i*5] = pdN[(dim,algo,'bias')]
      matN[:,2+i*5] = pdN[(dim,algo,'bias_err')]
    matN[:,3+i*5] = pdN[(dim,algo,'var')]
    matN[:,4+i*5] = pdN[(dim,algo,'var_err')]

    matM[:,0+i*5] = pdM[(dim,algo,'values')]
    if algo != 'exact':
      matM[:,1+i*5] = pdM[(dim,algo,'bias')]
      matM[:,2+i*5] = pdM[(dim,algo,'bias_err')]
    matM[:,3+i*5] = pdM[(dim,algo,'var')]
    matM[:,4+i*5] = pdM[(dim,algo,'var_err')]
  np.savetxt("plot_data/conv%dd_N.csv" % dim, matN, delimiter=",", header=','.join(columns),comments='')
  np.savetxt("plot_data/conv%dd_M.csv" % dim, matM, delimiter=",", header=','.join(columns),comments='')