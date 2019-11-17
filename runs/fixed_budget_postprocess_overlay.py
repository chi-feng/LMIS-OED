import numpy as np
from numpy import *
import os, sys
from common import *
from plotting import *

pd = { }

dim, sigma = 4, 0.4
design = 0.8
suffix = '_3'

budgets = np.logspace(3, 6.5, 8)
budgets = budgets[:-1]

budgets4d = budgets

base_ratio = 100

def get_samples(budget, ratio):
  N = 0.25 * (-ratio + np.sqrt(ratio) * np.sqrt(ratio + 8 * budget))
  M = (-ratio + np.sqrt(ratio) * np.sqrt(ratio + 8 * budget)) / (4 * ratio)
  if (budget < 10**3):
    return int(N+0.5), int(M+0.5)
  else:
    return int(N), int(M)

xmin, xmax = 10**(2.75), 10**(6.25)


algos =  ['prior', 'mis','exact']

trends = {
  'prior':{'name':'Prior'},
  'mis':{'name':'LMIS'},
  'exact':{'name':'Exact'}
}

colors = {'prior':'#e41a1c', 'is':'#377eb8', 'mis':'#4daf4a', 'exact':'#984ea3'} # from colorbrewer

exact = eig_marginal(dim, design, sigma)

for i, algo in enumerate(algos):

  trends[algo] = {'name':trends[algo]['name'],
                  'mean': np.zeros(len(budgets)), 'mean_err': np.zeros(len(budgets)),
                  'var': np.zeros(len(budgets)), 'var_err': np.zeros(len(budgets)),
                  'mse': np.zeros(len(budgets)),  'mse_err': np.zeros(len(budgets))}
  for j, budget in enumerate(budgets):
    N, M = get_samples(budget, base_ratio if not algo == 'prior' else 1.0 / base_ratio)
    filename = 'out/fixed_budget/%dd_%s_N%04d_M%04d%s' % (dim, algo, N, M, suffix)
    data = get_eig(filename)
    trends[algo]['mean'][j], trends[algo]['mean_err'][j] = mean(data)
    trends[algo]['var'][j], trends[algo]['var_err'][j] = var(data)
    trends[algo]['mse'][j], trends[algo]['mse_err'][j] = mse(data, exact)
    if algo == 'exact':
      trends[algo]['mse'][j] = trends[algo]['var'][j]
      trends[algo]['mse_err'][j] = trends[algo]['var_err'][j]

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

fig = plt.figure(1)
fig.set_size_inches(3.5, 4.5)
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log', nonposy='clip')
ymin, ymax = 1e8, 0
for i, algo in enumerate(algos):
  x_values = []
  for budget in budgets:
    N, M = get_samples(budget, base_ratio if not algo == 'prior' else 1.0 / base_ratio)
    x_values.append(N + N * (M + M));
  mmse = trends[algo]['mse']
  mse_err = trends[algo]['mse_err']
  # plt.errorbar(x_values, mmse, mse_err, lw=1, fmt='-', color=colors[algo], elinewidth=0.4, capthick=0.4, capsize=2, barsabove=True)
  pd[(dim, algo, 'evals')] = x_values;
  pd[(dim, algo, 'mse')] = mmse;
  pd[(dim, algo, 'mse_err')] = mse_err;

  ymin = min(ymin, np.min(mmse))
  ymax = max(ymax, np.max(mmse))
ymin = 10**np.floor(np.log10(ymin) - 0.1)
ymax = 10**np.ceil(np.log10(ymax) + 0.1)
plot_loglines(log10(xmin), log10(xmax), log10(ymin), log10(ymax))
plt.xlim([xmin, xmax])
plt.ylim([ymin, ymax])
plt.ylabel(r'Estimator MSE')
plt.xlabel(r'Model evaluations')
make_legend()
#save_fig('plot/fixed_budget_%dd_mse%s.pdf' % (dim, suffix))

fig = plt.figure(2)
fig.set_size_inches(3.5, 4.5)
plt.subplot(211)
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log', nonposy='clip')
ymin, ymax = 1e8, 0
for i, algo in enumerate(algos):
  if algo == 'exact':
    continue
  x_values = []
  for budget in budgets:
    N, M = get_samples(budget, base_ratio if not algo == 'prior' else 1.0 / base_ratio)
    x_values.append(N + N * (M + M));
  bias = (trends[algo]['mean'] - exact)
  print x_values
  print bias
  print  2 * (np.array(trends[algo]['mean_err']))
  # plt.errorbar(x_values, bias, 2 * (np.array(trends[algo]['mean_err'])), lw=1, fmt='-', color=colors[algo], elinewidth=0.4, capthick=0.4, capsize=2, barsabove=True)
  pd[(dim, algo, 'bias')] = bias;
  pd[(dim, algo, 'bias_err')] = 2 * (np.array(trends[algo]['mean_err']));
  ymin = min(ymin, np.min(bias))
  ymax = max(ymax, np.max(bias))
ymin = 10**np.floor(np.log10(ymin) - 0.1)
ymax = 10**np.ceil(np.log10(ymax) + 0.1)
plot_loglines(log10(xmin), log10(xmax), log10(ymin), log10(ymax))
plt.xlim([xmin, xmax])
plt.ylim([ymin, ymax])
plt.ylabel(r'Estimator bias')
# plt.xlabel(r'Model evaluations')
# make_legend()
plt.subplot(212)
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log', nonposy='clip')
ymin, ymax = 1e8, 0
for i, algo in enumerate(algos):
  variance = trends[algo]['var']
  # plt.errorbar(x_values, variance, 2 * (np.array(trends[algo]['var_err'])), lw=1, fmt='-', color=colors[algo], elinewidth=0.4, capthick=0.4, capsize=2, barsabove=True)
  pd[(dim, algo, 'var')] = variance;
  pd[(dim, algo, 'var_err')] = 2 * (np.array(trends[algo]['var_err']));
  ymin = min(ymin, np.min(variance))
  ymax = max(ymax, np.max(variance))
ymin = 10**np.floor(np.log10(ymin) - 0.1)
ymax = 10**np.ceil(np.log10(ymax) + 0.1)
plot_loglines(log10(xmin), log10(xmax), log10(ymin), log10(ymax))
plt.xlim([xmin, xmax])
plt.ylim([ymin, ymax])
plt.ylabel(r'Estimator variance')
plt.xlabel(r'Model evaluations')
#save_fig('plot/fixed_budget_%dd_error%s.pdf' % (dim, suffix))



budgets = np.array([
  1000.            ,
  3162.27766017   ,
  4500,
  7000,
  10000.         ,
  12000.         ,
  15000.         ,
  20000,
  31622.77660168,
  100000.           ,
  316227.76601684  ,
  1000000.        ])

budgets8d = budgets

suffix = '_3'

dim = 8

design = 0.8
exact = eig_marginal(dim, design, sigma)

for i, algo in enumerate(algos):

  trends[algo] = {'name':trends[algo]['name'],
                  'mean': np.zeros(len(budgets)), 'mean_err': np.zeros(len(budgets)),
                  'var': np.zeros(len(budgets)), 'var_err': np.zeros(len(budgets)),
                  'mse': np.zeros(len(budgets)),  'mse_err': np.zeros(len(budgets))}
  for j, budget in enumerate(budgets):
    N, M = get_samples(budget, base_ratio if not algo == 'prior' else 1.0 / base_ratio)
    filename = 'out/fixed_budget/%dd_%s_N%04d_M%04d%s' % (dim, algo, N, M, suffix)
    data = get_eig(filename)
    trends[algo]['mean'][j], trends[algo]['mean_err'][j] = mean(data)
    trends[algo]['var'][j], trends[algo]['var_err'][j] = var(data)
    trends[algo]['mse'][j], trends[algo]['mse_err'][j] = mse(data, exact)
    if algo == 'exact':
      trends[algo]['mse'][j] = trends[algo]['var'][j]
      trends[algo]['mse_err'][j] = trends[algo]['var_err'][j]

fig = plt.figure(1)
fig.set_size_inches(3.5, 4.5)
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log', nonposy='clip')
ymin, ymax = 1e8, 0
for i, algo in enumerate(algos):
  x_values = []
  for budget in budgets:
    N, M = get_samples(budget, base_ratio if not algo == 'prior' else 1.0 / base_ratio)
    x_values.append(N + N * (M + M));
  mmse = trends[algo]['mse']
  mse_err = trends[algo]['mse_err']
  # plt.errorbar(x_values, mmse, mse_err, lw=1, fmt='--', color=colors[algo], elinewidth=0.4, capthick=0.4, capsize=2, barsabove=True)
  pd[(dim, algo, 'evals')] = x_values;
  pd[(dim, algo, 'mse')] = mmse;
  pd[(dim, algo, 'mse_err')] = mse_err;

  ymin = min(ymin, np.min(mmse))
  ymax = max(ymax, np.max(mmse))
ymin = 10**np.floor(np.log10(ymin) - 0.1)
ymax = 10**np.ceil(np.log10(ymax) + 0.1)
# plot_loglines(log10(xmin), log10(xmax), log10(ymin), log10(ymax))
plt.xlim([xmin, xmax])
plt.ylim([ymin, ymax])
plt.ylabel(r'Estimator MSE')
plt.xlabel(r'Model evaluations')
make_legend()
#save_fig('plot/fixed_budget_mse.pdf')

fig = plt.figure(2)
fig.set_size_inches(3.5, 4.5)
plt.subplot(211)
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log', nonposy='clip')
ymin, ymax = 1e8, 0
for i, algo in enumerate(algos):
  if algo == 'exact':
    continue
  x_values = []
  for budget in budgets:
    N, M = get_samples(budget, base_ratio if not algo == 'prior' else 1.0 / base_ratio)
    x_values.append(N + N * (M + M));
  bias = (trends[algo]['mean'] - exact)
  print x_values
  print bias
  print  2 * (np.array(trends[algo]['mean_err']))
  # plt.errorbar(x_values, bias, 2 * (np.array(trends[algo]['mean_err'])), lw=1, fmt='--', color=colors[algo], elinewidth=0.4, capthick=0.4, capsize=2, barsabove=True)
  pd[(dim, algo, 'bias')] = bias;
  pd[(dim, algo, 'bias_err')] = 2 * (np.array(trends[algo]['mean_err']));
  ymin = min(ymin, np.min(bias))
  ymax = max(ymax, np.max(bias))
ymin = 10**np.floor(np.log10(ymin) - 0.1)
ymax = 10**np.ceil(np.log10(ymax) + 0.1)
# plot_loglines(log10(xmin), log10(xmax), log10(ymin), log10(ymax))
plt.xlim([xmin, xmax])
plt.ylim([ymin, ymax])
plt.ylabel(r'Estimator bias')
# plt.xlabel(r'Model evaluations')
# make_legend()

plt.subplot(212)
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log', nonposy='clip')
ymin, ymax = 1e8, 0
for i, algo in enumerate(algos):
  variance = trends[algo]['var']
  # plt.errorbar(x_values, variance, 2 * (np.array(trends[algo]['var_err'])), lw=1, fmt='--', color=colors[algo], elinewidth=0.4, capthick=0.4, capsize=2, barsabove=True)
  pd[(dim, algo, 'var')] = variance;
  pd[(dim, algo, 'var_err')] = 2 * (np.array(trends[algo]['var_err']));
  ymin = min(ymin, np.min(variance))
  ymax = max(ymax, np.max(variance))
ymin = 10**np.floor(np.log10(ymin) - 0.1)
ymax = 10**np.ceil(np.log10(ymax) + 0.1)
# plot_loglines(log10(xmin), log10(xmax), log10(ymin), log10(ymax))
plt.xlim([xmin, xmax])
plt.ylim([ymin, ymax])
plt.ylabel(r'Estimator variance')
plt.xlabel(r'Model evaluations')
#save_fig('plot/fixed_budget_error.pdf')


for dim in [4,8]:
  columns = []
  mat = np.zeros((len(pd[(dim,'prior','evals')]),7*len(algos)))
  for i, algo in enumerate(algos):
    columns = columns + [
      '%s%ddevals' % (algo,dim),
      '%s%ddmse' % (algo,dim),
      '%s%ddmseerr' % (algo,dim),
      '%s%ddbias' % (algo,dim),
      '%s%ddbiaserr' % (algo,dim),
      '%s%ddvar' % (algo,dim),
      '%s%ddvarerr' % (algo,dim),
    ]
    mat[:,0+i*7] = pd[(dim,algo,'evals')]
    mat[:,1+i*7] = pd[(dim,algo,'mse')]
    mat[:,2+i*7] = pd[(dim,algo,'mse_err')]
    if algo != 'exact':
      mat[:,3+i*7] = np.abs(pd[(dim,algo,'bias')])
      mat[:,4+i*7] = pd[(dim,algo,'bias_err')]
    mat[:,5+i*7] = pd[(dim,algo,'var')]
    mat[:,6+i*7] = pd[(dim,algo,'var_err')]
  np.savetxt("plot_data/fixed_budget_%dd.csv" % dim, mat, delimiter=",", header=','.join(columns),comments='')