from multiprocessing.pool import Pool
import subprocess
import numpy as np
import os, sys
from common import *
from plotting import *

#if sys.argv[1] == '2':
#  dim, sigma = 2, 0.1

#if sys.argv[1] == '4':
dim, sigma = 4, 0.4

npts = 51
designs = np.linspace(0, 1, npts)
trials = 2000

N = 1000
M1 = 100
M2 = 100

algos = ['prior', 'mis']

if 'noprior' in sys.argv:
  algos = algos[1:]

trends = {
  'prior':{'name':'Prior'},
  'is':{'name':'IS'},
  'mis':{'name':'MIS'},
  'post':{'name':'Posterior'}
}

isjoint = '_joint' if 'joint' in sys.argv else ''
for i, algo in enumerate(algos):
  print i, algo
  trends[algo]['mean'], trends[algo]['mean_err'] = np.zeros(npts), np.zeros(npts)
  trends[algo]['var'], trends[algo]['var_err'] = np.zeros(npts), np.zeros(npts)
  trends[algo]['lo'], trends[algo]['hi'] = np.zeros(npts), np.zeros(npts)
  for j in xrange(npts):
    data = get_eig('out/sweep/sweep%dd_%s%s_%03d_of_%03d' % (dim, isjoint, algo, j, npts))
    trends[algo]['mean'][j], trends[algo]['mean_err'][j] = mean(np.array(data))
    trends[algo]['var'][j], trends[algo]['var_err'][j] = var(np.array(data))
    trends[algo]['lo'][j] = np.percentile(np.array(data), 2.5)
    trends[algo]['hi'][j] = np.percentile(np.array(data), 97.5)

colors = {'prior':'#e41a1c', 'is':'#377eb8', 'mis':'#4daf4a', 'post':'#984ea3', 'exact':'#000000'} # from colorbrewer

styles = {}
for algo in algos:
  styles[algo] = plt.Line2D((0,1),(0,0), color=colors[algo], linestyle='-', linewidth=1)
styles['exact'] = plt.Line2D((0,1),(0,0), color=colors['exact'], linestyle='-', linewidth=0.5)

def make_legend(exact=False):
  labels = [r'%s' % trends[algo]['name'] for algo in algos];
  if exact: labels.append('Exact')
  label_styles = [styles[algo] for algo in algos];
  if exact: label_styles.append(styles['exact'])
  plt.legend(label_styles, labels,
    bbox_to_anchor=(0, 1.05, 1, 1.02), loc=3, ncol=3, mode="expand", borderaxespad=0., numpoints=1, prop={'size':10})
pd = {}
fig = plt.figure()
fig.set_size_inches(5, 3)
for i, algo in enumerate(algos):
  plt.errorbar(designs, trends[algo]['mean'], 2 * np.array(trends[algo]['mean_err']), lw=1, fmt='-', color=colors[algo], elinewidth=0.4, capthick=0.4, capsize=2, barsabove=True)
  y1 = trends[algo]['lo']
  y2 = trends[algo]['hi']
  plt.fill_between(designs,y1,y2,where=y2>=y1,facecolor=colors[algo],lw=0,alpha=0.25,interpolate=True)
  pd[(dim,algo,'u')] = trends[algo]['mean'];
  pd[(dim,algo,'u_err')] = 2 * np.array(trends[algo]['mean_err']);
  pd[(dim,algo,'u_lo')] = y1
  pd[(dim,algo,'u_hi')] = y2
exact = eig_marginal(dim, designs, sigma)
if 'joint' in sys.argv:
  exact = eig_joint(designs)
plt.plot(designs, exact, lw=0.5, color=colors['exact'])
plt.ylabel(r'Estimator mean')
plt.xlabel(r'Design $d$')
make_legend(exact=True)
plt.xlim([0,1])
# plt.ylim([0,4])
save_fig('plot/sweep%dd%s_mean_%s.pdf' % (dim, isjoint, 'noprior' if 'noprior' in sys.argv else 'prior'))


columns = ['designs','exact']
mat = np.zeros((len(designs),4*len(algos)+2))
mat[:,0] = designs;
mat[:,1] = eig_marginal(dim, designs, sigma);
if 'joint' in sys.argv:
  mat[:,1] = eig_joint(designs);
for i, algo in enumerate(algos):
  columns = columns + [
    '%s%dd' % (algo,dim),
    '%s%dderr' % (algo,dim),
    '%s%ddlo' % (algo,dim),
    '%s%ddhi' % (algo,dim),
  ]
  mat[:,0+i*4+2] = pd[(dim,algo,'u')]
  mat[:,1+i*4+2] = pd[(dim,algo,'u_err')]
  mat[:,2+i*4+2] = pd[(dim,algo,'u_lo')]
  mat[:,3+i*4+2] = pd[(dim,algo,'u_hi')]
np.savetxt("plot_data/sweep%dd%s.csv" % (dim, isjoint), mat, delimiter=",", header=','.join(columns),comments='')
