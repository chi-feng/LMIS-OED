from mossbauer import *
from plotting import *

suffix = 'center';
options['-poi'] = 1
options['-index'] = 0
if 'joint' in sys.argv:
  suffix = 'joint'
  options['-poi'] = 4
if 'height' in sys.argv:
  suffix = 'height'
  options['-index'] = 3

evaluations = np.logspace(3, 5, 5)
M_values = map(int, np.floor(np.sqrt(evaluations / 10) / 2 + 0.5))
N_values = map(int, np.floor(np.sqrt(evaluations / 10) * 10 + 0.5))

xmin, xmax = 10**(2.75), 10**(5.25)

algos =  ['prior', 'mis']

trends = {
  'prior':{'name':'Prior'},
  'mis':{'name':'MIS'},
}

colors = {'prior':'#e41a1c', 'mis':'#4daf4a'}

exact = 1.54

for i, algo in enumerate(algos):

  trends[algo]['N'] = {'mean': np.zeros(len(N_values)), 'mean_err': np.zeros(len(N_values)),
                       'var': np.zeros(len(N_values)),  'var_err': np.zeros(len(N_values)),
                       'mse': np.zeros(len(N_values)),  'mse_err': np.zeros(len(N_values))}
  for j, N in enumerate(N_values):
    data = get_eig('expdesign/out/mossbauer_convergence_%s_%s_%05d_%05d' % (suffix, algo, N, M))
    trends[algo]['N']['mean'][j], trends[algo]['N']['mean_err'][j] = mean(data)
    trends[algo]['N']['var'][j], trends[algo]['N']['var_err'][j] = var(data)
    trends[algo]['N']['mse'][j], trends[algo]['N']['mse_err'][j] = mse(data, exact)

styles = {}
for algo in algos:
  styles[algo] = plt.Line2D((0,1),(0,0), color=colors[algo], linestyle='-', linewidth=1)

def make_legend():
  labels = [r'%s' % trends[algo]['name'] for algo in algos];
  label_styles = [styles[algo] for algo in algos];
  plt.legend(label_styles, labels,
    bbox_to_anchor=(0, 1.05, 1, 1.02), loc=3, ncol=4, mode="expand", borderaxespad=0., numpoints=1, prop={'size':10})

for index in ['N']:

  x_values = evaluations
  fig = plt.figure()
  fig.set_size_inches(3, 4)
  ax = plt.gca()
  ax.set_xscale('log')
  ax.set_yscale('log', nonposy='clip')
  ymin, ymax = 1e8, 0
  for i, algo in enumerate(algos):
    mse = trends[algo][index]['mse']
    mse_err = trends[algo][index]['mse_err']
    rmse = np.sqrt(trends[algo][index]['mse'])
    rmse_err = 0.5 * trends[algo][index]['mse_err'] / rmse
    plt.errorbar(x_values, mse, mse_err, lw=1, fmt='-', color=colors[algo], elinewidth=0.4, capthick=0.4, capsize=2, barsabove=True)
    ymin = min(ymin, np.min(mse))
    ymax = max(ymax, np.max(mse))
  ymin = 10**np.floor(np.log10(ymin) - 0.1)
  ymax = 10**np.ceil(np.log10(ymax) + 0.1)
  plot_loglines(log10(xmin), log10(xmax), log10(ymin), log10(ymax))
  plt.xlim([xmin, xmax])
  plt.ylim([ymin, ymax])
  plt.ylabel(r'Estimator MSE')
  plt.xlabel(r'Model evaluations')
  make_legend()
  save_fig('expdesign/plots/mossbauer_convergence_mse.pdf')

  x_values = evaluations
  fig = plt.figure()
  fig.set_size_inches(3, 4)
  plt.subplot(211)
  ax = plt.gca()
  ax.set_xscale('log')
  ax.set_yscale('log', nonposy='clip')
  ymin, ymax = 1e8, 0
  for i, algo in enumerate(algos):
    if algo == 'exact':
      continue
    bias = np.abs(trends[algo][index]['mean'] - exact)
    print algo, trends[algo][index]['mean'], np.mean(trends[algo][index]['mean'])
    plt.errorbar(x_values, bias, 2 * (np.array(trends[algo][index]['mean_err'])), lw=1, fmt='-', color=colors[algo], elinewidth=0.4, capthick=0.4, capsize=2, barsabove=True)
    ymin = min(ymin, np.min(bias))
    ymax = max(ymax, np.max(bias))
  ymin = 10**np.floor(np.log10(ymin) - 0.1)
  ymax = 10**np.ceil(np.log10(ymax) + 0.1)
  plot_loglines(log10(xmin), log10(xmax), log10(ymin), log10(ymax))
  plt.xlim([xmin, xmax])
  plt.ylim([ymin, ymax])
  plt.ylabel(r'Estimator bias')
  plt.xlabel(r'Model evaluations')
  make_legend()
  plt.subplot(212)
  ax = plt.gca()
  ax.set_xscale('log')
  ax.set_yscale('log', nonposy='clip')
  ymin, ymax = 1e8, 0
  for i, algo in enumerate(algos):
    variance = trends[algo][index]['var']
    plt.errorbar(x_values, variance, 2 * (np.array(trends[algo][index]['var_err'])), lw=1, fmt='-', color=colors[algo], elinewidth=0.4, capthick=0.4, capsize=2, barsabove=True)
    ymin = min(ymin, np.min(variance))
    ymax = max(ymax, np.max(variance))
  ymin = 10**np.floor(np.log10(ymin) - 0.1)
  ymax = 10**np.ceil(np.log10(ymax) + 0.1)
  plot_loglines(log10(xmin), log10(xmax), log10(ymin), log10(ymax))
  plt.xlim([xmin, xmax])
  plt.ylim([ymin, ymax])
  plt.ylabel(r'Estimator variance')
  plt.xlabel(r'Model evaluations')
  save_fig('expdesign/plots/mossbauer_convergence_error.pdf')

