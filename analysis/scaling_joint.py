from common import *

options = linear_defaults.copy()
options['-design'] = '0.8'

if 'dim' in sys.argv:
  options['-dim'] = int(sys.argv[sys.argv.index('dim') + 1])

options['-poi'] = '1'
exact = eig_marginal(options['-dim'], float(options['-design']), float(options['-sigeps']))
if 'joint' in sys.argv:
  options['-poi'] = '4'
  exact = 4.05967

algorithms = ['prior']
W_values = np.array([ 1e3, 1e4, 1e5, 1e6, 1e7 ]);
xmin, xmax = (10**2.75, 10**7.25)
base_ratios = np.array([ 10, 1, 0.1]);

trials = 1
if 'trials' in sys.argv:
  trials = int(sys.argv[sys.argv.index('trials') + 1])

if 'marginal' in sys.argv:
  options['-useMarginal'] = 1

fixed = False
if 'fixed' in sys.argv:
  fixed = True
  options['-tag'] = 'fixed';

runs = []

def get_run(N, M, algo):
  options['-N'] = N;
  options['-M1'] = M;
  options['-M2'] = M;
  options['-maxComponents'] = min(N, M);
  options['-useExactPosterior'] = 0
  options['-biasingDistributionType'] = 'MVT'
  if algo is 'prior':
    options['-useIS'], options['-useMIS'] = 0, 0
  if algo is 'mis':
    options['-useIS'], options['-useMIS'] = 0, 1
  if algo is 'exact':
    options['-useIS'], options['-useMIS'] = 1, 0
    options['-useExactPosterior'] = 1
    options['-biasingDistributionType'] = 'MVN'
  return options.copy();

for algo in algorithms:
  for i, base_ratio in enumerate(base_ratios):
    constant = base_ratio * W_values[0]**(1./3);
    ratios = constant * W_values**(-1./3);
    if fixed:
      ratios = base_ratio * np.ones(len(W_values));
    alpha_values = np.sqrt(2 * ratios);
    N_values = np.sqrt(W_values) / alpha_values;
    M_values = N_values * alpha_values**2 / 2;

    N_values = map(int, np.ceil(N_values));
    M_values = map(int, np.ceil(M_values));
    print N_values, M_values
    for j in range(len(N_values)):
      for trial in xrange(trials):
        run = get_run(N_values[j], M_values[j], algo)
        runs.append(run)

random.shuffle(runs)

if 'execute' in sys.argv:
  execute_all(runs)

def print_row(label, vals, fmt):
  print '%12s' % label,
  for i, val in enumerate(vals):
    print fmt % val,
  print '\n',

if 'postprocess' in sys.argv:
  from plotting import *
  conn = sqlite3.connect(database)
  conn.row_factory = sqlite3.Row
  cursor = conn.cursor()
  properties = { }
  trends = {
    'prior': {'name':'Prior'}
  }
  for algo in algorithms:
    trends[algo] = { };
    for fixed in [0, 1]:
      for base_ratio in base_ratios:
        entries = len(W_values);
        i = (base_ratio, fixed)
        trends[algo][i] = {
          'mean': np.zeros(entries), 'mean_err': np.zeros(entries),
          'var': np.zeros(entries),  'var_err': np.zeros(entries),
          'mse': np.zeros(entries),  'mse_err': np.zeros(entries)}
        constant = base_ratio * W_values[0]**(1./3);
        ratios = constant * W_values**(-1./3);
        if fixed:
          ratios = base_ratio * np.ones(len(W_values));
        alpha_values = np.sqrt(2 * ratios);
        N_values = np.sqrt(W_values) / alpha_values;
        M_values = N_values * alpha_values**2 / 2;
        N_values = map(int, np.ceil(N_values));
        M_values = map(int, np.ceil(M_values));
        for j in range(len(N_values)):
          properties['N'] = N_values[j];
          properties['M1'] = M_values[j];
          properties['M2'] = M_values[j];
          properties['useMarginal'] = 0
          if algo is 'prior':
            properties['useMIS'] = 0
          if algo is 'mis':
            properties['useMIS'] = 1
          if fixed:
            properties['tag'] = 'fixed'
          data = get_eig(cursor, properties)
          trends[algo][i]['mean'][j], trends[algo][i]['mean_err'][j] = mean(data)
          trends[algo][i]['var'][j], trends[algo][i]['var_err'][j] = var(data)
          trends[algo][i]['mse'][j], trends[algo][i]['mse_err'][j] = mse(data, exact)
        trends[algo][i]['N'] = N_values;
        trends[algo][i]['M'] = M_values;
        trends[algo][i]['W'] = map(int, 2 * np.array(N_values) * np.array(M_values));
  pp.pprint(trends)

  for algo in algorithms:
    for fixed in [0, 1]:
      for base_ratio in base_ratios:
        i = (base_ratio, fixed)
        N = trends[algo][i]['N'];
        M = trends[algo][i]['M'];
        W = trends[algo][i]['W'];
        ratio = np.array(map(float,M))/np.array(map(float,N));
        print '\nBase ratio: %0.2f %s\n' % (base_ratio, '(fixed)' if fixed else '')
        print_row('goal ratio', trends[algo][i]['target_ratios'], '%8.2f');
        print_row('real ratio', ratio, '%8.2f');
        print_row('N', N, '%8d');
        print_row('M', M, '%8d');
        print_row('W', W, '%8d');
    print '\n'

  colors = ['#e41a1c', '#ff7f00', '#4daf4a', '#377eb8', '#984ea3', '#666666']
  fig = plt.figure()
  fig.set_size_inches(6, 4)

  gs = gridspec.GridSpec(2,2)
  gs.update(wspace = 0.33, hspace = 0.33)

  # gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])
  plt.subplot(gs[:,0])
  ax = plt.gca()
  ax.set_xscale('log')
  ax.set_yscale('log', nonposy='clip')
  ymin, ymax = 1e8, 0
  for i, algo in enumerate(algorithms):
    for index, base_ratio in enumerate(base_ratios):
      key = (base_ratio, 0)
      mse = trends[algo][key]['mse']
      mse_err = trends[algo][key]['mse_err']
      plt.errorbar(trends[algo][key]['W'], mse, mse_err, lw=1, fmt='-', color=colors[index], elinewidth=0.4, capthick=0.4, capsize=2, barsabove=True)
      ymin = min(ymin, np.min(mse))
      ymax = max(ymax, np.max(mse))
      key = (base_ratio, 1)
      mse = trends[algo][key]['mse']
      mse_err = trends[algo][key]['mse_err']
      plt.errorbar(trends[algo][key]['W'], mse, mse_err, lw=1, fmt='--', color=colors[index], elinewidth=0.4, capthick=0.4, capsize=2, barsabove=True)
      ymin = min(ymin, np.min(mse))
      ymax = max(ymax, np.max(mse))
  ymin = 10**np.floor(np.log10(ymin) - 0.1)
  ymax = 10**np.ceil(np.log10(ymax) + 0.1)
  plot_loglines(np.log10(xmin), np.log10(xmax), np.log10(ymin), np.log10(ymax))
  plt.xlim([xmin, xmax])
  plt.ylim([ymin, ymax])
  plt.ylabel(r'Estimator MSE')
  plt.xlabel(r'Model evaluations')

  labels1 = [r'%0.1f' % base_ratio for base_ratio in base_ratios];
  labels2 = [r'fixed' for base_ratio in base_ratios];
  styles1 = [plt.Line2D((0,1),(0,0), color=colors[i], linestyle='-', linewidth=1) for i in range(len(base_ratios))];
  styles2 = [plt.Line2D((0,1),(0,0), color=colors[i], linestyle='--', linewidth=1) for i in range(len(base_ratios))];
  labels = labels1+labels2;
  styles = styles1+styles2;
  labels[::2] = labels1;
  labels[1::2] = labels2;
  styles[::2] = styles1;
  styles[1::2] = styles2;

  plt.legend(styles, labels,
    bbox_to_anchor=(0, 1.05, 1, 1.02), loc=3, ncol=3, mode="expand", borderaxespad=0., numpoints=1, prop={'size':8})

  #plt.subplot2grid((2,2),(0,1),rowspan=1,colspan=1)
  plt.subplot(gs[0,1])
  ax = plt.gca()
  ax.set_xscale('log')
  ax.set_yscale('log', nonposy='clip')
  ymin, ymax = 1e8, 0
  for i, algo in enumerate(algorithms):
    for index, base_ratio in enumerate(base_ratios):
      key = (base_ratio, 0)
      mean = trends[algo][key]['mean']
      mean_err = trends[algo][key]['mean_err']
      plt.errorbar(trends[algo][key]['W'], np.abs(mean-exact), mean_err, lw=1, fmt='-', color=colors[index], elinewidth=0.4, capthick=0.4, capsize=2, barsabove=True)
      ymin = min(ymin, np.min(mean))
      ymax = max(ymax, np.max(mean))
      key = (base_ratio, 1)
      mean = trends[algo][key]['mean']
      mean_err = trends[algo][key]['mean_err']
      plt.errorbar(trends[algo][key]['W'], np.abs(mean-exact), mean_err, lw=1, fmt='--', color=colors[index], elinewidth=0.4, capthick=0.4, capsize=2, barsabove=True)
      ymin = min(ymin, np.min(mean-exact))
      ymax = max(ymax, np.max(mean-exact))
  ymin = 10**np.floor(np.log10(ymin) - 0.1)
  ymax = 10**np.ceil(np.log10(ymax) + 0.1)
  plot_loglines(np.log10(xmin), np.log10(xmax), np.log10(ymin), np.log10(ymax))
  plt.xlim([xmin, xmax])
  plt.ylim([ymin, ymax])
  plt.ylabel(r'Estimator Bias')
  plt.xlabel(r'Model evaluations')

  #plt.subplot2grid((2,2),(1,1),rowspan=1,colspan=1)
  plt.subplot(gs[1,1])
  ax = plt.gca()
  ax.set_xscale('log')
  ax.set_yscale('log', nonposy='clip')
  ymin, ymax = 1e8, 0
  for i, algo in enumerate(algorithms):
    for index, base_ratio in enumerate(base_ratios):
      key = (base_ratio, 0)
      var = trends[algo][key]['var']
      var_err = trends[algo][key]['var_err']
      plt.errorbar(trends[algo][key]['W'], np.abs(var), var_err, lw=1, fmt='-', color=colors[index], elinewidth=0.4, capthick=0.4, capsize=2, barsabove=True)
      ymin = min(ymin, np.min(var))
      ymax = max(ymax, np.max(var))
      key = (base_ratio, 1)
      var = trends[algo][key]['var']
      var_err = trends[algo][key]['var_err']
      plt.errorbar(trends[algo][key]['W'], np.abs(var), mean_err, lw=1, fmt='--', color=colors[index], elinewidth=0.4, capthick=0.4, capsize=2, barsabove=True)
      ymin = min(ymin, np.min(var))
      ymax = max(ymax, np.max(var))
  ymin = 10**np.floor(np.log10(ymin) - 0.1)
  ymax = 10**np.ceil(np.log10(ymax) + 0.1)
  plot_loglines(np.log10(xmin), np.log10(xmax), np.log10(ymin), np.log10(ymax))
  plt.xlim([xmin, xmax])
  plt.ylim([ymin, ymax])
  plt.ylabel(r'Estimator Variance')
  plt.xlabel(r'Model evaluations')

  save_fig('analysis/plots/scaling.pdf')
