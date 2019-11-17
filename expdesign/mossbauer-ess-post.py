from expdesign import *

algorithms = ['prior', 'mis']

trials = 1
for arg in sys.argv:
  if '--trials=' in arg:
    trials = int(arg[9:])

N = 1000
M = 100

algorithms = ['prior', 'mis']
names = {'prior':'Prior', 'is':'IS', 'mis':'LMIS'}
trends = {}

suffix = '';
if 'joint' in sys.argv:
  suffix = '_joint'

if 'height' in sys.argv:
  suffix = '_height'

def ess(w):
  return np.sum(w)**2 / np.sum(w**2)

def custom_ess(f, w):
  wt = np.abs(f) * w
  wt = wt / np.sum(wt)
  return 1.0 / np.sum(wt**2)

def read_dump(filename):
  print 'Reading %s' % filename
  data = {}
  data['logf_ML'] = read_binary(filename + 'logf_ML')
  data['logw_ML'] = read_binary(filename + 'logw_ML')
  data['logf_CL'] = read_binary(filename + 'logf_CL')
  data['logw_CL'] = read_binary(filename + 'logw_CL')
  data['logML'] = np.genfromtxt(filename + 'logML')
  data['logCL'] = np.genfromtxt(filename + 'logCL')
  data['ESS'] = np.genfromtxt(filename + 'ESS')
  return data

for algo in algorithms:
  print algo
  trends[algo] = {
    'ESS':[],
    'ML_ess':[],
    'ML_custom_ess':[],
    'CL_ess':[],
    'CL_custom_ess':[]}
  for trial in range(trials):
    filename = 'expdesign/out/ess/%s%s_%03d' % (suffix, algo, trial)
    dump = read_dump(filename)
    for i in xrange(N):
      if algo == 'prior':
        trends[algo]['ML_ess'].append(ess(np.ones(M)))
        trends[algo]['ML_custom_ess'].append(custom_ess(np.exp(dump['logf_ML'][i,:]),np.exp(np.zeros(M))))
        trends[algo]['CL_ess'].append(ess(np.ones(M)))
        trends[algo]['CL_custom_ess'].append(custom_ess(np.exp(dump['logf_CL'][i,:]),np.exp(np.zeros(M))))
      else:
        trends[algo]['ESS'].append(dump['ESS'][i])
        trends[algo]['ML_ess'].append(ess(np.exp(dump['logw_ML'][i,:])))
        trends[algo]['ML_custom_ess'].append(custom_ess(np.exp(dump['logf_ML'][i,:]),np.exp(dump['logw_ML'][i,:])))
        trends[algo]['CL_ess'].append(ess(np.exp(dump['logw_CL'][i,:])))
        trends[algo]['CL_custom_ess'].append(custom_ess(np.exp(dump['logf_CL'][i,:]),np.exp(dump['logw_CL'][i,:])))


from plotting import *

colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3']

def setcolors(bp, num):
  for i in range(num):
    plt.setp(bp['boxes'][i], color=colors[i])
    plt.setp(bp['caps'][i*2], color=colors[i])
    plt.setp(bp['caps'][i*2+1], color=colors[i])
    plt.setp(bp['whiskers'][i*2], color=colors[i])
    plt.setp(bp['whiskers'][i*2+1], color=colors[i])
    #plt.setp(bp['fliers'][i*2], color=colors[i])
    #plt.setp(bp['fliers'][i*2+1], color=colors[i])
    plt.setp(bp['medians'][i], color=colors[i])

fig = plt.figure()
fig.set_size_inches(3.5, 2.5)

algo_data = {}
tick_names = []

plot_algs = ['prior', 'mis']

for algo in plot_algs:

  tick_names.append(names[algo])

  ML_ess = trends[algo]['ML_ess']
  CL_ess = trends[algo]['CL_ess']

  ML_custom_ess = trends[algo]['ML_custom_ess']
  CL_custom_ess = trends[algo]['CL_custom_ess']

  algo_data[algo] = [ML_custom_ess, CL_custom_ess]

for i, algo in enumerate(plot_algs):
  bp = plt.boxplot(algo_data[algo], positions=[i*3+1,i*3+2], widths=0.8, sym='')
  setcolors(bp, 2)

hB, = plt.plot([1,1],color=colors[0])
hR, = plt.plot([1,1],color=colors[1])
plt.legend((hB, hR),(r'ESS($w_\mathrm{ML}$)', r'ESS($w_\mathrm{CL}$)'),loc=2,prop={'size':10})
hB.set_visible(False)
hR.set_visible(False)

plt.xlim([0,len(plot_algs)*3])
plt.ylim([0,M])
plt.ylabel('Customized ESS')

ax = plt.gca()
ax.yaxis.grid(False, linestyle='-', which='major', color='lightgrey', alpha=0)
ax.set_xticklabels(tick_names, rotation=0, fontsize=10)
ax.set_xticks([1.5 + 3 * i for i in xrange(len(plot_algs))])

save_fig('expdesign/plots/custom_ess%s.pdf' % suffix);
