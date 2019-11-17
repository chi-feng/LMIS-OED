from expdesign import *
from plotting import *

exact = eig_marginal(4, 0.5, 0.4)

neval = 8
nratio = 5

algos  = ['prior', 'mis', 'exact']
lws = {'prior':2,'mis':1}
labels = ['1:100', '1:10', '1:1', '10:1', '100:1']
colors = ['#e41a1c', '#ff7f00', '#4daf4a', '#377eb8', '#984ea3', '#666666']

M = np.genfromtxt('expdesign/out/mse_fraction_M')
N = np.genfromtxt('expdesign/out/mse_fraction_N')
evals = 2 * M * N

mat = {}

fig = plt.figure()
fig.set_size_inches(5, 5)
for algo in ['prior','mis']:
  mat[algo] = np.zeros((neval-2,nratio*2))
  columns=[]
  ax = plt.gca()
  ax.set_xscale('log')
  ax.set_yscale('log', nonposy='clip')
  for i in range(nratio):
    x = np.zeros(neval)
    y = np.zeros(neval)
    y_err = np.zeros(neval)
    for j in range(neval):
      x[j] = evals[i, j]
      data = get_eig('expdesign/out/mse_ratio_%s_%02d_%02d' % (algo, i, j))
      y[j], y_err[j] = mse(data, exact)
    plt.errorbar(x[2:], y[2:], y_err[2:], color=colors[i], label=labels[i], lw=lws[algo], fmt='-', elinewidth=0.5, capthick=0.5, capsize=2, barsabove=True)
    mat[algo][:,i*2] = x[2:]
    mat[algo][:,i*2+1] = y[2:]
    columns += ['x%d'%i,'y%d'%i]
  np.savetxt("plot_data/ratios_%s.csv" % (algo), mat[algo], delimiter=",", header=','.join(columns),comments='')

label_styles = [plt.Line2D((0,1),(0,0), color='#ffffff', linestyle='-', linewidth=1)]
label_styles += [plt.Line2D((0,1),(0,0), color=colors[i], linestyle='-', linewidth=2) for i in range(nratio)]
label_styles += [plt.Line2D((0,1),(0,0), color='#ffffff', linestyle='-', linewidth=1)]
label_styles += [plt.Line2D((0,1),(0,0), color='#ffffff', linestyle='-', linewidth=1)]
label_styles += [plt.Line2D((0,1),(0,0), color=colors[i], linestyle='-', linewidth=1) for i in range(nratio)]
labels = ['Prior', '1:100', '1:10', '1:1', '10:1', '100:1']
labels += ['', 'MIS', '1:100', '1:10', '1:1', '10:1', '100:1']
plt.legend(label_styles, labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., numpoints=1, prop={'size':10})


# plot_loglines(3, 6, -3, 2)
plt.xlim([10**3.8, 10**6.05])
plt.ylim([10**-4.0, 300])

ax = plt.gca()
ax.yaxis.grid(False, linestyle='-', which='major', color='#333333', alpha=0.2)
ax.xaxis.grid(False, linestyle='-', which='major', color='#333333', alpha=0.2)

plt.ylabel(r'Estimator MSE')
plt.xlabel(r'Model evaluations')

plt.savefig('expdesign/plots/mse_ratio.pdf', bbox_inches='tight')


