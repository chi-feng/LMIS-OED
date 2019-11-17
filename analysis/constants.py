import numpy as np
import sys

raw = np.genfromtxt(sys.argv[1])
data = {}
data['C1'] = raw[:,0]
data['C2'] = raw[:,1]
data['D1'] = raw[:,2]
data['D2'] = raw[:,3]
data['D3'] = raw[:,4]
data['CESS1'] = raw[:,5]
data['CESS2'] = raw[:,6]


labels = ['C1','C2','D1','D2','D3','CESS1','CESS2']
means = np.mean(raw,axis=0)
errors = np.std(raw,axis=0)/np.sqrt(np.size(raw,axis=0))
for i, label in enumerate(labels):
	if label == 'D3':
		means[i] = np.var(raw[:,i])
		errors[i] = np.var(raw[:,i]) * np.sqrt(2.0 / (np.size(raw,axis=0)-1))
	print '%s = %f +- %0f' % (labels[i], means[i], errors[i])

from plotting import *

colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3']
fig = plt.figure(figsize=(6, 2))
plt.subplot(1,2,1)
bins = 10**np.linspace(-3,3,80)
plt.hist(data['C1'], bins=bins, histtype='step',color=colors[1])
plt.axvline(x=np.mean(data['C1']),color='k')
plt.xscale('log')
plt.xlabel('$C_1$')
plt.subplot(1,2,2)
plt.hist(data['C2'], bins=bins, histtype='step',color=colors[2])
plt.axvline(x=np.mean(data['C2']),color='k')
plt.xscale('log')
plt.xlabel('$C_2$')
plt.savefig('C.pdf', bbox_inches='tight')


fig = plt.figure(figsize=(3, 2))
plt.hist(data['D3']**2-np.mean(data['D3']), np.linspace(-5,20,100), histtype='step',color=colors[1])
plt.axvline(x=np.var(data['D3']),color='k')
plt.xlabel('$D_3$')
plt.savefig('D3.pdf', bbox_inches='tight')


fig = plt.figure(figsize=(3, 2))
bins = np.linspace(70000,90000,100)
plt.hist(data['CESS1'], bins=bins, histtype='step',color=colors[1])
plt.axvline(x=np.mean(data['CESS1']),color=colors[1])
plt.hist(data['CESS2'], bins=bins, histtype='step',color=colors[2])
plt.axvline(x=np.mean(data['CESS2']),color=colors[2])
plt.xlabel(r'$\text{cESS}$ (max $10^4$)')
plt.legend(['cond','marg'],loc='upper center',fontsize='small')
plt.savefig('CESS.pdf', bbox_inches='tight')












'''
fig = plt.figure()
fig.set_size_inches(3.5, 2.5)



bp = plt.boxplot(data['C1']

n = 100

plt.subplot(2,2,1)
get_bins = lambda x,s: np.linspace(np.min(x), np.mean(x) + s * np.std(x), n)

plt.hist(data[:,0], get_bins(data[:,0],1))
plt.xlabel('C1')

plt.subplot(2,2,2)
plt.hist(data[:,1], get_bins(data[:,1],1))
plt.xlabel('C2')

plt.subplot(2,2,3)
plt.hist(data[:,4], get_bins(data[:,4],2))
plt.xlabel('D3')

plt.subplot(2,2,4)
plt.hist(data[:,5], get_bins(data[:,5],4))
plt.hist(data[:,6], get_bins(data[:,6],4))
plt.xlabel('CESS')

plt.show()


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


tick_names = [r'$\text{CESS}_1$', r'$\text{CESS}_2$']
bp = plt.boxplot(data, positions=[i*3+1,i*3+2], widths=0.8, sym='')
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

save_fig('analysis/plots/custom_ess_%s.pdf' % suffix);
'''