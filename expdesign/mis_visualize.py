from expdesign import *
from plotting import *

N = 100
M = 20
inputDim = 2

teal = '#4daf4a'
gray = '#aaaaaa'
red = '#e41a1c'
blue = '#377eb8'

def ess(w):
  return np.sum(w)**2 / np.sum(w**2)

def custom_ess(f, w):
  wt = np.abs(f) * w
  wt = wt / np.sum(wt)
  return 1.0 / np.sum(wt**2)

def read_dump(filename):
  data = {}
  data['X_prior'] = read_binary(filename + 'X_prior')
  data['logf_ML'] = read_binary(filename + 'logf_ML')
  data['logw_ML'] = read_binary(filename + 'logw_ML')
  data['logf_CL'] = read_binary(filename + 'logf_CL')
  data['logw_CL'] = read_binary(filename + 'logw_CL')
  data['logML'] = np.genfromtxt(filename + 'logML')
  data['logCL'] = np.genfromtxt(filename + 'logCL')
  data['ESS'] = np.genfromtxt(filename + 'ESS')
  data['postMean'] = [None] * N
  data['postCov'] = [None] * N
  data['postMeanExact'] = [None] * N
  data['postCovExact'] = [None] * N
  data['X_ML'] = [None] * N
  data['indices'] = []
  for line in open(filename + 'mixtureIndices'):
    data['indices'].append(map(int, line.split()))
  postMeanAgg = np.genfromtxt(filename + 'postMean')
  postCovAgg = np.genfromtxt(filename + 'postCov')
  postMeanExactAgg = np.genfromtxt(filename + 'postMeanExact')
  postCovExactAgg = np.genfromtxt(filename + 'postCovExact')
  X_MLAgg = read_binary(filename + 'X_ML')
  for i in xrange(N):
    data['postMean'][i] = postMeanAgg[i,:]
    data['postCov'][i] = postCovAgg[inputDim*i:inputDim*i+inputDim,:]
    data['postMeanExact'][i] = postMeanExactAgg[i,:]
    data['postCovExact'][i] = postCovExactAgg[inputDim*i:inputDim*i+inputDim,:]
    data['X_ML'][i] = X_MLAgg[M*i:M*i+M,:]
  return data

data = read_dump('viz/simple')

print data

from matplotlib.patches import Ellipse
from scipy.stats import chi2
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import matplotlib.lines as mlines

priorsamp = mlines.Line2D([], [], marker='o', markersize=3, markeredgewidth=0, color=gray,linestyle='None' )
yi = mlines.Line2D([], [], marker='*', markersize=10, markeredgewidth=0, color=red,linestyle='None' )
postcont = mlines.Line2D([], [], marker='o', markersize=10, markeredgewidth=1, markeredgecolor=red, color='#ffffff',linestyle='None' )
mixcont = mlines.Line2D([], [], marker='o', markersize=10, markeredgewidth=1, markeredgecolor=teal, color='#ffffff',linestyle='None' )
mixsamp = mlines.Line2D([], [], marker='o', markersize=6, markeredgewidth=0, color=teal,linestyle='None' )
MLcont = mlines.Line2D([], [], marker='o', markersize=10, markeredgewidth=1, markeredgecolor=blue, color='#ffffff',linestyle='None' )
MLsamp = mlines.Line2D([], [], marker='+', markersize=8, markeredgewidth=1, markeredgecolor=blue,color='#ffffff',linestyle='None' )
unused_MLsamp = mlines.Line2D([], [], marker='o', markersize=3, markeredgewidth=0, color=gray,linestyle='None' )

def make_legend():
  plt.legend([yi,postcont,mixcont,mixsamp,MLcont,MLsamp,unused_MLsamp],
    [r'$(\theta^{(i)},\eta^{(i)})$',
      r'$p(\theta,\eta|y^{(i)})$',
      r'$q_\text{post}^{(i)}$',
      r'$(\theta,\eta)\sim q_\text{post}^{(i)}$',
      r'$q_\text{marg}^{(i)}$',
      r'$(\theta,\eta)\sim q_\text{marg}^{(i)}$'],
    handletextpad=0.1, bbox_to_anchor=(0, 1.05, 1, 1.02), loc=3, ncol=3, mode="expand", borderaxespad=0., numpoints=1, prop={'size':9})

def plot_cov_ellipse(cov, pos, volume=.5, ax=None, fc='none', ec=[0,0,0], a=1, lw=1):

  def eigsorted(cov):
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:,order]

  if ax is None:
    ax = plt.gca()

  vals, vecs = eigsorted(cov)
  theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))

  kwrg = {'facecolor':fc, 'edgecolor':ec, 'alpha':a, 'linewidth':lw, 'zorder':11}

  # Width and height are "full" widths, not radius
  width, height = 2 * np.sqrt(chi2.ppf(volume,2)) * np.sqrt(vals)
  ellip = Ellipse(xy=pos, width=width, height=height, angle=theta, **kwrg)

  ax.add_artist(ellip)
  return ellip

def execute(i):
  artists0 = []
  print i
  fig = plt.figure()
  fig.set_size_inches(4, 4)
  plt.xlim([-3,3])
  plt.ylim([-3,3])
  plt.xlabel(r'$\theta$')
  plt.ylabel(r'$\eta$')
  plt.xticks([-2,0,2])
  plt.yticks([-2,0,2])
  plt.plot(data['X_prior'][:,0], data['X_prior'][:,1], 'go', color=gray, ms=2, mew=0, alpha=0.8)
  artist, = plt.plot(data['X_prior'][i,0], data['X_prior'][i,1], 'r*', ms=10, mew=0, zorder=10)
  artists0.append(artist)
  for ii in range(i):
    if ii == i:
      continue
    plt.plot(data['X_ML'][ii][:,0], data['X_ML'][ii][:,1], 'ko', color=gray, ms=2, mew=0,alpha=0.8)
  for tup in [(1.0, 0.68)]:#, (0.4, 0.5), (0.2, 0.75)]:
    artist = plot_cov_ellipse(data['postCovExact'][i], data['postMeanExact'][i], volume=tup[1], ec=red, a=tup[0])
    artists0.append(artist)
  make_legend()
  # plt.savefig('expdesign/plots/mis/%03da.pdf' % (i+1), bbox_inches='tight')
  artists = []
  artists2 = []
  for mi in data['indices'][i]:
    if mi == i:
      continue
    artist = plot_cov_ellipse(data['postCov'][mi]*5, data['postMean'][mi], volume=0.68, ec=teal, a=0.8)
    artists.append(artist)
  for mi in data['indices'][i]:
    if mi == i:
      continue
    mixsamps = data['X_ML'][mi]
    artist, = plt.plot(mixsamps[:,0], mixsamps[:,1], 'go', color=teal, ms=4, mew=0, alpha=0.8)
    artists2.append(artist)
  artist, = plt.plot(data['X_prior'][:,0], data['X_prior'][:,1], 'go', color=teal, ms=4, mew=0, alpha=0.8)
  artists2.append(artist)
  plt.savefig('expdesign/plots/mis/%03db.pdf' % (i+1), bbox_inches='tight')
  for artist in artists:
    artist.remove()
  for tup in [(1.0, 0.68)]:#, (0.4, 0.5), (0.2, 0.75)]:
    artist = plot_cov_ellipse(data['postCov'][i] * 5, data['postMean'][i], volume=tup[1], ec=blue, a=tup[0])
    artists0.append(artist)
  # plt.savefig('expdesign/plots/mis/%03dc.pdf' % (i+1), bbox_inches='tight')
  for artist in artists2:
    artist.remove()
  artist, = plt.plot(data['X_ML'][i][:,0], data['X_ML'][i][:,1], 'g+', ms=8, mew=1, color=blue, mec=blue)
  artists0.append(artist)
  plt.savefig('expdesign/plots/mis/%03dd.pdf' % (i+1), bbox_inches='tight')
  for artist in artists0:
    artist.remove()
  plt.plot(data['X_ML'][i][:,0], data['X_ML'][i][:,1], 'ko', color=blue, ms=4, mew=0,alpha=0.8)
  # plt.savefig('expdesign/plots/mis/%03de.pdf' % (i+1), bbox_inches='tight')
  plt.close()

runs = list(range(10))

pool = Pool(processes=multiprocessing.cpu_count())
result = pool.map(execute, runs)
# while not result.ready():
time.sleep(0.2)
#pool.close()
