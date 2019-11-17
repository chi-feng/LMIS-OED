from expdesign import *

algorithms = ['prior', 'mis']

evaluations = np.logspace(3, 6, 8)
ratios = np.array([0.01, 0.1, 1, 10, 100])

M = np.array([np.sqrt(evaluations / (2 * ratio)) for ratio in ratios])
N = np.array([M[i,:] * ratios[i] for i in range(len(ratios))])

M = np.round(M)
N = np.round(N)

np.savetxt('expdesign/out/mse_fraction_M', M)
np.savetxt('expdesign/out/mse_fraction_N', N)

options = default.copy()
options['--design'] = 0.5

trials = 1
for arg in sys.argv:
  if '--trials=' in arg:
    trials = int(arg[9:])

runs = []

def get_run(i, j, algo):
  N_ = int(N[i, j])
  M_ = int(M[i, j])
  options['--N'] = N_;
  options['--M1'] = M_;
  options['--M2'] = M_;
  options['--maxComponents'] = N_;
  options['--useExactPosterior'] = 0
  options['--biasingDistributionType'] = 'MVT'
  if algo is 'prior':
    options['--useIS'], options['--useMIS'] = 0, 0
  if algo is 'is':
    options['--useIS'], options['--useMIS'] = 1, 0
  if algo is 'mis':
    options['--useIS'], options['--useMIS'] = 0, 1
  if algo is 'exact':
    options['--useIS'], options['--useMIS'] = 1, 0
    options['--useExactPosterior'] = 1
    options['--biasingDistributionType'] = 'MVN'
  filename = 'expdesign/out/mse_ratio_%s_%02d_%02d' % (algo, i, j)
  return make_run(filename, options)

for algo in algorithms:
  for i in range(len(ratios)):
    for j in range(len(evaluations)):
      run = get_run(i, j, algo)
      for trial in xrange(trials):
        runs.append(run)

random.shuffle(runs)

execute_all(runs)
