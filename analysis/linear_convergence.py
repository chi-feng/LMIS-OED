from common import *

options = linear_defaults.copy()
options['-design'] = '0.8'

if 'dim' in sys.argv:
  options['-dim'] = int(sys.argv[sys.argv.index('dim') + 1])

evaluations = np.logspace(3, 6, 7)

if 'prior' in sys.argv:
  algorithms = ['prior']
  N_values, M_values = get_linear_NM('prior', evaluations);
elif 'mis' in sys.argv:
  algorithms = ['mis']
  N_values, M_values = get_linear_NM('mis', evaluations);
elif 'exact' in sys.argv:
  algorithms = ['exact']
  N_values, M_values = get_linear_NM('mis', evaluations);
else:
  print 'no algo'
  exit()

print N_values, M_values

trials = 1
if 'trials' in sys.argv:
  trials = int(sys.argv[sys.argv.index('trials') + 1])

if 'simulate' in sys.argv:
  options['-simulate'] = 1

if 'marginal' in sys.argv:
  options['-useMarginal'] = 1

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
  for trial in xrange(trials):
    for i in range(len(evaluations)):
      run = get_run(N_values[i], M_values[i], algo)
      runs.append(run)

random.shuffle(runs)

execute_all(runs)
