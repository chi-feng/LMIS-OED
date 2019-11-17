from mossbauer import *

options = defaults.copy()

evaluations = np.logspace(3, 6, 7)

options['-design'] = '-1.3,0.0,1.3'

algorithms = []
if 'prior' in sys.argv:
  algorithms.append('prior')
  N_values = map(int, np.floor(np.sqrt(evaluations / 20) + 0.5))
  M_values = map(int, 10 * np.sqrt(evaluations / 20) + 0.5);
elif 'mis' in sys.argv:
  algorithms.append('mis')
  M_values = map(int, np.floor(np.sqrt(evaluations / 10) / 2 + 0.5))
  N_values = map(int, np.floor(np.sqrt(evaluations / 10) * 10 + 0.5))

print('N = ' + str(N_values))
print('M = ' + str(M_values))
print('2NM = ' + str(np.array(N_values) * 2 * np.array(M_values)))

trials = 1
if 'trials' in sys.argv:
  trials = int(sys.argv[sys.argv.index('trials') + 1])

suffix = 'center';
options['-poi'] = 1
options['-index'] = 0
if 'joint' in sys.argv:
  suffix = 'joint'
  options['-poi'] = 4
if 'height' in sys.argv:
  suffix = 'height'
  options['-index'] = 3

runs = []

def get_run(N, M, algo):
  options['-N'] = N;
  options['-M1'] = M;
  options['-M2'] = M;
  if algo is 'prior':
    options['-useIS'], options['-useMIS'] = 0, 0
  if algo is 'mis':
    options['-useIS'], options['-useMIS'] = 0, 1
  filename = 'expdesign/out/mossbauer_convergence_%s_%s_%05d_%05d' % (suffix, algo, N, M)
  return make_run(filename, options)

for algo in algorithms:
  for trial in xrange(trials):
    for i in range(len(evaluations)):
      run = get_run(N_values[i], M_values[i], algo)
      runs.append(run)

random.shuffle(runs)

execute_all(runs)
