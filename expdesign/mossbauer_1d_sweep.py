# salloc -n1 -c16 srun python expdesign/mossbauer_1d_sweep.py mis baseline trials 1
from mossbauer import *

options = defaults.copy()

algorithms = []
if 'prior' in sys.argv:
  algorithms.append('prior')
if 'mis' in sys.argv:
  algorithms.append('mis')
if 'baseline' in sys.argv:
  algorithms.append('baseline')

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

designs = []
for d in np.linspace(0,2,21):
    designs.append('%0.1f,0.0,%0.1f' % (-d, d))

runs = []

def get_run(i, algo):
  options['-design'] = designs[i]
  options['-N'] = defaults['-N'];
  options['-M1'] = defaults['-M1'];
  options['-M2'] = defaults['-M2'];
  if algo is 'prior':
    options['-N'] = 100;
    options['-M1'] = 1000;
    options['-M2'] = 1000;
    options['-useIS'], options['-useMIS'] = 0, 0
  if algo is 'mis':
    options['-N'] = 1000;
    options['-M1'] = 100;
    options['-M2'] = 100;
    options['-useIS'], options['-useMIS'] = 0, 1
  if algo is 'baseline':
    options['-useIS'], options['-useMIS'] = 0, 0
    options['-N'] = 10;
    options['-M1'] = 1000000;
    options['-M2'] = 1000000;
  filename = 'expdesign/out/mossbauer_1d_sweep_%s_%s_%0.4d' % (suffix, algo, i)
  return make_run(filename, options)

for algo in algorithms:
  for i in range(len(designs)):
    run = get_run(i, algo)
    for j in range(trials):
      runs.append(run)

random.shuffle(runs)

execute_all(runs)
