from mossbauer import *

options = defaults.copy()

algorithms = ['prior']
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
for d1 in np.linspace(-2,2,21):
  for d2 in np.linspace(-2,2,21):
    designs.append('0.0,%0.1f,%0.1f' % (d1, d2))

runs = []

def get_run(i, algo):
  options['-design'] = designs[i]
  if algo is 'prior':
    options['-N'] = 100;
    options['-M1'] = 500;
    options['-M2'] = 500;
    options['-useIS'], options['-useMIS'] = 0, 0
  if algo is 'mis':
    options['-N'] = 1000;
    options['-M1'] = 50;
    options['-M2'] = 50;
    options['-useIS'], options['-useMIS'] = 0, 1
  filename = 'expdesign/out/mossbauer_3d_sweep_%s_%s_%0.4d' % (suffix, algo, i)
  return make_run(filename, options)

for algo in algorithms:
  for i in range(len(designs)):
    run = get_run(i, algo)
    for j in range(trials):
      runs.append(run)

random.shuffle(runs)

execute_all(runs)
