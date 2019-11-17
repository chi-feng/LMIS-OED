from common import *

options = mossbauer_defaults.copy()

if 'prior' in sys.argv:
  algorithms = ['prior']
elif 'mis' in sys.argv:
  algorithms = ['mis']
else:
  print 'no algo'
  exit()

trials = 1
if 'trials' in sys.argv:
  trials = int(sys.argv[sys.argv.index('trials') + 1])

if 'simulate' in sys.argv:
  options['-simulate'] = 1

suffix = 'center';
options['-poi'] = 1
options['-index'] = 0
if 'joint' in sys.argv:
  suffix = 'joint'
  options['-poi'] = 4
if 'height' in sys.argv:
  suffix = 'height'
  options['-index'] = 3

if 'marginal' in sys.argv:
  options['-useMarginal'] = 1

designs = []
for d1 in np.linspace(-2,2,21):
  for d2 in np.linspace(-2,2,21):
    designs.append('0.0,%0.1f,%0.1f' % (d1, d2))

runs = []

def get_run(i, design, algo):
  options['-design'] = design
  options['-tag'] = 'sweep3d_%03d' % i;
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
  return options.copy();

for algo in algorithms:
  for i, design in enumerate(designs):
    run = get_run(i, design, algo)
    for trial in xrange(trials):
      runs.append(run)

# random.shuffle(runs)

execute_all(runs)
