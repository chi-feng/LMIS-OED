from expdesign import *

algorithms = ['prior', 'mis']

trials = 1
for arg in sys.argv:
  if '--trials=' in arg:
    trials = int(arg[9:])

options = default.copy()
options['--experiment'] = 'mossbauer'
options['--type'] = '0'
options['--dim'] = '3'
options['--design'] = '-2.0,0,1.0'
options['--dof'] = '2.5'
options['--N'] = '1000';
options['--M1'] = '100';
options['--M2'] = '100';
options['--maxComponents'] = '1000';

options['--sigeps'] = '0.1'

suffix = '';
if 'joint' in sys.argv:
  suffix = '_joint'
  options['--poi'] = 4

if 'height' in sys.argv:
  suffix = '_height'
  options['--index'] = 3

runs = []

def get_run(i, algo):
  options['--dumpfile'] = 'expdesign/out/ess/%s%s_%03d' % (suffix, algo, i)
  if algo is 'prior':
    options['--useIS'], options['--useMIS'] = 0, 0
  if algo is 'mis':
    options['--useIS'], options['--useMIS'] = 0, 1
  filename = 'expdesign/out/ess/%s%s_%03d.run' % (suffix, algo, i)
  return make_run(filename, options)

for algo in algorithms:
  for i in range(trials):
    run = get_run(i, algo)
    run['mode'] = 'w'
    runs.append(run)

random.shuffle(runs)

print '%d runs' % len(runs)

execute_all(runs)
