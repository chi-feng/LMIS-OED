from expdesign import *

algorithms = ['prior']

designs = []
for d1 in np.linspace(-2, 2, 21):
  for d2 in np.linspace(-2, 2, 21):
    designs.append('%0.1f,%0.1f' % (d1, d2));

f = open('expdesign/out/mossbauer_designs', 'w')
for design in designs:
  f.write('%s\n' % design)
f.close()

options = default.copy()
options['--experiment'] = 'mossbauer'
options['--type'] = '1'
options['--dim'] = '2'
options['--dof'] = '6'
options['--N'] = '1000';
options['--M1'] = '200';
options['--M2'] = '200';
options['--maxComponents'] = '200';

runs = []

def get_run(i, algo):
  options['--design'] = designs[i]
  if algo is 'prior':
    options['--useIS'], options['--useMIS'] = 0, 0
  if algo is 'mis':
    options['--useIS'], options['--useMIS'] = 0, 1
  filename = 'expdesign/out/mossbauer_sweep_gamma_%s_%04d' % (algo, i)
  return make_run(filename, options)

for algo in algorithms:
  for i in range(len(designs)):
    run = get_run(i, algo)
    run['mode'] = 'a'
    runs.append(run)

random.shuffle(runs)

print '%d runs' % len(runs)

execute_all(runs)
