import sys, os, random
from common2 import linear_defaults, execute_all

# usage: python strange.py database db.sqlite dim 4 trials 10 [execute|postprocess]

options = linear_defaults.copy()
options['-design'] = '0.8'
options['-poi'] = '1'

if 'dim' in sys.argv:
    options['-dim'] = int(sys.argv[sys.argv.index('dim') + 1])

trials = 100
if 'trials' in sys.argv:
    trials = int(sys.argv[sys.argv.index('trials') + 1])

N = int(sys.argv[sys.argv.index('N') + 1])
M = int(sys.argv[sys.argv.index('M') + 1])

def get_run(N, M, algo, trial):
    options['-N'] = N
    options['-M1'] = M
    options['-M2'] = M
    options['-maxComponents'] = N
    options['-useExactPosterior'] = 0
    options['-biasingDistributionType'] = 'MVT'
    directory = 'ess/%s_%dd' % (algo, options['-dim'])
    if not os.path.exists(directory):
        os.makedirs(directory)
    options['-dumpfile'] = '%s/N%04d_M%04d_%04d_' % (directory, N, M, trial)
    if algo is 'prior':
        options['-useIS'], options['-useMIS'] = 0, 0
    if algo is 'mis':
        options['-useIS'], options['-useMIS'] = 0, 1
    return options.copy()

runs = []
for algo in ['prior', 'mis']:
    for trial in range(trials):
        run = get_run(N, M, algo, trial)
        runs.append(run)

random.shuffle(runs)

if 'execute' in sys.argv:
    execute_all(runs)
