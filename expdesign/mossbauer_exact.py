from mossbauer import *

options = defaults.copy()

options["-design"] = "-1.3,0.0,1.3"

algorithms = ["prior"]
N = 10
M = 1000000

trials = 1
if "trials" in sys.argv:
    trials = int(sys.argv[sys.argv.index("trials") + 1])

suffix = "center"

runs = []


def get_run(algo):
    options["-N"] = N
    options["-M1"] = M
    options["-M2"] = M
    if algo is "prior":
        options["-useIS"], options["-useMIS"] = 0, 0
    filename = "expdesign/out/mossbauer_baseline_%s_%s_%05d_%05d" % (suffix, algo, N, M)
    return make_run(filename, options)


for algo in algorithms:
    for trial in range(trials):
        run = get_run(algo)
        runs.append(run)

random.shuffle(runs)

execute_all(runs)
