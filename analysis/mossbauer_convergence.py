from common import *

# init_database('mossbauer_convergence');

options = mossbauer_defaults.copy()
options["-design"] = "-1.3,0.0,1.3"
options["-nugget"] = "0.001"

evaluations = np.logspace(3, 6, 7)
# evaluations = np.logspace(5, 6, 3)

if "prior" in sys.argv:
    algorithms = ["prior"]
    N_values, M_values = getMossbauerNM("prior", evaluations)
elif "mis" in sys.argv:
    algorithms = ["mis"]
    N_values, M_values = getMossbauerNM("mis", evaluations)
else:
    print "no algo"
    exit()

print N_values, M_values

trials = 1
if "trials" in sys.argv:
    trials = int(sys.argv[sys.argv.index("trials") + 1])
if "simulate" in sys.argv:
    options["-simulate"] = 1


suffix = "center"
options["-poi"] = 1
options["-index"] = 0
if "joint" in sys.argv:
    suffix = "joint"
    options["-poi"] = 4
if "height" in sys.argv:
    suffix = "height"
    options["-index"] = 3

if "marginal" in sys.argv:
    options["-useMarginal"] = 1

runs = []


def get_run(N, M, algo):
    options["-N"] = N
    options["-M1"] = M
    options["-M2"] = M
    if algo is "prior":
        options["-useIS"], options["-useMIS"] = 0, 0
    if algo is "mis":
        options["-useIS"], options["-useMIS"] = 0, 1
    return options.copy()


for algo in algorithms:
    for trial in xrange(trials):
        for i in range(len(evaluations)):
            run = get_run(N_values[i], M_values[i], algo)
            runs.append(run)

random.shuffle(runs)

execute_all(runs)
