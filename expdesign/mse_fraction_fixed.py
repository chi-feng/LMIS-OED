from expdesign import *

algorithms = ["prior", "is", "mis"]

evaluations = np.logspace(3, 6, 4)

options = default.copy()
options["--design"] = 0.5
trials = 1
if len(sys.argv) > 1:
    trials = int(sys.argv[1])

runs = []


def get_run(i, j, algo):
    M_ = 200
    N_ = int(evaluations[j] / 400)
    options["--N"] = N_
    options["--M1"] = M_
    options["--M2"] = M_
    options["--maxComponents"] = 20
    options["--useExactPosterior"] = 0
    options["--biasingDistributionType"] = "MVT"
    if algo is "prior":
        options["--useIS"], options["--useMIS"] = 0, 0
    if algo is "is":
        options["--useIS"], options["--useMIS"] = 1, 0
    if algo is "mis":
        options["--useIS"], options["--useMIS"] = 0, 1
    if algo is "exact":
        options["--useIS"], options["--useMIS"] = 1, 0
        options["--useExactPosterior"] = 1
        options["--biasingDistributionType"] = "MVN"
    filename = "expdesign/out/mse_fraction_%s_%02d_%02d" % (algo, i, j)
    return make_run(filename, options)


for algo in algorithms:
    for trial in range(trials):
        i = 6
        for j in range(4):
            run = get_run(i, j, algo)
            runs.append(run)

random.shuffle(runs)

print("%d runs" % len(runs))

execute_all(runs)
