from expdesign import *

algorithms = ["prior", "is", "mis"]

evaluations = np.logspace(3, 6, 4)
ratios = np.array([0.1, 0.5, 1, 2, 10, 20])

M = np.array([np.sqrt(evaluations / (2 * ratio)) for ratio in ratios])
N = np.array([M[i, :] * ratios[i] for i in range(len(ratios))])

M = np.round(M)
N = np.round(N)

np.savetxt("expdesign/out/mse_fraction_M", M)
np.savetxt("expdesign/out/mse_fraction_N", N)

options = default.copy()
options["--design"] = 0.5
trials = 1
if len(sys.argv) > 1:
    trials = int(sys.argv[1])

runs = []


def get_run(i, j, algo):
    N_ = int(N[i, j])
    M_ = int(M[i, j])
    options["--N"] = N_
    options["--M1"] = M_
    options["--M2"] = M_
    options["--maxComponents"] = N_
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
        # for i in range(len(ratios)):
        i = 5
        for j in range(len(evaluations)):
            run = get_run(i, j, algo)
            runs.append(run)

random.shuffle(runs)

print("%d runs" % len(runs))

execute_all(runs)
