import multiprocessing as mp
from multiprocessing.pool import Pool
import subprocess
import random
import numpy as np
import os, sys, time, datetime
from common import *

executable = "./build/Driver"
if not os.path.isfile(executable[2:]):
    print "could not find executable"
    sys.exit()

if sys.argv[1] == "2":
    dim, sigeps = 2, 0.1

if sys.argv[1] == "4":
    dim, sigeps = 4, 0.4

if sys.argv[1] == "8":
    dim, sigeps = 8, 0.4

print dim, sigeps

algorithms = ["prior", "mis", "exact"]

algorithms = ["mis"]

# budgets = np.logspace(3, 6, 7)
# budgets = np.logspace(3, 7, 9)
# budgets = np.logspace(6.5, 7, 2)
# budgets = np.logspace(2,2.5,2)
# budgets = [10**6.5]
# budgets = np.array([10**2,10**2.5,10**2.5,10**2.75])
budgets = np.logspace(3, 6.5, 8)
budgets = budgets[:-1]


budgets = np.array(
    [
        1000.0,
        3162.27766017,
        4500,
        7000,
        10000.0,
        12000.0,
        15000.0,
        20000,
        31622.77660168,
        100000.0,
        316227.76601684,
        1000000.0,
    ]
)


# budgets = np.array([12000])

print budgets

suffix = "_3"


base_ratio = 100


def get_samples(budget, ratio):
    N = 0.25 * (-ratio + np.sqrt(ratio) * np.sqrt(ratio + 8 * budget))
    M = (-ratio + np.sqrt(ratio) * np.sqrt(ratio + 8 * budget)) / (4 * ratio)
    if budget < 10 ** 3:
        return int(N + 0.5), int(M + 0.5)
    else:
        return int(N), int(M)


for budget in budgets:
    N, M = get_samples(budget, base_ratio)
    print N, M, N * (M + M) + N
    N, M = get_samples(budget, 1.0 / base_ratio)
    print N, M, N * (M + M) + N

design = 0.8

trials = 1000

options = {
    "-dim": dim,
    "-sigeps": sigeps,
    "-poi": 1,
    "-design": design,
    "-seed": 0,
    "-dumpfile": "none",
    "-outfile": "none",
    "-N": 0,
    "-M1": 0,
    "-M2": 0,
    "-useIS": 0,
    "-useMIS": 0,
    "-useReverseLikelihood": 1,
    "-useMinSampleDistance": 0,
    "-sort": 1,
    "-useExactPosterior": 0,
    "-maxComponents": 0,
    "-biasingDistributionType": "MVT",
    "-dof": 2.5,
    "-nugget": 1e-3,
}

runs = []


def get_run(budget, algo):

    N, M = get_samples(budget, base_ratio if not algo == "prior" else 1.0 / base_ratio)

    options["-N"] = N
    options["-M1"] = M
    options["-M2"] = M
    options["-maxComponents"] = min(N, M)
    options["-useExactPosterior"] = 0
    options["-biasingDistributionType"] = "MVT"
    if algo is "prior":
        options["-useIS"], options["-useMIS"] = 0, 0
    if algo is "is":
        options["-useIS"], options["-useMIS"] = 1, 0
    if algo is "mis":
        options["-useIS"], options["-useMIS"] = 0, 1
    if algo is "exact":
        options["-useIS"], options["-useMIS"] = 1, 0
        options["-useExactPosterior"] = 1
        options["-biasingDistributionType"] = "MVN"
        options["-M1"] = 2
        options["-M2"] = 2

    filename = "out/fixed_budget/%dd_%s_N%04d_M%04d%s" % (dim, algo, N, M, suffix)
    arguments = [executable]
    for k, v in options.iteritems():
        arguments.append("%s" % k)
        arguments.append("%s" % v)
    return {"filename": filename, "arguments": arguments}


for algo in algorithms:
    for trial in xrange(trials):
        for budget in budgets:
            run = get_run(budget, algo)
            runs.append(run)

random.shuffle(runs)


def execute(run):
    # print run
    output = subprocess.check_output(run["arguments"])
    f = open(run["filename"], "a")
    f.write(output)
    f.close()


pool = Pool(processes=mp.cpu_count())
result = pool.map_async(execute, runs, chunksize=1)
start_time = time.time()
while not result.ready():
    completed = len(runs) - result._number_left
    if completed > 0:
        eta = (time.time() - start_time) / completed * result._number_left / 60
        print "%6d left, ETA: %4.1f minutes" % (result._number_left, eta)
    sys.stdout.flush()
    time.sleep(1)

pool.close()

# nice parallel -j 6 --"command1""command2"
