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
    exit()

if sys.argv[1] == "2":
    dim, sigeps = 2, 0.1

if sys.argv[1] == "4":
    dim, sigeps = 4, 0.4

algorithms = ["is", "mis"]
N, M1, M2 = 500, 50, 50
design = 0.8
trials = 10000
nuggets = np.logspace(-4, 0, 11)

default_options = {
    "-dim": dim,
    "-sigeps": sigeps,
    "-poi": 1,
    "-design": design,
    "-seed": 0,
    "-dumpfile": "none",
    "-outfile": "none",
    "-N": N,
    "-M1": M1,
    "-M2": M2,
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

for algo in algorithms:
    for i, nugget in enumerate(nuggets):
        for trial in xrange(trials):
            options = get_options(default_options, algo)
            options["-nugget"] = nugget
            filename = "out/nugget/nugget%dd_%s_%d" % (dim, algo, i)
            arguments = [executable]
            for k, v in options.iteritems():
                arguments.append("%s" % k)
                arguments.append("%s" % v)
            runs.append({"filename": filename, "arguments": arguments})

random.shuffle(runs)


def execute(run):
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
