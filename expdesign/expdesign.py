import multiprocessing
from multiprocessing.pool import Pool
import subprocess
import struct
import numpy as np
import random
import os, sys, time

executable = "./build/Driver"

if not os.path.isfile(executable[2:]):
    print("could not find %s" % executable)
    # exit()

default = {
    "--experiment": "linear",
    "--dim": 4,
    "--poi": 1,
    "--sigeps": 0.4,
    "--N": 100,
    "--M1": 100,
    "--M2": 100,
    "--maxComponents": 100,
    "--useIS": 0,
    "--useMIS": 0,
    "--useExactPosterior": 0,
    "--useSortHeuristic": 1,
    "--useReverseLikelihood": 1,
    "--useMinSampleDistance": 0,
    "--biasingDistributionType": "MVT",
    "--dof": 2.5,
    "--nugget": 0.001,
    "--seed": 0,
    "--design": "0.8",
    "--dumpfile": "none",
    "--outfile": "none",
}


def read_binary(filename):
    f = open(filename, "rb")
    rows = struct.unpack("i", f.read(4))[0]
    cols = struct.unpack("i", f.read(4))[0]
    m = np.empty([rows, cols])
    for i in range(rows):
        for j in range(cols):
            m[i, j] = struct.unpack("d", f.read(8))[0]
    return m


def make_run(filename, options):
    arguments = [executable]
    for k, v in options.items():
        arguments.append("%s" % k)
        arguments.append("%s" % v)
    return {"filename": filename, "arguments": arguments}


def execute(run):
    # print ' '.join(run['arguments'])
    output = subprocess.check_output(run["arguments"])
    mode = run.get("mode", "a")
    f = open(run["filename"], mode)
    f.write(output)
    f.close()
    for line in output.split("\n"):
        if line[0:3] == "EIG":
            return float(line[6:])


def progressbar(fraction):
    barwidth = 70
    pos = int(barwidth * fraction)
    bar = "["
    for i in range(barwidth):
        if i < pos:
            bar += "="
        elif i == pos:
            bar += ">"
        else:
            bar += " "
    bar += "] "
    return bar


def execute_all(runs, mode="a"):
    pool = Pool(processes=multiprocessing.cpu_count())
    result = pool.map_async(execute, runs, chunksize=1)

    start_time = time.time()
    while not result.ready():
        completed = len(runs) - result._number_left
        if completed > 0:
            rate = 60.0 * completed / (time.time() - start_time)
            percent = 100.0 * completed / len(runs)
            print(
                "%04d/%04d %4.1f%% %s %4d/min  ETA: %3.1f\r"
                % (
                    completed,
                    len(runs),
                    percent,
                    progressbar(percent / 100),
                    rate,
                    1.0 / rate * result._number_left,
                ),
                end=" ",
            )
        sys.stdout.flush()
        time.sleep(0.2)
    print("\n")
    pool.close()


def mean(array):
    return np.mean(array), np.std(array) / np.sqrt(array.size)


def var(array):
    return np.var(array), np.var(array) * np.sqrt(2.0 / (array.size - 1))


def mse(array, exact):
    bias, bias_err = mean(array)
    bias = bias - exact
    vari, vari_err = var(array)
    return (
        vari + bias ** 2,
        np.sqrt((vari_err / vari) ** 2 + (2 * bias_err / bias) ** 2)
        * (vari + bias ** 2),
    )


def get_eig(filename, length=0):
    if not os.path.isfile(filename):
        print("ERROR: Could not find file " + filename)
        return np.array([float("nan")])
    print("Reading %s" % filename, end=" ")
    data = []
    for line in open(filename, "r"):
        if line[0:3] == "EIG":
            data.append(float(line[6:]))
    print("got %d entries" % len(data))
    return np.array(data)


def eig_marginal(dim, d, s):
    # print 'dim = %d, design = %f, sigma=%f' % (dim, d, s)
    if dim == 2:
        return -0.5 * np.log(
            (s ** 2 * (26.0 + 25.0 * (-2.0 + d) * d + s ** 2))
            / (
                (1.0 + 25.0 * (-1.0 + d) * d) ** 2
                + (27.0 + 50.0 * (-1.0 + d) * d) * s ** 2
                + s ** 4
            )
        )
    if dim == 4:
        return -0.5 * np.log(
            (s ** 2 * (52.0 + 5.0 * (-14.0 + 5.0 * d) * d + s ** 2))
            / (
                (3.0 + 5.0 * (-7.0 + 5.0 * d) * d) ** 2
                + 5.0 * (11.0 + 2.0 * d * (-7.0 + 5.0 * d)) * s ** 2
                + s ** 4
            )
        )
