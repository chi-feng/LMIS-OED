import numpy as np
import struct


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


def eig_joint(d, k, s):
    return -0.5 * np.log(
        s ** 8
        / (
            ((1.0 + (-1.0 + d) * k) ** 2 + s ** 2) ** 2
            * (
                (3.0 + d * k * (-2.0 + (-1.0 + d) * k)) ** 2
                + (10.0 + k * (4.0 + k + 2.0 * d * (-2.0 + (-1.0 + d) * k))) * s ** 2
                + s ** 4
            )
        )
    )


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
    print("Reading %s" % filename, end=" ")
    data = []
    for line in open(filename, "r"):
        if line[0:3] == "EIG":
            data.append(float(line[6:]))
    print("got %d entries" % len(data))
    """
  if length > 0 and len(data) > length:
    print 'Note: truncated to %d entries' % length
    return np.array(data[0:length])
  """
    return np.array(data)


def read_binary(filename):
    f = open(filename, "rb")
    rows = struct.unpack("i", f.read(4))[0]
    cols = struct.unpack("i", f.read(4))[0]
    m = np.empty([rows, cols])
    for i in range(rows):
        for j in range(cols):
            m[i, j] = struct.unpack("d", f.read(8))[0]
    return m


def get_options(default, algo):
    options = default.copy()
    options["-useExactPosterior"] = 0
    options["-biasingDistributionType"] = "MVT"
    if algo == "prior":
        options["-useIS"], options["-useMIS"] = 0, 0
    if algo == "is":
        options["-useIS"], options["-useMIS"] = 1, 0
    if algo == "mis":
        options["-useIS"], options["-useMIS"] = 0, 1
    if algo == "exact":
        options["-useIS"], options["-useMIS"] = 1, 0
        options["-useExactPosterior"] = 1
        options["-biasingDistributionType"] = "MVN"
    return options
