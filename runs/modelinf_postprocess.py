import numpy as np
import sys, os
import operator
from plotting import *

dim = int(sys.argv[1])
xdata = np.genfromtxt(sys.argv[2])
ydata = np.genfromtxt(sys.argv[3])
prefix = sys.argv[4]
acf = np.genfromtxt(prefix + ".acf")
chain = np.genfromtxt(prefix + ".chain")
mean = np.mean(chain, axis=0)  # colwise
dens = np.genfromtxt(prefix + ".logdens")
max_index, max_value = max(enumerate(dens), key=operator.itemgetter(1))
MAP = chain[max_index, :]

print "mean ", mean
print "MAP  ", MAP


def realization(x, dim, alpha, xi):

    lam = np.zeros(len(x))
    i = 0
    for d in range(dim):
        lam[d] = alpha[i]
        i += 1
        for j in range(0, d + 1):
            lam[d] += alpha[i] * xi[j]
            i += 1

    y = np.ones(len(x)) * lam[0]
    for i in range(1, dim):
        y += lam[i] * x ** i

    return y


def format_bp(bp, n, color):
    for i in range(n):
        for component in ["boxes", "whiskers", "fliers", "medians", "caps"]:
            for j in range(len(bp[component])):
                plt.setp(bp[component][j], color=color, linewidth=1, linestyle="-")
            for j in range(len(bp["whiskers"])):
                plt.setp(bp["whiskers"][j], alpha=0.5)
            for j in range(len(bp["caps"])):
                plt.setp(bp["caps"][j], alpha=0.5)


nsamp = 10000
samples = np.zeros((nsamp, len(xdata)))
for i in range(nsamp):
    xi = np.random.normal(0, 0.4, dim)
    samples[i, :] = realization(xdata, dim, mean, xi)

fig = plt.figure(figsize=(5, 6))

mean_handle = plt.Line2D(
    (0, 1), (0, 0), color=colorbrewer(1), linestyle="-", linewidth=1, label="Mean"
)
MAP_handle = plt.Line2D(
    (0, 1), (0, 0), color=colorbrewer(0), linestyle="-", linewidth=1, label="MAP"
)
data_handle = plt.Line2D(
    [], [], color="k", linewidth=0, marker="o", markersize=4, label="data"
)
plt.subplot(211)

plt.legend(
    handles=[mean_handle, MAP_handle, data_handle],
    bbox_to_anchor=(0, 1.05, 1, 1.02),
    loc=3,
    ncol=3,
    mode="expand",
    borderaxespad=0.0,
    numpoints=1,
    prop={"size": 10},
)

statistics = [list(samples[:, i]) for i in range(len(xdata))]
bp = plt.boxplot(statistics, positions=xdata, widths=0.035, sym="")
format_bp(bp, len(xdata), colorbrewer(1))

x = np.linspace(-0.1, 1.1, 100)
plt.plot(x, realization(x, dim, mean, np.zeros(dim)), color=colorbrewer(1))

plt.plot(xdata, ydata, "ko", ms=4)
plt.ylim([-0.2, 0.6])
plt.xlim([-0.1, 1.1])

plt.subplot(212)
samples = np.zeros((nsamp, len(xdata)))
for i in range(nsamp):
    xi = np.random.normal(0, 0.4, dim)
    samples[i, :] = realization(xdata, dim, MAP, xi)
statistics = [list(samples[:, i]) for i in range(len(xdata))]
bp = plt.boxplot(statistics, positions=xdata, widths=0.035, sym="")
format_bp(bp, len(xdata), colorbrewer(0))

x = np.linspace(-0.1, 1.1, 100)
plt.plot(x, realization(x, dim, MAP, np.zeros(dim)), color=colorbrewer(0))

plt.plot(xdata, ydata, "ko", ms=4)
plt.ylim([-0.2, 0.6])
plt.xlim([-0.1, 1.1])

"""
plt.subplot(212)
plt.plot(range(1, len(acf)+1), acf)
"""

plt.savefig(prefix + ".pdf", bbox_inches="tight")
os.system("open %s.pdf" % prefix)
