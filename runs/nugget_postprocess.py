import numpy as np
from numpy import *
import os, sys
from common import *
from plotting import *

if sys.argv[1] == "2":
    dim, sigma = 2, 0.1

if sys.argv[1] == "4":
    dim, sigma = 4, 0.4

design = 0.8
algorithms = ["is", "mis"]
nuggets = np.logspace(-4, 0, 11)
xmax = 10 ** (0.5)
xmin = 10 ** (-4.5)
exact = eig_marginal(dim, design, sigma)

for algo in algorithms:
    trends[algo] = {}
    for trend in [
        "mean",
        "mean_err",
        "bias",
        "bias_err",
        "var",
        "var_err",
        "mse",
        "mse_err",
    ]:
        trends[algo][trend] = zeros(len(nuggets))
    for i, nugget in enumerate(nuggets):
        filename = "out/nugget/nugget%dd_%s_%d" % (dim, algo, i)
        data = get_eig(filename)
        trends[algo]["mean"][i], trends[algo]["mean_err"][i] = mean(data)
        trends[algo]["var"][i], trends[algo]["var_err"][i] = var(data)
        trends[algo]["mse"][i], trends[algo]["mse_err"][i] = mse(data, exact)
    trends[algo]["bias"] = trends[algo]["mean"] - exact
    trends[algo]["bias_err"] = trends[algo]["mean_err"]


def plot_trend(trend, ylabel):
    fig = plt.figure()
    fig.set_size_inches(5, 3)
    ax = plt.gca()
    ax.set_xscale("log")
    ax.set_yscale("log", nonposy="clip")
    ymin, ymax = 1e8, 0
    for i, algo in enumerate(algorithms):
        values = trends[algo][trend]
        errors = trends[algo][trend + "_err"]
        plt.errorbar(
            nuggets,
            values,
            errors,
            lw=1,
            fmt="-",
            color=colors[algo],
            elinewidth=0.4,
            capthick=0.4,
            capsize=2,
            barsabove=True,
        )
        ymin = min(ymin, np.min(values))
        ymax = max(ymax, np.max(values))
    ymin = 10 ** np.floor(np.log10(ymin) - 0.5)
    ymax = 10 ** np.ceil(np.log10(ymax) + 0.1)
    plt.legend(["IS", "MIS"], loc=2)
    plot_loglines(log10(xmin), log10(xmax), log10(ymin), log10(ymax))
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.ylabel(ylabel)
    plt.xlabel(r"Nugget")

    save_fig("plot/nugget%dd_%s.pdf" % (dim, trend))


plot_trend("bias", "Bias")
plot_trend("var", "Variance")
plot_trend("mse", "MSE")
