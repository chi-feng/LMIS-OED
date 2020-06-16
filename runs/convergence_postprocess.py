import numpy as np
from numpy import *
import os, sys
from common import *
from plotting import *

if sys.argv[1] == "2":
    dim, sigma = 2, 0.1

if sys.argv[1] == "4":
    dim, sigma = 4, 0.4

if sys.argv[1] == "8":
    dim, sigma = 8, 0.4

N_values = [50, 100, 200, 400, 800, 1600]
M_values = [50, 100, 200, 400, 800, 1600]

N_values = [50, 100, 200, 400, 800]
M_values = [50, 100, 200, 400, 800]

N_fixed = 800
M_fixed = 800

xmin, xmax = 40, 2000

design = 0.8

trials = 1000

algos = ["prior", "mis", "exact"]

trends = {
    "prior": {"name": "Prior"},
    # 'is':{'name':'LIS'},
    "mis": {"name": "LMIS"},
    "exact": {"name": "Exact"},
}

colors = {
    "prior": "#e41a1c",
    "is": "#377eb8",
    "mis": "#4daf4a",
    "exact": "#984ea3",
}  # from colorbrewer

exact = eig_marginal(dim, design, sigma)

for i, algo in enumerate(algos):

    trends[algo]["N"] = {
        "mean": np.zeros(len(N_values)),
        "mean_err": np.zeros(len(N_values)),
        "var": np.zeros(len(N_values)),
        "var_err": np.zeros(len(N_values)),
    }
    for j, N in enumerate(N_values):
        data = get_eig("out/convergence%dd_%s_N%04d_M%04d_3" % (dim, algo, N, M_fixed))
        trends[algo]["N"]["mean"][j], trends[algo]["N"]["mean_err"][j] = mean(data)
        trends[algo]["N"]["var"][j], trends[algo]["N"]["var_err"][j] = var(data)

    trends[algo]["M"] = {
        "mean": np.zeros(len(M_values)),
        "mean_err": np.zeros(len(M_values)),
        "var": np.zeros(len(M_values)),
        "var_err": np.zeros(len(M_values)),
    }
    for j, M in enumerate(M_values):
        data = get_eig("out/convergence%dd_%s_N%04d_M%04d_3" % (dim, algo, N_fixed, M))
        trends[algo]["M"]["mean"][j], trends[algo]["M"]["mean_err"][j] = mean(data)
        trends[algo]["M"]["var"][j], trends[algo]["M"]["var_err"][j] = var(data)

styles = {}
for algo in algos:
    styles[algo] = plt.Line2D(
        (0, 1), (0, 0), color=colors[algo], linestyle="-", linewidth=1
    )


def make_legend():
    labels = [r"%s" % trends[algo]["name"] for algo in algos]
    label_styles = [styles[algo] for algo in algos]
    plt.legend(
        label_styles,
        labels,
        bbox_to_anchor=(0, 1.05, 1, 1.02),
        loc=3,
        ncol=4,
        mode="expand",
        borderaxespad=0.0,
        numpoints=1,
        prop={"size": 10},
    )


for index in ["N", "M"]:

    x_values = N_values if index is "N" else M_values
    label = "Outer" if index is "N" else "Inner"
    fig = plt.figure()
    fig.set_size_inches(4, 5)
    plt.subplot(211)
    ax = plt.gca()
    ax.set_xscale("log")
    ax.set_yscale("log", nonposy="clip")
    ymin, ymax = 1e8, 0
    for i, algo in enumerate(algos):
        if algo == "exact":
            continue
        bias = trends[algo][index]["mean"] - exact
        print algo, trends[algo][index]["mean"], np.mean(trends[algo][index]["mean"])
        plt.errorbar(
            x_values,
            np.abs(bias),
            2 * (np.array(trends[algo][index]["mean_err"])),
            lw=1,
            fmt="-",
            color=colors[algo],
            elinewidth=0.4,
            capthick=0.4,
            capsize=2,
            barsabove=True,
        )
        ymin = min(ymin, np.min(np.abs(bias)))
        ymax = max(ymax, np.max(np.abs(bias)))
        for j, value in enumerate(bias):
            if value < 0:
                plt.scatter(
                    (x_values[j]),
                    (np.abs(value)),
                    s=80,
                    facecolors="none",
                    edgecolors="r",
                )
    ymin = 10 ** np.floor(np.log10(ymin) - 0.5)
    ymax = 10 ** np.ceil(np.log10(ymax) + 0.5)
    plot_loglines(log10(xmin), log10(xmax), log10(ymin), log10(ymax))
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.ylabel(r"Estimator bias")
    # plt.xlabel(r'%s Monte Carlo samples $%s$' % (label, index))
    make_legend()
    plt.subplot(212)
    ax = plt.gca()
    ax.set_xscale("log")
    ax.set_yscale("log", nonposy="clip")
    ymin, ymax = 1e8, 0
    for i, algo in enumerate(algos):
        variance = trends[algo][index]["var"]
        plt.errorbar(
            x_values,
            variance,
            2 * (np.array(trends[algo][index]["var_err"])),
            lw=1,
            fmt="-",
            color=colors[algo],
            elinewidth=0.4,
            capthick=0.4,
            capsize=2,
            barsabove=True,
        )
        ymin = min(ymin, np.min(variance))
        ymax = max(ymax, np.max(variance))
    ymin = 10 ** np.floor(np.log10(ymin) - 0.5)
    ymax = 10 ** np.ceil(np.log10(ymax) + 0.5)
    plot_loglines(log10(xmin), log10(xmax), log10(ymin), log10(ymax))
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.ylabel(r"Estimator variance")
    plt.xlabel(r"%s Monte Carlo samples $%s$" % (label, index))
    filename = "plot/convergence%dd_%s_sweep.pdf" % (dim, index)
    save_fig(filename)
    os.system("gnome-open %s" % filename)
