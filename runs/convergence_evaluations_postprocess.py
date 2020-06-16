import numpy as np
from numpy import *
import os, sys
from common import *
from plotting import *

if sys.argv[1] == "2":
    dim, sigma = 2, 0.1

if sys.argv[1] == "4":
    dim, sigma = 4, 0.4

evaluations = np.logspace(3, 6, 11)
M_values = map(int, np.floor(np.sqrt(evaluations / 10) / 2 + 0.5))
N_values = map(int, np.floor(np.sqrt(evaluations / 10) * 10 + 0.5))

xmin, xmax = 10 ** (2.75), 10 ** (6.25)

design = 0.5

trials = 10000

algos = ["prior", "is", "mis", "exact"]

trends = {
    "prior": {"name": "Prior"},
    "is": {"name": "IS"},
    "mis": {"name": "MIS"},
    "exact": {"name": "Posterior"},
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
        "mse": np.zeros(len(N_values)),
        "mse_err": np.zeros(len(N_values)),
    }
    for j, N in enumerate(N_values):
        data = get_eig(
            "out/convergence_eval/eval%dd_%s_N%04d_M%04d" % (dim, algo, N, M_values[j]),
            trials,
        )
        trends[algo]["N"]["mean"][j], trends[algo]["N"]["mean_err"][j] = mean(data)
        trends[algo]["N"]["var"][j], trends[algo]["N"]["var_err"][j] = var(data)
        trends[algo]["N"]["mse"][j], trends[algo]["N"]["mse_err"][j] = mse(data, exact)

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


for index in ["N"]:

    x_values = evaluations
    fig = plt.figure()
    fig.set_size_inches(3, 4)
    ax = plt.gca()
    ax.set_xscale("log")
    ax.set_yscale("log", nonposy="clip")
    ymin, ymax = 1e8, 0
    for i, algo in enumerate(algos):
        mse = trends[algo][index]["mse"]
        mse_err = trends[algo][index]["mse_err"]
        rmse = np.sqrt(trends[algo][index]["mse"])
        rmse_err = 0.5 * trends[algo][index]["mse_err"] / rmse
        plt.errorbar(
            x_values,
            mse,
            mse_err,
            lw=1,
            fmt="-",
            color=colors[algo],
            elinewidth=0.4,
            capthick=0.4,
            capsize=2,
            barsabove=True,
        )
        ymin = min(ymin, np.min(mse))
        ymax = max(ymax, np.max(mse))
    ymin = 10 ** np.floor(np.log10(ymin) - 0.1)
    ymax = 10 ** np.ceil(np.log10(ymax) + 0.1)
    plot_loglines(log10(xmin), log10(xmax), log10(ymin), log10(ymax))
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.ylabel(r"Estimator MSE")
    plt.xlabel(r"Model evaluations")
    make_legend()
    save_fig("plot/eval%dd_%s_mse.pdf" % (dim, index))
    os.system("gnome-open plot/eval%dd_%s_mse.pdf" % (dim, index))

    x_values = evaluations
    fig = plt.figure()
    fig.set_size_inches(3, 4)
    plt.subplot(211)
    ax = plt.gca()
    ax.set_xscale("log")
    ax.set_yscale("log", nonposy="clip")
    ymin, ymax = 1e8, 0
    for i, algo in enumerate(algos):
        if algo == "exact":
            continue
        bias = np.abs(trends[algo][index]["mean"] - exact)
        print algo, trends[algo][index]["mean"], np.mean(trends[algo][index]["mean"])
        plt.errorbar(
            x_values,
            bias,
            2 * (np.array(trends[algo][index]["mean_err"])),
            lw=1,
            fmt="-",
            color=colors[algo],
            elinewidth=0.4,
            capthick=0.4,
            capsize=2,
            barsabove=True,
        )
        ymin = min(ymin, np.min(bias))
        ymax = max(ymax, np.max(bias))
    ymin = 10 ** np.floor(np.log10(ymin) - 0.1)
    ymax = 10 ** np.ceil(np.log10(ymax) + 0.1)
    plot_loglines(log10(xmin), log10(xmax), log10(ymin), log10(ymax))
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.ylabel(r"Estimator bias")
    plt.xlabel(r"Model evaluations")
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
    ymin = 10 ** np.floor(np.log10(ymin) - 0.1)
    ymax = 10 ** np.ceil(np.log10(ymax) + 0.1)
    plot_loglines(log10(xmin), log10(xmax), log10(ymin), log10(ymax))
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.ylabel(r"Estimator variance")
    plt.xlabel(r"Model evaluations")
    save_fig("plot/eval%dd_%s_error.pdf" % (dim, index))
    os.system("gnome-open plot/eval%dd_%s_error.pdf" % (dim, index))
