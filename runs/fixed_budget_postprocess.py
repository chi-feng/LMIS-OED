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

suffix = "_3"

print budgets
base_ratio = 100


def get_samples(budget, ratio):
    N = 0.25 * (-ratio + np.sqrt(ratio) * np.sqrt(ratio + 8 * budget))
    M = (-ratio + np.sqrt(ratio) * np.sqrt(ratio + 8 * budget)) / (4 * ratio)
    if budget < 10 ** 3:
        return int(N + 0.5), int(M + 0.5)
    else:
        return int(N), int(M)


xmin, xmax = 10 ** (2.75), 10 ** (6.75)

design = 0.8

# algos =  ['prior', 'is', 'mis','exact']
algos = ["prior", "mis", "exact"]
# algos =  ['prior', 'is','exact']

trends = {
    "prior": {"name": "Prior"},
    "is": {"name": "LIS"},
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

    trends[algo] = {
        "name": trends[algo]["name"],
        "mean": np.zeros(len(budgets)),
        "mean_err": np.zeros(len(budgets)),
        "var": np.zeros(len(budgets)),
        "var_err": np.zeros(len(budgets)),
        "mse": np.zeros(len(budgets)),
        "mse_err": np.zeros(len(budgets)),
    }
    for j, budget in enumerate(budgets):
        N, M = get_samples(
            budget, base_ratio if not algo == "prior" else 1.0 / base_ratio
        )
        filename = "out/fixed_budget/%dd_%s_N%04d_M%04d%s" % (dim, algo, N, M, suffix)
        data = get_eig(filename)
        # fig = plt.figure()
        # fig.set_size_inches(4,3)
        # plt.hist(data, 100)
        # save_fig('plot/%dd_%s_N%04d_M%04d_3_hist.pdf' % (dim, algo, N, M))

        trends[algo]["mean"][j], trends[algo]["mean_err"][j] = mean(data)
        trends[algo]["var"][j], trends[algo]["var_err"][j] = var(data)
        trends[algo]["mse"][j], trends[algo]["mse_err"][j] = mse(data, exact)
        if algo == "exact":
            trends[algo]["mse"][j] = trends[algo]["var"][j]
            trends[algo]["mse_err"][j] = trends[algo]["var_err"][j]

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


fig = plt.figure()
fig.set_size_inches(3.5, 4.5)
ax = plt.gca()
ax.set_xscale("log")
ax.set_yscale("log", nonposy="clip")
ymin, ymax = 1e8, 0
for i, algo in enumerate(algos):
    x_values = []
    for budget in budgets:
        N, M = get_samples(
            budget, base_ratio if not algo == "prior" else 1.0 / base_ratio
        )
        x_values.append(N + N * (M + M))
    mse = trends[algo]["mse"]
    mse_err = trends[algo]["mse_err"]
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
save_fig("plot/fixed_budget_%dd_mse%s.pdf" % (dim, suffix))

fig = plt.figure()
fig.set_size_inches(3.5, 4.5)
plt.subplot(211)
ax = plt.gca()
ax.set_xscale("log")
ax.set_yscale("log", nonposy="clip")
ymin, ymax = 1e8, 0
for i, algo in enumerate(algos):
    if algo == "exact":
        continue
    x_values = []
    for budget in budgets:
        N, M = get_samples(
            budget, base_ratio if not algo == "prior" else 1.0 / base_ratio
        )
        x_values.append(N + N * (M + M))
    bias = np.abs(trends[algo]["mean"] - exact)
    plt.errorbar(
        x_values,
        bias,
        2 * (np.array(trends[algo]["mean_err"])),
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
# plt.xlabel(r'Model evaluations')
# make_legend()
plt.subplot(212)
ax = plt.gca()
ax.set_xscale("log")
ax.set_yscale("log", nonposy="clip")
ymin, ymax = 1e8, 0
for i, algo in enumerate(algos):
    variance = trends[algo]["var"]
    plt.errorbar(
        x_values,
        variance,
        2 * (np.array(trends[algo]["var_err"])),
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
save_fig("plot/fixed_budget_%dd_error%s.pdf" % (dim, suffix))
