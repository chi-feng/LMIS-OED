from expdesign import *
from plotting import *

exact = eig_marginal(4, 0.5, 0.4)

algos = ["prior", "mis"]
labels = ["1:100", "1:10", "1:1", "10:1", "100:1"]
colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#666666"]

M = np.genfromtxt("expdesign/out/mse_fraction_M")
N = np.genfromtxt("expdesign/out/mse_fraction_N")
evals = 2 * M * N

for algo in algos:
    print(algo)
    fig = plt.figure()
    fig.set_size_inches(6, 4)
    ax = plt.gca()
    ax.set_xscale("log")
    ax.set_yscale("log", nonposy="clip")
    handles = [None] * 6
    for i in range(5):
        x = np.zeros(4)
        y = np.zeros(4)
        y_err = np.zeros(4)
        for j in range(4):
            x[j] = evals[i, j]
            data = get_eig("expdesign/out/mse_ratio_%s_%02d_%02d" % (algo, i, j))
            y[j], y_err[j] = mse(data, exact)
        handles[i] = plt.errorbar(
            x,
            y,
            y_err,
            color=colors[i],
            label=labels[i],
            lw=1,
            fmt="-",
            elinewidth=0.5,
            capthick=0.5,
            capsize=2,
            barsabove=True,
        )

    plt.xlim([1e3, 1e6])
    plt.ylim([2e-4, 400])
    plot_loglines(3, 6, -3, 2)
    plt.ylabel(r"Estimator MSE")
    plt.xlabel(r"Model evaluations")
    plt.legend(
        bbox_to_anchor=(0, 1.05, 1, 1.02),
        loc=3,
        ncol=6,
        mode="expand",
        borderaxespad=0.0,
        numpoints=1,
        prop={"size": 8},
    )
    plt.savefig("expdesign/plots/mse_ratio_%s.pdf" % (algo), bbox_inches="tight")
