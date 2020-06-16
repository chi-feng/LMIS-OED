from common import *
from plotting import *

# init_database('mossbauer_convergence');

prefix = ""

if "exact" in sys.argv:
    exact = float(sys.argv[sys.argv.index("exact") + 1])

properties = {
    "experiment": "mossbauer",
    "design": "-1.3,0.0,1.3",
    "sigeps": "0.1",
    "useMarginal": "0",
    "nugget": "0.001",
}

properties["poi"] = 1
properties["index"] = 0
if "joint" in sys.argv:
    properties["poi"] = 4
if "height" in sys.argv:
    properties["index"] = 3

algos = ["prior", "mis"]

evaluations = np.logspace(3, 6, 7)
N_values = {algo: getMossbauerNM(algo, evaluations)[0] for algo in algos}
M_values = {algo: getMossbauerNM(algo, evaluations)[1] for algo in algos}
# print N_values
# print M_values

xmin, xmax = (10 ** 2.75, 10 ** 6.25)

trends = {"prior": {"name": "Prior"}, "mis": {"name": "LMIS"}}

colors = {"prior": "#e41a1c", "mis": "#4daf4a", "mis_marg": "#377eb8"}

conn = sqlite3.connect("mossbauer_convergence.sqlite")
conn.row_factory = sqlite3.Row
cursor = conn.cursor()

for i, algo in enumerate(algos):

    trends[algo]["N"] = {
        "mean": np.zeros(len(evaluations)),
        "mean_err": np.zeros(len(evaluations)),
        "var": np.zeros(len(evaluations)),
        "var_err": np.zeros(len(evaluations)),
        "mse": np.zeros(len(evaluations)),
        "mse_err": np.zeros(len(evaluations)),
    }

    total = np.array(N_values[algo]) * 2 * np.array(M_values[algo])
    # print total;

    for j in range(len(evaluations)):
        N = N_values[algo][j]
        M = M_values[algo][j]

        properties["N"] = N
        properties["M1"] = M
        properties["M2"] = M
        properties["useMarginal"] = 0
        if algo is "prior":
            properties["useMIS"] = 0
        if algo is "mis":
            properties["useMIS"] = 1
        if algo is "mis_marg":
            properties["useMIS"] = 1
            properties["useMarginal"] = 1

        # data = get_eig('expdesign/out/mossbauer%s_convergence_%s_%s_%05d_%05d' % (prefix, suffix, algo, N, M))
        data = get_eig(cursor, properties)
        trends[algo]["N"]["mean"][j], trends[algo]["N"]["mean_err"][j] = mean(data)
        trends[algo]["N"]["var"][j], trends[algo]["N"]["var_err"][j] = var(data)
        trends[algo]["N"]["mse"][j], trends[algo]["N"]["mse_err"][j] = mse(data, exact)

    print trends[algo]["N"]["mean"]

print trends

conn.close()


styles = {}
for algo in algos:
    styles[algo] = plt.Line2D(
        (0, 1), (0, 0), color=colors[algo], linestyle="-", linewidth=1
    )

# print trends


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


index = "N"

fig = plt.figure()
fig.set_size_inches(3, 4)
ax = plt.gca()
ax.set_xscale("log")
ax.set_yscale("log", nonposy="clip")
ymin, ymax = 1e8, 0
pd = {}
for i, algo in enumerate(algos):
    pd[algo] = {}
    mse = trends[algo][index]["mse"]
    mse_err = trends[algo][index]["mse_err"]
    rmse = np.sqrt(trends[algo][index]["mse"])
    rmse_err = 0.5 * trends[algo][index]["mse_err"] / rmse
    x_values = np.array(N_values[algo]) * 2 * np.array(M_values[algo])
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
    pd[algo]["W"] = x_values
    pd[algo]["mse"] = mse
    pd[algo]["mse_err"] = mse_err

    ymin = min(ymin, np.min(mse))
    ymax = max(ymax, np.max(mse))
ymin = 10 ** np.floor(np.log10(ymin) - 0.1)
ymax = 10 ** np.ceil(np.log10(ymax) + 0.1)
plot_loglines(np.log10(xmin), np.log10(xmax), np.log10(ymin), np.log10(ymax))
plt.xlim([xmin, xmax])
plt.ylim([ymin, ymax])
plt.ylabel(r"Estimator MSE")
plt.xlabel(r"Model evaluations")
make_legend()
save_fig("analysis/plots/mossbauer%s_convergence_mse.pdf" % (prefix))

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
    # print algo, trends[algo][index]['mean'], np.mean(trends[algo][index]['mean'])
    x_values = np.array(N_values[algo]) * 2 * np.array(M_values[algo])
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
    pd[algo]["bias"] = bias
    pd[algo]["bias_err"] = 2 * (np.array(trends[algo][index]["mean_err"]))
    ymin = min(ymin, np.min(bias))
    ymax = max(ymax, np.max(bias))
ymin = 10 ** np.floor(np.log10(ymin) - 0.1)
ymax = 10 ** np.ceil(np.log10(ymax) + 0.1)
plot_loglines(np.log10(xmin), np.log10(xmax), np.log10(ymin), np.log10(ymax))
plt.xlim([xmin, xmax])
plt.ylim([ymin, ymax])
plt.ylabel(r"Estimator bias")
plt.xlabel(r"")
# make_legend()
plt.subplot(212)
ax = plt.gca()
ax.set_xscale("log")
ax.set_yscale("log", nonposy="clip")
ymin, ymax = 1e8, 0
for i, algo in enumerate(algos):
    variance = trends[algo][index]["var"]
    x_values = np.array(N_values[algo]) * 2 * np.array(M_values[algo])
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
    pd[algo]["var"] = variance
    pd[algo]["var_err"] = 2 * (np.array(trends[algo][index]["var_err"]))
    ymin = min(ymin, np.min(variance))
    ymax = max(ymax, np.max(variance))
ymin = 10 ** np.floor(np.log10(ymin) - 0.1)
ymax = 10 ** np.ceil(np.log10(ymax) + 0.1)
plot_loglines(np.log10(xmin), np.log10(xmax), np.log10(ymin), np.log10(ymax))
plt.xlim([xmin, xmax])
plt.ylim([ymin, ymax])
plt.ylabel(r"Estimator variance")
plt.xlabel(r"Model evaluations")
save_fig("analysis/plots/mossbauer%s_convergence_error.pdf" % (prefix))

columns = [
    "prior_W",
    "prior_mse",
    "prior_mse_err",
    "prior_bias",
    "prior_bias_err",
    "prior_var",
    "prior_var_err",
]
columns += [
    "mis_W",
    "mis_mse",
    "mis_mse_err",
    "mis_bias",
    "mis_bias_err",
    "mis_var",
    "mis_var_err",
]
mat = np.zeros((len(pd["prior"]["W"]), len(columns)))
mat[:, 0] = pd["prior"]["W"]
mat[:, 1] = pd["prior"]["mse"]
mat[:, 2] = pd["prior"]["mse_err"]
mat[:, 3] = pd["prior"]["bias"]
mat[:, 4] = pd["prior"]["bias_err"]
mat[:, 5] = pd["prior"]["var"]
mat[:, 6] = pd["prior"]["var_err"]

mat[:, 7] = pd["mis"]["W"]
mat[:, 8] = pd["mis"]["mse"]
mat[:, 9] = pd["mis"]["mse_err"]
mat[:, 10] = pd["mis"]["bias"]
mat[:, 11] = pd["mis"]["bias_err"]
mat[:, 12] = pd["mis"]["var"]
mat[:, 13] = pd["mis"]["var_err"]

np.savetxt(
    "plot_data/mossbauer_convergence.csv",
    mat,
    delimiter=",",
    header=",".join(columns),
    comments="",
)
