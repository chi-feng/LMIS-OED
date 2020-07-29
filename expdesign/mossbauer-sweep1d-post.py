from expdesign import *
from plotting import *

npts = 51
xmin, xmax = 0, 2
# algos = ['prior', 'mis', 'baseline']
algos = ["prior", "mis"]
labels = {"prior": "Prior", "mis": "LMIS", "baseline": "Prior baseline"}
trends = {"prior": {}, "mis": {}, "baseline": {}, "mis_baseline": {}}

suffix = ""
if "height" in sys.argv:
    suffix = "_height"
if "joint" in sys.argv:
    suffix = "_joint"

d = np.linspace(0, 2, npts)

for algo in algos:
    trends[algo]["mean"] = np.zeros(npts)
    trends[algo]["mean_err"] = np.zeros(npts)
    trends[algo]["var"] = np.zeros(npts)
    trends[algo]["var_err"] = np.zeros(npts)
    for i in range(npts):
        data = get_eig(
            "expdesign/out/mossbauer3d_sweep1d%s_%s_%04d" % (suffix, algo, i)
        )
        print(data)
        if len(data) > 0:
            trends[algo]["mean"][i], trends[algo]["mean_err"][i] = mean(data)
            trends[algo]["var"][i], trends[algo]["var_err"][i] = var(data)


fig = plt.figure()
fig.set_size_inches(3.5, 2.5)
for i, algo in enumerate(algos):
    plt.errorbar(
        d,
        trends[algo]["mean"],
        trends[algo]["mean_err"],
        color=colorbrewer(i),
        lw=1,
        fmt="-",
        elinewidth=0.5,
        capthick=0.5,
        capsize=2,
        barsabove=True,
        label=labels[algo],
    )
plt.grid(b=True, which="major", color="black", alpha=0.25, linestyle="-")
plt.xlabel(r"$d_1 = -d_2$")
plt.ylabel(r"Expected Information Gain")
plt.legend(numpoints=1, prop={"size": 10}, loc=4)
save_fig("%s.pdf" % "expdesign/plots/mossbauer3d_sweep1d%s_%s" % (suffix, algo))


exit()
