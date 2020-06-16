from common import *

options = linear_defaults.copy()
options["-design"] = "0.8"

if "dim" in sys.argv:
    options["-dim"] = int(sys.argv[sys.argv.index("dim") + 1])

options["-poi"] = "1"
exact = eig_marginal(
    options["-dim"], float(options["-design"]), float(options["-sigeps"])
)
if "joint" in sys.argv:
    options["-poi"] = "4"
    exact = 4.05967
print exact

algorithms = ["prior"]
W_values = np.array([1e4, 1e5, 1e6, 1e7])
xmin, xmax = (10 ** 3.75, 10 ** 7.25)
# base_ratios = np.array([ 1000, 100, 10, 1]);
base_ratios = np.array([int(sys.argv[sys.argv.index("ratio") + 1])])


trials = 1
if "trials" in sys.argv:
    trials = int(sys.argv[sys.argv.index("trials") + 1])

if "marginal" in sys.argv:
    options["-useMarginal"] = 1

fixed = False
if "fixed" in sys.argv:
    fixed = True
    options["-tag"] = "fixed"

runs = []


def get_run(N, M, algo):
    options["-N"] = N
    options["-M1"] = M
    options["-M2"] = M
    options["-maxComponents"] = min(N, M)
    options["-useExactPosterior"] = 0
    options["-biasingDistributionType"] = "MVT"
    if algo is "prior":
        options["-useIS"], options["-useMIS"] = 0, 0
    if algo is "mis":
        options["-useIS"], options["-useMIS"] = 0, 1
    if algo is "exact":
        options["-useIS"], options["-useMIS"] = 1, 0
        options["-useExactPosterior"] = 1
        options["-biasingDistributionType"] = "MVN"
    return options.copy()


for algo in algorithms:
    for i, base_ratio in enumerate(base_ratios):
        constant = base_ratio * W_values[0] ** (1.0 / 3)
        ratios = constant * W_values ** (-1.0 / 3)
        if fixed:
            ratios = base_ratio * np.ones(len(W_values))
        alpha_values = np.sqrt(2 * ratios)
        N_values = np.sqrt(W_values) / alpha_values
        N_values = map(int, np.floor(N_values + 0.5))
        M_values = np.array(N_values) * alpha_values ** 2 / 2
        M_values = map(int, np.floor(M_values + 0.5))
        print N_values, M_values
        for j in range(len(N_values)):
            for trial in xrange(trials):
                run = get_run(N_values[j], M_values[j], algo)
                runs.append(run)

random.shuffle(runs)

if "execute" in sys.argv:
    execute_all(runs)


def print_row(label, vals, fmt):
    print "%12s" % label,
    for i, val in enumerate(vals):
        print fmt % val,
    print "\n",


if "postprocess" in sys.argv:
    from plotting import *

    pd = {}

    br = base_ratios[0]

    db1 = sys.argv[sys.argv.index("database") + 1]
    db2 = sys.argv[sys.argv.index("database") + 2]
    print db1, db2
    conn1 = sqlite3.connect(db1)
    conn1.row_factory = sqlite3.Row
    cursor1 = conn1.cursor()
    conn2 = sqlite3.connect(db2)
    conn2.row_factory = sqlite3.Row
    cursor2 = conn2.cursor()

    properties = {}
    trends = {"prior": {"name": "Prior"}}
    for algo in algorithms:
        trends[algo] = {}
        for fixed in [0, 1]:
            for base_ratio in base_ratios:
                entries = len(W_values)
                i = (base_ratio, fixed)
                trends[algo][i] = {
                    "mean": np.zeros(entries),
                    "mean_err": np.zeros(entries),
                    "var": np.zeros(entries),
                    "var_err": np.zeros(entries),
                    "mse": np.zeros(entries),
                    "mse_err": np.zeros(entries),
                }
                constant = base_ratio * W_values[0] ** (1.0 / 3)
                ratios = constant * W_values ** (-1.0 / 3)
                if fixed:
                    ratios = base_ratio * np.ones(len(W_values))
                alpha_values = np.sqrt(2 * ratios)
                N_values = np.sqrt(W_values) / alpha_values
                N_values = map(int, np.floor(N_values + 0.5))
                M_values = np.array(N_values) * alpha_values ** 2 / 2
                M_values = map(int, np.floor(M_values + 0.5))
                for j in range(len(N_values)):
                    properties["N"] = N_values[j]
                    properties["M1"] = M_values[j]
                    properties["M2"] = M_values[j]
                    """
          properties['useMarginal'] = 0
          if algo is 'prior':
            properties['useMIS'] = 0
          if algo is 'mis':
            properties['useMIS'] = 1
          if fixed:
            properties['tag'] = 'fixed'
            data = get_eig(cursor2, properties)
          else:
            properties.pop('tag', None)
            data = get_eig(cursor1, properties)
          """
                    data = np.append(
                        get_eig(cursor1, properties), get_eig(cursor2, properties)
                    )
                    trends[algo][i]["mean"][j], trends[algo][i]["mean_err"][j] = mean(
                        data
                    )
                    trends[algo][i]["var"][j], trends[algo][i]["var_err"][j] = var(data)
                    trends[algo][i]["mse"][j], trends[algo][i]["mse_err"][j] = mse(
                        data, exact
                    )
                trends[algo][i]["target_ratios"] = ratios
                trends[algo][i]["N"] = N_values
                trends[algo][i]["M"] = M_values
                trends[algo][i]["W"] = map(
                    int, 2 * np.array(N_values) * np.array(M_values)
                )

                pd[(fixed, "N")] = trends[algo][i]["N"]
                pd[(fixed, "M")] = trends[algo][i]["M"]
                pd[(fixed, "W")] = trends[algo][i]["W"]
                pd[(fixed, "bias")] = np.abs(trends[algo][i]["mean"] - exact)
                pd[(fixed, "var")] = trends[algo][i]["var"]
                pd[(fixed, "mse")] = trends[algo][i]["mse"]
                pd[(fixed, "bias_err")] = trends[algo][i]["mean_err"]
                pd[(fixed, "var_err")] = trends[algo][i]["var_err"]
                pd[(fixed, "mse_err")] = trends[algo][i]["mse_err"]

    for fixed in [0, 1]:
        label = "optimal" if fixed == 0 else "fixed"
        columns = [
            "N",
            "M",
            "W",
            "bias",
            "biase_err",
            "var",
            "var_err",
            "mse",
            "mse_err",
        ]
        mat = np.zeros((len(pd[(fixed, "N")]), 9))
        mat[:, 0] = pd[(fixed, "N")]
        mat[:, 1] = pd[(fixed, "M")]
        mat[:, 2] = pd[(fixed, "W")]
        mat[:, 3] = pd[(fixed, "bias")]
        mat[:, 4] = pd[(fixed, "bias_err")]
        mat[:, 5] = pd[(fixed, "var")]
        mat[:, 6] = pd[(fixed, "var_err")]
        mat[:, 7] = pd[(fixed, "mse")]
        mat[:, 8] = pd[(fixed, "mse_err")]
        np.savetxt(
            "plot_data/scaling_%d_%s.csv" % (br, label),
            mat,
            delimiter=",",
            header=",".join(columns),
            comments="",
        )

    for algo in algorithms:
        for fixed in [0, 1]:
            for base_ratio in base_ratios:
                i = (base_ratio, fixed)
                N = trends[algo][i]["N"]
                M = trends[algo][i]["M"]
                W = trends[algo][i]["W"]
                ratio = np.array(map(float, M)) / np.array(map(float, N))
                print "\nBase ratio: %0.2f %s\n" % (
                    base_ratio,
                    "(fixed)" if fixed else "",
                )
                print_row("goal ratio", trends[algo][i]["target_ratios"], "%12.2f")
                print_row("real ratio", ratio, "%12.2f")
                print_row("N", N, "%12d")
                print_row("M", M, "%12d")
                print_row("W", W, "%12d")
                print_row("bias", trends[algo][i]["mean"] - exact, "%12.4f")
                print_row("var", trends[algo][i]["var"], "%12.4f")
                print_row("mse", trends[algo][i]["mse"], "%12.4f")
        print "\n"

    colors = ["#e41a1c", "#ff7f00", "#4daf4a", "#377eb8", "#984ea3", "#666666"]
    fig = plt.figure()
    fig.set_size_inches(6, 4)

    gs = gridspec.GridSpec(2, 2)
    gs.update(wspace=0.33, hspace=0.33)

    # gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])
    plt.subplot(gs[:, 0])
    ax = plt.gca()
    ax.set_xscale("log")
    ax.set_yscale("log", nonposy="clip")
    ymin, ymax = 1e8, 1e-2
    for i, algo in enumerate(algorithms):
        for index, base_ratio in enumerate(base_ratios):
            key = (base_ratio, 0)
            mse = trends[algo][key]["mse"]
            mse_err = trends[algo][key]["mse_err"]
            plt.errorbar(
                trends[algo][key]["W"],
                mse,
                mse_err,
                lw=1,
                fmt="-",
                color=colors[index],
                elinewidth=0.4,
                capthick=0.4,
                capsize=2,
                barsabove=True,
            )
            ymin = min(ymin, np.min(mse))
            ymax = max(ymax, np.max(mse))

            plt.annotate(
                "$N=%d$\n$M=%dN=%d$"
                % (
                    trends[algo][key]["N"][0],
                    trends[algo][key]["target_ratios"][0] + 0.5,
                    trends[algo][key]["M"][0],
                ),
                xy=(trends[algo][key]["W"][0], mse[0]),
                arrowprops=dict(arrowstyle="->"),
                size=9,
                xytext=(trends[algo][key]["W"][0] * 5, mse[0] * 1.5),
            )

            plt.annotate(
                "$N=%d$\n$M=%dN$"
                % (trends[algo][key]["N"][-1], trends[algo][key]["target_ratios"][-1]),
                xy=(trends[algo][key]["W"][-1], mse[-1]),
                arrowprops=dict(arrowstyle="->"),
                size=9,
                xytext=(trends[algo][key]["W"][-1] / 100, mse[-1] / 1.5),
            )

            key = (base_ratio, 1)
            mse = trends[algo][key]["mse"]
            mse_err = trends[algo][key]["mse_err"]
            plt.errorbar(
                trends[algo][key]["W"],
                mse,
                mse_err,
                lw=1,
                fmt="--",
                color=colors[index],
                elinewidth=0.4,
                capthick=0.4,
                capsize=2,
                barsabove=True,
            )
            ymin = min(ymin, np.min(mse))
            ymax = max(ymax, np.max(mse))

            plt.annotate(
                "$N=%d$\n$M=%dN$"
                % (trends[algo][key]["N"][-1], trends[algo][key]["target_ratios"][-1]),
                xy=(trends[algo][key]["W"][-1], mse[-1]),
                arrowprops=dict(arrowstyle="->"),
                size=9,
                xytext=(trends[algo][key]["W"][-1] / 10, mse[-1] * 10),
            )

    ymin = 10 ** np.floor(np.log10(ymin) - 0.1)
    ymax = 10 ** np.ceil(np.log10(ymax) + 0.1)
    # plot_loglines(np.log10(xmin), np.log10(xmax), np.log10(ymin), np.log10(ymax))
    plt.grid(linestyle="-", alpha=0.3)
    # plt.grid()
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.ylabel(r"MSE")
    plt.xlabel(r"Model evaluations ($W=2MN$)")

    labels1 = [r"%0.1f" % base_ratio for base_ratio in base_ratios]
    labels2 = [r"fixed" for base_ratio in base_ratios]
    styles1 = [
        plt.Line2D((0, 1), (0, 0), color=colors[i], linestyle="-", linewidth=1)
        for i in range(len(base_ratios))
    ]
    styles2 = [
        plt.Line2D((0, 1), (0, 0), color=colors[i], linestyle="--", linewidth=1)
        for i in range(len(base_ratios))
    ]
    labels = labels1 + labels2
    styles = styles1 + styles2
    labels[::2] = labels1
    labels[1::2] = labels2
    styles[::2] = styles1
    styles[1::2] = styles2

    styles = [
        plt.Line2D((0, 1), (0, 0), color=colors[i], linestyle="-", linewidth=1),
        plt.Line2D((0, 1), (0, 0), color=colors[i], linestyle="--", linewidth=1),
    ]

    labels = [r"$\alpha^2\propto W^{-1/3}$", r"$\alpha^2=\text{const}$"]

    plt.legend(
        styles,
        labels,
        bbox_to_anchor=(0, 1.05, 1, 1.02),
        loc=3,
        ncol=3,
        mode="expand",
        borderaxespad=0.0,
        numpoints=1,
        prop={"size": 8},
    )

    # plt.subplot2grid((2,2),(0,1),rowspan=1,colspan=1)
    plt.subplot(gs[0, 1])
    ax = plt.gca()
    ax.set_xscale("log")
    ax.set_yscale("log", nonposy="clip")
    ymin, ymax = 1e8, 1e-2
    for i, algo in enumerate(algorithms):
        for index, base_ratio in enumerate(base_ratios):
            key = (base_ratio, 0)
            mean = trends[algo][key]["mean"]
            mean_err = trends[algo][key]["mean_err"]
            plt.errorbar(
                trends[algo][key]["W"],
                np.abs(mean - exact),
                mean_err,
                lw=1,
                fmt="-",
                color=colors[index],
                elinewidth=0.4,
                capthick=0.4,
                capsize=2,
                barsabove=True,
            )
            ymin = min(ymin, np.min(mean))
            ymax = max(ymax, np.max(mean))
            key = (base_ratio, 1)
            mean = trends[algo][key]["mean"]
            mean_err = trends[algo][key]["mean_err"]
            plt.errorbar(
                trends[algo][key]["W"],
                np.abs(mean - exact),
                mean_err,
                lw=1,
                fmt="--",
                color=colors[index],
                elinewidth=0.4,
                capthick=0.4,
                capsize=2,
                barsabove=True,
            )
            ymin = min(ymin, np.min(np.abs(mean - exact)))
            ymax = max(ymax, np.max(np.abs(mean - exact)))
    ymin = 10 ** np.floor(np.log10(ymin) - 0.1)
    ymax = 10 ** np.ceil(np.log10(ymax) + 0.1)
    # plot_loglines(np.log10(xmin), np.log10(xmax), np.log10(ymin), np.log10(ymax))
    plt.grid(linestyle="-", alpha=0.3)
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    # plt.grid()
    plt.ylabel(r"Bias")
    plt.xlabel(r"Model evaluations")

    # plt.subplot2grid((2,2),(1,1),rowspan=1,colspan=1)
    plt.subplot(gs[1, 1])
    ax = plt.gca()
    ax.set_xscale("log")
    ax.set_yscale("log", nonposy="clip")
    ymin, ymax = 1e8, 1e-2
    for i, algo in enumerate(algorithms):
        for index, base_ratio in enumerate(base_ratios):
            key = (base_ratio, 0)
            var = trends[algo][key]["var"]
            var_err = trends[algo][key]["var_err"]
            plt.errorbar(
                trends[algo][key]["W"],
                np.abs(var),
                var_err,
                lw=1,
                fmt="-",
                color=colors[index],
                elinewidth=0.4,
                capthick=0.4,
                capsize=2,
                barsabove=True,
            )
            ymin = min(ymin, np.min(var))
            ymax = max(ymax, np.max(var))
            key = (base_ratio, 1)
            var = trends[algo][key]["var"]
            var_err = trends[algo][key]["var_err"]
            plt.errorbar(
                trends[algo][key]["W"],
                np.abs(var),
                mean_err,
                lw=1,
                fmt="--",
                color=colors[index],
                elinewidth=0.4,
                capthick=0.4,
                capsize=2,
                barsabove=True,
            )
            ymin = min(ymin, np.min(var))
            ymax = max(ymax, np.max(var))
    ymin = 10 ** np.floor(np.log10(ymin) - 0.1)
    ymax = 10 ** np.ceil(np.log10(ymax) + 0.1)
    # plot_loglines(np.log10(xmin), np.log10(xmax), np.log10(ymin), np.log10(ymax))
    plt.grid(linestyle="-", alpha=0.3)

    # plt.grid()
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.ylabel(r"Variance")
    plt.xlabel(r"Model evaluations")

    save_fig("analysis/plots/scaling_%d.pdf" % (base_ratios[0]))
