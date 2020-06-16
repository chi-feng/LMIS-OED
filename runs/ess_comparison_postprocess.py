import numpy as np
import os, sys, time, datetime
from common import *

if sys.argv[1] == "2":
    dim, sigeps = 2, 0.1

if sys.argv[1] == "4":
    dim, sigeps = 4, 0.4

N_values = [100, 1000]

samples = 10000

algorithms = ["prior", "mis"]
names = {"prior": "Prior", "mis": "LMIS", "mis_nosort": "LMIS (no sort)"}

design = 0.8


def read_dump(filename):
    print "Reading %s" % filename,
    data = {}
    data["logf_ML"] = read_binary(filename + "logf_ML")
    data["logw_ML"] = read_binary(filename + "logw_ML")
    data["logf_CL"] = read_binary(filename + "logf_CL")
    data["logw_CL"] = read_binary(filename + "logw_CL")
    data["logML"] = np.genfromtxt(filename + "logML")
    data["logCL"] = np.genfromtxt(filename + "logCL")
    # data['priorLogDensities'] = np.genfromtxt(filename + 'priorLogDensities')
    # data['priorLogLikelihoods'] = read_binary(filename + 'priorLogLikelihoods')
    data["ESS"] = np.genfromtxt(filename + "ESS")
    print "done"
    return data


trends = {}


def ess(w):
    return np.sum(w) ** 2 / np.sum(w ** 2)


def custom_ess(f, w):
    wt = np.abs(f) * w
    wt = wt / np.sum(wt)
    return 1.0 / np.sum(wt ** 2)


for algo in algorithms:
    trends[algo] = {}
    for N in N_values:
        trends[algo][str(N)] = {
            "ESS": [],
            "ML_ess": [],
            "ML_custom_ess": [],
            "CL_ess": [],
            "CL_custom_ess": [],
        }
        print "processing trends[%s][%s]" % (algo, str(N))
        trials = samples / N
        for trial in xrange(trials):
            filename = "out/ess/%dd_%s_N%04d_%04d" % (dim, algo, N, trial)
            dump = read_dump(filename)
            for i in xrange(N):
                if algo == "prior":
                    trends[algo][str(N)]["ML_ess"].append(ess(np.ones(N)))
                    trends[algo][str(N)]["ML_custom_ess"].append(
                        custom_ess(np.exp(dump["logf_ML"][i, :]), np.exp(np.zeros(N)))
                    )
                    trends[algo][str(N)]["CL_ess"].append(ess(np.ones(N)))
                    trends[algo][str(N)]["CL_custom_ess"].append(
                        custom_ess(np.exp(dump["logf_CL"][i, :]), np.exp(np.zeros(N)))
                    )
                else:
                    trends[algo][str(N)]["ESS"].append(dump["ESS"][i])
                    trends[algo][str(N)]["ML_ess"].append(
                        ess(np.exp(dump["logw_ML"][i, :]))
                    )
                    trends[algo][str(N)]["ML_custom_ess"].append(
                        custom_ess(
                            np.exp(dump["logf_ML"][i, :]), np.exp(dump["logw_ML"][i, :])
                        )
                    )
                    trends[algo][str(N)]["CL_ess"].append(
                        ess(np.exp(dump["logw_CL"][i, :]))
                    )
                    trends[algo][str(N)]["CL_custom_ess"].append(
                        custom_ess(
                            np.exp(dump["logf_CL"][i, :]), np.exp(dump["logw_CL"][i, :])
                        )
                    )

from plotting import *

colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3"]


def setcolors(bp, num):
    for i in range(num):
        plt.setp(bp["boxes"][i], color=colors[i])
        plt.setp(bp["caps"][i * 2], color=colors[i])
        plt.setp(bp["caps"][i * 2 + 1], color=colors[i])
        plt.setp(bp["whiskers"][i * 2], color=colors[i])
        plt.setp(bp["whiskers"][i * 2 + 1], color=colors[i])
        plt.setp(bp["fliers"][i * 2], color=colors[i])
        plt.setp(bp["fliers"][i * 2 + 1], color=colors[i])
        plt.setp(bp["medians"][i], color=colors[i])


for N in N_values:

    fig = plt.figure()
    fig.set_size_inches(3.5, 2.5)

    algo_data = {}
    tick_names = []

    plot_algs = ["prior", "mis"]

    for algo in plot_algs:

        tick_names.append(names[algo])

        ML_ess = trends[algo][str(N)]["ML_ess"]
        CL_ess = trends[algo][str(N)]["CL_ess"]

        ML_custom_ess = trends[algo][str(N)]["ML_custom_ess"]
        CL_custom_ess = trends[algo][str(N)]["CL_custom_ess"]

        algo_data[algo] = [ML_custom_ess, CL_custom_ess]

    for i, algo in enumerate(plot_algs):
        bp = plt.boxplot(
            algo_data[algo], positions=[i * 3 + 1, i * 3 + 2], widths=0.8, sym=""
        )
        print bp
        setcolors(bp, 2)

    (hB,) = plt.plot([1, 1], color=colors[0])
    (hR,) = plt.plot([1, 1], color=colors[1])
    plt.legend(
        (hB, hR),
        (r"ESS($w_\mathrm{marg}$)", r"ESS($w_\mathrm{cond}$)"),
        loc=2,
        prop={"size": 10},
    )
    hB.set_visible(False)
    hR.set_visible(False)

    plt.xlim([0, len(plot_algs) * 3])
    plt.ylim([0, N])

    ax = plt.gca()
    ax.yaxis.grid(False, linestyle="-", which="major", color="lightgrey", alpha=0)
    ax.set_xticklabels(tick_names, rotation=0, fontsize=10)
    ax.set_xticks([1.5 + 3 * i for i in xrange(len(plot_algs))])

    save_fig("plot/custom_ess%dd_N%04d.pdf" % (dim, N))

    """
  fig = plt.figure()
  fig.set_size_inches(4, 3)

  algo_data = {}
  tick_names = []

  plot_algs = ['prior','mis']

  for algo in plot_algs:

    tick_names.append(names[algo])

    ML_custom_ess = trends[algo][str(N)]['ML_custom_ess']
    CL_custom_ess = trends[algo][str(N)]['CL_custom_ess']
    ESS = trends[algo][str(N)]['ESS']

    algo_data[algo] = [ML_custom_ess, CL_custom_ess, ESS]

  for i, algo in enumerate(plot_algs):
    bp = plt.boxplot(algo_data[algo], positions=[i*4+1,i*4+2,i*4+3], widths=0.8, sym='')
    setcolors(bp, 3)

  h1, = plt.plot([1,1],color=colors[0])
  h2, = plt.plot([1,1],color=colors[1])
  h3, = plt.plot([1,1],color=colors[2])
  plt.legend((h1, h2, h3),(r'ESS($w_\mathrm{marg}$)', r'ESS($w_\mathrm{cond}$)', r'ESS($w_\mathrm{post}$)'),loc=2,prop={'size':10})
  h1.set_visible(False)
  h2.set_visible(False)
  h3.set_visible(False)

  plt.xlim([0,len(plot_algs)*4])

  ax = plt.gca()
  ax.yaxis.grid(False, linestyle='-', which='major', color='lightgrey', alpha=0)
  ax.set_xticklabels(tick_names, rotation=0, fontsize=10)
  ax.set_xticks([2 + 4 * i for i in xrange(len(plot_algs))])

  save_fig('plot/ess%dd_N%04d.pdf' % (dim, N));
  """
