import numpy as np
import os, sys
import matplotlib
import math

matplotlib.use("PDF")
matplotlib.rc("font", size="10")

if True:
    matplotlib.rc("text", usetex=True)
    # matplotlib.rc('font', family='sans-serif')
    matplotlib.rc("font", family="serif")
    # matplotlib.rcParams['font.family'] = 'Arial';
    # matplotlib.rc('font' ,serif='Computer Modern Bright')
    matplotlib.rc("font", serif="Computer Modern")
    matplotlib.rcParams[
        "text.latex.preamble"
    ] = r"\usepackage{amsmath}\usepackage{amsfonts}\usepackage{amssymb}"
    # matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}\usepackage{amsfonts}\usepackage{amssymb}\usepackage{helvet} \renewcommand{\familydefault}{\sfdefault}\usepackage{sfmath}'
    matplotlib.rc("font", size="10")
    # matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{cmbright}"]

from matplotlib import cm, pyplot as plt
from matplotlib.mlab import griddata
from matplotlib import gridspec

trendlabels = {
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


def plot_loglines(ymin, ymax, xmin, xmax):
    xmin = int(np.floor(xmin - 0.5))
    xmax = int(np.ceil(xmax + 0.5))
    ymin = int(np.floor(ymin - 0.5))
    ymax = int(np.ceil(ymax + 0.5))
    print xmin, xmax, ymin, ymax
    ax = plt.gca()
    ax.xaxis.set_ticks_position("none")
    ax.yaxis.set_ticks_position("none")
    for power in range(ymin, ymax):
        for i in range(1, 9):
            x = 10 ** power + 10 ** power * i
            plt.axvline(x=x, color="k", alpha=0.1, linewidth=0.3)

    for power in range(xmin, xmax):
        for i in range(1, 9):
            x = 10 ** power + 10 ** power * i
            plt.axhline(y=x, color="k", alpha=0.1, linewidth=0.3)

    for power in range(xmin, xmax):
        plt.axhline(y=10 ** power, color="k", alpha=0.3, linewidth=0.5)
    for power in range(ymin, ymax):
        plt.axvline(x=10 ** power, color="k", alpha=0.3, linewidth=0.5)


def save_fig(filename):
    plt.savefig(filename, bbox_inches="tight", dpi=300)


def colorbrewer(i):
    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3"]  # from colorbrewer
    return colors[i]
