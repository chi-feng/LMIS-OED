from common import *
from plotting import *

import scipy.ndimage as ndimage

from matplotlib import gridspec

conn = sqlite3.connect("mossbauer_sweep3d.sqlite")
conn.row_factory = sqlite3.Row
cursor = conn.cursor()

properties = {"experiment": "mossbauer"}

npts = 21
xmin, xmax = -2, 2
ymin, ymax = -2, 2
algos = ["mis"]
vmax = 2

suffix = "center"
properties["poi"] = 1
properties["index"] = 0
if "joint" in sys.argv:
    suffix = "joint"
    properties["poi"] = 4
if "height" in sys.argv:
    suffix = "height"
    properties["index"] = 3

d0 = np.zeros(npts)
d1 = np.linspace(-2, 2, npts)
d2 = np.linspace(-2, 2, npts)

eig = np.zeros(npts * npts)

for algo in algos:
    eig = np.zeros(npts * npts)

    if "cached" in sys.argv:
        eig = np.loadtxt("analysis/sweep3d_eig_%s_%s.txt.gz" % (algo, suffix))
    else:
        for i in xrange(npts * npts):
            properties["tag"] = "sweep3d_%03d" % i
            if algo is "prior":
                properties["useMIS"] = 0
            if algo is "mis":
                properties["useMIS"] = 1
            data = get_eig(cursor, properties)
            eig[i] = np.mean(data)
            np.savetxt("analysis/sweep3d_eig_%s_%s.txt.gz" % (algo, suffix), eig)

    eig2d = eig.reshape((npts, npts)).T
    eig = np.zeros((npts, npts))
    for i in xrange(npts):
        for j in xrange(npts):
            eig[i, j] = (eig2d[i, j] + eig2d[j, i]) / 2
    # eig = ndimage.gaussian_filter(eig, sigma=0.5, order=0)
    X, Y = np.meshgrid(d1, d2)

    fig = plt.figure()
    fig.set_size_inches(3.5, 2.5)
    plt.pcolor(
        X, Y, np.clip(eig, 0, np.max(eig)), cmap=matplotlib.cm.viridis, vmin=0, vmax=1.6
    )
    plt.colorbar()
    plt.xticks(np.linspace(xmin, xmax, 5))
    plt.yticks(np.linspace(ymin, ymax, 5))
    plt.grid(b=True, which="major", color="black", alpha=0.15, linestyle="-")
    plt.xlabel(r"$d_1$")
    plt.ylabel(r"$d_2$")
    if suffix == "center":
        plt.plot(1.3, -1.3, "k*", ms=8)
        plt.plot(-1.3, 1.3, "k*", ms=8)
    if suffix == "height":
        plt.plot(1.95, -1.95, "k*", ms=8)
        plt.plot(-1.95, 1.95, "k*", ms=8)

    save_fig("%s.pdf" % "analysis/plots/mossbauer_sweep3d_%s_%s" % (suffix, algo))
    os.system("open %s.pdf" % "analysis/plots/mossbauer_sweep3d_%s_%s" % (suffix, algo))

exit()
