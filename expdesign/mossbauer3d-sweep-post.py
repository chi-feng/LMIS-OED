from mossbauer import *
from plotting import *
import scipy.ndimage as ndimage

npts = 21
xmin, xmax = -2, 2
ymin, ymax = -2, 2
algos = ["prior", "mis"]

suffix = ""
if "height" in sys.argv:
    suffix = "_height"
if "joint" in sys.argv:
    suffix = "_joint"

d0 = np.zeros(npts)
d1 = np.linspace(-2, 2, npts)
d2 = np.linspace(-2, 2, npts)

eig = np.zeros(npts * npts)

for algo in algos:
    eig = np.zeros(npts * npts)
    for i in range(npts * npts):
        data = get_eig("expdesign/out/mossbauer_3d_sweep%s_%s_%04d" % (suffix, algo, i))
        eig[i] = np.mean(data)
    eig2 = eig.reshape((npts, npts)).T
    eig = np.zeros((npts, npts))
    for i in range(npts):
        for j in range(npts):
            eig[i, j] = (eig2[i, j] + eig2[j, i]) / 2
    eig = ndimage.gaussian_filter(eig, sigma=0.5, order=0)
    X, Y = np.meshgrid(d1, d2)

    fig = plt.figure()
    fig.set_size_inches(3.5, 2.5)
    plt.contourf(X, Y, eig, 40, cmap=matplotlib.cm.CMRmap)
    plt.colorbar()
    plt.xticks(np.linspace(xmin, xmax, 5))
    plt.yticks(np.linspace(ymin, ymax, 5))
    plt.grid(b=True, which="major", color="black", alpha=0.15, linestyle="-")
    plt.xlabel(r"$d_1$")
    plt.ylabel(r"$d_2$")
    save_fig("%s.pdf" % "expdesign/plots/mossbauer_3d_sweep%s_%s" % (suffix, algo))
    os.system(
        "open %s.pdf" % "expdesign/plots/mossbauer_3d_sweep%s_%s" % (suffix, algo)
    )

exit()
