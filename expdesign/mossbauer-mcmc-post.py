from expdesign import *
from plotting import *
from scipy import stats


def plot_trace(chain, burn=0):
    N, dim = np.shape(chain)
    mean = [np.mean(chain[burn:, i]) for i in range(dim)]
    stdev = [np.std(chain[burn:, i]) for i in range(dim)]
    fig = plt.figure(figsize=(7, 2 * dim))
    for i in range(dim):
        print("plotting trace %d of %d" % (i + 1, dim))
        plt.subplot(dim, 2, i * 2 + 1)
        plt.plot(chain[0:1000, i], "k-", lw=0.5)
        plt.ylabel(r"$\theta_%d$" % (i + 1), rotation=0)
        plt.ylim([mean[i] - 3.0 * stdev[i], mean[i] + 3.0 * stdev[i]])
        plt.grid()
        plt.subplot(dim, 2, i * 2 + 2)
        plt.plot(list(range(N - 1000, N)), chain[N - 1000 :, i], "k-", lw=0.5)
        plt.ylabel(r"$\theta_%d$" % (i + 1), rotation=0)
        plt.ylim([mean[i] - 3.0 * stdev[i], mean[i] + 3.0 * stdev[i]])
        plt.grid()


def plotmatrix(chain, priors=None, names=None, burn=1000, levels=10, colormap=cm.jet):
    N, dim = np.shape(chain)
    mean = [np.mean(chain[burn:, i]) for i in range(dim)]
    stdev = [np.std(chain[burn:, i]) for i in range(dim)]
    fig = plt.figure(figsize=(2 * dim, 2 * dim))
    for i in range(dim):
        for j in range(dim):
            print(
                "plotting matrix subplot %d of %d (%d, %d)"
                % (i * dim + j + 1, dim * dim, i, j)
            )
            xmin, xmax = priors[i][0][0], priors[i][0][-1]
            ymin, ymax = priors[j][0][0], priors[j][0][-1]
            X, Y = np.mgrid[xmin:xmax:51j, ymin:ymax:51j]
            positions = np.vstack([X.ravel(), Y.ravel()])
            if i > j:
                ax = plt.subplot(dim, dim, i * dim + j + 1)
                values = np.vstack([chain[burn:, i], chain[burn:, j]])
                kernel = stats.gaussian_kde(values)
                Z = np.rot90(np.reshape(kernel(positions).T, X.shape))
                plt.contour(X, Y, Z, 20, cmap=plt.cm.Blues, lw=0.5)
                plt.xticks([], [])
                plt.yticks([], [])
                plt.grid()
            if i == j:
                ax = plt.subplot(dim, dim, i * dim + j + 1)
                kernel = stats.gaussian_kde(chain[burn:, i])
                x = np.linspace(xmin, xmax, 101)
                plt.plot(x, kernel(x), "k-", lw=2)
                plt.xlim([xmin, xmax])
                plt.yticks([], [])
                # plt.title(r'$p(\theta_{%d}|y)$' % (i+1))
                plt.title(r"$p(\text{%s}|y)$" % names[i])
                if priors is not None:
                    plt.plot(priors[i][0], priors[i][1], "k--", lw=1)
                xloc = plt.MaxNLocator(5)
                ax.xaxis.set_major_locator(xloc)
                plt.grid()


def gaussian_prior(mu=0, sigma=1):
    xmin, xmax = mu - 1.5 * sigma, mu + 1.5 * sigma
    x = np.linspace(xmin, xmax, 101)
    y = np.exp(-(x - mu) ** 2 / (2.0 * sigma * sigma)) / (np.sqrt(2.0 * np.pi) * sigma)
    return x, y


def main():
    base = sys.argv[1]
    chain = np.genfromtxt(base + ".chain")

    N, dim = np.shape(chain)

    priors = [
        gaussian_prior(0, 1),
        gaussian_prior(0, 0.3),
        gaussian_prior(0, 0.3),
        gaussian_prior(1, 0.2),
    ]
    names = ["center", "width", "height", "offset"]
    # plot_trace(chain)
    # plt.savefig('%s.trace.pdf' % filename)

    plotmatrix(chain, priors=priors, names=names)
    plt.savefig("%s_matrix.pdf" % sys.argv[2], bbox_inches="tight")


if __name__ == "__main__":
    main()
