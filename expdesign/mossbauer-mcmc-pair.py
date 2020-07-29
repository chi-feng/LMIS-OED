from expdesign import *
from plotting import *
from scipy import stats


def plot_trace(chain, burn=0):
    N, dim = np.shape(chain)
    mean = [np.mean(chain[burn:, i]) for i in range(dim)]
    stdev = [np.std(chain[burn:, i]) for i in range(dim)]
    fig = plt.figure(figsize=(10, 3 * dim))
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
                plt.contour(X, Y, Z, 10, cmap=plt.cm.Blues, lw=0.5)
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


def savepair(prefix, chain1, chain2, priors, burn):
    N, dim = np.shape(chain1)
    for i in range(dim):
        for j in range(dim):
            print(
                "computing matrix subplot %d of %d (%d, %d)"
                % (i * dim + j + 1, dim * dim, i, j)
            )
            xmin, xmax = priors[i][0][0], priors[i][0][-1]
            ymin, ymax = priors[j][0][0], priors[j][0][-1]
            X, Y = np.mgrid[xmin:xmax:73j, ymin:ymax:73j]
            positions = np.vstack([X.ravel(), Y.ravel()])
            if i > j:
                values = np.vstack([chain1[burn:, i], chain1[burn:, j]])
                kernel = stats.gaussian_kde(values)
                Z = np.rot90(np.reshape(kernel(positions).T, X.shape))
                np.savetxt(prefix + "1_%d_%d.out" % (i, j), Z)
            if i < j:
                values = np.vstack([chain2[burn:, j], chain2[burn:, i]])
                kernel = stats.gaussian_kde(values)
                positions = np.vstack([Y.ravel(), X.ravel()])
                Z = np.rot90(np.reshape(kernel(positions).T, X.shape))
                np.savetxt(prefix + "2_%d_%d.out" % (i, j), Z)
            if i == j:
                x = np.linspace(xmin, xmax, 101)
                kernel = stats.gaussian_kde(chain1[burn:, i])
                np.savetxt(prefix + "1_%d_%d.out" % (i, j), kernel(x))
                kernel = stats.gaussian_kde(chain2[burn:, i])
                np.savetxt(prefix + "2_%d_%d.out" % (i, j), kernel(x))


border_color = "#999999"


def plotpair(dim, prefix="", priors=None, names=None, burn=1000, levels=[5, 11]):
    fig = plt.figure(figsize=(2 * dim, 2 * dim))
    for i in range(dim):
        for j in range(dim):
            print(
                "plotting matrix subplot %d of %d (%d, %d)"
                % (i * dim + j + 1, dim * dim, i, j)
            )
            xmin, xmax = priors[i][0][0], priors[i][0][-1]
            ymin, ymax = priors[j][0][0], priors[j][0][-1]
            X, Y = np.mgrid[xmin:xmax:73j, ymin:ymax:73j]
            positions = np.vstack([X.ravel(), Y.ravel()])

            if i > j:
                ax = plt.subplot(dim, dim, i * dim + j + 1)
                Z = np.loadtxt(prefix + "1_%d_%d.out" % (i, j))
                cl = np.linspace(np.min(Z), np.max(Z), levels[1])
                plt.contourf(X, Y, Z, cl, cmap=plt.cm.Blues)

                if j == 0:
                    ax = plt.gca()
                    ax.yaxis.tick_left()
                    ax.yaxis.set_label_position("left")
                    ax.set_ylabel(names[i])
                else:
                    plt.yticks([], [])

                if i == dim - 1:
                    ax = plt.gca()
                    ax.xaxis.tick_bottom()
                    ax.xaxis.set_label_position("bottom")
                    xloc = plt.MaxNLocator(5)
                    ax.xaxis.set_major_locator(xloc)
                    plt.xlabel(names[j])
                else:
                    plt.xticks([], [])

                plt.setp(list(ax.spines.values()), color=border_color)
                plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=border_color)
                plt.setp(list(ax.spines.values()), linewidth=0.5)

                ax = plt.subplot(dim, dim, j * dim + i + 1)
                Z = np.loadtxt(prefix + "2_%d_%d.out" % (j, i)).T
                Z = np.rot90(np.rot90(np.rot90(Z)))
                cl = np.linspace(np.min(Z), np.max(Z), levels[1])
                plt.contourf(X, Y, Z, cl, cmap=plt.cm.Oranges)

                if j == 0:
                    ax = plt.gca()
                    ax.xaxis.tick_top()
                    ax.xaxis.set_label_position("top")
                    ax.set_xlabel(names[i])
                    xloc = plt.MaxNLocator(5)
                    ax.xaxis.set_major_locator(xloc)
                else:
                    plt.xticks([], [])

                if i == dim - 1:
                    ax = plt.gca()
                    ax.yaxis.tick_right()
                    ax.yaxis.set_label_position("right")
                    ax.set_ylabel(names[j])
                    yloc = plt.MaxNLocator(5)
                    ax.yaxis.set_major_locator(yloc)
                else:
                    plt.yticks([], [])

                plt.setp(list(ax.spines.values()), color=border_color)
                plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=border_color)
                plt.setp(list(ax.spines.values()), linewidth=0.5)

            if i == j:
                ax = plt.subplot(dim, dim, i * dim + j + 1)
                x = np.linspace(xmin, xmax, 101)
                if priors is not None:
                    plt.plot(priors[i][0], priors[i][1], "k-", alpha=0.5, lw=0.5)
                pdf = np.loadtxt(prefix + "2_%d_%d.out" % (i, j))
                plt.plot(x, pdf, "-", lw=1.5, color="#ff7b31")
                pdf = np.loadtxt(prefix + "1_%d_%d.out" % (i, j))
                plt.plot(x, pdf, "-", lw=1.5, color="#0082bb")
                # plt.title(r'$p(\theta_{%d}|y)$' % (i+1))
                # plt.title(r'%s' % names[i])
                plt.setp(list(ax.spines.values()), color=border_color)
                plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=border_color)
                plt.setp(list(ax.spines.values()), linewidth=0.5)
                plt.xlim([xmin, xmax])

                if j == 0:
                    if i == 0:
                        plt.yticks([], [])
                        ax = plt.gca()
                        ax.xaxis.tick_top()
                        ax.xaxis.set_label_position("top")
                        ax.set_xlabel(names[i])
                        xloc = plt.MaxNLocator(5)
                        ax.xaxis.set_major_locator(xloc)
                else:
                    plt.yticks([], [])

                if i == dim - 1:
                    ax = plt.gca()
                    ax.xaxis.tick_bottom()
                    ax.xaxis.set_label_position("bottom")
                    xloc = plt.MaxNLocator(5)
                    ax.xaxis.set_major_locator(xloc)
                    plt.xlabel(names[j])
                else:
                    if not (i == 0 and j == 0):
                        plt.xticks([], [])


def plotover(dim, prefix="", priors=None, names=None, burn=1000, levels=[5, 9]):
    fig = plt.figure(figsize=(2 * dim, 2 * dim))
    for i in range(dim):
        for j in range(dim):
            print(
                "plotting matrix subplot %d of %d (%d, %d)"
                % (i * dim + j + 1, dim * dim, i, j)
            )
            xmin, xmax = priors[j][0][0], priors[j][0][-1]
            ymin, ymax = priors[i][0][0], priors[i][0][-1]
            X, Y = np.mgrid[xmin:xmax:73j, ymin:ymax:73j]
            positions = np.vstack([X.ravel(), Y.ravel()])

            if i > j:
                ax = plt.subplot(dim, dim, i * dim + j + 1)

                Z = np.loadtxt(prefix + "1_%d_%d.out" % (i, j))
                cl = np.linspace(np.min(Z), np.max(Z), levels[1])
                plt.contour(X, Y, Z, cl, cmap=plt.cm.Blues, alpha=1)

                Z = np.loadtxt(prefix + "2_%d_%d.out" % (j, i)).T
                Z = np.rot90(np.rot90(Z))
                cl = np.linspace(np.min(Z), np.max(Z), levels[1])
                plt.contour(X, Y, Z, cl, cmap=plt.cm.Oranges, alpha=0.8)

                if j == 0:
                    ax = plt.gca()
                    ax.yaxis.tick_left()
                    ax.yaxis.set_label_position("left")
                    ax.set_ylabel(names[i])
                else:
                    plt.yticks([], [])

                if i == dim - 1:
                    ax = plt.gca()
                    ax.xaxis.tick_bottom()
                    ax.xaxis.set_label_position("bottom")
                    xloc = plt.MaxNLocator(5)
                    ax.xaxis.set_major_locator(xloc)
                    plt.xlabel(names[j])
                else:
                    plt.xticks([], [])

                plt.setp(list(ax.spines.values()), color=border_color)
                plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=border_color)
                plt.setp(list(ax.spines.values()), linewidth=0.5)

            if i == j:
                ax = plt.subplot(dim, dim, i * dim + j + 1)
                x = np.linspace(xmin, xmax, 101)
                if priors is not None:
                    plt.plot(priors[i][0], priors[i][1], "k-", alpha=0.5, lw=0.5)
                pdf = np.loadtxt(prefix + "2_%d_%d.out" % (i, j))
                plt.plot(x, pdf, "-", lw=1.5, color="#ff7b31")
                pdf = np.loadtxt(prefix + "1_%d_%d.out" % (i, j))
                plt.plot(x, pdf, "-", lw=1.5, color="#0082bb")
                plt.xlim([xmin, xmax])

                ax = plt.gca()
                ax.xaxis.tick_top()
                ax.xaxis.set_label_position("top")
                ax.set_xlabel(names[i])
                xloc = plt.MaxNLocator(5)
                ax.xaxis.set_major_locator(xloc)
                # plt.xticks([], [])

                plt.yticks([], [])

                if i == dim - 1:
                    ax = plt.gca()
                    ax.xaxis.tick_bottom()
                    ax.xaxis.set_label_position("bottom")
                    xloc = plt.MaxNLocator(5)
                    ax.xaxis.set_major_locator(xloc)
                    plt.xlabel(names[j])
                else:
                    pass

                # plt.title(r'$p(\theta_{%d}|y)$' % (i+1))
                # plt.title(r'%s' % names[i])
                plt.setp(list(ax.spines.values()), color=border_color)
                plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=border_color)
                plt.setp(list(ax.spines.values()), linewidth=0.5)
                # plt.grid()


def gaussian_prior(mu=0, sigma=1, window=2):
    xmin, xmax = mu - window * sigma, mu + window * sigma
    x = np.linspace(xmin, xmax, 101)
    y = np.exp(-(x - mu) ** 2 / (2.0 * sigma * sigma)) / (np.sqrt(2.0 * np.pi) * sigma)
    return x, y


def main():
    # base = sys.argv[1]
    # chain = np.genfromtxt(base + '.chain')

    priors = [
        gaussian_prior(0, 1, window=1),
        gaussian_prior(0, 0.3),
        gaussian_prior(0, 0.3),
        gaussian_prior(1, 0.2, window=1.5),
    ]
    names = [
        r"$\text{center}$",
        r"$\exp(\text{width})$",
        r"$\exp(\text{height})$",
        r"$\exp(\text{offset})$",
    ]

    if len(sys.argv) > 1:
        base1 = "expdesign/out/mcmc/design1"
        base2 = "expdesign/out/mcmc/design2"
        chain1 = np.genfromtxt(base1 + ".chain")
        chain2 = np.genfromtxt(base2 + ".chain")
        savepair("mossbauer", chain1, chain2, priors, burn=int(sys.argv[1]))
    # plot_trace(chain)
    # plt.savefig('%s.trace.pdf' % filename)

    # plotmatrix(chain, priors=priors, names=names)
    # plt.savefig('%s_matrix.pdf' % sys.argv[2], bbox_inches='tight')

    # plotpair(4, prefix='mossbauer', priors=priors, names=names)
    # plt.savefig('pair.pdf', bbox_inches='tight')
    # plt.clf()
    plotover(4, prefix="mossbauer", priors=priors, names=names)
    plt.savefig("overlay.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
