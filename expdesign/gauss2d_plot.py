from plotting import *

s = 0.1
d = np.linspace(0, 1, 101)
Ujoint = 0.5 * np.log(((1 - d) ** 2 + s ** 2) * (d ** 2 + s ** 2) / (s ** 4))
Umarg = 0.5 * np.log(1 + d ** 2 / s ** 2)

fig = plt.figure()
fig.set_size_inches(3.5, 2.5)
plt.plot(d, Ujoint, "k-", color="#e41a1c", label="Joint")
plt.plot(d, Umarg, "k-", color="#377eb8", label="Marginal")

plt.xlabel(r"Design $d$")
plt.ylabel(r"Expected Information Gain")
plt.legend(numpoints=1, prop={"size": 10}, loc=4)

save_fig("plots/gauss2d_plot.pdf")
