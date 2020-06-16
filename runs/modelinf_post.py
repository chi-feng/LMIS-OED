import numpy as np
import sys, os
from plotting import *

# import seaborn as sns

data = np.genfromtxt(sys.argv[1], skip_header=0)

gridpts = 201

logPosterior = data[:, -1].reshape((gridpts, gridpts)).T

x = np.linspace(0, 0.5, gridpts)
y = np.linspace(0, 0.5, gridpts)
X, Y = np.meshgrid(x, y)

print "plotting"
"""
print x[275]

fig = plt.figure()
fig.set_size_inches(5, 4)
plt.plot(y, logPosterior[:,275])

save_fig('%s-conditional.pdf' % sys.argv[1])
os.system('open %s-conditional.pdf' % sys.argv[1])
exit()
"""

fig = plt.figure()
fig.set_size_inches(5, 4)
plt.contourf(X, Y, np.exp(logPosterior), 100)  #
plt.colorbar()
plt.xticks(np.linspace(np.min(x), np.max(x), 11))
plt.grid(b=True, which="major", color="black", alpha=0.25, linestyle="-")
plt.xlabel(r"$\alpha_{10}$")
plt.ylabel(r"$\alpha_{11}$")
plt.title(r"Unnormalized posterior density $p(\alpha|y)$")
save_fig("%s.pdf" % sys.argv[1])
os.system("open %s.pdf" % sys.argv[1])
exit()

fig = plt.figure()
fig.set_size_inches(5, 3)
plt.plot((0, 0.5, 1.0), (0.4, 0.6, 0.4), "ko")
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.yticks(np.linspace(0, 1, 11))
plt.grid(b=True, which="major", color="black", alpha=0.25, linestyle="-")
plt.xlim([-0.1, 1.1])
plt.ylim([0, 1])
plt.title(r"Data")
save_fig("data.pdf")
