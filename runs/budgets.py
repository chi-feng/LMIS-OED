import numpy as np

budgets = np.logspace(3, 6.5, 8)
budgets = budgets[:-1]
budgets = np.array(
    [
        1000.0,
        3162.27766017,
        4500,
        7000,
        10000.0,
        12000.0,
        15000.0,
        20000,
        31622.77660168,
        100000.0,
        316227.76601684,
        1000000.0,
    ]
)

base_ratio = 100


def get_samples(budget, ratio):
    N = 0.25 * (-ratio + np.sqrt(ratio) * np.sqrt(ratio + 8 * budget))
    M = (-ratio + np.sqrt(ratio) * np.sqrt(ratio + 8 * budget)) / (4 * ratio)
    if budget < 10 ** 3:
        return int(N + 0.5), int(M + 0.5)
    else:
        return int(N), int(M)


Ns = []
Ms = []

for budget in budgets:
    N, M = get_samples(budget, base_ratio)
    Ns.append(N)
    Ms.append(M)

from matplotlib import pyplot as plt

plt.plot(Ns, Ms, "o")
plt.xlabel("N")
plt.ylabel("M")
plt.grid()
plt.show()
