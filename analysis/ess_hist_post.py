import pickle
from common2 import *

options = linear_defaults.copy()
options["-design"] = "0.8"
options["-dim"] = 4
options["-poi"] = "1"

trials = 800

settings = [{"N": 320, "M": 100}, {"N": 320, "M": 1000}, {"N": 1000, "M": 1000}]

if True:
    data = {}
    columns = [
        "cESS_ML",
        "cESS_CL",
        "logML",
        "logCL",
        "logMLexact",
        "logCLexact",
        "priorLogDensities",
    ]
    for algo in ["prior", "mis"]:
        data[algo] = []
        for setting in settings:
            info = {}
            for trial in range(trials):
                prefix = "dump/%s_%dd/N%04d_M%04d_%04d_" % (
                    algo,
                    options["-dim"],
                    setting["N"],
                    setting["M"],
                    trial,
                )
                print(prefix)
                if trial == 0:
                    for column in columns:
                        info[column] = np.genfromtxt(prefix + column)
                else:
                    for column in columns:
                        info[column] = np.concatenate(
                            (info[column], np.genfromtxt(prefix + column))
                        )
            data[algo].append(info)
    pickle.dump(data, open("dump/ess.p", "wb"))

data = pickle.load(open("dump/ess.p", "rb"))
