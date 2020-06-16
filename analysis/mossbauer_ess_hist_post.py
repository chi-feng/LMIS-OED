import pickle
from common2 import *

options = mossbauer_defaults.copy()
options["-design"] = "-1.3,0.0,1.3"
options["-nugget"] = "0.001"

trials = 800

settings = [{"N": 1000, "M": 100}]

if True:
    data = {}
    columns = ["cESS_ML", "cESS_CL", "logML", "logCL", "priorLogDensities"]
    for algo in ["prior", "mis"]:
        data[algo] = []
        for setting in settings:
            info = {}
            for trial in range(trials):
                prefix = "mossbauer_ess/%s_%dd/N%04d_M%04d_%04d_" % (
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
    pickle.dump(data, open("mossbauer_ess/ess.p", "wb"))

data = pickle.load(open("mossbauer_ess/ess.p", "rb"))
