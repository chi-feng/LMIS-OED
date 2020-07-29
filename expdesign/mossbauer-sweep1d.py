from expdesign import *

algorithms = ["prior", "mis"]

designs = []
for d in np.linspace(0, 2, 51):
    designs.append("%f,0.0,%f" % (d, -d))

f = open("expdesign/out/mossbauer3d_designs1d", "w")
for design in designs:
    f.write("%s\n" % design)
f.close()

trials = 1
for arg in sys.argv:
    if "--trials=" in arg:
        trials = int(arg[9:])

options = default.copy()
options["-experiment"] = "mossbauer"
options["-dim"] = "3"
options["-dof"] = "4"
options["-N"] = "1000"
options["-M1"] = "100"
options["-M2"] = "100"
options["-maxComponents"] = "100"
options["-sigeps"] = "0.2"

suffix = ""
if "joint" in sys.argv:
    suffix = "_joint"
    options["-poi"] = 4

if "height" in sys.argv:
    suffix = "_height"
    options["-index"] = 3

runs = []


def get_run(i, algo):
    options["-design"] = designs[i]
    if algo is "prior":
        options["-useIS"], options["-useMIS"] = 0, 0
    if algo is "mis":
        options["-useIS"], options["-useMIS"] = 0, 1
    filename = "expdesign/out/mossbauer3d_sweep1d%s_%s_%04d" % (suffix, algo, i)
    return make_run(filename, options)


for algo in algorithms:
    for i in range(len(designs)):
        run = get_run(i, algo)
        run["mode"] = "a"
        for j in range(trials):
            runs.append(run)

random.shuffle(runs)

print("%d runs" % len(runs))

execute_all(runs)
