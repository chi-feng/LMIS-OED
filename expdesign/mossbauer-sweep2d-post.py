from expdesign import *

algorithms = ["prior", "mis"]

designs = []
for d1 in np.linspace(-2, 2, 21):
    for d2 in np.linspace(-2, 2, 21):
        designs.append("%0.1f,0.0,%0.1f" % (d1, d2))

f = open("expdesign/out/mossbauer3d_designs", "w")
for design in designs:
    f.write("%s\n" % design)
f.close()

options = default.copy()
options["--experiment"] = "mossbauer"
options["--index"] = "0"
options["--dim"] = "3"
options["--dof"] = "4"
options["--N"] = "1000"
options["--M1"] = "200"
options["--M2"] = "200"
options["--maxComponents"] = "100"
options["--sigeps"] = "0.2"

trials = 1
for arg in sys.argv:
    if "--trials=" in arg:
        trials = int(arg[9:])

suffix = ""
if "joint" in sys.argv:
    suffix = "_joint"
    options["--poi"] = "4"

if "height" in sys.argv:
    suffix = "_height"
    options["--index"] = "3"

runs = []


def get_run(i, algo):
    options["--design"] = designs[i]
    if algo is "prior":
        options["--useIS"], options["--useMIS"] = 0, 0
    if algo is "mis":
        options["--useIS"], options["--useMIS"] = 0, 1
    filename = "expdesign/out/mossbauer3d_sweep%s_%s_%04d" % (suffix, algo, i)
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
