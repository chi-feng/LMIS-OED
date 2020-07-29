from mossbauer import *

N = 10
M = 1000000
data = get_eig(
    "expdesign/out/mossbauer_baseline_%s_%s_%05d_%05d" % ("center", "prior", N, M)
)

print(mean(data), var(data))
