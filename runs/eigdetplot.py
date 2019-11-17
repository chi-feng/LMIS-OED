import numpy as np
import os, sys
from common import *

dim, sigma = 4, 0.4

npts = 201
designs = np.linspace(0, 1, npts)

columns = ['designs','Umarg','Ujoint','detGpost','detg11']
mat = np.zeros((len(designs),5))
mat[:,0] = designs;
mat[:,1] = eig_marginal(dim, designs, sigma);
mat[:,2] = eig_joint(designs);
mat[:,3] = np.exp(-2*mat[:,2])
mat[:,4] = np.exp(-2*mat[:,1])

np.savetxt("plot_data/eigdet.csv", mat, delimiter=",", header=','.join(columns),comments='')