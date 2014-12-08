
import numpy
import matplotlib.pyplot as plt

infile = 'data/fiducial/cov_8192_fiducial_fiducial_64.npz'

data = numpy.load(infile)

x = data['x']
mean = data['mean']
cov = data['cov']

std = numpy.zeros(len(mean))

for i in range(len(mean)):
    for j in range(len(mean)):
        cov[i][j] -= mean[i] * mean[j]

for i in range(len(mean)):
    std[i] = numpy.sqrt(cov[i][i])

plt.errorbar(x, mean, yerr=std, fmt='.')

plt.yscale('log')

plt.show()
