
import numpy
import matplotlib.pyplot as plt

infile = 'data/fiducial/cov_8192_fiducial_fiducial_fullcov.npz'

data = numpy.load(infile)

x = data['x']
mean = data['mean']
cov = data['cov']

r = data['r']



std = numpy.zeros(len(mean), dtype=numpy.float64)

cov = cov - numpy.outer(mean, mean)

for i in range(len(std)):
    std[i] = numpy.sqrt(cov[i][i])

for i in r:
    r0 = i[0]
    r1 = i[1]
    
    filter = range(r0,r1)
    
    plt.figure()

    plt.errorbar(x[filter], mean[filter], yerr=std[filter], fmt='.')

    plt.yscale('log')

plt.show()
