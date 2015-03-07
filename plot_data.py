
import numpy
import matplotlib.pyplot as plt

infile = 'data/fiducial/cov_8192_fiducial_fiducial_fullcov.npz'

data = numpy.load(infile)

x = data['x']
mean = data['mean']
cov = data['cov']

r = data['r']
"""
corrcoef = numpy.zeros((len(mean),len(mean)), dtype=numpy.float64)

for i in range(len(mean)):
    for j in range(len(mean)):
        if cov[i][j] != 0.0:
            corrcoef[i][j] = cov[i][j] / numpy.sqrt(cov[i][i] * cov[j][j])
        else:
            corrcoef[i][j] = 0.0



plt.pcolor(corrcoef)
plt.colorbar()
"""
std = numpy.zeros(len(mean), dtype=numpy.float64)

for i in range(len(std)):
    std[i] = numpy.sqrt(cov[i][i])

for i in range(len(r)):
    r0 = r[i][0]
    r1 = r[i][1]
    
    filter = range(r0,r1)
    
    plt.figure()

    if i == 0:
        plt.title("Power Spectrum")
        plt.xlabel("l")
        plt.ylabel("C_l")
    else:
        plt.title("PDF")
        plt.xlabel("Projected Mass (M_solar)")
        plt.ylabel("N")

    plt.errorbar(x[filter], mean[filter], yerr=std[filter], fmt='.')

    if i == 0:
        plt.xscale('log')

    plt.yscale('log')

    

plt.show()
