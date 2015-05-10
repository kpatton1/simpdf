
import pysimpdf.cosmology

import healpy

import numpy

import matplotlib.pyplot as plt

c = pysimpdf.cosmology.Cosmology()

nside = 256

h = 0.69

x = []

y = []

for zi in numpy.arange(0.05,0.8,0.05):
    da = c.Da(zi)

    cell_size = 4.0 * numpy.pi / (12.0 * nside * nside)

    ang_size = numpy.sqrt(cell_size)

    dist = ang_size * da

    inv_dist = 1.0 / dist / h

    print zi, dist, inv_dist
    
    x.append(zi)
    y.append(inv_dist)


plt.plot(x, y)

plt.xlabel('z')
plt.ylabel('1/(cell width) [h/Mpc]')

plt.title('nside 256 redshift scaling')

plt.show()
