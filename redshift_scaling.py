
import pysimpdf.cosmology

import healpy

import numpy

c = pysimpdf.cosmology.Cosmology()

nside = 256

h = 0.69

for zi in numpy.arange(0.0,0.8,0.05):
    da = c.Da(zi)

    cell_size = 4.0 * numpy.pi / (12.0 * nside * nside)

    ang_size = numpy.sqrt(cell_size)

    dist = ang_size * da

    inv_dist = 1.0 / dist / h

    print zi, dist, inv_dist

