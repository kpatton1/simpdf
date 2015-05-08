
import numpy
import matplotlib.pyplot as plt
import healpy

divs=384

nside=256

map_data_raw = healpy.fitsfunc.read_map('data_real/fiducial/r1/healpix_1024_fiducial_r1.fits')
map_data_noise = healpy.fitsfunc.read_map('data_real/fiducial/r1/noise_1024_fiducial_fiducial_r1n1.fits')

if nside != 1024:
    map_data_raw = healpy.pixelfunc.ud_grade(map_data_raw, nside, power=0, order_in='NESTED', order_out='NESTED', dtype=numpy.float32)
    map_data_noise = healpy.pixelfunc.ud_grade(map_data_noise, nside, power=0, order_in='NESTED', order_out='NESTED', dtype=numpy.float32)

map_data_meansubtract = numpy.zeros(len(map_data_noise), dtype=numpy.float32)

npix = len(map_data_noise)

for i in range(divs):
    map_data_meansubtract[npix/divs*i:npix/divs*(i+1)] = map_data_noise[npix/divs*i:npix/divs*(i+1)] - numpy.mean(map_data_noise[npix/divs*i:npix/divs*(i+1)])

std = numpy.std(map_data_meansubtract)

print std

map_data_lograw = numpy.log(map_data_raw)

map_data_log = numpy.log(map_data_meansubtract)

#map_data_log2 = numpy.log(1.0 + map_data_meansubtract*10)/10

map_data_log2 = numpy.log(-map_data_meansubtract)


b = 800
r = (-20.0,20.0)

hist_raw,bin_edges = numpy.histogram(map_data_raw, bins=b, range=r)
hist_noise,bin_edges = numpy.histogram(map_data_noise, bins=b, range=r)
hist_meansubtract,bin_edges = numpy.histogram(map_data_meansubtract, bins=b, range=r)
hist_log,bin_edges = numpy.histogram(map_data_log, bins=b, range=r)
hist_log2,bin_edges = numpy.histogram(map_data_log2, bins=b, range=r)
hist_lograw,bin_edges = numpy.histogram(map_data_lograw, bins=b, range=r)

x=0.5*(bin_edges[:-1] + bin_edges[1:])



plt.plot(x,hist_raw, label='raw')
plt.plot(x,hist_lograw, label='log(raw)')
plt.plot(x,hist_noise, label='noise')
plt.plot(x,hist_meansubtract, label='meansubtract')
plt.plot(x,hist_log, label='log(K)')
plt.plot(x,hist_log2, label='log(-K)')

plt.yscale('log')

plt.xlabel("Convergence")
plt.ylabel("N")

plt.title("Convergence PDF " + str(nside))

plt.legend()

plt.show()


