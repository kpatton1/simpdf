

import numpy

ranges = [5]

data = numpy.load('data/fiducial/cov_8192_fiducial_fiducial_fullcovmoments.npz')

dom_data = numpy.load('data/dcov_8192_delta_om_0.32_fiducial_fullcovmoments.npz')
ds8_data = numpy.load('data/dcov_8192_delta_s8_0.9_fiducial_fullcovmoments.npz')
dh_data = numpy.load('data/dcov_8192_delta_h_0.76_fiducial_fullcovmoments.npz')
dok_data = numpy.load('data/dcov_8192_delta_ok_0.1_fiducial_fullcovmoments.npz')
dw_data = numpy.load('data/dcov_8192_delta_w_-0.9_fiducial_fullcovmoments.npz')

dg_data = numpy.load('data/dcov_8192_fiducial_delta_g_1.1_fullcovmoments.npz')
dq_data = numpy.load('data/dcov_8192_fiducial_delta_q_1.1_fullcovmoments.npz')

mean = data['mean']
cov = data['cov']
x = data['x']

r = data['r']

dom_mean = dom_data['mean']
ds8_mean = ds8_data['mean']
dh_mean = dh_data['mean']
dok_mean = dok_data['mean']
dw_mean = dw_data['mean']

dg_mean = dg_data['mean']
dq_mean = dq_data['mean']

filter1 = []

for n in ranges:
    r1 = r[n][0]
    r2 = r[n][1]

    filter1.extend(range(r1,r2))

cov = cov[filter1][:,filter1]
mean = mean[filter1]
x = x[filter1]
dom_mean = dom_mean[filter1]
ds8_mean = ds8_mean[filter1]
dh_mean = dh_mean[filter1]
dok_mean = dok_mean[filter1]
dw_mean = dw_mean[filter1]
dg_mean = dg_mean[filter1]
dq_mean = dq_mean[filter1]

numpy.savez('cov_pdf256.npz', cov=cov, mean=mean, x=x, dmean_dom=dom_mean, dmean_ds8=ds8_mean, dmean_dh=dh_mean, dmean_dok=dok_mean, dmean_dw=dw_mean, dmean_dg=dg_mean, dmean_dq=dq_mean)
