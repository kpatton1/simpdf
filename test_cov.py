
import numpy
import matplotlib.pyplot as plt

ranges = [5]
#ranges = [0]

data = numpy.load('data/fiducial/cov_8192_fiducial_fiducial_fullcov.npz')

dom_data = numpy.load('data/dcov_8192_delta_om_0.32_fiducial_fullcov.npz')
ds8_data = numpy.load('data/dcov_8192_delta_s8_0.9_fiducial_fullcov.npz')
dg_data = numpy.load('data/dcov_8192_fiducial_delta_g_1.1_fullcov.npz')
dq_data = numpy.load('data/dcov_8192_fiducial_delta_q_1.1_fullcov.npz')

mean = data['mean']
cov = data['cov']
x = data['x']

r = data['r']

#cov = cov - numpy.outer(mean, mean)

filter1 = []

dom_mean = dom_data['mean']
ds8_mean = ds8_data['mean']
dg_mean = dg_data['mean']
dq_mean = dq_data['mean']

s = numpy.shape(cov)

prefilter = []
'''
for i in range(len(dmean)):
    if cov[i][i] != 0.0:
        prefilter.append(i)

precov = cov[prefilter][:,prefilter]
predmean = dmean[prefilter]

for i in range(len(precov)):
    precov[i][i] = (1.0 + 0.00001) * precov[i][i]

e,v = numpy.linalg.eigh(precov)

for i in range(len(e)):
    temp = numpy.dot(predmean,v[i])
    e[i] = temp * temp / e[i]

plt.plot(range(len(e)), e, color='blue')
plt.plot(range(len(e)), -e, color='red')
#plt.plot(range(len(e)),numpy.ones(len(e))*1.0e25)
plt.yscale('log')
'''
for n in ranges:
    r1 = r[n][0]
    r2 = r[n][1]

    filter1.extend(range(r1,r2))

cov = cov[filter1][:,filter1]
mean = mean[filter1]
x=x[filter1]
dom_mean = dom_mean[filter1]
ds8_mean = ds8_mean[filter1]
dg_mean = dg_mean[filter1]
dq_mean = dq_mean[filter1]
'''
s = numpy.shape(cov)

e,v = numpy.linalg.eigh(cov)

for i in range(len(e)):
    temp = numpy.dot(dmean,v[i])
    e[i] = temp * temp / e[i]

plt.figure()

plt.plot(range(len(e)), e, color='blue')
plt.plot(range(len(e)), -e, color='red')
#plt.plot(range(len(e)),numpy.ones(len(e))*1.0e25)
plt.yscale('log')
'''
filter2 = []

for i in range(len(dom_mean)):
    if cov[i][i] != 0.0 and mean[i] >= 50.0/(200.0*384):
        filter2.append(i)

cov = cov[filter2][:,filter2]
mean = mean[filter2]
x=x[filter2]
dom_mean = dom_mean[filter2]
ds8_mean = ds8_mean[filter2]
dg_mean = dg_mean[filter2]
dq_mean = dq_mean[filter2]

sd = numpy.zeros(len(cov),dtype=numpy.float64)

for i in range(len(cov)):
    sd[i] = numpy.sqrt(cov[i][i])
    dom_mean[i] = dom_mean[i] / sd[i]
    ds8_mean[i] = ds8_mean[i] / sd[i]
    dq_mean[i] = dq_mean[i] / sd[i]
    dg_mean[i] = dg_mean[i] / sd[i]

#for i in range(len(cov)):
#    dom_mean[i] = dom_mean[i] / mean[i]
#    ds8_mean[i] = ds8_mean[i] / mean[i]
#    dq_mean[i] = dq_mean[i] / mean[i]
#    dg_mean[i] = dg_mean[i] / mean[i]

#dom_tot = numpy.sqrt(numpy.inner(dom_mean,dom_mean))
#ds8_tot = numpy.sqrt(numpy.inner(ds8_mean,ds8_mean))
#dq_tot = numpy.sqrt(numpy.inner(dq_mean,dq_mean))
#dg_tot = numpy.sqrt(numpy.inner(dg_mean,dg_mean))

#for i in range(len(cov)):
#    dom_mean[i] = dom_mean[i] / dom_tot
#    ds8_mean[i] = ds8_mean[i] / ds8_tot
#    dq_mean[i] = dq_mean[i] / dq_tot
#    dg_mean[i] = dg_mean[i] / dg_tot

#plt.figure()
#plt.errorbar(x, numpy.zeros(len(mean)), yerr=sd, color='black',label='mean sd')
#plt.plot(x, dom_mean, color='blue',label='delta om')
#plt.plot(x, ds8_mean, color='green',label='delta s8')
#plt.plot(x, dq_mean, color='red',label='delta q')
#plt.plot(x, dg_mean, color='orange',label='delta g')

#plt.legend(loc=0)

#plt.ylabel('dData/dParam / PDF bin sd')
#plt.xlabel('PDF bin')

#plt.xscale('log')

#plt.figure()
#plt.errorbar(x, mean, yerr=sd)
#plt.xscale('log')
#plt.yscale('log')

#plt.show()

for i in range(len(cov)):
    for j in range(len(cov)):
        cov[i][j] = cov[i][j] / sd[i] / sd[j]

trace = numpy.trace(cov)

shrink = 1e-13
#shrink = 0.0

cov = cov + shrink * numpy.identity(len(cov),dtype=numpy.float64)

#cov = (1.0 - shrink) * cov + shrink * numpy.identity(len(cov),dtype=numpy.float64)

s = numpy.shape(cov)

cov_inv = numpy.linalg.inv(cov)

e,v = numpy.linalg.eigh(cov)

plt.figure()

plt.plot(range(len(e)), e, color='blue')
plt.plot(range(len(e)), -e, color='red')
#plt.plot(range(len(e)),numpy.ones(len(e))*1.0e25)
plt.yscale('log')
plt.ylabel('magnitude')
plt.xlabel('eigenvalue')
plt.title('covariance eigenvalues')

for i in range(len(e)):
    temp = numpy.dot(ds8_mean,v[i])
    e[i] = temp * temp / e[i]

tot = numpy.sum(e[10:])

print 'sd eigen: ' + str(numpy.sqrt(1.0/tot))
print 'sd inv: ' + str(numpy.sqrt(1.0/numpy.dot(ds8_mean, numpy.dot(cov_inv, ds8_mean))))

plt.show()

for i in range(0,10): #range(len(v)):
    print i, e[i]
    plt.plot(range(len(v[i])), v[i]*ds8_mean)
plt.show()
