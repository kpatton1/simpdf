
import numpy
import matplotlib.pyplot as plt

#ranges = [0]
#title = 'PS (nside 256)'

#ranges = [5]
#title = 'PDF (nside 256)'

#ranges = [3]
#title = 'PDF (nside 1024)'
#ranges = [10]
#title = 'Moments (nside 1024)'

#ranges = [10]
#title = 'Moments (nside 1024)'

ranges = [3]
title = 'PDF (nside 256)'

data = numpy.load('data_real/fiducial/cov_1024_fiducial_fiducial_test_pdf256.npz')

dom_data = numpy.load('data_real/delta_om_0.32/cov_1024_delta_om_0.32_fiducial_test_pdf256.npz')
ds8_data = numpy.load('data_real/delta_s8_0.9/cov_1024_delta_s8_0.9_fiducial_test_pdf256.npz')
#dg_data = numpy.load('data_real/fiducial/cov_1024_fiducial_fiducial_test_pdf256.npz')
#dq_data = numpy.load('data_real/fiducial/cov_1024_fiducial_fiducial_test_pdf256.npz')

#delta_om = 0.32 - 0.29
#delta_s8 = 0.90 - 0.82
#delta_g = 1.1 - 1.0
#delta_q = 1.1 - 1.0

mean = data['mean']
cov = data['cov']
x = data['x']

r = data['r']

#cov = cov - numpy.outer(mean, mean)

filter1 = []

dom_mean = dom_data['mean']
ds8_mean = ds8_data['mean']
#dg_mean = dg_data['mean']
#dq_mean = dq_data['mean']

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
#dg_mean = dg_mean[filter1]
#dq_mean = dq_mean[filter1]
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
    if cov[i][i] != 0.0 and mean[i] >= 100.0/(20.0*384):
        filter2.append(i)

cov = cov[filter2][:,filter2]
mean = mean[filter2]
x=x[filter2]
dom_mean = dom_mean[filter2]
ds8_mean = ds8_mean[filter2]
#dg_mean = dg_mean[filter2]
#dq_mean = dq_mean[filter2]


#dom_mean = dom_mean / mean * 0.32 
#ds8_mean = ds8_mean / mean * 0.82
#dg_mean = dg_mean / mean * 1.0
#dq_mean = dq_mean / mean * 1.0

sd = numpy.zeros(len(cov),dtype=numpy.float64)

for i in range(len(cov)):
    sd[i] = numpy.sqrt(cov[i][i])
    #dom_mean[i] = dom_mean[i] 
    #ds8_mean[i] = ds8_mean[i] 
    #dq_mean[i] = dq_mean[i] 
    #dg_mean[i] = dg_mean[i]

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


plt.figure()
plt.errorbar(x, mean, yerr=sd, color='black',label='fiducial')
plt.plot(x, dom_mean, color='green',label='delta om')
plt.plot(x, ds8_mean, color='blue',label='delta s8')
#plt.scatter(x, dq_mean, color='red',label='delta q')
#plt.scatter(x, -dq_mean, color='red',label='-delta q', marker='x')
#plt.scatter(x, dg_mean, color='orange',label='delta g')

plt.legend()

plt.xlabel('Convergence')
plt.ylabel('N')
plt.title(title)

plt.yscale('log')
#plt.ylim([1e-2,1e2])
plt.xlim(-0.03,0.13)
plt.ylim([1e-2,5.0e2])

#plt.xscale('log')

#plt.figure()
#plt.errorbar(x, mean, yerr=sd)
#plt.xscale('log')
#plt.yscale('log')

plt.show()

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
    temp = numpy.dot(dq_mean,v[i])
    e[i] = temp * temp / e[i]

tot = numpy.sum(e)

print 'sd eigen: ' + str(numpy.sqrt(1.0/tot))
print 'sd inv: ' + str(numpy.sqrt(1.0/numpy.dot(dq_mean, numpy.dot(cov_inv, dq_mean))))

plt.show()

for i in range(0,10): #range(len(v)):
    print i, e[i]
    plt.plot(range(len(v[i])), v[i]*dq_mean)
plt.show()
