
import numpy
import matplotlib.pyplot as plt

#ranges = [0,3]
ranges = [0]

data = numpy.load('data/fiducial/cov_8192_fiducial_fiducial_fullcov.npz')

dom_data = numpy.load('data/delta_om_0.32/dcov_8192_delta_om_0.32_fiducial_fullcov.npz')
ds8_data = numpy.load('data/delta_s8_0.9/dcov_8192_delta_s8_0.9_fiducial_fullcov.npz')

mean = data['mean']
cov = data['cov']

r = data['r']

#cov = cov - numpy.outer(mean, mean)

filter1 = []

dmean = dom_data['mean']

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
dmean = dmean[filter1]
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

for i in range(len(dmean)):
    if cov[i][i] != 0.0:
        filter2.append(i)

cov = cov[filter2][:,filter2]
dmean = dmean[filter2]

sd = numpy.zeros(len(cov),dtype=numpy.float64)

for i in range(len(cov)):
    sd[i] = numpy.sqrt(cov[i][i])
    dmean[i] = dmean[i] / sd[i]

for i in range(len(cov)):
    for j in range(len(cov)):
        cov[i][j] = cov[i][j] / sd[i] / sd[j]

trace = numpy.trace(cov)

shrink = 1.0e-11
#shrink = 0.0

cov = cov + shrink * numpy.identity(len(cov),dtype=numpy.float64)

#cov = (1.0 - shrink) * cov + shrink * numpy.identity(len(cov),dtype=numpy.float64)

s = numpy.shape(cov)

e,v = numpy.linalg.eigh(cov)

for i in range(len(e)):
    temp = numpy.dot(dmean,v[i])
    e[i] = temp * temp / e[i]

for i in e:
    print i

plt.figure()

plt.plot(range(len(e)), e, color='blue')
plt.plot(range(len(e)), -e, color='red')
#plt.plot(range(len(e)),numpy.ones(len(e))*1.0e25)
plt.yscale('log')
plt.ylabel('magnitude')
plt.xlabel('eigenvalue')
plt.title('covariance eigenvalues')

tot = numpy.sum(e)

print 'sd sigma_m: ' + str(numpy.sqrt(1.0/tot))

plt.show()
