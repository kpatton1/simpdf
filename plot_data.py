
import numpy
import matplotlib.pyplot as plt

infile = 'data/fiducial/cov_8192_fiducial_fiducial_fullcovmoments.npz'

data = numpy.load(infile)

x = data['x']
mean = data['mean']
cov = data['cov']

r = data['r']

info  = data['i']


#corrcoef = numpy.zeros((len(mean),len(mean)), dtype=numpy.float64)

#for i in range(len(mean)):
#    for j in range(len(mean)):
#        if cov[i][j] != 0.0 and mean[i] >= 0.0/(200.0*384) and mean[j] >= 0.0/(200.0*384):
#            corrcoef[i][j] = cov[i][j] / numpy.sqrt(cov[i][i] * cov[j][j])
#        else:
#            corrcoef[i][j] = 0.0


#filter = (mean >= 10.0/(200.0*384))
#corrcoef = corrcoef[1168:1368,1568:1768]

#plt.pcolor(corrcoef)
#plt.colorbar()
#plt.show()

#exit(0)

std = numpy.zeros(len(mean), dtype=numpy.float64)

for i in range(len(std)):
    std[i] = numpy.sqrt(cov[i][i])

for i in range(len(r)):
    r0 = r[i][0]
    r1 = r[i][1]

    t = info[i][1]
    nside = info[i][0]
    
    filter = range(r0,r1)
    
    plt.figure()

    if t == 'ps':
        plt.title("Power Spectrum " + str(nside))
        plt.xlabel("l")
        plt.ylabel("C_l")
        plt.ylim(1e20,1e23)
        plt.xscale('log')
        plt.yscale('log')
    elif t == 'pdf':
        plt.title("PDF " + str(nside))
        plt.xlabel("Projected Mass (M_solar)")
        plt.ylabel("N")
        x[filter] = x[filter] - 1.0e16
        plt.xlim(-1.0e16,1.0e16)
        plt.yscale('log')
    elif t == 'moments':
        plt.title("Moments " + str(nside))
        plt.xlabel("n")
        plt.ylabel("<x^n>")
        plt.xlim(1,10)
        plt.yscale('log')

    filter2 = (mean[filter] >= 50.0/(200.0*384.0))

    plt.errorbar(x[filter][filter2], mean[filter][filter2], yerr=std[filter][filter2], fmt='.')
    plt.errorbar(x[filter][~filter2], mean[filter][~filter2], yerr=std[filter][~filter2], fmt='.',color='red')
    

plt.show()
