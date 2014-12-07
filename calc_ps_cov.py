
import healpy
import numpy
import matplotlib.pyplot as plt

infile = 'data/fiducial/r1/noise_8192_fiducial_fiducial_r1n1.fits'
outfile = 'test_cls.npz'

nside = 1024

divs = 384

map_data_orig = healpy.fitsfunc.read_map(infile,dtype=numpy.float32)

print 'read in map!'

map_data = healpy.pixelfunc.ud_grade(map_data_orig, nside, order_in='NESTED', order_out='NESTED')

mdata_1024 = healpy.pixelfunc.ud_grade(map_data_orig, 1024, order_in='NESTED', order_out='NESTED')
mdata_512 = healpy.pixelfunc.ud_grade(map_data_orig, 512, order_in='NESTED', order_out='NESTED')
mdata_256 = healpy.pixelfunc.ud_grade(map_data_orig, 256, order_in='NESTED', order_out='NESTED')
mdata_128 = healpy.pixelfunc.ud_grade(map_data_orig, 128, order_in='NESTED', order_out='NESTED')

npix_1024 = healpy.pixelfunc.nside2npix(1024)
npix_512 = healpy.pixelfunc.nside2npix(512)
npix_256 = healpy.pixelfunc.nside2npix(256)
npix_128 = healpy.pixelfunc.nside2npix(128)

print 'downgraded map!'

npix = healpy.pixelfunc.nside2npix(nside)
r = numpy.arange(npix)

nd = npix / divs

ncls = 3*nside

mean = numpy.zeros(ncls, dtype=numpy.float32)
cov = numpy.zeros((ncls,ncls), dtype=numpy.float32)

for n in range(divs):

    print 'div ' + str(n)

    mdata_survey_1024 = numpy.zeros(npix_1024, dtype=numpy.float32)
#    mdata_survey_512 = numpy.zeros(npix_512, dtype=numpy.float32)
#    mdata_survey_256 = numpy.zeros(npix_256, dtype=numpy.float32)
#    mdata_survey_128 = numpy.zeros(npix_128, dtype=numpy.float32)

#    mdata_survey_1024 = numpy.full(npix_1024, healpy.pixelfunc.UNSEEN, dtype=numpy.float32)
#    mdata_survey_512 = numpy.full(npix_512, healpy.pixelfunc.UNSEEN, dtype=numpy.float32)
#    mdata_survey_256 = numpy.full(npix_256, healpy.pixelfunc.UNSEEN, dtype=numpy.float32)
#    mdata_survey_128 = numpy.full(npix_128, healpy.pixelfunc.UNSEEN, dtype=numpy.float32)

#    mdata_survey_1024[npix_1024/divs*n:npix_1024/divs*(n+1)] = mdata_1024[npix_1024/divs*n:npix_1024/divs*(n+1)] - numpy.mean(mdata_1024[npix_1024/divs*n:npix_1024/divs*(n+1)])
#    mdata_survey_512[npix_512/divs*n:npix_512/divs*(n+1)] = mdata_512[npix_512/divs*n:npix_512/divs*(n+1)] - numpy.mean(mdata_512[npix_512/divs*n:npix_512/divs*(n+1)])
#    mdata_survey_256[npix_256/divs*n:npix_256/divs*(n+1)] = mdata_256[npix_256/divs*n:npix_256/divs*(n+1)] - numpy.mean(mdata_256[npix_256/divs*n:npix_256/divs*(n+1)])
#    mdata_survey_128[npix_128/divs*n:npix_128/divs*(n+1)] = mdata_128[npix_128/divs*n:npix_128/divs*(n+1)] - numpy.mean(mdata_128[npix_128/divs*n:npix_128/divs*(n+1)])

    mdata_survey_1024[0:npix_1024/divs] = mdata_1024[npix_1024/divs*n:npix_1024/divs*(n+1)] - numpy.mean(mdata_1024[npix_1024/divs*n:npix_1024/divs*(n+1)])
#    mdata_survey_512[0:npix_512/divs] = mdata_512[npix_512/divs*n:npix_512/divs*(n+1)] - numpy.mean(mdata_512[npix_512/divs*n:npix_512/divs*(n+1)])
#    mdata_survey_256[0:npix_256/divs] = mdata_256[npix_256/divs*n:npix_256/divs*(n+1)] - numpy.mean(mdata_256[npix_256/divs*n:npix_256/divs*(n+1)])
#    mdata_survey_128[0:npix_128/divs] = mdata_128[npix_128/divs*n:npix_128/divs*(n+1)] - numpy.mean(mdata_128[npix_128/divs*n:npix_128/divs*(n+1)])

#    mdata_survey_1024 = healpy.ma(mdata_survey_1024)
#    mdata_survey_512 = healpy.ma(mdata_survey_512)
#    mdata_survey_256 = healpy.ma(mdata_survey_256)
#    mdata_survey_128 = healpy.ma(mdata_survey_128)

    mdata_survey_1024 = healpy.pixelfunc.reorder(mdata_survey_1024, n2r = True)
#    mdata_survey_512 = healpy.pixelfunc.reorder(mdata_survey_512, n2r = True)
#    mdata_survey_256 = healpy.pixelfunc.reorder(mdata_survey_256, n2r = True)
#    mdata_survey_128 = healpy.pixelfunc.reorder(mdata_survey_128, n2r = True)

    cls_1024 = healpy.sphtfunc.anafast(mdata_survey_1024)
#    cls_512 = healpy.sphtfunc.anafast(mdata_survey_512)
#    cls_256 = healpy.sphtfunc.anafast(mdata_survey_256)
#    cls_128 = healpy.sphtfunc.anafast(mdata_survey_128)

    cls_1024[0] = cls_1024[1]
#    cls_512[0] = cls_512[1]
#    cls_256[0] = cls_256[1]
#    cls_128[0] = cls_128[1]

    for i in range(ncls):
        mean[i] = mean[i] + cls_1024[i]

    for i in range(ncls):
        for j in range(ncls):
            cov[i][j] = cov[i][j] + cls_1024[i] * cls_1024[j]


#    l_1024 = numpy.arange(len(cls_1024))
#    l_512 = numpy.arange(len(cls_512))
#    l_256 = numpy.arange(len(cls_256))
#    l_128 = numpy.arange(len(cls_128))

#    plt.plot(l_1024, l_1024*(l_1024+1)*cls_1024, label='1024')
#    plt.plot(l_512, l_512*(l_512+1)*cls_512, label='512')
#    plt.plot(l_256, l_256*(l_256+1)*cls_256, label='256')
#    plt.plot(l_128, l_128*(l_128+1)*cls_128, label='128')

#    plt.plot(l_1024, cls_1024, label=str(n))
#    plt.plot(l_512, cls_512, label='512')
#    plt.plot(l_256, cls_256, label='256')
#    plt.plot(l_128, cls_128, label='128')

#plt.yscale('log')
#plt.xscale('log')

#plt.ylim(1.0e14,1.0e16)

#plt.legend(loc='upper left')

#plt.show()

#    print 'calculated cls!'

#    mean += cls
#    for x in range(ncls):
#        for y in range(ncls):
#            cov[x][y] += cls[x] * cls[y]

mean = mean / float(divs)
cov = cov / float(divs)

std = numpy.zeros(ncls, dtype=numpy.float32)

for i in range(ncls):
    std[i] = numpy.sqrt(cov[i][i] - mean[i] * mean[i])

numpy.savez(outfile,mean=mean,cov=cov,std=std)
    
