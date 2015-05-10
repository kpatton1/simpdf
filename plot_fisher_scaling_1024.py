
import numpy
import matplotlib.pyplot as plt

N = 20.0 * 384.0

def calc_fisher_fom(F):
    
    i = 2
    j = 3

    F_inv = numpy.linalg.pinv(F)
    
    ixgrid = numpy.ix_([j,i],[j,i])
    
    subM = F_inv[ixgrid]
    
    #print subM
    
    #sigma1 = numpy.sqrt(subM[0,0])
    #sigma2 = numpy.sqrt(subM[1,1])
    
    #subF = numpy.linalg.pinv(subM)
    
    l, v = numpy.linalg.eig(subM)
    l = numpy.sqrt(l)
    
    area = l[0]*l[1]*numpy.pi
    
    #print names[j] + ' sd: ' + str(sigma1)
    #print names[i] + ' sd: ' + str(sigma2)
    #print 'area: ' + str(area)
    #print 'fom: ' + str(1.0/area)

    return 1.0 / area

#file_list = []
#F_list = []
x = []
y = []

params = ['h', 'ok', 'om', 's8', 'w', 'q', 'g']

prior = numpy.zeros([7,7])

prior[5,5] = 1/(0.1*0.1)
prior[6,6] = 1/(0.05*0.05)

for i in range(0,16):
    data = numpy.load('data_real/F_1024_test_pdf1024_' + str(i) + '.npz')
    bins = 50.0 * (i+1)
    F = data['F']
    #F = (N-bins-2.0)/(N-2.0) * F
    F += prior
    fom = calc_fisher_fom(F)
    y.append(fom)
    x.append(bins)

for i in range(0,13):
    data = numpy.load('data_real/F_1024_test_pdf1024_2_' + str(i) + '.npz')
    bins = 800.0 + 100.0 * i
    F = data['F']
    #F = (N-bins-2.0)/(N-2.0) * F
    F += prior
    fom = calc_fisher_fom(F)
    y.append(fom)
    x.append(bins)

for i in range(0,5):
    data = numpy.load('data_real/F_1024_test_pdf1024_3_' + str(i) + '.npz')
    bins = 2000.0 + 200.0 * i
    if i == 4:
        bins = 3000
    F = data['F']
    #F = (N-bins-2.0)/(N-2.0) * F
    F += prior
    fom = calc_fisher_fom(F)
    y.append(fom)
    x.append(bins)

for i in range(0,2):
    data = numpy.load('data_real/F_1024_test_pdf1024_4_' + str(i) + '.npz')
    bins = 4000.0 + 1000.0 * i
    F = data['F']
    #F = (N-bins-2.0)/(N-2.0) * F
    F += prior
    fom = calc_fisher_fom(F)
    y.append(fom)
    x.append(bins)

'''
for file in file_list:
    data = numpy.load(file)
    F = data['F']
    F += prior
    F_list.append(F)

for F in F_list:
    area = calc_fisher_area(F)
    y_list.append(area)
'''

plt.xlim([0,5100])

plt.xlabel('bins')
plt.ylabel('fom')

plt.title('PDF 1024 resolution scaling')

plt.scatter(x,y)
plt.show()


