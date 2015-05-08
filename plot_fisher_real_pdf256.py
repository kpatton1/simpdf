

import numpy
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.patches import Ellipse

N = 20.0 * 384.0

fisher = 'data_real/F_1024_test_simple_5.npz'
basetitle = 'PDF 256 scaling'

fisher2 = 'data_real/F_1024_test_pdf256_2.npz'

fisher3 = 'data_real/F_1024_test_pdf256_3.npz'

fisher4 = 'data_real/F_1024_test_pdf256_2_3.npz'

fisher5 = 'data_real/F_1024_test_pdf256_4.npz'

#fisher = 'data/F_8192_fullcov_7.npz'
#basetitle = 'PDF (nside 1024)'

#fisher = 'data/F_8192_fullcov_3.npz'
#basetitle = 'PDF (nside 4096)'

data = numpy.load(fisher)
F = data['F'] #* (N - 400.0 - 2.0) / (N - 2.0)
params = data['params']

data = numpy.load(fisher2)
F2 = data['F'] #* (N - 800.0 - 2.0) / (N - 2.0)

data = numpy.load(fisher3)
F3 = data['F'] #* (N - 1600.0 - 2.0) / (N - 2.0)

data = numpy.load(fisher4)
F4 = data['F'] #* (N - 2400.0 - 2.0) / (N - 2.0)

data = numpy.load(fisher5)
F5 = data['F'] #* (N - 3200.0 - 2.0) / (N - 2.0)


n = len(params)

names = params

print names

centers = {'h':0.69,'ok':0.0, 'w':-1.0, 's8':0.82, 'om':0.29, 'q':1.0}

#print F4

F[-2,-2] += 1/(0.1*0.1)
F2[-2,-2] += 1/(0.1*0.1)
F3[-2,-2] += 1/(0.1*0.1)
F4[-2,-2] += 1/(0.1*0.1)
F5[-2,-2] += 1/(0.1*0.1)
F[-1,-1] += 1/(0.05*0.05)
F2[-1,-1] += 1/(0.05*0.05)
F3[-1,-1] += 1/(0.05*0.05)
F4[-1,-1] += 1/(0.05*0.05)
F5[-1,-1] += 1/(0.05*0.05)
#F[4,4] += 1/(0.0001*0.0001)
#F2[4,4] += 1/(0.0001*0.0001)
#F3[4,4] += 1/(0.0001*0.0001)
#F4[4,4] += 1/(0.0001*0.0001)
#F5[4,4] += 1/(0.0001*0.0001)
#F[1,1] += 1/(0.0001*0.0001)
#F2[1,1] += 1/(0.0001*0.0001)
#F3[1,1] += 1/(0.0001*0.0001)
#F4[1,1] += 1/(0.0001*0.0001)
#F5[1,1] += 1/(0.0001*0.0001)
#F[0,0] += 1/(0.1*0.1)
#F2[0,0] += 1/(0.1*0.1)
#F3[0,0] += 1/(0.1*0.1)
#F4[0,0] += 1/(0.1*0.1)
#F5[0,0] += 1/(0.1*0.1)




def plot_fisher(F, basetitle, i,j, color, rescale=False):
    
    F_inv = numpy.linalg.pinv(F)

    ixgrid = numpy.ix_([j,i],[j,i])

    subM = F_inv[ixgrid]

    print subM
    
    sigma1 = numpy.sqrt(subM[0,0])
    sigma2 = numpy.sqrt(subM[1,1])

    subF = numpy.linalg.pinv(subM)
    
    l, v = numpy.linalg.eig(subM)
    l = numpy.sqrt(l)
    
    area = l[0]*l[1]*numpy.pi
    
    print names[j] + ' sd: ' + str(sigma1)
    print names[i] + ' sd: ' + str(sigma2)
    print 'area: ' + str(area)
    print 'fom: ' + str(1.0/area)

    ax = plt.subplot(111)
    ax.set_title(basetitle + ': ' + names[j] + ' vs ' + names[i])
    
    ax.set_xlabel(names[j])
    ax.set_ylabel(names[i])
    
    center_i = centers[params[i]]
    center_j = centers[params[j]]
    
    ell = Ellipse(xy=(center_j,center_i),width=l[0]*2, height=l[1]*2, angle=numpy.rad2deg(-numpy.arccos(v[0,0])))
    ell.set_facecolor('none')
    ell.set_edgecolor(color)
    ax.add_artist(ell)

    if rescale:
        plt.xlim([center_j-numpy.sqrt(subM[0,0])*1.2,center_j+numpy.sqrt(subM[0,0])*1.2])
        plt.ylim([center_i-numpy.sqrt(subM[1,1])*1.2,center_i+numpy.sqrt(subM[1,1])*1.2])
        plt.scatter([center_j],[center_i])

for i in range(n):
    for j in range(i):

        if names[i] != 's8':
            continue
        if names[j] != 'om':
            continue

        handles = []

        plot_fisher(F, basetitle, i,j,'blue',False)
        plt.plot([],color='blue',label='400 bins')

        plot_fisher(F2, basetitle, i, j, 'green',True)
        plt.plot([],color='green',label='800 bins')

        plot_fisher(F3, basetitle, i, j, 'red', False)
        plt.plot([],color='red',label='1600 bins')
        
        plot_fisher(F4,basetitle,i,j, 'black',False)
        plt.plot([],color='black',label='2400 bins')

        plot_fisher(F5, basetitle, i, j, 'purple', False)
        plt.plot([],color='purple',label='3200 bins')
        

        
        plt.legend()
        
        plt.xlim([0.15,0.43])
        plt.ylim([0.61,1.03])

        plt.show()
