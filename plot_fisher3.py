

import numpy
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.patches import Ellipse

fisher = 'data/F_8192_fullcov_0.npz'
basetitle = 'PS vs PDF'

fisher2 = 'data/F_8192_fullcov_7.npz'

fisher3 = 'data/F_8192_fullcov_07.npz'

fisher4 = 'data/F_8192_fullcov_0372.npz'




#fisher = 'data/F_8192_fullcov_7.npz'
#basetitle = 'PDF (nside 1024)'

#fisher = 'data/F_8192_fullcov_3.npz'
#basetitle = 'PDF (nside 4096)'

data = numpy.load(fisher)
F = data['F']
params = data['params']

data = numpy.load(fisher2)
F2 = data['F']

data = numpy.load(fisher3)
F3 = data['F']

data = numpy.load(fisher4)
F4 = data['F']

n = len(params)

names = params

centers = {'h':0.69,'ok':0.0, 'w':-1.0, 's8':0.82, 'om':0.29, 'q':1.0}

#print F4

F[-1,-1] += 1/(0.1*0.1)
F2[-1,-1] += 1/(0.1*0.1)
F3[-1,-1] += 1/(0.1*0.1)
F4[-1,-1] += 1/(0.1*0.1)

def plot_fisher(F, basetitle, i,j, color, rescale=False):
    
    F_inv = numpy.linalg.inv(F)

    ixgrid = numpy.ix_([j,i],[j,i])

    subM = F_inv[ixgrid]

    print subM

    subF = numpy.linalg.inv(subM)
    
    l, v = numpy.linalg.eig(subM)
    l = numpy.sqrt(l)

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

        plt.plot([],color='blue',label='PS 256')

        plot_fisher(F2, basetitle, i, j, 'green',False)

        plt.plot([],color='green',label='PDF 1024')

        plot_fisher(F3, basetitle, i, j, 'red', False)
            
        plt.plot([],color='red',label='PS 256 + PDF 1024')

        plot_fisher(F4, basetitle, i, j, 'purple', False)

        plt.plot([],color='purple',label='PS 256 + PDF 4096,1024,256')
        
        plt.legend()
        
        plt.xlim([0.21,0.37])
        plt.ylim([0.75,0.89])

        plt.show()
