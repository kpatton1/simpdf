

import numpy
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.patches import Ellipse

fisher0 = 'data/F_8192_fullcov_1357.npz'

fisher1 = 'data/F_8192_fullcov_35.npz'
basetitle = 'PDF resolution dependence'

fisher2 = 'data/F_8192_fullcov_1.npz'

fisher3 = 'data/F_8192_fullcov_3.npz'

fisher4 = 'data/F_8192_fullcov_5.npz'

fisher5 = 'data/F_8192_fullcov_7.npz'



#fisher = 'data/F_8192_fullcov_7.npz'
#basetitle = 'PDF (nside 1024)'

#fisher = 'data/F_8192_fullcov_3.npz'
#basetitle = 'PDF (nside 4096)'

data = numpy.load(fisher0)
F0 = data['F']

data = numpy.load(fisher1)
F1 = data['F']
params = data['params']

data = numpy.load(fisher2)
F2 = data['F']

data = numpy.load(fisher3)
F3 = data['F']

data = numpy.load(fisher4)
F4 = data['F']

data = numpy.load(fisher5)
F5 = data['F']

n = len(params)

names = params

centers = {'h':0.69,'ok':0.0, 'w':-1.0, 's8':0.82, 'om':0.29, 'q':1.0}

#print F4

F0[-1,-1] += 1/(0.1*0.1)
F1[-1,-1] += 1/(0.1*0.1)
F2[-1,-1] += 1/(0.1*0.1)
F3[-1,-1] += 1/(0.1*0.1)
F4[-1,-1] += 1/(0.1*0.1)
F5[-1,-1] += 1/(0.1*0.1)

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
        
        plot_fisher(F0, basetitle, i, j, 'black', False)
        plt.plot([],color='black',label='PDF 4096+1024+256+64')

        plot_fisher(F1, basetitle, i,j,'blue',False)

        plt.plot([],color='blue',label='PDF 1024+256')

        plot_fisher(F2, basetitle, i, j, 'green',False)

        plt.plot([],color='green',label='PDF 4096')

        plot_fisher(F3, basetitle, i, j, 'red', False)
            
        plt.plot([],color='red',label='PDF 1024')

        plot_fisher(F4, basetitle, i, j, 'purple', False)

        plt.plot([],color='purple',label='PDF 256')

        plot_fisher(F5, basetitle, i, j, 'orange', False)

        plt.plot([],color='orange',label='PDF 64')
        
        plt.legend()
        
        plt.xlim([0.17,0.42])
        plt.ylim([0.70,0.92])

        plt.show()
