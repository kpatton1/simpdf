

import numpy
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.patches import Ellipse

fisher = 'data/F_8192_fullcov_0.npz'
basetitle = 'PS (nside 256)'

#fisher = 'data/F_8192_fullcov_7.npz'
#basetitle = 'PDF (nside 1024)'

#fisher = 'data/F_8192_fullcov_3.npz'
#basetitle = 'PDF (nside 4096)'

data = numpy.load(fisher)

F = data['F']
F_inv = data['F_inv']
params = data['params']

n = len(params)

names = params

centers = {'h':0.69,'ok':0.0, 'w':-1.0, 's8':0.82, 'om':0.29, 'q':1.0}

print F

F_noq = F[:-1,:-1]

F_priorq = F.copy()
F_priorq[-1,-1] += 1/(0.1*0.1)

for p in params:
    print p

def plot_fisher(F, basetitle, color, rescale=False):
    
    F_inv = numpy.linalg.inv(F)

    ixgrid = numpy.ix_([j,i],[j,i])

    subM = F_inv[ixgrid]

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

        plot_fisher(F, basetitle, 'blue',True)

        plt.plot([],color='blue',label='Marginalize over q')

        plot_fisher(F_priorq, basetitle, 'green',False)

        plt.plot([],color='green',label='Prior on q')

        if i < n-1 and j < n-1:

            plot_fisher(F_noq, basetitle, 'red', False)
            
            plt.plot([],color='red',label='Exactly known q')
        
        plt.legend()

        plt.show()