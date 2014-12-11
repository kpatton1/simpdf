

import numpy
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

#fisher = 'data/F_8192_4096.npz'
fisher = 'data/F_8192_fullcov_7.npz'

data = numpy.load(fisher)

F = data['F']
F_inv = data['F_inv']
params = data['params']

n = len(params)

names = params

centers = {'h':0.69,'ok':0.0, 'w':-1.0, 's8':0.82, 'om':0.29, 'q':1.0}

print F

F_noq = F[:-1,:-1]
F_inv_noq = numpy.linalg.inv(F_noq)

F_priorq = F
F_priorq[-1,-1] += 1/(0.1*0.1)
F_inv_priorq = numpy.linalg.inv(F_priorq)

for i in range(n):
    for j in range(i):
        
        ixgrid = numpy.ix_([j,i],[j,i])

        subM = F_inv[ixgrid]

        subF = numpy.linalg.inv(subM)

        l, v = numpy.linalg.eig(subM)
        l = numpy.sqrt(l)

        ax = plt.subplot(111)
        ax.set_title(names[j] + ' vs ' + names[i])

        ax.set_xlabel(names[j])
        ax.set_ylabel(names[i])

        center_i = centers[params[i]]
        center_j = centers[params[j]]

        ell = Ellipse(xy=(center_j,center_i),width=l[0]*2, height=l[1]*2, angle=numpy.rad2deg(-numpy.arccos(v[0,0])))
        ell.set_facecolor('none')
        ell.set_edgecolor('blue')
        ax.add_artist(ell)

        plt.xlim([center_j-numpy.sqrt(subM[0,0])*1.2,center_j+numpy.sqrt(subM[0,0])*1.2])
        plt.ylim([center_i-numpy.sqrt(subM[1,1])*1.2,center_i+numpy.sqrt(subM[1,1])*1.2])
        plt.scatter([center_j],[center_i])

        subM = F_inv_priorq[ixgrid]

        subF = numpy.linalg.inv(subM)

        l, v = numpy.linalg.eig(subM)
        l = numpy.sqrt(l)

        center_i = centers[params[i]]
        center_j = centers[params[j]]

        ell3 = Ellipse(xy=(center_j,center_i),width=l[0]*2, height=l[1]*2, angle=numpy.rad2deg(-numpy.arccos(v[0,0])))
        ell3.set_facecolor('none')
        ell3.set_edgecolor('green')
        ax.add_artist(ell3)

        if i < n-1 and j < n-1:
            subM = F_inv_noq[ixgrid]

            subF = numpy.linalg.inv(subM)

            l, v = numpy.linalg.eig(subM)
            l = numpy.sqrt(l)

            center_i = centers[params[i]]
            center_j = centers[params[j]]

            ell2 = Ellipse(xy=(center_j,center_i),width=l[0]*2, height=l[1]*2, angle=numpy.rad2deg(-numpy.arccos(v[0,0])))
            ell2.set_facecolor('none')
            ell2.set_edgecolor('red')
            ax.add_artist(ell2)



        plt.show()
