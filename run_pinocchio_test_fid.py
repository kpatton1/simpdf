
import os
from pysimpdf.core import Simulation

nsim_fid = 1
nsim_delta = 1
nnoise = 1

f = ['f', 0.0, nsim_fid, nnoise]

t = [['h', 0.76, nsim_delta, nnoise],
     ['ok', 0.10, nsim_delta, nnoise],
     ['om', 0.32, nsim_delta, nnoise],
     ['s8', 0.90, nsim_delta, nnoise],
     ['w', -0.90, nsim_delta, nnoise],
     ['q', 1.1, nsim_delta, nnoise]]

#l = [f] + t

l = [f]

print l

adict={(4096,'pdf'):(200,-10.0e15,10.0e15), (2048,'pdf'):(200,-10.0e15,10.0e15), (1024,'pdf'):(200,-10.0e15,10.0e15), (512,'pdf'):(200,-10.0e15,10.0e15), (256,'pdf'):(200,-10.0e15,10.0e15), (128,'pdf'):(200,-10.0e15,10.0e15), (64,'pdf'):(200,-10.0e15,10.0e15), (256,'ps'):(768)}

aname='fullcov'

basedir='data'

ls = []

for i in l:
    ls.append(Simulation(basedir,i,adict,aname))

for s in ls:
    s.run_pinocchio()

for s in ls:
    s.convert_healpix()
    s.add_noise()
    
    s.clean_measures()
    s.clean_covariances()
    
    s.measure_data()
    s.generate_covariances()

