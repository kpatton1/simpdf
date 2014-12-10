
import os
from pysimpdf.core import Simulation

nsim_fid = 20
nsim_delta = 5
nnoise = 1

f = ['f', 0.0, nsim_fid, nnoise]

t = [['h', 0.76, nsim_delta, nnoise],
     ['ok', 0.10, nsim_delta, nnoise],
     ['om', 0.32, nsim_delta, nnoise],
     ['s8', 0.90, nsim_delta, nnoise],
     ['w', -0.90, nsim_delta, nnoise],
     ['q', 1.1, nsim_delta, nnoise]]


l = [f] + t

#l = t + [f]

print l

adict={(4096,'pdf'):(200,-10.0e15,10.0e15), (2048,'pdf'):(200,-10.0e15,10.0e15), (1024,'pdf'):(200,-10.0e15,10.0e15), (512,'pdf'):(200,-10.0e15,10.0e15), (256,'pdf'):(200,-10.0e15,10.0e15), (128,'pdf'):(200,-10.0e15,10.0e15), (64,'pdf'):(200,-10.0e15,10.0e15), (256,'ps'):(768)}

for key in adict.keys():
    print key

aname='fullcov'

basedir='data'

fs = Simulation(basedir,f,adict,aname)

ts = []

for i in t:
    ts.append(Simulation(basedir,i,adict,aname))

ls = ts + [fs]

#for s in ls:
#    s.run_pinocchio()

#for s in ls:
#    s.convert_healpix()
#    s.add_noise()
#    s.measure_data()
#    s.generate_covariances()
#    s.combine_covariances()


#for s in ts:
#    s.diff_covariances(fs)

fs.calc_fisher(ts,[0])

fs.calc_fisher(ts,[3])
fs.calc_fisher(ts,[0,3])

fs.calc_fisher(ts,[7])
fs.calc_fisher(ts,[0,7])

fs.calc_fisher(ts,[0,3,7,2])
#fs.calc_fisher(ts,[3,1,7,4,2,6,5])
