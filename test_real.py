
import os
from pysimpdf.core import Simulation

nsim_fid = 1
nsim_delta = 1
nnoise = 1

f = ['f', 0.0, nsim_fid, nnoise]

t_h = ['h', 0.76, nsim_delta, nnoise] 
t_ok = ['ok', 0.10, nsim_delta, nnoise]
t_om = ['om', 0.32, nsim_delta, nnoise]
t_s8 = ['s8', 0.90, nsim_delta, nnoise]
t_w = ['w', -0.90, nsim_delta, nnoise]
t_q = ['q', 1.1, nsim_delta, nnoise]
t_g = ['g', 1.1, nsim_delta, nnoise]

t = [t_h, t_ok, t_om, t_s8, t_w, t_q, t_g]

#l = [f] + t

l = [f]

print l

alist=[(256,'ps',768), (1024,'pdf',200,-1.0,1.0), (512,'pdf',200,-1.0,1.0), (256,'pdf',200,-1.0,1.0), (128,'pdf',200,-1.0,1.0), (64,'pdf',200,-1.0,1.0), (1024, 'moments', 1.0, 2, 10), (512, 'moments', 1.0, 2, 10), (256, 'moments', 1.0, 2, 10), (128, 'moments', 1.0, 2, 10), (64, 'moments', 1.0, 2, 10), (1024,'logpdf',200,-1.0,1.0), (512,'logpdf',200,-1.0,1.0), (256,'logpdf',200,-1.0,1.0), (128,'logpdf',200,-1.0,1.0), (64,'logpdf',200,-1.0,1.0), (1024, 'logmoments', 2, 10), (512, 'logmoments', 2, 10), (256, 'logmoments', 2, 10), (128, 'logmoments', 2, 10), (64, 'logmoments', 2, 10)]

aname='real_test'

basedir='data_real'

ls = []

for i in l:
    ls.append(Simulation(basedir,i,alist,aname))

for s in ls:
    s.run_pinocchio()

for s in ls:
    s.convert_healpix()
    s.add_noise()
#    s.measure_data()
#    s.combine_measures()

#for s in ls:
#    s.generate_covariances()

#for s in ts:
#    s.diff_covariances()

#fs.calc_fisher(ts,[0])
#fs.calc_fisher(ts,[1])
#fs.calc_fisher(ts,[3])
