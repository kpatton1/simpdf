
import os
from pysimpdf.core import Simulation

nsim_fid = 200
nsim_delta = 50
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

l = [f] + t

#l = t + [f]

print l

alist=[(256,'ps',768), (4096,'pdf',200,-10.0e15,10.0e15), (2048,'pdf',200,-10.0e15,10.0e15), (1024,'pdf',200,-10.0e15,10.0e15), (512,'pdf',200,-10.0e15,10.0e15), (256,'pdf',200,-10.0e15,10.0e15), (128,'pdf',200,-10.0e15,10.0e15), (64,'pdf',200,-10.0e15,10.0e15), (4096, 'moments', 1.0e14, 2, 10), (2048, 'moments', 1.0e14, 2, 10), (1024, 'moments', 1.0e14, 2, 10), (512, 'moments', 1.0e14, 2, 10), (256, 'moments', 1.0e14, 2, 10), (128, 'moments', 1.0e14, 2, 10), (64, 'moments', 1.0e14, 2, 10)]

for i in range(len(alist)):
    print i, '->', alist[i]

aname='fullcovmoments'

basedir='data'

fs = Simulation(basedir,f,alist,aname)

ts = []

for i in t:
    ts.append(Simulation(basedir,i,alist,aname))

ls = ts + [fs]

dostuff=False

if dostuff:

    for s in ls:
        s.run_pinocchio()

    for s in ls:
        s.convert_healpix()
        s.add_noise()
        s.measure_data()
        s.combine_measures()

#fs.clean_covariances()

for s in ls:
    #s.clean_covariances()
    s.generate_covariances()

for s in ts:
    #s.clean_diff_covariances()
    s.diff_covariances(fs)

fs.calc_fisher(ts,[0])
fs.calc_fisher(ts,[3])
fs.calc_fisher(ts,[0,3])

fs.calc_fisher(ts,[1])
fs.calc_fisher(ts,[3])
fs.calc_fisher(ts,[5])
fs.calc_fisher(ts,[7])

fs.calc_fisher(ts, [0,5])

fs.calc_fisher(ts, [3,5])
fs.calc_fisher(ts, [0,3,5])

fs.calc_fisher(ts, [3,5,7])
fs.calc_fisher(ts, [0,3,5,7])

fs.calc_fisher(ts,[1,3,5,7])

fs.calc_fisher(ts,[0,1,3,5,7])
