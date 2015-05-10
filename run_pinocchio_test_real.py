
import os
from pysimpdf.core import Simulation

nsim_fid = 20
nsim_delta = 5
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

#alist=[(256,'ps',768), (1024,'pdf',200,-1.0,1.0), (512,'pdf',200,-1.0,1.0), (256,'pdf',200,-1.0,1.0), (128,'pdf',200,-1.0,1.0), (64,'pdf',200,-1.0,1.0), (1024, 'moments', 1.0, 2, 10), (512, 'moments', 1.0, 2, 10), (256, 'moments', 1.0, 2, 10), (128, 'moments', 1.0, 2, 10), (64, 'moments', 1.0, 2, 10), (256,'logps',768), (1024,'logpdf',200,-1.0,1.0), (512,'logpdf',200,-1.0,1.0), (256,'logpdf',200,-1.0,1.0), (128,'logpdf',200,-1.0,1.0), (64,'logpdf',200,-1.0,1.0), (1024, 'logmoments', 2, 10), (512, 'logmoments', 2, 10), (256, 'logmoments', 2, 10), (128, 'logmoments', 2, 10), (64, 'logmoments', 2, 10)]
#aname='real'

#alist=[(1024,'pdf',200,-1.0,1.0), (256,'pdf',200,-1.0,1.0), (1024, 'moments', 1.0, 2, 10), (256, 'moments', 1.0, 2, 10), (1024,'logpdf',800,-20.0,20.0), (256,'logpdf',800,-20.0,20.0), (1024, 'logmoments', 1, 10), (256, 'logmoments', 1, 10)]
#aname='test_simple'

#alist=[(256,'ps',768), (256,'logps',768)]
#aname='test_ps'

#alist=[(1024,'pdf',50, -1.0,1.0),(1024,'pdf',100, -1.0,1.0),(1024,'pdf',150, -1.0,1.0),(1024,'pdf',200, -1.0,1.0), (1024,'pdf',250, -1.0,1.0), (1024,'pdf',300, -1.0,1.0), (1024,'pdf',350, -1.0,1.0), (1024,'pdf',400, -1.0,1.0), (1024,'pdf',450, -1.0,1.0), (1024,'pdf',500, -1.0,1.0), (1024,'pdf',550, -1.0,1.0), (1024,'pdf',600, -1.0,1.0), (1024,'pdf',650, -1.0,1.0), (1024,'pdf',700, -1.0,1.0), (1024,'pdf',750, -1.0,1.0), (1024,'pdf',800, -1.0,1.0)]
#aname='test_pdf1024'

#alist=[(1024,'pdf',800, -1.0,1.0),(1024,'pdf',900, -1.0,1.0),(1024,'pdf',1000, -1.0,1.0),(1024,'pdf',1100, -1.0,1.0), (1024,'pdf',1200, -1.0,1.0), (1024,'pdf',1300, -1.0,1.0), (1024,'pdf',1400, -1.0,1.0), (1024,'pdf',1500, -1.0,1.0), (1024,'pdf',1600, -1.0,1.0), (1024,'pdf',1700, -1.0,1.0), (1024,'pdf',1800, -1.0,1.0), (1024,'pdf',1900, -1.0,1.0), (1024,'pdf',2000, -1.0,1.0)]
#aname='test_pdf1024_2'

#alist=[(1024,'pdf',2000, -1.0,1.0),(1024,'pdf',2200, -1.0,1.0),(1024,'pdf',2400, -1.0,1.0),(1024,'pdf',2600, -1.0,1.0), (1024,'pdf',3000, -1.0,1.0)]
#aname='test_pdf1024_3'

#alist=[(1024,'pdf',4000, -1.0,1.0), (1024,'pdf',5000, -1.0,1.0)]
#aname='test_pdf1024_4'

alist=[(256,'pdf',200, -1.0,1.0),(256,'pdf',400, -1.0,1.0),(256,'pdf',800, -1.0,1.0),(256,'pdf',1600, -1.0,1.0), (256,'pdf',3200, -1.0,1.0)]
aname='test_pdf256'

#alist=[(256,'pdf',300, -1.0,1.0),(256,'pdf',600, -1.0,1.0),(256,'pdf',1200, -1.0,1.0),(256,'pdf',2400, -1.0,1.0)]
#aname='test_pdf256_2'


for i in range(len(alist)):
    print i, '->', alist[i]

basedir='data_real'

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
        #s.add_noise()
        #s.measure_data()

#fs.clean_covariances()

for s in ls:
    s.add_noise()
    s.measure_data()
    s.combine_measures()
    #s.clean_covariances()
    s.generate_covariances()

for s in ts:
    #s.clean_diff_covariances()
    s.diff_covariances(fs)

fs.calc_fisher(ts,[0])
fs.calc_fisher(ts,[1])
fs.calc_fisher(ts,[2])
fs.calc_fisher(ts,[3])
fs.calc_fisher(ts,[4])
fs.calc_fisher(ts,[5])
fs.calc_fisher(ts,[6])
fs.calc_fisher(ts,[7])
fs.calc_fisher(ts,[8])
fs.calc_fisher(ts,[9])
fs.calc_fisher(ts,[10])
fs.calc_fisher(ts,[11])
fs.calc_fisher(ts,[12])
fs.calc_fisher(ts,[13])
fs.calc_fisher(ts,[14])
fs.calc_fisher(ts,[15])

#fs.calc_fisher(ts,[0])
#fs.calc_fisher(ts,[3])
#fs.calc_fisher(ts,[0,3])

#fs.calc_fisher(ts,[1])
#fs.calc_fisher(ts,[3])
#fs.calc_fisher(ts,[5])
#fs.calc_fisher(ts,[7])

#fs.calc_fisher(ts, [0,5])

#fs.calc_fisher(ts, [3,5])
#fs.calc_fisher(ts, [0,3,5])

#fs.calc_fisher(ts, [3,5,7])
#fs.calc_fisher(ts, [0,3,5,7])

#fs.calc_fisher(ts,[1,3,5,7])

#fs.calc_fisher(ts,[0,1,3,5,7])
