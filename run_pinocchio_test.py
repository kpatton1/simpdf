
import os
from pysimpdf.core import Simulation

f = ['f', 0.0, 20, 5]

t = [['h', 0.76, 5, 5],
     ['ok', 0.10, 5, 5],
     ['om', 0.32, 5, 5],
     ['s8', 0.90, 5, 5],
     ['w', -0.90, 5, 5],
     ['q', 1.1, 5, 5]]

#l = [f] + t

l = t + [f]

print l

maps={4096:(4096,384,110,(-1.0e15,10.0e15))}
nside=8192

basedir='data'

ls = []

for i in l:
    ls.append(Simulation(basedir,i,nside,maps))


for s in ls:
    s.run_pinocchio()

for s in ls:
    s.convert_healpix()
    s.add_noise()

