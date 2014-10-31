
import os
from pysimpdf.core import Simulation

f = ['f', 0.0, 1, 1]

t = [['h', 0.76, 5, 5],
     ['ok', 0.10, 5, 5],
     ['om', 0.32, 5, 5],
     ['s8', 0.90, 5, 5],
     ['w', -0.90, 5, 5]]

maps={4096:(4096,384,110,(-1.0e15,10.0e15))}
nside=8192

basedir='data'

s = Simulation(basedir,f,nside,maps)

s.run_pinocchio()

#s.convert_healpix(8192)

#s.add_noise()

