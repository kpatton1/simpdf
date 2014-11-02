


import os
import os.path

from pysimpdf.core import Simulation

def covert_files(s):
    
    cwd = os.getcwd()
    
    if not os.path.exists(s.basedir):
        print 'Error! basedir does not exist: ' + s.basedir
        return
    
    os.chdir(s.basedir)
    
    if not os.path.exists(s.cosmo.name):
        print 'Error! simulation dir does not exist: ' + s.cosmo.name
        return
    
    os.chdir(s.cosmo.name)
    
    workingdir = os.getcwd()
    
    for i in range(1,s.nsim+1):
        os.chdir(workingdir)
        
        rundir = 'r' + str(i)
        
        if not os.path.exists(rundir):
            print 'Error! rundir does not exist: ' + rundir
            return
        
        os.chdir(rundir)

        pinocchio_output = 'pinocchio.' + self.cosmo.name + '_' + rundir + '.plc.out'

        if not os.path.exists(pinocchio_output):
            print 'Error! pinocchio output does not exist: ' + pinocchio_output
            return

        healpix_output = 'healpix_' + str(s.analysis.nside) + '_' + s.cosmo.name + '_' + rundir + '.fits'
        
        if not os.path.exists(healpix_output):
            old_healpix_output = 'healpix_' + s.cosmo.name + '_' + rundir + '_' + str(s.analysis.nside) + '.fits'
            
            if os.path.exists(old_healpix_output):
                os.rename(old_healpix_output, healpix_output)
            else:
                print 'Error! healpix output does not exist: ' + healpix_output
                print 'Error! old healpix output does not exist: ' + old_healpix_output
                return
        
        for j in range(1,s.nnoise+1):
            noise_output = 'noise_' + str(s.analysis.nside) + '_' + s.cosmo.name + '_' + s.survey.name + '_' + rundir + 'n' + str(j) + '.fits'
            
            if not os.path.exists(noise_output):
                old_noise_output = 'noise_' + str(s.cosmo.name) + '_' + rundir + 'n' + str(j) + '_' + str(s.analysis.nside) + '.fits'
                
                if os.path.exists(old_noise_output):
                    os.rename(old_noise_output, noise_output)
                else:
                    print 'Error! noise output does not exist: ' + noise_output
                    print 'Error! old noise output does not exist: ' + old_noise_output
                
    
    os.chdir(cwd)

maps={4096:(4096,384,110,(-1.0e15,10.0e15))}
nside=8192

basedir = 'data'

f = ['f', 0.0, 20, 5]

t = [['h', 0.76, 5, 5],
     ['ok', 0.10, 5, 5],
     ['om', 0.32, 5, 5],
     ['s8', 0.90, 5, 5],
     ['w', -0.90, 5, 5]]

l = t + [f]

ls = []

for i in l:
    ls.append(Simulation(basedir,i,nside,maps))

for s in ls:
    convert_filenames(s)