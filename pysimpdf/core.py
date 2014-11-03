
PINOCCHIO_CORE_SURVEY_DEFAULT_NAME='test'
PINOCCHIO_CORE_SURVEY_DEFAULT_SEED=0
PINOCCHIO_CORE_SURVEY_DEFAULT_E=0.30
PINOCCHIO_CORE_SURVEY_DEFAULT_NG=10
PINOCCHIO_CORE_SURVEY_DEFAULT_ZL=0.4
PINOCCHIO_CORE_SURVEY_DEFAULT_ZS=0.8
PINOCCHIO_CORE_SURVEY_DEFAULT_Q=1.0

PINOCCHIO_CORE_ANALYSIS_DEFAULT_NSIDE=8192
PINOCCHIO_CORE_ANALYSIS_DEFAULT_PDFS={'test':(4096,384,110,(-1.0e15,10.0e15))}

PINOCCHIO_CORE_SIMPARAMS_DEFAULT_ZPLC=0.6
PINOCCHIO_CORE_SIMPARAMS_DEFAULT_BOXSIZE=1600.0
PINOCCHIO_CORE_SIMPARAMS_DEFAULT_GRIDSIZE=800

import os
import sys
import numpy
import subprocess
import healpy

from cosmology import Cosmology

class SimParams:

    def __init__(self, zplc=PINOCCHIO_CORE_SIMPARAMS_DEFAULT_ZPLC, boxsize=PINOCCHIO_CORE_SIMPARAMS_DEFAULT_BOXSIZE, gridsize=PINOCCHIO_CORE_SIMPARAMS_DEFAULT_GRIDSIZE):
        self.zplc = zplc
        self.boxsize = boxsize
        self.gridsize = gridsize

class Survey:

    def __init__(self, name=PINOCCHIO_CORE_SURVEY_DEFAULT_NAME, seed=PINOCCHIO_CORE_SURVEY_DEFAULT_SEED, e=PINOCCHIO_CORE_SURVEY_DEFAULT_E, ng=PINOCCHIO_CORE_SURVEY_DEFAULT_NG, zl=PINOCCHIO_CORE_SURVEY_DEFAULT_ZL, zs=PINOCCHIO_CORE_SURVEY_DEFAULT_ZS, q=PINOCCHIO_CORE_SURVEY_DEFAULT_Q):
        
        self.name = name
        self.seed = seed
        self.e = e
        self.ng = ng

        self.zl = zl
        self.zs = zs

        self.q = q

class Analysis:
    
    def __init__(self, nside=PINOCCHIO_CORE_ANALYSIS_DEFAULT_NSIDE, pdfs=PINOCCHIO_CORE_ANALYSIS_DEFAULT_PDFS):
        self.nside = nside
        self.pdfs = pdfs

class SimulationSet:

    def __init__(self, basedir, fiducial, sim_types, map_sizes, bins, mass_range):
        

        self.sims = sims
        self.basedir = basedir
        self.bins = bins
        self.mass_range = mass_range

        self.sims = {}
        #self.fiducial =
        
        type_key = fiducial[0]
        type_delta = fiducial[1]
        type_num_sims = fiducial[2]
        type_num_noise = fiducial[3]
        
        type_name = type_key + '_' + str(type_delta)
            

        for t in self.types:
            type_key = t[0]
            type_delta = t[1]
            type_num_sims = t[2]
            type_num_noise = t[3]

            type_name = type_key + '_' + str(type_delta)
            
            if type_name in self.sims:
                print 'Error! Already simulation type: ' + type_name
                continue
            
            temp_cosmo = Cosmology(seed=i)
            temp_sim = Simulation(temp_cosmo)
        
            self.sims[type_key].append(temp_sim)
            



class Simulation:

    

    def __init__(self,basedir,t,nside,pdfs):
        
        self.basedir = os.path.abspath(basedir)
        
        self.t = t

        param = t[0]
        val = t[1]
        self.nsim = t[2]
        self.nnoise = t[3]

        if param == 'h' or param == 'w' or param == 'om' or param == 'ok' or param == 's8':
            cname = 'delta_' + param + '_' + str(val)
        else:
            cname = 'fiducial'

        if param == 'e' or param == 'q':
            sname = 'delta_' + param + '_' + str(val)
        else:
            sname = 'fiducial'

        self.cosmo = Cosmology(name=cname)
        self.survey = Survey(name=sname)
        self.analysis = Analysis(nside=nside,pdfs=pdfs)
        self.simparams = SimParams()

        if param == 'h':
            self.cosmo.h = val
            
        elif param == 'w':
            self.cosmo.w0 = val

        elif param == 'om':
            self.cosmo.om = val

        elif param == 's8':
            self.cosmo.s8 = val

        elif param == 'ok':
            self.cosmo.ok = val

        elif param == 'e':
            self.survey.e = val

        elif param == 'q':
            self.survey.q = val

        else:
            if param != 'f':
                print 'Simulation error! invalid parameter type: ' + param
            
    def run_pinocchio(self):
        
        cwd = os.getcwd()
        
        if not os.path.exists(self.basedir):
            print 'Pinocchio error! basedir does not exist: ' + self.basedir
            return

        os.chdir(self.basedir)

        if not os.path.exists(self.cosmo.name):
            os.mkdir(self.cosmo.name)

        if not os.path.exists(self.cosmo.name):
            print 'Pinocchio error! could not create simulation dir: ' + self.cosmo.name
            return
        
        os.chdir(self.cosmo.name)

        workingdir = os.getcwd()
        
        for i in range(1,self.nsim+1):
            os.chdir(workingdir)
            
            rundir = 'r' + str(i)
            
            if not os.path.exists(rundir):
                os.mkdir(rundir)
                
            if not os.path.exists(rundir):
                print 'Pinocchio error! could not create rundir: ' + rundir
                return
            
            os.chdir(rundir)
                    
            if not os.path.exists('parameters') or not os.path.exists('outputs'):
                self.cosmo.seed = i
                generate_parameter_files(self.cosmo, self.cosmo.name + '_' + rundir,self.simparams)

            if not os.path.exists('parameters') or not os.path.exists('outputs'):
                print 'Pinocchio error! could not create parameter files'
                return

            pinocchio_output = 'pinocchio.' + self.cosmo.name + '_' + rundir + '.plc.out'

            if not os.path.exists(pinocchio_output):
                run_pinocchio_simulation()

            if not os.path.exists(pinocchio_output):
                print 'Pinocchio error! could not generate output: ' + str(pinocchio_output)
                return

        os.chdir(cwd)

    def convert_healpix(self):
        
        cwd = os.getcwd()

        if not os.path.exists(self.basedir):
            print 'Healpix error! basedir does not exist: ' + self.basedir
            return

        os.chdir(self.basedir)

        if not os.path.exists(self.cosmo.name):
            print 'Healpix error! simulation dir does not exist: ' + self.cosmo.name
            return
        
        os.chdir(self.cosmo.name)

        workingdir = os.getcwd()
        
        for i in range(1,self.nsim+1):
            os.chdir(workingdir)
            
            rundir = 'r' + str(i)
                
            if not os.path.exists(rundir):
                print 'Healpix error! rundir does not exist: ' + rundir
                return
            
            os.chdir(rundir)
            
            pinocchio_output = 'pinocchio.' + self.cosmo.name + '_' + rundir + '.plc.out'

            if not os.path.exists(pinocchio_output):
                print 'Healpix error! pinocchio output does not exist: ' + pinocchio_output
                return

            healpix_output = 'healpix_' + str(self.analysis.nside) + '_' + self.cosmo.name + '_' + rundir + '.fits'

            if not os.path.exists(healpix_output):
                healpixify(pinocchio_output,healpix_output,self.analysis.nside)
                
            if not os.path.exists(healpix_output):
                print 'Healpix error! could not generate: ' + healpix_output
                return

        os.chdir(cwd)

    def add_noise(self):
        
        cwd = os.getcwd()

        if not os.path.exists(self.basedir):
            print 'Noise error! basedir does not exist: ' + self.basedir
            return

        os.chdir(self.basedir)

        if not os.path.exists(self.cosmo.name):
            print 'Noise error! simulation dir does not exist: ' + self.cosmo.name
            return
        
        os.chdir(self.cosmo.name)

        workingdir = os.getcwd()
        
        for i in range(1,self.nsim+1):
            os.chdir(workingdir)
            
            rundir = 'r' + str(i)
                
            if not os.path.exists(rundir):
                print 'Noise error! rundir does not exist: ' + rundir
                return
            
            os.chdir(rundir)

            healpix_output = 'healpix_' + str(self.analysis.nside) + '_' + self.cosmo.name + '_' + rundir + '.fits'
                
            if not os.path.exists(healpix_output):
                print 'Noise error! healpix output does not exist: ' + healpix_output
                return

            for j in range(1,self.nnoise+1):
                noise_output = 'noise_' + str(self.analysis.nside) + '_' + self.cosmo.name + '_' + self.survey.name + '_' + rundir + 'n' + str(j) + '.fits'
                
                if not os.path.exists(noise_output):
                    add_noise(healpix_output,noise_output,self.cosmo,self.survey,self.analysis)
                    
                if not os.path.exists(noise_output):
                    print 'Noise error! could not generate: ' + noise_output
                    continue
        
        os.chdir(cwd)

    def generate_covariances(self):
        
        cwd = os.getcwd()
        
        if not os.path.exists(self.basedir):
            print 'Covariance error! basedir does not exist: ' + self.basedir
            return
        
        os.chdir(self.basedir)
        
        if not os.path.exists(self.cosmo.name):
            print 'Covariance error! simulation dir does not exist: ' + self.cosmo.name
            return
        
        os.chdir(self.cosmo.name)
        
        workingdir = os.getcwd()
        
        for i in range(1,self.nsim+1):
            os.chdir(workingdir)
            
            rundir = 'r' + str(i)
            
            if not os.path.exists(rundir):
                print 'Covariance error! rundir does not exist: ' + rundir
                return
            
            os.chdir(rundir)

            for j in range(1,self.nnoise+1):
                noise_output = 'noise_' + str(self.analysis.nside) + '_' + self.cosmo.name + '_' + self.survey.name + '_' + rundir + 'n' + str(j) + '.fits'
                
                if not os.path.exists(noise_output):
                    print 'Covariance error! noise output does not exist: ' + noise_output
                    continue
                
                cov_outputs = {}

                for key in self.analysis.pdfs.keys():
                    cov_output = 'cov_' + str(self.analysis.nside) + '_' + self.cosmo.name + '_' + self.survey.name + '_' + str(key) + '_' + rundir + 'n' + str(j) + '.npz'
                    
                    if not os.path.exists(cov_output):
                        cov_outputs[key] = cov_output
                
                calc_covariances(noise_output, cov_outputs, self.analysis.pdfs)
                
        os.chdir(cwd)
        

    def plot_analysis(self):
        pass

def generate_parameter_files(cosmo,name,simparams):

    print 'Generating config file, name: ' + name + '  seed: ' + str(cosmo.seed)

    parameters_file = open('parameters','w')

    output= 'RunFlag ' + name + '\n'
    
    output+= 'Omega0 ' + str(cosmo.om) + '\n'
    output+= 'OmegaLambda ' + str(cosmo.ol) + '\n'
    output+= 'OmegaBaryon ' + str(cosmo.ob) + '\n'
    output+= 'Hubble100 ' + str(cosmo.h) + '\n'
    output+= 'Sigma8 ' + str(cosmo.s8) + '\n'
    output+= 'PowerSpectrum ' + str(cosmo.ns) + '\n'
    output+= 'DEw0 ' + str(cosmo.w0) + '\n'
    output+= 'DEw1 ' + str(cosmo.w1) + '\n'

    output+= 'RandomSeed ' + str(cosmo.seed) + '\n'
    
    output+= 'BoxSize ' + str(simparams.boxsize) + '\n'
    output+= 'GridSize ' + str(simparams.gridsize) + '\n'
    output+= 'StartingzForPLC ' + str(simparams.zplc) + '\n'

    extra_output = '''OutputList outputs
TabulatedEoSfile no
CatalogInAscii
FileWithInputSpectrum no
InputSpectrum_UnitLength_in_cm 0
WDM_PartMass_in_kev 0.0
'''

    output += extra_output

    parameters_file.write(output)
    
    parameters_file.close()

    outputs_file = file('outputs','w')
    
    outputs_file.write('0.0\n')
    outputs_file.close()
    

def run_pinocchio_simulation():

    print 'Running pinocchio'

    stdout_file = open('runtime-output.txt','w')
    ret = subprocess.call(["mpirun","-np","8","/n/des/patton.161/pinocchio/pinocchio-3.0.x","parameters"],stdout=stdout_file,stderr=subprocess.STDOUT)
    stdout_file.close()

    return ret

def healpixify(infile,outfile,nside):

    print 'Converting pinocchio output to healpix: ' + infile + ' -> ' + outfile

    npix = healpy.pixelfunc.nside2npix(nside)

    resol = healpy.pixelfunc.nside2resol(nside) * 180.0 / numpy.pi * 60.0
    pixelarea = healpy.pixelfunc.nside2pixarea(nside)  * 180.0 / numpy.pi * 180.0 / numpy.pi * 60.0 * 60.0

    print 'nside:',nside,'npix:',npix,'resol[arcmin]:',resol,'pixelarea[arcmin^2]:',pixelarea

    print 'Infile:',infile,'outfile:',outfile

    print 'Loading data...'

    data = numpy.loadtxt(infile)

    print 'Separating data...'

    mass = data[:,8]
    theta = data[:,9]
    phi = data[:,10]
    z = data[:,1]

    print 'Allocating map memory...'
    
    map_data = numpy.zeros(npix,dtype=numpy.float32)

    print 'Converting angles to radians...'

    theta_rad = -theta * (numpy.pi / 180.0)
    phi_rad = phi * (numpy.pi / 180.0)

    print 'Converting angles to pixels...'
    
    pix = healpy.pixelfunc.ang2pix(nside, theta_rad, phi_rad, nest=True)

    for i in range(len(pix)):
        p = pix[i]
        m = mass[i]
        map_data[p] += m

    healpy.fitsfunc.write_map(outfile, map_data)

def add_noise(infile,outfile,cosmo,survey,analysis):

    print 'Adding noise to healpix map: ' + infile + ' -> ' + outfile
    
    sigma_e = survey.e

    zl = survey.zl
    zs = survey.zs

    dl = cosmo.Da(zl)
    ds = cosmo.Da(zs)
    dls = cosmo.Da_rel(zl,zs)

    e_crit = cosmo.c * cosmo.c * ds / (4.0 * numpy.pi * cosmo.G * dls * dl) # [Msolar / Mpc^2]

    e_crit = e_crit * dl * dl * numpy.pi * numpy.pi / 60.0 / 60.0 / 180.0 / 180.0 # [Msolar / arcmin^2]

    print '  e_crit: ' + str(e_crit)
    
#    e_crit = 5.49308e+14

    n_gal = survey.ng

#    nside = analysis.nside

    map_data = healpy.fitsfunc.read_map(infile, dtype=numpy.float32)

    nside = healpy.pixelfunc.get_nside(map_data)

    pixelarea = healpy.pixelfunc.nside2pixarea(nside)  * 180.0 / numpy.pi * 180.0 / numpy.pi * 60.0 * 60.0
    
    sigma = survey.q * survey.e * e_crit * pixelarea / numpy.sqrt(8.0 * numpy.pi * numpy.pi * n_gal * pixelarea)
        
    map_data = map_data + sigma * numpy.random.randn(12*nside*nside)
    
    healpy.fitsfunc.write_map(outfile, map_data)

def calc_covariances(infile,outfiles,pdfs):
    
    if len(outfiles) < 1:
        return
    
    map_data = healpy.fitsfunc.read_map(infile)
    
    for key in outfiles.keys():
        
        outfile = outfiles[key]
        pdf = pdfs[key]
        
        print 'Calculating covariances: ' + infile + ' -> ' + outfile
    
        mean_dict = {}
        cov_dict = {}

        map_size = pdf[0]
        divs = pdf[1]
        bins = pdf[2]
        mass_range = pdf[3]
        
        mean = numpy.zeros(bins,dtype=numpy.float32)
        cov = numpy.zeros((bins,bins),dtype=numpy.float32)
        
        map_temp = healpy.pixelfunc.ud_grade(map_data,map_size,power=-2,order_in='NESTED',order_out='NESTED')
        
        npix = healpy.pixelfunc.nside2npix(map_size)
        
        for n in range(divs):
            temp = map_temp[npix/divs*n:npix/divs*(n+1)]
            temp = temp - numpy.mean(temp)
            hist,bin_edges = numpy.histogram(temp, bins=bins, range=mass_range)
            mean += hist
            for x in range(bins):
                for y in range(bins):
                    cov[x][y] += hist[x] * hist[y]
        
        mean = mean / float(divs)
        cov = cov / float(divs)
    
        x = 0.5*(bin_edges[1:]+bin_edges[:-1])
    
        numpy.savez(outfile,x=x,mean=mean,cov=cov)
