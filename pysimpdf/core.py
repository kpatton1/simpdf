
PINOCCHIO_CORE_SURVEY_DEFAULT_NAME='test'
PINOCCHIO_CORE_SURVEY_DEFAULT_NSIDE=8192
PINOCCHIO_CORE_SURVEY_DEFAULT_SEED=0
PINOCCHIO_CORE_SURVEY_DEFAULT_E=0.30
PINOCCHIO_CORE_SURVEY_DEFAULT_NG=10
PINOCCHIO_CORE_SURVEY_DEFAULT_ZL=0.4
PINOCCHIO_CORE_SURVEY_DEFAULT_ZS=0.8
PINOCCHIO_CORE_SURVEY_DEFAULT_Q=1.0

PINOCCHIO_CORE_ANALYSIS_DEFAULT_NAME='test'
PINOCCHIO_CORE_ANALYSIS_DEFAULT_ADICT={}
PINOCCHIO_CORE_ANALYSIS_DEFAULT_DIVS=384

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

    def __init__(self, name=PINOCCHIO_CORE_SURVEY_DEFAULT_NAME, nside=PINOCCHIO_CORE_SURVEY_DEFAULT_NSIDE, seed=PINOCCHIO_CORE_SURVEY_DEFAULT_SEED, e=PINOCCHIO_CORE_SURVEY_DEFAULT_E, ng=PINOCCHIO_CORE_SURVEY_DEFAULT_NG, zl=PINOCCHIO_CORE_SURVEY_DEFAULT_ZL, zs=PINOCCHIO_CORE_SURVEY_DEFAULT_ZS, q=PINOCCHIO_CORE_SURVEY_DEFAULT_Q):
        
        self.name = name
        self.nside = nside
        self.seed = seed
        self.e = e
        self.ng = ng

        self.zl = zl
        self.zs = zs

        self.q = q

class Analysis:
    
    def __init__(self, name=PINOCCHIO_CORE_ANALYSIS_DEFAULT_NAME, adict=PINOCCHIO_CORE_ANALYSIS_DEFAULT_ADICT, divs=PINOCCHIO_CORE_ANALYSIS_DEFAULT_DIVS):

        self.name = name
        self.adict = adict
        self.divs = divs

        self.alist = []

        for n,t in adict.keys():
            params = adict[(n,t)]
            
            if t == 'pdf':
                map_size = n
                
                bins = params[0]
                bin_min = params[1]
                bin_max = params[2]

                a = AnalysisPDF(n, self.divs, bins, bin_min, bin_max)
                a.info = ((n,t),params)

                self.alist.append(a)

            elif t == 'ps':
                map_size = n
                
                lmin = 0
                lmax = 3 * map_size
                
                a = AnalysisPS(n, self.divs, lmin, lmax)
                a.info = ((n,t),params)

                self.alist.append(a)

            else:
                print 'Error! Unknown analysis type: ' + str(t)
                print str((n,t)) + ' -> ' + str(params)

class AnalysisPS:

    def __init__(self, map_size, divs, lmin=0, lmax=-1):
        
        if lmax < 0 or lmax > 3 * map_size:
            lmax = 3 * map_size

        if lmin > lmax:
            lmin = lmax

        self.map_size = map_size
        self.divs = divs
        self.lmin = lmin
        self.lmax = lmax

        self.bins = lmax - lmin

    def x(self):

        x = numpy.arange(self.lmin,self.lmax)

        return x

    def process_data(self, data):

        npix = healpy.pixelfunc.nside2npix(self.map_size)

        ps_data = numpy.zeros(npix, dtype=numpy.float32)

        ps_data[0:npix/self.divs] = data

        ps_data_n2r = healpy.pixelfunc.reorder(ps_data, n2r = True)

        cls = healpy.sphtfunc.anafast(ps_data_n2r)

        return cls[self.lmin:self.lmax]

class AnalysisPDF:

    def __init__(self, map_size, divs, bins, bin_min, bin_max):

        self.map_size = map_size
        self.divs = divs
        self.bins = bins
        self.bin_min = bin_min
        self.bin_max = bin_max

    def x(self):
        step = numpy.float32(self.bin_max - self.bin_min) / self.bins

        x = numpy.arange(self.bins, dtype=numpy.float32)

        x = (x + 0.5) * step

        return x

    def process_data(self, data):
        
        hist,bin_edges = numpy.histogram(data, bins=self.bins, range=(self.bin_min,self.bin_max))

        return hist
        

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

    

    def __init__(self,basedir,t,adict,aname='fiducial'):
        
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

        #aname = 'fiducial'

        self.cosmo = Cosmology(name=cname)
        self.survey = Survey(name=sname)
        self.analysis = Analysis(name=aname,adict=adict)
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
        
        print 'Running pinocchio for ' + str(self.cosmo.name)
        
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
        
        print 'Converting to healpix for ' + str(self.cosmo.name)
        
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

            healpix_output = 'healpix_' + str(self.survey.nside) + '_' + self.cosmo.name + '_' + rundir + '.fits'

            if not os.path.exists(healpix_output):
                healpixify(pinocchio_output,healpix_output,self.survey.nside)
                
            if not os.path.exists(healpix_output):
                print 'Healpix error! could not generate: ' + healpix_output
                return

        os.chdir(cwd)

    def add_noise(self):
        
        print 'Adding noise to ' + str(self.cosmo.name) + ' with survey ' + str(self.survey.name)
        
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

            healpix_output = 'healpix_' + str(self.survey.nside) + '_' + self.cosmo.name + '_' + rundir + '.fits'
                
            if not os.path.exists(healpix_output):
                print 'Noise error! healpix output does not exist: ' + healpix_output
                return

            for j in range(1,self.nnoise+1):
                noise_output = 'noise_' + str(self.survey.nside) + '_' + self.cosmo.name + '_' + self.survey.name + '_' + rundir + 'n' + str(j) + '.fits'
                
                if not os.path.exists(noise_output):
                    add_noise(healpix_output,noise_output,self.cosmo,self.survey)
                    
                if not os.path.exists(noise_output):
                    print 'Noise error! could not generate: ' + noise_output
                    continue
        
        os.chdir(cwd)

    def measure_data(self):
        
        print 'Measuring data for ' + str(self.cosmo.name) + ' with survey ' + str(self.survey.name) + ' and analysis ' + str(self.analysis.name)
    
        cwd = os.getcwd()
        
        if not os.path.exists(self.basedir):
            print 'Measurement error! basedir does not exist: ' + self.basedir
            return

        os.chdir(self.basedir)
        
        if not os.path.exists(self.cosmo.name):
            print 'Measurement error! simulation dir does not exist: ' + self.cosmo.name
            return

        os.chdir(self.cosmo.name)
    
        workingdir = os.getcwd()
        
        for i in range(1,self.nsim+1):
            os.chdir(workingdir)
            
            rundir = 'r' + str(i)
            
            if not os.path.exists(rundir):
                print 'Measurement error! rundir does not exist: ' + rundir
                return
        
            os.chdir(rundir)
            
            for j in range(1,self.nnoise+1):
                noise_output = 'noise_' + str(self.survey.nside) + '_' + self.cosmo.name + '_' + self.survey.name + '_' + rundir + 'n' + str(j) + '.fits'
                
                if not os.path.exists(noise_output):
                    print 'Measurement error! noise output does not exist: ' + noise_output
                    continue

                measure_output = 'measure_' + str(self.survey.nside) + '_' + self.cosmo.name + '_' + self.survey.name + '_' + self.analysis.name + '_' + rundir + 'n' + str(j) + '.npz'
                
                if not os.path.exists(measure_output):
                    measure_data(noise_output, measure_output, self.analysis)
                
                if not os.path.exists(measure_output):
                    print 'Measurement error! could not generate: ' + measure_output
    
        os.chdir(cwd)

    def clean_measures(self):
    
        print 'Cleaning measurements for ' + str(self.cosmo.name) + ' with survey ' + str(self.survey.name) + ' and analysis ' + str(self.analysis.name)
        
        cwd = os.getcwd()
        
        if not os.path.exists(self.basedir):
            print 'Measurement error! basedir does not exist: ' + self.basedir
            return
        
        os.chdir(self.basedir)
        
        if not os.path.exists(self.cosmo.name):
            print 'Measurement error! simulation dir does not exist: ' + self.cosmo.name
            return
    
        os.chdir(self.cosmo.name)
        
        workingdir = os.getcwd()
        
        for i in range(1,self.nsim+1):
            os.chdir(workingdir)
            
            rundir = 'r' + str(i)
            
            if not os.path.exists(rundir):
                print 'Measurement error! rundir does not exist: ' + rundir
                return
        
            os.chdir(rundir)
            
            for j in range(1,self.nnoise+1):
                measure_output = 'measure_' + str(self.survey.nside) + '_' + self.cosmo.name + '_' + self.survey.name + '_' + self.analysis.name + '_' + rundir + 'n' + str(j) + '.npz'
                
                if os.path.exists(measure_output):
                    os.remove(measure_output)

        os.chdir(cwd)

    def generate_covariances(self):
        
        print 'Generating covariances for ' + str(self.cosmo.name) + ' with survey ' + str(self.survey.name) + ' and analysis ' + str(self.analysis.name)
        
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

                measure_output = 'measure_' + str(self.survey.nside) + '_' + self.cosmo.name + '_' + self.survey.name + '_' + self.analysis.name + '_' + rundir + 'n' + str(j) + '.npz'
                
                if not os.path.exists(measure_output):
                    print 'Covariance error! measurement output does not exist: ' + measure_output
                    continue
                
                cov_output = 'cov_' + str(self.survey.nside) + '_' + self.cosmo.name + '_' + self.survey.name + '_' + self.analysis.name + '_' + rundir + 'n' + str(j) + '.npz'

                if not os.path.exists(cov_output):
                    calc_covariance(measure_output, cov_output)
                        
                if not os.path.exists(cov_output):
                    print 'Covariance error! could not generate: ' + cov_output
        
        os.chdir(cwd)

    def clean_covariances(self):
    
        print 'Cleaning covariances for ' + str(self.cosmo.name) + ' with survey ' + str(self.survey.name) + ' and analysis ' + str(self.analysis.name)
        
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
            
                cov_output = 'cov_' + str(self.survey.nside) + '_' + self.cosmo.name + '_' + self.survey.name + '_' + self.analysis.name + '_' + rundir + 'n' + str(j) + '.npz'
                
                if os.path.exists(cov_output):
                    os.remove(cov_output)
        
        os.chdir(cwd)


    def combine_covariances(self):

        cwd = os.getcwd()
        
        if not os.path.exists(self.basedir):
            print 'Covariance averaging error! basedir does not exist: ' + self.basedir
            return
        
        os.chdir(self.basedir)
        
        if not os.path.exists(self.cosmo.name):
            print 'Covariance averaging error! simulation dir does not exist: ' + self.cosmo.name
            return
        
        os.chdir(self.cosmo.name)
        
        workingdir = os.getcwd()

        cov_outputs = []
    
        for i in range(1,self.nsim+1):
            os.chdir(workingdir)
        
            rundir = 'r' + str(i)
        
            if not os.path.exists(rundir):
                print 'Covariance averaging error! rundir does not exist: ' + rundir
                return
        
            os.chdir(rundir)

            for j in range(1,self.nnoise+1):
                cov_output = 'cov_' + str(self.survey.nside) + '_' + self.cosmo.name + '_' + self.survey.name + '_' + self.analysis.name + '_' + rundir + 'n' + str(j) + '.npz'
                
                cov_output = os.path.abspath(cov_output)

                if not os.path.exists(cov_output):
                    print 'Covariance averaging error! covariance does not exist: ' + cov_output
                    return
                
                cov_outputs.append(cov_output)

        os.chdir(workingdir)
        
        avg_cov_output = 'cov_' + str(self.survey.nside) + '_' + self.cosmo.name + '_' + self.survey.name + '_' + self.analysis.name + '.npz'

        avg_covariances(cov_outputs, avg_cov_output)
        
        os.chdir(cwd)

    def diff_covariances(self, fid):
        
        cwd = os.getcwd()
        
        if not os.path.exists(self.basedir):
            print 'Covariance differencing error! basedir does not exist: ' + self.basedir
            return
        
        os.chdir(self.basedir)

        param = self.t[0]

        if param == 'h':
            delta = self.cosmo.h - fid.cosmo.h
            
        elif param == 'w':
            delta = self.cosmo.w0 - fid.cosmo.w0

        elif param == 'om':
            delta = self.cosmo.om - fid.cosmo.om

        elif param == 's8':
            delta = self.cosmo.s8 - fid.cosmo.s8

        elif param == 'ok':
            delta = self.cosmo.ok - fid.cosmo.ok

        elif param == 'e':
            delta = self.survey.e - fid.survey.e

        elif param == 'q':
            delta = self.survey.q - fid.survey.q

        else:
            print 'Covariance differencing error! invalid parameter type: ' + param
            return


        delta_cov_output = 'cov_' + str(self.survey.nside) + '_' + self.cosmo.name + '_' + self.survey.name + '_' + self.analysis.name + '.npz'
        fid_cov_output = 'cov_' + str(fid.survey.nside) + '_' + fid.cosmo.name + '_' + fid.survey.name + '_' + fid.analysis.name + '.npz'

        delta_cov_output = self.basedir + '/' + self.cosmo.name + '/' + delta_cov_output
        fid_cov_output = fid.basedir + '/' + fid.cosmo.name + '/' + fid_cov_output

        if not os.path.exists(delta_cov_output):
            print 'Covariance differencing error! delta cov output does not exist: ' + delta_cov_output
            return

        if not os.path.exists(fid_cov_output):
            print 'Covariance differencing error! fiducial cov output does not exist: ' + fid_cov_output
            return

        dcov_out = 'dcov_' + str(self.survey.nside) + '_' + self.cosmo.name + '_' + self.survey.name + '_' + self.analysis.name + '.npz'

        diff_covariances(delta_cov_output, fid_cov_output, dcov_out, delta)
                
        os.chdir(cwd)

    def calc_fisher(self, deltas, ranges):
        
        cwd = os.getcwd()
        
        if not os.path.exists(self.basedir):
            print 'Fisher error! basedir does not exist: ' + self.basedir
            return
        
        os.chdir(self.basedir)

        fid_cov_output = 'cov_' + str(self.survey.nside) + '_' + self.cosmo.name + '_' + self.survey.name + '_' + self.analysis.name + '.npz'
        fid_cov_output = self.basedir + '/' + self.cosmo.name + '/' + fid_cov_output

        if not os.path.exists(fid_cov_output):
            print 'Fisher error! fiducial cov output does not exist: ' + fid_cov_output
            return

        dcov_outs = []
        params = []

        for delta in deltas:
            dcov_out = 'dcov_' + str(delta.survey.nside) + '_' + delta.cosmo.name + '_' + delta.survey.name + '_' + delta.analysis.name + '.npz'
            dcov_out = delta.basedir + '/' + dcov_out
            
            if not os.path.exists(dcov_out):
                print 'Fisher error! delta cov output does not exist: ' + dcov_out
                return

            param = delta.t[0]

            dcov_outs.append(dcov_out)
            params.append(param)

        code = ''

        for r in ranges:
            code = code + str(r)

        fisher_out = 'F_' + str(delta.survey.nside) + '_' + self.analysis.name + '_' + code + '.npz'

        calc_fisher(fid_cov_output, dcov_outs, params, fisher_out, ranges)
                
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

def add_noise(infile,outfile,cosmo,survey):

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
    
    n_gal = survey.ng

    map_data = healpy.fitsfunc.read_map(infile, dtype=numpy.float32)

    nside = healpy.pixelfunc.get_nside(map_data)

    pixelarea = healpy.pixelfunc.nside2pixarea(nside)  * 180.0 / numpy.pi * 180.0 / numpy.pi * 60.0 * 60.0
    
    sigma = survey.q * survey.e * e_crit * pixelarea / numpy.sqrt(8.0 * numpy.pi * numpy.pi * n_gal * pixelarea)
        
    map_data = map_data + sigma * numpy.random.randn(12*nside*nside)
    
    healpy.fitsfunc.write_map(outfile, map_data)

def measure_data(infile,outfile,analysis):

    map_data_raw = healpy.fitsfunc.read_map(infile, dtype=numpy.float32)

    mdata = {}

    divs = analysis.divs
    alist = analysis.alist

    measurements = []
    
    x = []
    
    ranges = []
    
    info = []
    
    xlen = 0

    for a in alist:
        map_size = a.map_size
        
        a_x = a.x()
        
        a_range = (xlen,xlen+len(a_x))
        
        ranges.append(a_range)
        
        xlen = xlen+len(a_x)
        
        info.append(a.info)
        
        x.extend(a_x)

        if map_size not in mdata.keys():
            mdata[map_size] = healpy.pixelfunc.ud_grade(map_data_raw,map_size,power=-2,order_in='NESTED',order_out='NESTED')


    for n in range(divs):
        
        print 'div ' + str(n)
        
        measurement = []
        
        acount = 0
        
        for a in alist:
            print 'analysis ' + str(acount)
            map_size = a.map_size
            npix = healpy.pixelfunc.nside2npix(map_size)
            data = mdata[map_size][npix/divs*n:npix/divs*(n+1)]
            data = data - numpy.mean(data)

            a_m = a.process_data(data)

            measurement.extend(a_m)
        
            acount += 1

        measurements.append(measurement)

    numpy.savez(outfile,x=x,m=measurements,r=ranges,info=info)

def calc_covariance(infile,outfile):
    
    data = numpy.load(infile)
    
    x = numpy.array(data['x'],dtype=numpy.float64)
    m = numpy.array(data['m'],dtype=numpy.float64)
    r = data['r']
    info = data['info']
    
    n = len(x)
    
    mean = numpy.zeros(n,dtype=numpy.float64)
    cov = numpy.zeros((n,n),dtype=numpy.float64)
    
    count = 0
    
    for i in m:
        mean += i
        count += 1
        cov += numpy.outer(i,i)

    if count > 0:
        mean = mean / float(count)
        cov = cov / float(count)
        
    numpy.savez(outfile,x=x,mean=mean,cov=cov,r=r,info=info)

def avg_covariances(infiles, outfile):
    count = 0.0
    
    first = True
    
    for f in infiles:
        
        data = numpy.load(f)
        
        if first:
            x = data['x']
            info = data['info']
            r = data['r']
            bins = len(x)
            mean_avg = numpy.zeros(bins,dtype=numpy.float64)
            cov_avg = numpy.zeros((bins,bins),dtype=numpy.float64)
            first = False

        cov = data['cov']
        mean = data['mean']

        mean_avg += mean
        cov_avg += cov
        count += 1.0

    mean_avg = mean_avg / count
    cov_avg = cov_avg / count

    numpy.savez(outfile, x=x, mean=mean_avg, cov=cov_avg, r=r, info=info)


def diff_covariances(infile, fidfile, outfile, delta):

    data = numpy.load(infile)
    fid = numpy.load(fidfile)

    data_x = data['x']

    data_mean = data['mean']
    data_cov = data['cov']
    
    data_r = data['r']
    data_info = data['info']

    fid_mean = fid['mean']
    fid_cov = fid['cov']

    diff_mean = (data_mean - fid_mean) / delta
    diff_cov = (data_cov - fid_cov) / delta

    numpy.savez(outfile, x=data_x, mean=diff_mean, cov=diff_cov, r=data_r, info=data_info)


def calc_fisher(fid_cov_output, dcov_outs, params, fisher_out, ranges):
    
    fid = numpy.load(fid_cov_output)

    mean = fid['mean']
    cov = fid['cov']
    
    r = fid['r']
    
    filter1 = []
    
    for n in ranges:
        r1 = r[n][0]
        r2 = r[n][1]
        filter1.extend(range(r1,r2))

    cov = cov - numpy.outer(mean, mean)

    mean = mean[filter1]
    cov = cov[filter1][:,filter1]

    filter2 = []

    for i in range(len(mean)):
        if cov[i][i] != 0.0:
            filter2.append(i)

    cov = cov[filter2][:,filter2]

    cov_inv = numpy.linalg.inv(cov)

    n = len(dcov_outs)

    diffs = []

    for dcov_out in dcov_outs:

        #print dcov_out
        data = numpy.load(dcov_out)
        diff = data['mean']
        #print numpy.shape(diff)
        diff = diff[filter1]
        diff = diff[filter2]
        #print diff
        #print numpy.shape(diff)
        diffs.append(diff)
        
    F = numpy.zeros([n,n])

    for i in range(n):
        for j in range(n):
            #print numpy.shape(diffs[i])
            #print numpy.shape(cov_inv)
            #print numpy.shape(diffs[j])

            #print diffs[i]
            #print cov_inv
            #print diffs[j]

            d = numpy.dot(diffs[i], numpy.dot(cov_inv, diffs[j]))
            #print d
            F[i][j] = d

    F_inv = numpy.linalg.inv(F)

    numpy.savez(fisher_out, F=F, F_inv=F_inv, params=params)
