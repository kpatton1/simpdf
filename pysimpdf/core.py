
PINOCCHIO_CORE_SURVEY_DEFAULT_NAME='test'
PINOCCHIO_CORE_SURVEY_DEFAULT_NSIDE=8192
PINOCCHIO_CORE_SURVEY_DEFAULT_SEED=0
PINOCCHIO_CORE_SURVEY_DEFAULT_E=0.30
PINOCCHIO_CORE_SURVEY_DEFAULT_NG=10
PINOCCHIO_CORE_SURVEY_DEFAULT_ZL=0.4
PINOCCHIO_CORE_SURVEY_DEFAULT_ZS=0.8
PINOCCHIO_CORE_SURVEY_DEFAULT_Q=1.0
PINOCCHIO_CORE_SURVEY_DEFAULT_G=1.0

PINOCCHIO_CORE_ANALYSIS_DEFAULT_NAME='test'
PINOCCHIO_CORE_ANALYSIS_DEFAULT_ALIST=[]
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

    def __init__(self, name=PINOCCHIO_CORE_SURVEY_DEFAULT_NAME, nside=PINOCCHIO_CORE_SURVEY_DEFAULT_NSIDE, seed=PINOCCHIO_CORE_SURVEY_DEFAULT_SEED, e=PINOCCHIO_CORE_SURVEY_DEFAULT_E, ng=PINOCCHIO_CORE_SURVEY_DEFAULT_NG, zl=PINOCCHIO_CORE_SURVEY_DEFAULT_ZL, zs=PINOCCHIO_CORE_SURVEY_DEFAULT_ZS, q=PINOCCHIO_CORE_SURVEY_DEFAULT_Q, g=PINOCCHIO_CORE_SURVEY_DEFAULT_G):
        
        self.name = name
        self.nside = nside
        self.seed = seed
        self.e = e
        self.ng = ng

        self.zl = zl
        self.zs = zs

        self.q = q
        self.g = g

class Analysis:
    
    def __init__(self, name=PINOCCHIO_CORE_ANALYSIS_DEFAULT_NAME, alist=PINOCCHIO_CORE_ANALYSIS_DEFAULT_ALIST, divs=PINOCCHIO_CORE_ANALYSIS_DEFAULT_DIVS):

        self.name = name
        self.alist_orig = alist
        self.divs = divs

        self.alist = []

        for a in alist:

            n = a[0]
            t = a[1]
            
            if t == 'pdf':
                map_size = n
                
                bins = a[2]
                bin_min = a[3]
                bin_max = a[4]

                a_pdf = AnalysisPDF(n, self.divs, bins, bin_min, bin_max)
                a_pdf.info = a

                self.alist.append(a_pdf)

            elif t == 'ps':
                map_size = n
                
                lmin = 0
                lmax = 3 * map_size
                
                a_ps = AnalysisPS(n, self.divs, lmin, lmax)
                a_ps.info = a

                self.alist.append(a_ps)
            
            elif t == 'moments':
                map_size = n
                
                scale = a[2]
                mmin = a[3]
                mmax = a[4]
            
                a_moments = AnalysisMoments(n, self.divs, scale, mmin, mmax)
                a_moments.info = a
            
                self.alist.append(a_moments)

            else:
                print 'Error! Unknown analysis type: ' + str(t)
                print a

    def process(self, map_data_raw):

        scaled_maps = {}

        divs = self.divs
        alist = self.alist

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
            
            if map_size not in scaled_maps.keys():
                scaled_maps[map_size] = healpy.pixelfunc.ud_grade(map_data_raw,map_size,power=-2,order_in='NESTED',order_out='NESTED',dtype=numpy.float64)

        for n in range(divs):
        
            print 'div ' + str(n)
        
            measurement = []
        
            acount = 0
            
            for a in alist:
                print 'analysis ' + str(acount)
                map_size = a.map_size
                npix = healpy.pixelfunc.nside2npix(map_size)
                chunk = scaled_maps[map_size][npix/divs*n:npix/divs*(n+1)]
                chunk = chunk - numpy.mean(chunk)

                a_m = a.process_chunk(chunk)

                measurement.extend(a_m)
                
                acount += 1

            measurements.append(measurement)

        return measurements,x,ranges,info


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

    def process_chunk(self, data):

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

    def process_chunk(self, data):
        
        hist,bin_edges = numpy.histogram(data, bins=self.bins, range=(self.bin_min,self.bin_max))

        return hist

class AnalysisMoments:

    def __init__(self, map_size, divs, scale, mmin=2, mmax=10):
        self.map_size = map_size
        self.divs = divs
        self.scale = scale
        self.mmin = mmin
        self.mmax = mmax

    def x(self):
        x = numpy.arange(self.mmin, self.mmax)

        return x

    def process_chunk(self, data):

        moment_data = data / self.scale
        
        moments = numpy.zeros(self.mmax - self.mmin, dtype=numpy.float32)

        if self.mmin == 2:
            temp = moment_data
        else:
            temp = numpy.pow(moment_data, self.mmin-1)

        for i in range(0, self.mmax-self.mmin):
            temp = temp * moment_data
            moments[i] = numpy.mean(temp)

        return moments

class AnalysisLogPDF:

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

    def process_chunk(self, data):
        
        hist,bin_edges = numpy.histogram(data, bins=self.bins, range=(self.bin_min,self.bin_max))

        return hist

class AnalysisLogMoments:

    def __init__(self, map_size, divs, scale, mmin=2, mmax=10):
        self.map_size = map_size
        self.divs = divs
        self.scale = scale
        self.mmin = mmin
        self.mmax = mmax

    def x(self):
        x = numpy.arange(self.mmin, self.mmax)

        return x

    def process_chunk(self, data):

        moment_data = data / self.scale
        
        moments = numpy.zeros(self.mmax - self.mmin, dtype=numpy.float32)

        if self.mmin == 2:
            temp = moment_data
        else:
            temp = numpy.pow(moment_data, self.mmin-1)

        for i in range(0, self.mmax-self.mmin):
            temp = temp * moment_data
            moments[i] = numpy.mean(temp)

        return moments

class Simulation:

    
    def pinocchio_output(self,i):
        pinocchio_output = self.basedir + '/' + self.cosmo.name + '/r' + str(i) + '/pinocchio.' + self.cosmo.name + '_r' + str(i) + '.plc.out'
        return pinocchio_output

    def healpix_output(self,i):
        healpix_output = self.basedir + '/' + self.cosmo.name + '/r' + str(i) + '/healpix_' + str(self.survey.nside) + '_' + self.cosmo.name + '_r' + str(i) + '.fits'
        return healpix_output

    def noise_output(self,i,j):
        noise_output = self.basedir + '/' + self.cosmo.name + '/r' + str(i) + '/noise_' + str(self.survey.nside) + '_' + self.cosmo.name + '_' + self.survey.name + '_r' + str(i) + 'n' + str(j) + '.fits'
        return noise_output

    def measure_output(self,i,j):
        measure_output = self.basedir + '/' + self.cosmo.name + '/r' + str(i) + '/measure_' + str(self.survey.nside) + '_' + self.cosmo.name + '_' + self.survey.name + '_' + self.analysis.name + '_r' + str(i) + 'n' + str(j) + '.npz'
        return measure_output

    def combined_measure_output(self):
        combined_measure_output = self.basedir + '/' + self.cosmo.name + '/measure_' + str(self.survey.nside) + '_' + self.cosmo.name + '_' + self.survey.name + '_' + self.analysis.name + '.npz'
        return combined_measure_output

    def cov_output(self):
        cov_output = self.basedir + '/' + self.cosmo.name + '/cov_' + str(self.survey.nside) + '_' + self.cosmo.name + '_' + self.survey.name + '_' + self.analysis.name + '.npz'
        return cov_output

    def dcov_output(self):
        dcov_output = self.basedir + '/dcov_' + str(self.survey.nside) + '_' + self.cosmo.name + '_' + self.survey.name + '_' + self.analysis.name + '.npz'
        return dcov_output

    def fisher_output(self,code):
        fisher_output = self.basedir + '/F_' + str(self.survey.nside) + '_' + self.analysis.name + '_' + code + '.npz'
        return fisher_output

    def cosmodir(self):
        cosmodir = self.basedir + '/' + self.cosmo.name
        return cosmodir

    def rundir(self,i):
        rundir = self.basedir + '/' + self.cosmo.name + '/r' + str(i)
        return rundir

    def check_dir(self,d):
        if not os.path.exists(d):
            os.mkdir(d)
        elif not os.path.isdir(d):
            print 'Error! Invalid directory: ' + d
            return True

        if not os.path.exists(d):
            print 'Error! Could not create directory: ' + d
            return True
        
        return False

    def check_basedir(self):

        d = self.basedir

        return self.check_dir(d)

    def check_cosmodir(self):

        if self.check_basedir():
            return True

        d = self.cosmodir()

        return self.check_dir(d)

    def check_rundir(self,i):

        if self.check_cosmodir():
            return True

        d = self.rundir(i)

        return self.check_dir(d)
        

    def __init__(self,basedir,t,alist,aname='fiducial'):
        
        self.basedir = os.path.abspath(basedir)
        
        self.t = t

        param = t[0]
        val = t[1]
        self.nsim = t[2]
        self.nnoise = t[3]
        
        self.param = param

        if param == 'h' or param == 'w' or param == 'om' or param == 'ok' or param == 's8':
            cname = 'delta_' + param + '_' + str(val)
        else:
            cname = 'fiducial'

        if param == 'e' or param == 'q' or param == 'g':
            sname = 'delta_' + param + '_' + str(val)
        else:
            sname = 'fiducial'

        #aname = 'fiducial'

        self.cosmo = Cosmology(name=cname)
        self.survey = Survey(name=sname)
        self.analysis = Analysis(name=aname,alist=alist)
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
        
        elif param == 'g':
            self.survey.g = val

        else:
            if param != 'f':
                print 'Simulation error! invalid parameter type: ' + param
            
    def run_pinocchio(self):
        
        print 'Running pinocchio for ' + str(self.cosmo.name)
        
        cwd = os.getcwd()
        
        for i in range(1,self.nsim+1):

            if self.check_rundir(i):
                os.chdir(cwd)
                return

            os.chdir(self.rundir(i))
                    
            if not os.path.exists('parameters') or not os.path.exists('outputs'):
                self.cosmo.seed = i
                generate_parameter_files(self.cosmo, self.cosmo.name + '_r' + str(i),self.simparams)

            if not os.path.exists('parameters') or not os.path.exists('outputs'):
                print 'Pinocchio error! could not create parameter files'
                os.chdir(cwd)
                return

            pinocchio_output = self.pinocchio_output(i)

            if not os.path.exists(pinocchio_output):
                run_pinocchio_simulation()

            if not os.path.exists(pinocchio_output):
                print 'Pinocchio error! could not generate output: ' + str(pinocchio_output)
                os.chdir(cwd)
                return

        os.chdir(cwd)

    def convert_healpix(self):
        
        print 'Converting to healpix for ' + str(self.cosmo.name)
        
        for i in range(1,self.nsim+1):
            
            pinocchio_output = self.pinocchio_output(i)

            if not os.path.exists(pinocchio_output):
                print 'Healpix error! pinocchio output does not exist: ' + pinocchio_output
                return

            healpix_output = self.healpix_output(i)

            if not os.path.exists(healpix_output):
                healpixify(pinocchio_output,healpix_output,self.survey.nside)
                
            if not os.path.exists(healpix_output):
                print 'Healpix error! could not generate: ' + healpix_output
                return

    def add_noise(self):
        
        print 'Adding noise to ' + str(self.cosmo.name) + ' with survey ' + str(self.survey.name)
        
        for i in range(1,self.nsim+1):

            healpix_output = self.healpix_output(i)
                
            if not os.path.exists(healpix_output):
                print 'Noise error! healpix output does not exist: ' + healpix_output
                return

            for j in range(1,self.nnoise+1):
                noise_output = self.noise_output(i,j)
                
                if not os.path.exists(noise_output):
                    add_noise(healpix_output,noise_output,self.cosmo,self.survey)
                    
                if not os.path.exists(noise_output):
                    print 'Noise error! could not generate: ' + noise_output
                    return

    def measure_data(self):
        
        print 'Measuring data for ' + str(self.cosmo.name) + ' with survey ' + str(self.survey.name) + ' and analysis ' + str(self.analysis.name)
        
        for i in range(1,self.nsim+1):
            
            for j in range(1,self.nnoise+1):
                
                
                noise_output = self.noise_output(i,j)

                if not os.path.exists(noise_output):
                    print 'Measurement error! noise output does not exist: ' + noise_output
                    return

                measure_output = self.measure_output(i,j)
                
                if not os.path.exists(measure_output):
                    measure_data(noise_output, measure_output, self.analysis)
                
                if not os.path.exists(measure_output):
                    print 'Measurement error! could not generate: ' + measure_output
                    return

    def combine_measures(self):

        print 'Combining measures for ' + str(self.cosmo.name) + ' with survey ' + str(self.survey.name) + ' and analysis ' + str(self.analysis.name)

        measure_outputs = []

        for i in range(1,self.nsim+1):
            
            for j in range(1,self.nnoise+1):

                measure_output = self.measure_output(i,j)

                if not os.path.exists(measure_output):
                    print 'Combining error! measure output does not exist: ' + measure_output
                    return
                
                measure_outputs.append(measure_output)

        combined_measure_output = self.combined_measure_output()


        if not os.path.exists(combined_measure_output):
            combine_measures(measure_outputs,combined_measure_output)

        if not os.path.exists(combined_measure_output):
            print 'Combining error! could not generate: ' + combined_measure_output

    def generate_covariances(self):
        
        print 'Generating covariances for ' + str(self.cosmo.name) + ' with survey ' + str(self.survey.name) + ' and analysis ' + str(self.analysis.name)
        
        combined_measure_output = self.combined_measure_output()

        if not os.path.exists(combined_measure_output):
            print 'Covariance error! combined measurement output does not exist: ' + combined_measure_output
            return

        cov_output = self.cov_output()

        if not os.path.exists(cov_output):
            calc_covariance(combined_measure_output, cov_output, self.survey)

        if not os.path.exists(cov_output):
            print 'Covariance error! could not generate: ' + cov_output

    def diff_covariances(self, fid):
        
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
        
        elif param == 'g':
            delta = self.survey.g - fid.survey.g

        else:
            print 'Covariance differencing error! invalid parameter type: ' + param
            return

        delta_cov_output = self.cov_output()
        fid_cov_output = fid.cov_output()

        if not os.path.exists(delta_cov_output):
            print 'Covariance differencing error! delta cov output does not exist: ' + delta_cov_output
            return

        if not os.path.exists(fid_cov_output):
            print 'Covariance differencing error! fiducial cov output does not exist: ' + fid_cov_output
            return

        dcov_output = self.dcov_output()
        
        if not os.path.exists(dcov_output):
            diff_covariances(delta_cov_output, fid_cov_output, dcov_output, delta)
                
        if not os.path.exists(dcov_output):
            print 'Covariance differencing error! could not generate: ' + dcov_output

    def calc_fisher(self, deltas, ranges):
        
        fid_cov_output = self.cov_output()

        if not os.path.exists(fid_cov_output):
            print 'Fisher error! fiducial cov output does not exist: ' + fid_cov_output
            return

        dcov_outs = []
        params = []

        for delta in deltas:
            dcov_out = delta.dcov_output()
            
            if not os.path.exists(dcov_out):
                print 'Fisher error! delta cov output does not exist: ' + dcov_out
                return

            param = delta.t[0]

            dcov_outs.append(dcov_out)
            params.append(param)

        code = ''

        for r in ranges:
            code = code + str(r)

        fisher_out = self.fisher_output(code)

        calc_fisher(fid_cov_output, dcov_outs, params, fisher_out, ranges)

    def plot_analysis(self):
        pass

    def clean_measures(self):
    
        print 'Cleaning measurements for ' + str(self.cosmo.name) + ' with survey ' + str(self.survey.name) + ' and analysis ' + str(self.analysis.name)
        
        for i in range(1,self.nsim+1):
            
            for j in range(1,self.nnoise+1):

                measure_output = self.measure_output(i,j)
                
                if os.path.exists(measure_output):
                    os.remove(measure_output)

    def clean_covariances(self):
    
        print 'Cleaning covariances for ' + str(self.cosmo.name) + ' with survey ' + str(self.survey.name) + ' and analysis ' + str(self.analysis.name)
      
        cov_output = self.cov_output()

        if os.path.exists(cov_output):
            os.remove(cov_output)

    def clean_diff_covariances(self):

        print 'Cleaning diff covariances for ' + str(self.cosmo.name) + ' with survey ' + str(self.survey.name) + ' and analysis ' + str(self.analysis.name)
        
        dcov_output = self.dcov_output()

        if os.path.exists(dcov_output):
            os.remove(dcov_output)
        

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

def healpixify(infile,outfile,survey):

    print 'Converting pinocchio output to healpix: ' + infile + ' -> ' + outfile

    nside = survey.nside
    zs = survey.zs

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
        zi = z[i]
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

    pixelarea = healpy.pixelfunc.nside2pixarea(nside) * 180.0 / numpy.pi * 180.0 / numpy.pi * 60.0 * 60.0
    
    sigma = survey.q * survey.e * e_crit / 2.0 * pixelarea / numpy.sqrt(2.0 * n_gal * pixelarea)
        
    map_data = survey.g * (map_data + sigma * numpy.random.randn(12*nside*nside))
    
    healpy.fitsfunc.write_map(outfile, map_data)

def measure_data(infile,outfile,analysis):

    map_data_raw = healpy.fitsfunc.read_map(infile, dtype=numpy.float32)

    m,x,r,i = analysis.process(map_data_raw)

    numpy.savez(outfile,x=x,m=m,r=r,i=i)

def combine_measures(infiles,outfile):

    first = True

    m = []

    for f in infiles:
        
        data = numpy.load(f)

        temp_m = data['m']

        if first:
            x = data['x']
            r = data['r']
            i = data['i']

        m.extend(temp_m)

    numpy.savez(outfile,x=x,m=m,r=r,i=i)

def calc_covariance(infile,outfile,survey):
    
    data = numpy.load(infile)
    
    x = data['x']
    m = data['m']
    r = data['r']
    i = data['i']
    
    cov = numpy.cov(m,rowvar=0)
    mean = numpy.mean(m,axis=0)
        
    numpy.savez(outfile,x=x,mean=mean,cov=cov,r=r,i=i)

def diff_covariances(infile, fidfile, outfile, delta):

    data = numpy.load(infile)
    fid = numpy.load(fidfile)

    data_x = data['x']

    data_mean = data['mean']
    data_cov = data['cov']
    
    data_r = data['r']
    data_i = data['i']

    fid_mean = fid['mean']
    fid_cov = fid['cov']

    diff_mean = (data_mean - fid_mean) / delta
    diff_cov = (data_cov - fid_cov) / delta

    numpy.savez(outfile, x=data_x, mean=diff_mean, cov=diff_cov, r=data_r, i=data_i)


def calc_fisher(fid_cov_output, dcov_outs, params, fisher_out, ranges):
    
    fid = numpy.load(fid_cov_output)

    mean = fid['mean']
    cov = fid['cov']
    
    r = fid['r']

    info = fid['i']
    
    filter1 = []
    
    for n in ranges:
        r1 = r[n][0]
        r2 = r[n][1]

        t = info[n][1]

        if t == 'pdf':
            for i in range(r1,r2):
                if cov[i][i] != 0.0 and mean[i] >= 50.0/(200.0*384.0):
                    filter1.append(i)
        else:
            filter1.extend(range(r1,r2))

#    filter2 = []

#    for i in range(len(mean)):
#        if cov[i][i] != 0.0 and mean[i] >= 50.0/(200.0*384.0):
#            filter2.append(i)

#    cov = cov - numpy.outer(mean, mean)

    mean = mean[filter1]
    cov = cov[filter1][:,filter1]

#    cov = cov[filter2][:,filter2]

    #e,v = numpy.linalg.eigh(cov)

    sd = numpy.zeros(len(cov), dtype=numpy.float64)

    for i in range(len(cov)):
       sd[i] = numpy.sqrt(cov[i][i])

    for i in range(len(cov)):
       for j in range(len(cov)):
            cov[i][j] = cov[i][j] / sd[i] / sd[j]

    #shrink = 0.00000000001
    #shrink = 1e-16
    #shrink = 1e-28

    #shrink = 1e-13
#    shrink = 0.0
    #cov = cov + shrink * numpy.identity(len(cov),dtype=numpy.float64)

    cov_inv = numpy.linalg.pinv(cov)

    #cov_inv_2 = numpy.dot(cov_inv, cov_inv)

    #cov_inv_3 = numpy.dot(cov_inv_2, cov_inv)

    #cov_inv = cov_inv + shrink * cov_inv_2

    n = len(dcov_outs)

    diffs = []

    for dcov_out in dcov_outs:

        data = numpy.load(dcov_out)
        diff = data['mean']

        diff = diff[filter1]
#        diff = diff[filter2]
        
        for i in range(len(diff)):
            diff[i] = diff[i] / sd[i]

        diffs.append(diff)
        
    F = numpy.zeros([n,n])

    for i in range(n):
        for j in range(n):

            d = 0.0

            '''
            for l in range(len(e)):
                el = e[l]

                vl = v[l]

                d = d + numpy.dot(diffs[i], vl) * numpy.dot(diffs[j], vl) / el
            '''
            d = numpy.dot(diffs[i], numpy.dot(cov_inv, diffs[j]))

            F[i][j] = d

    F_inv = numpy.linalg.inv(F)

    numpy.savez(fisher_out, F=F, F_inv=F_inv, params=params)
