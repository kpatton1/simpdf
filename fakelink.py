
import os

sourcedir=os.path.abspath('data')
linkdir='data_real'

subdir_fid = 'fiducial'

subdir_deltas = ['delta_h_0.76','delta_ok_0.1', 'delta_om_0.32', 'delta_s8_0.9','delta_w_-0.9']

fiducial_sims = 200
delta_sims = 50

if not os.path.exists(linkdir):
    os.mkdir(linkdir)

if not os.path.exists(linkdir+'/'+subdir_fid):
    os.mkdir(linkdir+'/'+subdir_fid)

for i in range(1,fiducial_sims+1):
    rundir = subdir_fid + '/r' + str(i)

    if not os.path.exists(linkdir + '/' + rundir):
        os.mkdir(linkdir + '/' + rundir)

    parameters = 'parameters'
    outputs = 'outputs'
    pinocchio = 'pinocchio.' + subdir_fid + '_r' + str(i) + '.plc.out'
    healpix = 'healpix_8192_' + subdir_fid + '_r' + str(i) + '.fits'

    #files = [parameters, outputs, pinocchio, healpix]
    files = [parameters, outputs, pinocchio]
    
    for f in files:
        if not os.path.exists(linkdir + '/' + rundir + '/' + f):
            os.symlink(sourcedir + '/' + rundir + '/' + f, linkdir + '/' + rundir + '/' + f)

for subdir_delta in subdir_deltas:
    if not os.path.exists(linkdir + '/' + subdir_delta):
        os.mkdir(linkdir + '/' + subdir_delta)
    
for i in range(1,delta_sims+1):
    for subdir_delta in subdir_deltas:
        rundir = subdir_delta + '/r' + str(i)

        if not os.path.exists(linkdir + '/' + rundir):
            os.mkdir(linkdir + '/' + rundir)

        parameters = 'parameters'
        outputs = 'outputs'
        pinocchio = 'pinocchio.' + subdir_delta + '_r' + str(i) + '.plc.out'
        healpix = 'healpix_8192_' + subdir_delta + '_r' + str(i) + '.fits'
        
        #files = [parameters, outputs, pinocchio, healpix]
        files = [parameters, outputs, pinocchio]
    
        for f in files:
            if not os.path.exists(linkdir + '/' + rundir + '/' + f):
                os.symlink(sourcedir + '/' + rundir + '/' + f, linkdir + '/' + rundir + '/' + f)
