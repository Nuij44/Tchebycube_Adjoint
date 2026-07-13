import os
import glob
import argparse
import numpy as np
import h5py as h5
import spectral_tools as sp
#TODO: from tecio import *
#python3 h5conv_cyl.py -dir tchebycube/disk_accretion/outdir/snapshots/

parser = argparse.ArgumentParser()

parser.add_argument("h5db")
parser.add_argument("-dir", action="store_true")
#parser.add_argument("h5mean")
#parser.add_argument("-dir_mean", action="store_true")

args = parser.parse_args()

h5mesh = glob.glob(args.h5db+"/*_grid.h5")[0]
h5db = args.h5db
#h5mean = args.h5mean


def CheckDir(target_dir):
    if not os.path.isdir(target_dir):
        print("The target directory", target_dir, "doesn't exist")
        os.mkdir(target_dir)



# ---- Main ----

if __name__ == "__main__":

    CheckDir(h5db)
    
    #Read grid
    h5m = h5.File(h5mesh, 'r')
    
    nb_snap, = h5m.attrs['id'][:]
    tc, = h5m.attrs['tc'][:]

    x,y,z = np.transpose(h5m['/x1'][:,:,:]), \
            np.transpose(h5m['/x2'][:,:,:]), \
            np.transpose(h5m['/x3'][:,:,:])

    h5m.close()    

    
    ops = sp.SpectralDiscretization(
        xmin=[0.0, 0.0, -1.0],
        xmax=[4.35*np.pi,1.05*np.pi, 1.0],
        n=np.shape(x),
        bases=["fourier", "fourier","chebyshev"]
    )

    print(np.shape(x))

    
    # Store files in a list Using [0-9] pattern
    h5dlist = []
    for files in glob.glob(h5db + '/*_*[0-9].h5'):
        h5dlist.append(os.path.basename(files))
    h5dlist.sort()

    
    ope = sp.SpectralDiscretization(
        xmin=[0.0, 0.0, -1.0],
        xmax=[1.75*np.pi, 1.2*np.pi, 1.0],
        n=[129,97,257],
        bases=["fourier", "fourier","chebyshev"]
    )
    i=0
    print(h5db+h5dlist[i])

    h5d = h5.File(h5db+h5dlist[i], 'r')
    ua,uz,ur =     np.transpose(h5d['/u1'][:,:,:]), \
                   np.transpose(h5d['/u2'][:,:,:]), \
                   np.transpose(h5d['/u3'][:,:,:])
    h5d.close()

    interp = sp.SpectralInterpolate(ops, ope)

    ua_f = interp @ ua
    uz_f = interp @ uz
    ur_f = interp @ ur

    outfile=h5.File('init_interpol.dat','a')

    init_a = outfile.create_dataset(name='/u1', data=np.transpose(ua_f), dtype=np.float64)
    init_z = outfile.create_dataset(name='/u2', data=np.transpose(uz_f), dtype=np.float64)
    init_r = outfile.create_dataset(name='/u3', data=np.transpose(ur_f), dtype=np.float64)

    outfile.close()
