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


def Norme2(ops,h5snap,R):
    h5d = h5.File(h5db+h5snap, 'r')
    ua,uz,ur =     np.transpose(h5d['/u1'][:,:,:]), \
                   np.transpose(h5d['/u2'][:,:,:]), \
                   np.transpose(h5d['/u3'][:,:,:])
    h5d.close()

    AL2 = np.sqrt(ops.inner(ua,ua*R))
    ZL2 = np.sqrt(ops.inner(uz,uz*R))
    RL2 = np.sqrt(ops.inner(ur,ur*R))

    return AL2+ZL2+RL2
    
def write_nrj(dir,list,R,ops):
    with open("Nrj_L2.dat","w") as f:
        for i in range(1,np.size(list)):
            print(i,Norme2(ops,list[i],R))
            f.write(f"{i}    {Norme2(ops,list[i],R)}\n")

           
def mean(list,r):
    h5d = h5.File(h5db+list[1], 'r')
    ua,uz,ur =     np.transpose(h5d['/u1'][:,:,:]), \
                   np.transpose(h5d['/u2'][:,:,:]), \
                   np.transpose(h5d['/u3'][:,:,:])
    h5d.close()

    A = -1.0/3.0
    B = 4.0/3.0
        
    ua = ua 
    
    for i in range(2,np.size(list)):
        print("mean",i)
        h5d = h5.File(h5db+list[i], 'r')
        ua = ua + np.transpose(h5d['/u1'][:,:,:]) 
        uz = uz + np.transpose(h5d['/u2'][:,:,:])
        ur = ur + np.transpose(h5d['/u3'][:,:,:])

    return ua/np.size(list),uz/np.size(list),ur/np.size(list)


def RMS(h5list,r):

    mean_a,mean_z,mean_r = mean(h5list,r)
    
    #Computing fluctuation
    h5d = h5.File(h5db+h5list[1], 'r')
    ua,uz,ur =     np.transpose(h5d['/u1'][:,:,:]), \
                   np.transpose(h5d['/u2'][:,:,:]), \
                   np.transpose(h5d['/u3'][:,:,:])
    h5d.close()

    A = -2.0/3.0
    B = 4.0/3.0
    
    fa = np.square(ua - mean_a)
    fz = np.square(uz - mean_z)
    fr = np.square(ur - mean_r)
    
    for i in range(2,np.size(list)):
        h5d = h5.File(h5db+list[i], 'r')
        fa = fa + np.square(mean_a - np.transpose(h5d['/u1'][:,:,:]))
        fz = fz + np.square(mean_z - np.transpose(h5d['/u2'][:,:,:]))
        fr = fr + np.square(mean_r - np.transpose(h5d['/u3'][:,:,:]))

    return fa/np.size(list),fz/np.size(list),fr/np.size(list),mean_a,mean_z,mean_r
    

# ---- Main ----

if __name__ == "__main__":

    CheckDir(h5db)
    
    #Read grid
    h5m = h5.File(h5mesh, 'r')
    
    nb_snap, = h5m.attrs['id'][:]
    tc, = h5m.attrs['tc'][:]

    phi,z,r = np.transpose(h5m['/x1'][:,:,:]), \
            np.transpose(h5m['/x2'][:,:,:]), \
            np.transpose(h5m['/x3'][:,:,:])

    h5m.close()    

    
    ops = sp.SpectralDiscretization(
        xmin=[0.0, 0.0, 1.0],
        xmax=[2*np.pi, 2*np.pi, 2.0],
        n=np.shape(phi),
        bases=["fourier", "fourier","chebyshev"]
    )

    print(np.shape(phi))

    
    # Store files in a list Using [0-9] pattern
    h5dlist = []
    for files in glob.glob(h5db + '/*_*[0-9].h5'):
        h5dlist.append(os.path.basename(files))
    h5dlist.sort()

    
    ope = sp.SpectralDiscretization(
        xmin=[0.0, 0.0, 1.0],
        xmax=[2*np.pi, 2*np.pi, 2.0],
        n=[257,197,97],
        bases=["fourier", "fourier","chebyshev"]
    )
    i=40
    print(h5dlist[i])
    
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
