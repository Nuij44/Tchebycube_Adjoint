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
    h5d = h5.File(h5db+list[0], 'r')
    ua,uz,ur =     np.transpose(h5d['/u1'][:,:,:]), \
                   np.transpose(h5d['/u2'][:,:,:]), \
                   np.transpose(h5d['/u3'][:,:,:])
    h5d.close()

    for i in range(1,np.size(list)):
        print("mean",i)
        h5d = h5.File(h5db+list[i], 'r')
        ua = ua + np.transpose(h5d['/u1'][:,:,:])
        uz = uz + np.transpose(h5d['/u2'][:,:,:])
        ur = ur + np.transpose(h5d['/u3'][:,:,:])

    return ua/np.size(list),uz/np.size(list),ur/np.size(list)


def RMS(h5list,r):

    mean_a,mean_z,mean_r = mean(h5list,r)

    A = -1.0/3.0
    B = 4.0/3.0
        
    #Computing fluctuation
    h5d = h5.File(h5db+h5list[0], 'r')
    ua,uz,ur =     np.transpose(h5d['/u1'][:,:,:]), \
                   np.transpose(h5d['/u2'][:,:,:]), \
                   np.transpose(h5d['/u3'][:,:,:])
    h5d.close()

    
    fa = np.square(-ua  + mean_a)
    fz = np.square(-uz + mean_z)
    fr = np.square(-ur + mean_r)
    
    for i in range(1,np.size(list)):
        h5d = h5.File(h5db+list[i], 'r')
        fa = fa + np.square(mean_a - np.transpose(h5d['/u1'][:,:,:]))
        fz = fz + np.square(mean_z - np.transpose(h5d['/u2'][:,:,:]))
        fr = fr + np.square(mean_r - np.transpose(h5d['/u3'][:,:,:]))

    return fa/np.size(list),fz/np.size(list),fr/np.size(list),mean_a,mean_z,mean_r


def RMS_V2(h5list,r):
    h5d = h5.File(h5db+h5list[0], 'r')
    ua,uz,ur =     np.transpose(h5d['/u1'][:,:,:]), \
                   np.transpose(h5d['/u2'][:,:,:]), \
                   np.transpose(h5d['/u3'][:,:,:])
    h5d.close()


    UA = ua#[ua for i in range(np.size(h5list))]
    UZ = uz
    UR = ur
    
    for i in range(1,np.size(h5list)):
        print("mean",i)
        h5d = h5.File(h5db+h5list[i], 'r')
        UA = UA + np.transpose(h5d['/u1'][:,:,:])
        UZ = UZ + np.transpose(h5d['/u2'][:,:,:])
        UR = UR + np.transpose(h5d['/u3'][:,:,:])


    mean_a = UA/np.size(h5list)
    mean_z = UZ/np.size(h5list)
    mean_r = UR/np.size(h5list)

    FA = np.zeros_like(UA)
    FZ = FA
    FR = FA
    
    for i in range(np.size(h5list)):
        print("load snap",i)
        h5d = h5.File(h5db+h5list[i], 'r')
        UA = np.transpose(h5d['/u1'][:,:,:])
        UZ = np.transpose(h5d['/u2'][:,:,:])
        UR = np.transpose(h5d['/u3'][:,:,:])

        FA = FA + (UA - mean_a)**2
        FZ = FZ + (UZ - mean_z)**2
        FR = FR + (UR - mean_r)**2


    fa = FA/np.size(h5list)
    fz = FZ/np.size(h5list)
    fr = FR/np.size(h5list)
        
    rms_a = ops.integrate(fa, axis=(0,1))
    rms_z = ops.integrate(fz, axis=(0,1))
    rms_r = ops.integrate(fr, axis=(0,1))

    return rms_a,rms_z,rms_r,mean_a,mean_z,mean_r
            
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
        xmax=[2*np.pi, np.pi, 2.0],
        n=np.shape(phi),
        bases=["fourier", "fourier","chebyshev"]
    )


    
    # Store files in a list Using [0-9] pattern
    h5dlist = []
    for files in glob.glob(h5db + '/*_*[0-9].h5'):
        h5dlist.append(os.path.basename(files))
    h5dlist.sort()

    '''
    rms_a,rms_z,rms_r = ArtBilson(h5dlist[300:600],r)
    A,Z,R = ops.nodes

    with open("RMS_a.dat","w") as f:
        for i in range(0,np.size(rms_a)):
            print(R[i],rms_a[i])
            f.write(f"{R[i]}    {rms_a[i]}\n")

    with open("RMS_z.dat","w") as f:
        for i in range(0,np.size(rms_z)):
            print(R[i],rms_z[i])
            f.write(f"{R[i]}    {rms_z[i]}\n")

    with open("RMS_r.dat","w") as f:
        for i in range(0,np.size(rms_r)):
            print(R[i],rms_r[i])
            f.write(f"{R[i]}    {rms_r[i]}\n")

    exit()
    '''    
    #Kinetic energy of the flow
    write_nrj(h5db,h5dlist,r,ops)
    #    exit()
    #Compute RMS
    print("Taille list:",np.size(h5dlist))
    #    print(h5dlist)
    #    print(h5dlist[-2:np.size(h5dlist)])
    rms_a,rms_z,rms_r,mean_a,mean_z,mean_r = RMS_V2(h5dlist[-30:np.size(h5dlist)],r)
    A,Z,R = ops.nodes

    with open("RMS_a.dat","w") as f:
        for i in range(np.size(rms_a)):
            print(R[i],rms_a[i])
            f.write(f"{R[i]}    {rms_a[i]}\n")

    with open("RMS_z.dat","w") as f:
        for i in range(np.size(rms_z)):
            print(R[i],rms_z[i])
            f.write(f"{R[i]}    {rms_z[i]}\n")

    with open("RMS_r.dat","w") as f:
        for i in range(np.size(rms_r)):
            print(R[i],rms_r[i])
            f.write(f"{R[i]}    {rms_r[i]}\n")

    pro_a = ops.integrate(mean_a, axis=(0,1))
    with open("mean_a.dat","w") as f:
        for i in range(np.size(pro_a)):
            print(R[i],pro_a[i])
            f.write(f"{R[i]}    {pro_a[i]}\n")

            
    exit()

        
    #integration sans jaccobienne car on intégre seulement dans theta et z

#    moy_a = np.mean(rms_a,axis=[0,1])
#    int_azi = np.mean(moy_a,axis=0)


    
    int_a = ops.integrate_x(rms_a)
    int_z = ops.integrate_x(rms_z)
    int_r = ops.integrate_x(rms_r)
    
    ops_2d = sp.SpectralDiscretization(
        xmin=[0.0, 1.0],
        xmax=[np.pi, 2.0],
        n=np.shape(phi)[1:3],
        bases=[ "fourier","chebyshev"])
        
    int_azi = ops_2d.integrate_x(int_a)
    int_ver = ops_2d.integrate_x(int_z)
    int_rad = ops_2d.integrate_x(int_r)
    
    A,Z,R = ops.nodes

    with open("RMS_a.dat","w") as f:
        for i in range(0,np.size(int_azi)):
            print(R[i],int_azi[i])
            f.write(f"{R[i]}    {int_azi[i]}\n")

    with open("RMS_z.dat","w") as f:
        for i in range(0,np.size(int_ver)):
            print(R[i],int_ver[i])
            f.write(f"{R[i]}    {int_ver[i]}\n")

    with open("RMS_r.dat","w") as f:
        for i in range(0,np.size(int_rad)):
            print(R[i],int_rad[i])
            f.write(f"{R[i]}    {int_rad[i]}\n")

    int_a = ops.integrate_x(mean_a)
    int_azi = ops_2d.integrate_x(int_a)
    with open("mean_a.dat","w") as f:
        for i in range(0,np.size(int_azi)):
            print(R[i],int_azi[i])
            f.write(f"{R[i]}    {int_azi[i]}\n")

    A = -1.0/3.0
    B = 4.0/3.0
    phi = A*r + B/r
            
    int_a = ops.integrate_x(phi)
    int_az = ops_2d.integrate_x(int_a)
    with open("laminar_a.dat","w") as f:
        for i in range(0,np.size(int_az)):
            print(R[i],int_az[i])
            f.write(f"{R[i]}    {int_az[i]}\n")

    
