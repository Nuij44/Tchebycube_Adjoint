import os
import glob
import argparse
import numpy as np
import h5py as h5
from pyevtk.hl import gridToVTK
#TODO: from tecio import *
#python3 h5conv_cyl.py -vtk tchebycube/disk_accretion/outdir/snapshots/

parser = argparse.ArgumentParser()

parser.add_argument("h5db")
parser.add_argument("-vtk", action="store_true")
parser.add_argument("-plt", action="store_true")

args = parser.parse_args()

h5mesh = glob.glob(args.h5db+"/*_grid.h5")[0]
h5db = args.h5db


def CheckDir(target_dir):
    if not os.path.isdir(target_dir):
        print("The target directory", target_dir, "doesn't exist")
        os.mkdir(target_dir)

def h52vtkmesh(h5mesh):
# --- Read and Re-order h5 file h5datas ---                                                                   
    h5m = h5.File(h5mesh, 'r')
    
    nb_snap, = h5m.attrs['id'][:]
    tc, = h5m.attrs['tc'][:]

    phi,z,r = np.transpose(h5m['/x1'][:,:,:]), \
            np.transpose(h5m['/x2'][:,:,:]), \
            np.transpose(h5m['/x3'][:,:,:])

    h5m.close()

    Nx = len(phi[:,1,1])
    Ny = len(phi[1,:,1])
    Nz = len(phi[1,1,:])

    phi_aug = np.zeros((Nx+1,Ny,Nz),dtype=float)
    phi_aug[0:Nx,:,:] = phi[:,:,:]
    phi_aug[Nx,:,:] = 2.*np.pi*np.ones((Ny,Nz),dtype=float)

    z_aug = np.zeros((Nx+1,Ny,Nz),dtype=float)
    z_aug[0:Nx,:,:] = z[:,:,:]
    z_aug[Nx,:,:] = z[0,:,:]
    
    r_aug = np.zeros((Nx+1,Ny,Nz),dtype=float)
    r_aug[0:Nx,:,:] = r[:,:,:]
    r_aug[Nx,:,:] = r[0,:,:]
    
#    print("len Phi : ",len(phi[:,1,1]))

#    for i in range(len(phi[:,1,1])):
#        print("Phi(",i,") = ",phi[i,1,1])
    
    x = r_aug * np.cos(phi_aug)
    y = r_aug * np.sin(phi_aug)
    
    return x, y, z_aug

def h52vtk(dir, h5data, x, y, z):

    CheckDir(dir)
    stride_x,stride_y,stride_z = 1,1,1

    h5d = h5.File(h5data, 'r')
    T,p,u1,u2,u3 = np.transpose(h5d['/T'][:,:,:]), \
                   np.transpose(h5d['/p'][:,:,:]), \
                   np.transpose(h5d['/u1'][:,:,:]), \
                   np.transpose(h5d['/u2'][:,:,:]), \
                   np.transpose(h5d['/u3'][:,:,:])
    h5d.close()

    Nx = len(u1[:,1,1])
    Ny = len(u1[1,:,1])
    Nz = len(u1[1,1,:])

    u1_aug = np.zeros((Nx+1,Ny,Nz),dtype=float)
    u1_aug[0:Nx,:,:] = u1[:,:,:]
    u1_aug[Nx,:,:] = u1[0,:,:]

    u2_aug = np.zeros((Nx+1,Ny,Nz),dtype=float)
    u2_aug[0:Nx,:,:] = u2[:,:,:]
    u2_aug[Nx,:,:] = u2[0,:,:]

    u3_aug = np.zeros((Nx+1,Ny,Nz),dtype=float)
    u3_aug[0:Nx,:,:] = u3[:,:,:]
    u3_aug[Nx,:,:] = u3[0,:,:]

    p_aug = np.zeros((Nx+1,Ny,Nz),dtype=float)
    p_aug[0:Nx,:,:] = p[:,:,:]
    p_aug[Nx,:,:] = p[0,:,:]

    T_aug = np.zeros((Nx+1,Ny,Nz),dtype=float)
    T_aug[0:Nx,:,:] = T[:,:,:]
    T_aug[Nx,:,:] = T[0,:,:]

    
    x,y,z = x[::stride_x,::stride_y,::stride_z], \
            y[::stride_x,::stride_y,::stride_z], \
            z[::stride_x,::stride_y,::stride_z]
    u1_aug,u2_aug,u3_aug = u1_aug[::stride_x,::stride_y,::stride_z], \
        u2_aug[::stride_x,::stride_y,::stride_z], \
        u3_aug[::stride_x,::stride_y,::stride_z]
    
    nx,ny,nz = np.shape(x)
    vtkfile = os.path.basename(h5data)
    vtkfile, ext = os.path.splitext(vtkfile)
    gridToVTK(dir+"/"+vtkfile, x,y,z, pointData = {"T":T_aug,"p":p_aug,"u1":u1_aug,"u2":u2_aug,"u3":u3_aug})    


def h52plt(dir, h5mesh, h5data):
    CheckDir(dir)
    stride_x,stride_y,stride_z = 1,1,1

    # --- Read and Re-order h5 file h5datas ---
    h5db = h5.File(h5mesh, 'r')
    nb_snap, = h5db.attrs['id'][:]
    tc, = h5db.attrs['tc'][:]
    
    x,y,z = np.transpose(h5db['/x1'][:,:,:]), \
            np.transpose(h5db['/x2'][:,:,:]), \
            np.transpose(h5db['/x3'][:,:,:])

    h5db = h5.File(h5data, 'r')
    T,p,u1,u2,u3 = np.transpose(h5db['/T'][:,:,:]), \
                   np.transpose(h5db['/p'][:,:,:]), \
                   np.transpose(h5db['/u1'][:,:,:]), \
                   np.transpose(h5db['/u2'][:,:,:]), \
                   np.transpose(h5db['/u3'][:,:,:])
    h5db.close()
    
    x,y,z = x[::stride_x,::stride_y,::stride_z], \
            y[::stride_x,::stride_y,::stride_z], \
            z[::stride_x,::stride_y,::stride_z]
    u1,u2,u3 = u1[::stride_x,::stride_y,::stride_z], \
               u2[::stride_x,::stride_y,::stride_z], \
               u3[::stride_x,::stride_y,::stride_z]
    
    nx,ny,nz = np.shape(x)
    pltfile, ext = os.path.splitext(h5data)
    open_file("dir"+pltfile, 'Title', ['x','y','z','T','p','u1','u2','u3'])
    create_ordered_zone('Zone', (nz,ny,nx) )
    zone_write_values(x) 
    zone_write_values(y)
    zone_write_values(z) 
    zone_write_values(u1)
    zone_write_values(u2)
    close_file()




# ---- Main ----

# --- Read mesh and datas than convert to VTK or PLT format ---
# --- h5 dirpath

if args.plt:
    h52plt("./plt/", h5mesh, h5data)
elif args.vtk:
# Store files in a list Using [0-9] pattern
    h5dlist = []
    for files in glob.glob(h5db + '/*_*[0-9].h5'):
        h5dlist.append(os.path.basename(files))
    h5dlist.sort()

# ---- Read h5 Mesh ----
x, y, z = h52vtkmesh(h5mesh)
# ---- Convert h5 to vtk ----
for f in h5dlist:
    h52vtk("./vtk", h5db+"/"+f, x, y, z)
