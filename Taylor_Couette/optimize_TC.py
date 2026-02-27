#!/usr/bin/env python3
# Florence Marcotte & Yannick Ponty
# Mai 2025
# python3 optimize_SSP_VWI.py --np 16 --input run_SSP_VWI.in --iter 10 --output Result --exe_dir dev --E0 1e-6

import sys
import os
import os.path
import shutil
import time
from math import sqrt
import glob
import fileinput
from collections import OrderedDict
from multiprocessing import Pool
import multiprocessing
import argparse
import h5py
import numpy as np
import numpy.random as rd
import subprocess
from functools import partial
import lib_optimize
import Sphere_Grad_Descent as sgd

#####################################################################
class MyException(Exception):
    pass
#####################################################################
def OptionValue(pattern,**kwargs):
    try: 
        value=kwargs[pattern]
        return value
    except KeyError:
        raise MyException(pattern+' Option is not defined !')
###################################################################
def OptionValueD(pattern, default, **kwargs):
    try: 
        value=kwargs[pattern]
        return value
    except:
        print('Warning '+pattern+' option is not defined ! default='+str(default))
        return default
#########################################################################
def read_grad_strat(str_file,type='float64'):

    try:
        file = h5py.File('./'+str_file, 'r')
        grad_u1 = file['/save/u1'][()]
        grad_u2 = file['/save/u2'][()]
        grad_u3 = file['/save/u3'][()]
        
        #transposition des données pour avoir un ordre colonne major
        #data_read = [np.transpose(grad_u1, axes=(2,1,0)),np.transpose(grad_u2, axes=(2,1,0)),np.transpose(grad_u3, axes=(2,1,0))]
        
        data_read = [np.reshape(grad_u1, (N[0],N[1],N[2])),np.reshape(grad_u2, (N[0],N[1],N[2])),np.reshape(grad_u3, (N[0],N[1],N[2]))]
    except: 
        print("The current directory is "+os.getcwd())
        print("%s is not inside this directory ?" % str_file) 
        sys.exit(1)
       
    return data_read
####################################################################################
def write_init_h5(U,V,W,output_name):
    '''X :: condition initiale
       output_name :: nom du fichier h5 contenant la condition initiale
    '''
    outfile=h5py.File(output_name,'a')
    #Ajout des groupes pour l'écriture du fichier de condition initales
    groupe_dump_init = outfile.create_group('dump')

    dump = outfile['/dump']

    dump_u = dump.create_dataset(name='u1', data=U, dtype=np.float64)
    dump_v = dump.create_dataset(name='u2', data=V, dtype=np.float64)
    dump_w = dump.create_dataset(name='u3', data=W, dtype=np.float64)

    outfile.close()
####################################################################################
def write_init_strat_h5(U,V,W,T,output_name):
    '''X :: condition initiale
       output_name :: nome du fichier h5 contenant la condition initiale
    '''
    outfile=h5py.File(output_name,'a')
    #Ajout des groupes pour l'écriture du fichier de condition initales
    groupe_dump_init = outfile.create_group('dump')

    dump = outfile['/dump']

    dump_u = dump.create_dataset(name='u1', data=U, dtype=np.float64)
    dump_v = dump.create_dataset(name='u2', data=V, dtype=np.float64)
    dump_w = dump.create_dataset(name='u3', data=W, dtype=np.float64)
    dump_w = dump.create_dataset(name='T' , data=T, dtype=np.float64)

    outfile.close()
###############################################################################
def last_energy_file(energyfile):
    time, energy, enstrophy= np.loadtxt(energyfile,unpack=True,skiprows=1)
    return energy[-1]
####################################################################################
def sum_energy_file(energyfile):
    time, energy = np.loadtxt(energyfile,unpack=True,skiprows=0)
    sum=np.sum(energy)*time[0] 
    print(sum)
    return sum
####################################################################################
def fun_tchebycube(vector, **kwargs):
    """[Give back J with input vector and a DNS simulation of Tchebycube
        char : [char] character of vector 
        nprocs: [int] number of proc for running cubby
        sum_energy : [bool] J is computed with integral of the energies along time
        
    Args:
        vector ([float vector]): [input vector]
    Returns:
        [float]: [J cost function]
    """

    nprocs = options.np    
    input_dir = options.input_fun
    os.mkdir('data_fun')
    
    write_init_h5(vector[0],vector[1],vector[2],'data_fun/init.h5')
 
    command_base="./"+options.exe_dir+"/DAL_TC.x --input "+input_dir+" --output-dir data_fun --init-data data_fun/init.h5 > resu_fun.txt"

    command="srun -n "+str(nprocs)+" "+command_base                                                                
    print(command)   
    
    os.system(command)
    J2=sum_energy_file('data_fun/timevar') #calcul v^2 sans facteur 1/2
    #print("DNS J2 :",J2)    
    os.system("rm -fr data_fun resu_fun_*")
    J = np.ravel(J2)
    return J
######################################################################################
def static_vars(**kwargs):
    def decorate(func):
        for k in kwargs:
            setattr(func, k, kwargs[k])
        return func
    return decorate
################################################################################################    
@static_vars(counter=-1)
def count():
    count.counter += 1
    return count.counter
###############################################################################################
def jac_tchebycube(vector, **kwargs):
    """[function producing the gradient (adjoint) by DNS-Adjoint compute by Tchebycube
        nprocs: [int] number of proc for running cubby
        sum_energy : [bool] J is computed with integral of the energies along time
    
    Args:
        vector ([float array]): [input vector to optimize]
    Returns:
        [float array]: [return adjoint of vector (gradient )]
    """ 
    
    nprocs = options.np
    input_dir = options.input_jac
    os.mkdir('data_fun')

    write_init_h5(vector[0],vector[1],vector[2],'data_fun/init.h5')

    #Run avec looping adjoint  et sauvegarde des iterations
    command_base="./"+options.exe_dir+"/DAL_TC.x --input "+input_dir+" --output-dir data_fun --init-data data_fun/init.h5 --do-adjoint > resu_jac_dns.txt"

    command="srun -n "+str(nprocs)+" "+command_base
    
    print(command)   
    os.system(command)

    grad_adj=read_grad_strat("data_fun/grad.h5")

    J2=sum_energy_file('data_fun/timevar') #Sans le facteur 1/2
    J = np.ravel(J2)
    #print("Jacobienne J2 :",J2)    
    os.system("rm -rf data_fun")
    return J,grad_adj
    
####################################################################################################
def fun_tchebycubeS(Lvector, **kwargs):
    vector=Lvector[0]
    J=fun_tchebycube(vector, **kwargs)
    return J
####################################################################################################    
def jac_tchebycubeS(Lvector, **kwargs):
    vector=Lvector[0]
    J,data_adj=jac_tchebycube(vector, **kwargs)
    return J,[data_adj]

####################################################################################################
def prod_f90(x1,x2):
    """[Scalar Product of vector with Fortran90]
    
    Args:
        x1 ([float vector]): [vector of ncomp components]
        x2 ([float vector]): [vector of ncomp components]
    Returns:
        [float]: [scalar product]
    """    

    os.mkdir("prod_scal")
    write_init_h5(x1[0],x1[1],x1[2],'prod_scal/x1.h5')
    write_init_h5(x2[0],x2[1],x2[2],'prod_scal/x2.h5')

    input_dir = options.input
    nprocs = options.np

    command_base=" "+options.exe_dir+"/prod_scal_space.x --input "+input_dir+" --a-file prod_scal/x1.h5 --b-file prod_scal/x2.h5 > prod_scal.txt"

    command="srun -np "+str(nprocs)+command_base                                                                
    os.system(command)

    f = open("prod_scal/res",'r')
    res = float(f.read())
    f.close()
    os.system("rm -rf prod_scal")

    return res
    

####################################################################################################
def prodn(x1,x2,ncomp):
    """[Scalar Product of vector with ncomp components]
    
    Args:
        x1 ([float vector]): [vector of ncomp components]
        x2 ([float vector]): [vector of ncomp components]
    Returns:
        [float]: [scalar product]
    """    
    x1 = np.asarray(x1)
    x2 = np.asarray(x2)
    n=x1.size/ncomp
    prod_scalar=np.vdot(x1,x2)
    prod_scalar/=n
    return prod_scalar
####################################################################################################
def prod3(x1,x2):
    """[Scalar Product of vector with 3 components]
    
    Args:
        x1 ([float vector]): [vector of 3 components]
        x2 ([float vector]): [vector of 3 components]
    Returns:
        [float]: [scalar product]
    """    
    x1 = np.asarray(x1)
    x2 = np.asarray(x2)
    n=x1.size/3
    prod_scalar=np.vdot(x1,x2)
    prod_scalar/=n
    return prod_scalar
##################################################################################################
def prod6(x1,x2):
    """[Scalar Product of vector with 6 components]
    
    Args:
        x1 ([float vector]): [vector of 6 components]
        x2 ([float vector]): [vector of 6 components]
    Returns:
        [float]: [scalar product]
    """
    x1 = np.asarray(x1)
    x2 = np.asarray(x2)
    n=x1.size/6
    prod_scalar=np.vdot(x1,x2)
    prod_scalar/=n
    return prod_scalar
####################################################################################################
def normalization(vector, Energ_value, prod):
    """[Normalize a raw vector into Energ_value using prod function]
    
    Args:
        vector ([float]): [vector]
        Energ_value ([float]): [value to normalize]
        prod ([function]): [scalar product]
    Returns:
        [float]: [vector normalized]
    """    
    vector = np.asarray(vector)
    E=prod(vector,vector)
    vector*=sqrt(Energ_value)/sqrt(E)
    return vector
####################################################################################################
def options_parser():
    # Manage the options :
    parser = argparse.ArgumentParser(description="Add help")  # On peut possiblement rajouter le usage
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose mode")
    group_tchebycube=parser.add_argument_group('options tchebycube')
    group_tchebycube.add_argument("--input_fun", help="File name of the parmaters needed to run Tchebycube")                      
    group_tchebycube.add_argument("--input_jac", help="File name of the parmaters needed to run Tchebycube adjoint")                      
    group_tchebycube.add_argument("--exe_dir", default='build', help="Directory of tchebycube.x, by default it is build")                      
    group_tchebycube.add_argument("--np", type=int, default=1,
                         help="Process number for each run (default=1)")

    group_opt = parser.add_argument_group('Optimization_rotation options')
    group_opt.add_argument("--output", help="Output Directory")                      
    group_opt.add_argument("--tol", type=float, default=1e-6,
                        help="Error tolerance on optimized loop")
    group_opt.add_argument("--iter", type=int, default=10,
                        help="Max iteration of the optimized loop")
    group_opt.add_argument("--grad", action="store_true",default=False,
                        help="Enable the mode option grad, by default is conjgrad")
    group_opt.add_argument("--superls", action="store_true",
                        help="active the option tri for linesearch")
    group_opt.add_argument("--E0", type=float, default=1.0,
                        help="Energy Budget (normalization constraint) E0=1.0 default")
    group_opt.add_argument("--armijo", action="store_true",default=False,
                        help="Use armijo linesearch")
    group_opt.add_argument("--conjgrad", action="store_true",
                        help="Enable the mode option conjgrad, by default is grad")

    group_sph = parser.add_argument_group('Optimization multi-sphere options')
    group_sph.add_argument("--multi", action="store_true",default=False,
                        help="enable multi sphere optimization")
    group_sph.add_argument("--alpha",type=float, default=1.0,
                        help="Value for constant alpha_k [default=1.0]")
    group_sph.add_argument("--LS",default='LS_wolfe',
                        help="linesearch: default(LS_wolfe); or write --LS LS_armijo")
    return parser.parse_args()
############################################################################################
if __name__ == '__main__':


    beginning_time = time.time()
    options = options_parser()
    if (options.verbose) : print(vars(options))
    dict_options=vars(options)


    prod=lambda x, y: lib_optimize.prodn(x, y, 4)

    #Initialisation à partir d'un random
    
    N = [129, 33, 65]
    init_ua = 2.*np.random.random(N)-np.ones(N)
    init_uz = 2.*np.random.random(N)-np.ones(N)
    init_ur = 2.*np.random.random(N)-np.ones(N)

    os.system("mkdir -p "+options.output)    

    try:
        file = h5py.File('Data/1Z_RE_10000/run1/cond_init/init_6.h5', 'r')
        init_ua = file['/dump/u1'][()]
        init_uz = file['/dump/u2'][()]
        init_ur = file['/dump/u3'][()]
        
        #transposition des données pour avoir un ordre colonne major
        #data_read = [np.transpose(grad_u1, axes=(2,1,0)),np.transpose(grad_u2, axes=(2,1,0)),np.transpose(grad_u3, axes=(2,1,0))]
        
        data_dns_start = [np.reshape(init_ua, (N[0],N[1],N[2])),np.reshape(init_uz, (N[0],N[1],N[2])),np.reshape(init_ur, (N[0],N[1],N[2]))]
    except: 
        print("The current directory is "+os.getcwd())
        print("%s is not inside this directory ?" % str_file) 
        sys.exit(1)

    
    data_dns_start=[init_ua, init_uz, init_ur]
    vector_begin=normalization(data_dns_start, options.E0, prod)
    
    prod_value_start=prod(vector_begin,vector_begin)
    if (options.verbose) : print('Starting value of the scalar product='+str(prod_value_start))
    if options.conjgrad : 
        ConjGrad=True
        if (options.verbose) : print('option conjgrad descent actif')
    else : 
        ConjGrad=False

    if (options.multi):
        data_opt =sgd.Optimise_On_Multi_Sphere([vector_begin], [options.E0], fun_tchebycubeS, jac_tchebycubeS, prod,kwargs_f=dict_options, \
                  err_tol=options.tol, alpha_k =options.alpha, LS = options.LS, CG = ConjGrad)

    else:
        data_opt = lib_optimize.optimize_rotation(fun_tchebycube, vector_begin, jac_tchebycube, prod, 1, **dict_options)

    print(data_opt)

    print('Successfull optimization of in  %.2f sec.' % (time.time() - beginning_time))

    #Ecriture du resultat au format h5

    write_init_h5(data_opt.X_opt[0],data_opt.X_opt[1],data_opt.X_opt[2],options.output+"/init_opt.h5")

#########################################################################################################


