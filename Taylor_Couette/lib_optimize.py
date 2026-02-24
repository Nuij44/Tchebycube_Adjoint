#!/usr/bin/env python3
#
# lib_optimize.py
#
# Florence Marcotte & Yannick Ponty
# version 0.4   avril 2021 (pass willis test, stratification test)
###################################################################
###########   Library needed  #####################################
import numpy as np
import sys
import os
from math import sqrt, cos, sin
import numpy.matlib as mb
import argparse
import copy
import h5py
import matplotlib.pyplot as plt
import lib_linesearch
import lib_util as util
from functools import reduce
###################################################################
class result():
    """[class for result of optimize_rotation]
    """
    def __init__(self):
        self.X_opt=np.asarray([])
        self.Iterations=0
        self.counterloop_line_search=[]
        self.sum_line_search=0
        self.residuLoop=[]
        self.Residual=[]
        self.Alpha=[]
        self.J=[]
        
    def __str__(self):
        s= ( '# Optimize_rotation \n'
             +'Total_iterations     = '+str(self.Iterations)    +'\n'
             +'sum_of_line_search = '+str(self.sum_line_search)+'\n'
             +'iteration_line_search = '+str(self.counterloop_line_search)+'\n'
             +'Residual = '+str(self.residuLoop)+'\n'
             +'Residual_rescaled = '+str(self.Residual)+'\n'
             +'Alpha ='+str(self.Alpha)+'\n'
             +'J='+str(self.J)+'\n'
             #+'# x= array('+str(self.X_opt)+')\n'
             )
        return s
    
###################################################################
##  
def optimize_rotation(fun, X0, jac, innerprod, direction, DnsAdj=True,**kwargs):
    """[using rotation optimization method Douglas 98
        as options:
        tol : tolerance of residu 2 (1.0e-6)
        conjgrad : choose conjgrad computation instead of grad default
        verbose : print result as each time step
        superls : activate option for the basic linesearch superls         
    Args:
        fun ([function]): [get X0 and give back J]
        X0 ([float vector]): [input raw vector]
        jac ([function]): [take X0 give Adjoint gradient]
        innerprod ([function]): [take x1, x2 give dot product]
        direction ([int]): [1 for seek the maximun, -1 to seek the minimum]
    Returns:
        [class result]: [object result design for fancy output]
    """
    tol=util.OptionValue('tol',**kwargs)
    conjgrad=util.OptionValue('conjgrad',**kwargs)
    superls=util.OptionValue('superls',**kwargs)
    verbose=util.OptionValue('verbose',**kwargs)
    armijo=util.OptionValue('armijo',**kwargs)
    iter=util.OptionValueD('iter',0,warn=False,**kwargs)
    kwargs2=copy.deepcopy(kwargs)
    kwargs2['DnsAdj']= DnsAdj
    DnsAdj = True
    
    x0=copy.deepcopy(X0)
    x0=np.asarray(x0)
    R=result()
    
    file_exist=os.path.exists("opt_matlab.txt")
    f2= open("opt_matlab.txt", "a")
    if (not file_exist) :
        f2.write('% i  J  Res  Res_rescaled  alpha  nbr_ls \n')
    
    if (DnsAdj) :
        J,w = jac(x0,**kwargs2)
    else:
        w = jac(x0,**kwargs2)   # in case J can not be computed in jac 
        J = fun(x0,**kwargs)
    J_history=[]
    Count_ls=[]
    Mu_ls=[]
    J_history.append(J)  
    residuLoop=[]
    Residual=[]

    w0=w
    wp0 = w - innerprod(w,x0)/innerprod(x0,x0)*x0 # projected gradient 
    d0=direction*wp0
    mu=0.01
    x=lib_linesearch.update_vector_rotation(x0,mu,d0,innerprod, **kwargs)
    iota = 0
    res = 1.0

    os.system("mkdir -p cond_init")    
    
    while ( res > tol ):
        iota = iota+1
        write_init_h5(x[0],x[1],x[2],'cond_init/init_'+str(iota)+'.h5')
        if (DnsAdj) :
            J,w = jac(x,**kwargs2)
        else:
            w = jac(x,**kwargs2)
            J = fun(x,**kwargs)
        J_history.append(J) # 
        wp = w - innerprod(w,x)/innerprod(x,x)*x
        
        if (conjgrad) : 
            b = innerprod(wp,wp-wp0)/innerprod(wp0,wp0) #Polack-Ribiere
            if (b < 0.0 or b> 1.0): b=0.0 
            d0 = d0 - innerprod(d0,x)/innerprod(x,x)*x
        else :
            b = 0    
            
        d = direction*wp  + b*d0    
        d0=d
        if (iota> 16 and superls): mu0=0.1
        else : mu0=16 #### mu0=1
        
        if (armijo):
            mu, count_ls=lib_linesearch.line_search_armijo(fun, x, d, w, J, direction, innerprod, **kwargs)
        else:
            mu, count_ls=lib_linesearch.line_search(fun, x, d, direction, innerprod, superls, mu0, **kwargs)
            
        Mu_ls.append(mu)
        Count_ls.append(count_ls)
        wp0=wp
        x=lib_linesearch.update_vector_rotation(x,mu,d,innerprod, **kwargs)
        res0=sqrt(innerprod(wp,wp))
        res=sqrt(innerprod(wp,wp)/innerprod(w,w))
        residuLoop.append(res0)
        Residual.append(res)
        
        iteration=iota+iter
        f2.write(str(iteration)+'  '+str(J)+'  '+str(res0)+'  '+str(res)+'  '+str(mu)+'  '+str(count_ls)+'\n')
        f2.flush()
        
        R.X_opt=x
        R.J=J_history
        R.Iterations=iota+iter
        R.counterloop_line_search=Count_ls
        R.sum_line_search=sum(Count_ls)
        R.Alpha=Mu_ls
        R.residuLoop=residuLoop
        R.Residual=Residual

        if verbose : print(R,flush=True)
        f = open("opt_python.txt", "w")
        f.write('#'+str(kwargs)+'\n')
        f.write(str(R))
        f.flush()
        f.close()  
    f2.close()  
    return R   

####################################################################################################
def prodv(x1,x2):
    x1 = np.asarray(x1)
    x2 = np.asarray(x2)
    n=x1.size
    prod_scalar=np.vdot(x1,x2)
    prod_scalar/=n
    return prod_scalar
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
##################################################################
####  test function for the main
##################################################################
#################################################################
dim=1000
xmin=-0.5
xmax=0.5
ymin=-0.5
ymax=0.5
#################################################################
##################################################################
def fun_test(Z,**kwargs ):
    x = np.linspace(xmin, xmax, dim)
    y = np.linspace(ymin, ymax, dim)
    X, Y = np.meshgrid(x, y)
    #y=mb.repmat(mb.linspace(ymin,ymax,dim),dim,1)
    Ye=np.ravel(Y)
    result=prodv(Z,Ye)
    return result
##################################################################
def jac_test(Z, **kwargs ):
    DnsAdj=util.OptionValueD('DnsAdj',False,**kwargs)
    x = np.linspace(xmin, xmax, dim)
    y = np.linspace(ymin, ymax, dim)
    X, Y = np.meshgrid(x, y)
    Ye=np.ravel(Y)
    if (DnsAdj): #  option not available in this test 
        J=fun_test(Z,**kwargs)
        return J,Ye
    else:
        return Ye
##################################################################
def init(prod):
    norm=1
    np.random.seed(options.seed)
    Xold=np.random.rand(dim,dim)
    Xold=Xold-np.mean(np.mean(Xold))
    prod=lambda x, y: prodn(x, y, 2)
    E=prod(Xold,Xold)
    X=Xold/sqrt(E)*norm
    return  X
##################################################################
def one_to_two(x):
    X=copy.deepcopy(x)
    X2d=np.reshape(X,(dim,dim))
    return X2d
##################################################################
def sol(x):
    X=copy.deepcopy(x)
    X[::-1].sort()
    J=fun_test(X)
    X2d=np.reshape(X,((dim,dim)) )
    return J,X2d
##################################################################
def normdiff(x1,x2):
    X1=copy.deepcopy(x1)
    X2=copy.deepcopy(x2)
    X1=np.ravel(X1)
    X2=np.ravel(X2)
    diff=X1-X2
    diff=np.absolute(diff)
    max_diff=max(diff)
    return max_diff
#########################################################################
def plot(Z):
    x = np.linspace(xmin, xmax, dim)
    y = np.linspace(ymin, ymax, dim)
    X, Y = np.meshgrid(x, y)
    plt.figure()
    plt.imshow(Z, interpolation="bicubic", 
       origin="lower", extent=[xmin,xmax,ymin,ymax])
    plt.colorbar()
    plt.show()
#####################################################################################
def options_parser():
    parser = argparse.ArgumentParser(description="Add help")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose mode")
    parser.add_argument("--tol", type=float, default=1e-04,
                        help="Error tolerance on optimized loop")
    parser.add_argument("--conjgrad", action="store_true",
                        help="Enable the mode option conjgrad, by default is grad")
    parser.add_argument("--superls", action="store_true",
                        help="active the option tri for linesearch")
    parser.add_argument("--armijo", action="store_true",
                        help="Use armijo linesearch")                   
    parser.add_argument("--seed", type=int, default=42,
                        help="seed of the random 2D array")
    parser.add_argument("-p", "--plot", action="store_true",
                        help="Enable plot 2D array")
    parser.add_argument("--E0", type=float, default=1.0,
                        help="Energy Budget (normalization constraint)")
    return parser.parse_args()
#########################################################################################
#########################################################################################
if __name__ == '__main__':
    options = options_parser()
    if (options.verbose) : print(options)
    dict_options=vars(options)
    
    prod=lambda x, y: prodn(x, y, 2)
    
    Xo=init(prod)
    if (options.plot) : plot(Xo)
    xo=np.ravel(Xo)
    
    J_Xo=fun_test(xo)
    J_sol,Xo_sol=sol(xo)
    print('J_Xo='+str(J_Xo)+' J_sol='+str(J_sol))
    if (options.plot) : plot(Xo_sol)
    
    resu_rot = optimize_rotation(fun_test, xo, jac_test, prod, -1, DnsAdj=False, **dict_options)
    print(resu_rot)
    Xrot=one_to_two(resu_rot.X_opt)
    if (options.plot) : plot(Xrot)
    max_diff=normdiff(Xo_sol,Xrot)
    print('abs(Jsol - Jopt) = '+str(abs(J_sol-resu_rot.J)))
    print('max of abs(Xsol[i] - Xopt[i]) ='+str(max_diff))
    
    if (max_diff < 8.0e-3):
        print('Test lib_optimize.py succeed !')
        sys.exit(0)
    else:
        print('Test lib_optimize.py failed !')
        sys.exit(1)
########################################################################################        
