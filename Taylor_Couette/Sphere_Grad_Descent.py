#! /usr/bin/python3
import os
import numpy as np
import copy
from warnings import warn
from scipy.optimize import linesearch as ls
from scipy.optimize._linesearch import scalar_search_armijo
from scipy.optimize._linesearch import scalar_search_wolfe2
import lib_util as util
class LineSearchWarning(RuntimeWarning):
    pass;

# Paul Mannix 2022
# We acknowledge the scipy.optimize.linesearch library for the
# modified
# LS_armijo_multiple, LS_wolfe_multiple
# 
# which call scalar_search_armijo, scalar_search_wolfe2
# from scipy.optimize.linesearch


###################################################################
class result():
    
    """
    class for result of optimize_rotation
    
    Inputs:
    components - integer the number of norm constraints 

    Returns:
    None
    """
    
    def __init__(self,components):
        
        self.N=components;
        self.X_opt=np.asarray([])
        self.Iterations=0
        self.counterloop_line_search=[]
        self.Function_Evals=0
        self.Gradient_Evals=0
        self.Residual=[]
        self.Residual_rescaled=[]
        self.Step_Size=[]
        self.Function_Value=[]     

    def __str__(self):

        s= ( '#Optimize multi sphere  \n'
             +'Total_iterations     = '+str(self.Iterations)    +'\n'
             +'sum_of_line_search = '+str(self.Function_Evals)+'\n'
             +'iteration_line_search = '+str(self.counterloop_line_search)+'\n'
             +'Gradient_evaluations = '+str(self.Gradient_Evals)+'\n'
             +'Residual = '+str(self.Residual[0])              +'\n'
             +'Residual_rescaled= '+str(self.Residual_rescaled[0])              +'\n'
             +'Alpha='+str(self.Step_Size)     +'\n'
             +'J='+str(self.Function_Value)+'\n'
             #+'  X_opt              = array('+str(self.X_opt)+')\n'
             );
        return s
##################################################################

#------------------------------------------------------------------------------
# Armijo line and scalar searches
#------------------------------------------------------------------------------

def LS_armijo_multiple(f, inner_prod, M_0, X_k, g_k, d_k,  old_fval, args_f=(), args_IP = (), kwargs_f = {}, kwargs_IP={}, alpha0=1.0, c1=1e-4):

    """
    Minimize over alpha, the function ``phi(α) = f( R_xk(α_k*d_k) )``,
    where X_k+1 = R_xk =  √M_0*(X_k + α_k*d_k)/|| X_k + α_k*d_k ||_2

    is the rectraction based update see:
    Nicolas Boumal. An introduction to optimization on smooth manifolds, 2020
    
    Parameters
    ----------
    f : callable f(x,*args_f,**kwargs_f)
        Objective function to be minimised
    inner_prod : callable IP(x,y,*args_IP,**kwargs_IP)
        Inner product
    M_0: list of floats
         Spherical manifold radius <X_0,X_0> = M_0     
    X_k : list of array_like
        Current point.
    g_k: list of array_like
        tangent gradient of current point    
    d_k : list of array_like
        **like** the Search direction.  
    old_fval : float
        Value of `f` at point `X_k`.
    args_f : tuple, optional
        Optional arguments to pass to J(X_k).
    args_IP : tuple, optional
        Optional arguments to pass to inner_prod <F,G>.
    c1 : float, optional
        Value to control stopping criterion.
    alpha0 : scalar, optional
        Value of `alpha` at start of the optimization.
    Returns
    -------
    alpha
    f_count
    f_val_at_alpha
    Notes
    -----
    Uses the interpolation algorithm (Armijo backtracking) as suggested by
    Wright and Nocedal in 'Numerical Optimization', 1999, pp. 56-57
    """

    X_k = np.atleast_1d(X_k)
    fc = [0]

    def phi(alpha1):
        
        fc[0] += 1;
        
        #apply norm constraints
        X_new = copy.deepcopy(X_k);
        for index,c_i in enumerate(M_0):
            X_new[index]  = Update_vector(X_k[index],alpha1,d_k[index],c_i,inner_prod,args_IP,kwargs_IP);
        return f( X_new, *args_f, **kwargs_f)

    if old_fval is None:
        phi0 = phi(0.)
    else:
        phi0 = old_fval  # compute f(xk) -- done in past loop

    # Compute the derivative w.r.t alpha at alpha=0        
    derphi0=0.    
    for index,_ in enumerate(M_0):
        derphi0 += inner_prod(g_k[index],d_k[index],*args_IP,**kwargs_IP);
    
    alpha, phi1 = scalar_search_armijo(phi, phi0, derphi0, c1=c1, alpha0=alpha0)
    return alpha, fc[0], phi1;


#------------------------------------------------------------------------------
# Strong Wolfe line and scalar searches require
# 0 < c1 < c2 < 0.5 when using Fletcher-Reeves
# see H.Sato & T.Iwai, A New globally convergent Riemannian CG method, 2015
#------------------------------------------------------------------------------

def LS_wolfe_multiple(f, myfprime, inner_prod, M_0, X_k, g_k, d_k, old_fval=None, old_old_fval=None, args_f=(), args_IP = (), kwargs_f = {}, kwargs_IP={}, c1=1e-4, c2=0.4, amax=None, extra_condition=None, maxiter=50):
    
    """
    Find alpha that satisfies strong Wolfe conditions by ....

    Minimizing over alpha, the function ``phi(α) = f( R_xk(α_k*d_k) )``,
    
    where X_k+1 = R_xk =  √M_0*(X_k + α_k*d_k)/|| X_k + α_k*d_k ||_2

    is the rectraction based update see:
    Nicolas Boumal. An introduction to optimization on smooth manifolds, 2020

    Parameters
    ----------
    f : callable f(x,*args_f,**kwargs_f)
        Objective function.
    myfprime : callable f'(x,*args_f,**kwargs_f)
        Objective function gradient.
    inner_prod : callable IP(x,y,*args_IP,**kwargs_IP)
        Inner product
    M_0 : list of float(s)
        Constraint radii       
    xk : list of ndarray(s)
        Starting point.
    pk : list of ndarray(s)
        Search direction.
    gfk : list of ndarray(s), optional
        Gradient value for x=xk (xk being the current parameter
        estimate). Will be recomputed if omitted.
    old_fval : float, optional
        Function value for x=xk. Will be recomputed if omitted.
    old_old_fval : float, optional
        Function value for the point preceding x=xk.
    args : tuple, optional
        Additional arguments passed to objective function.
    c1 : float, optional
        Parameter for Armijo condition rule.
    c2 : float, optional
        Parameter for curvature condition rule.
    amax : float, optional
        Maximum step size
    extra_condition : callable, optional
        A callable of the form ``extra_condition(alpha, x, f, g)``
        returning a boolean. Arguments are the proposed step ``alpha``
        and the corresponding ``x``, ``f`` and ``g`` values. The line search
        accepts the value of ``alpha`` only if this
        callable returns ``True``. If the callable returns ``False``
        for the step length, the algorithm will continue with
        new iterates. The callable is only called for iterates
        satisfying the strong Wolfe conditions.
    maxiter : int, optional
        Maximum number of iterations to perform.
    Returns
    -------
    alpha : float or None
        Alpha for which ``x_new = x0 + alpha * pk``,
        or None if the line search algorithm did not converge.
    fc : int
        Number of function evaluations made.
    gc : int
        Number of gradient evaluations made.
    new_fval : float or None
        New function value ``f(x_new)=f(x0+alpha*pk)``,
        or None if the line search algorithm did not converge.
    old_fval : float
        Old function value ``f(x0)``.
    new_slope : float or None
        The local slope along the search direction at the
        new value ``<myfprime(x_new), pk>``,
        or None if the line search algorithm did not converge.
    Notes
    -----
    Uses the line search algorithm to enforce strong Wolfe
    conditions. See Wright and Nocedal, 'Numerical Optimization',
    1999, pp. 59-61.
    """

    fc = [0]
    gc = [0]
    gval = [None]
    gval_alpha = [None]

    def phi(alpha1):
        fc[0] += 1;

        #apply norm constraints
        X_new = copy.deepcopy(X_k);
        for index,c_i in enumerate(M_0):
            X_new[index]  = Update_vector(X_k[index],alpha1,d_k[index],c_i,inner_prod,args_IP,kwargs_IP);
        
        return f( X_new, *args_f,**kwargs_f)    


    fprime = myfprime
        
    def derphi(alpha1):
        gc[0] += 1

        X_kp1 = copy.deepcopy(X_k);
        g_kp1 = copy.deepcopy(g_k);
        Tdkm1 = copy.deepcopy(g_k);
        derphi1=0.
        
        #apply norm constraints
        for index,c_i in enumerate(M_0):
            X_kp1[index]  = Update_vector(X_k[index],alpha1,d_k[index],c_i,inner_prod,args_IP,kwargs_IP);
        
        # Calculate the Euclidean gradient
        Nab_Jkp1 = fprime(X_kp1,*args_f,**kwargs_f)
        
        # Compute the tangent gradient and perform vector transport of d_k 
        g_kp1[index] = tangent_vector(X_kp1[index],Nab_Jkp1[index],inner_prod,args_IP,kwargs_IP)
        Tdkm1[index] = transport_vector(X_kp1[index],d_k[index],inner_prod,args_IP,kwargs_IP)
        derphi1     += inner_prod(g_kp1[index],Tdkm1[index],*args_IP,**kwargs_IP); # Compute the derivate w.r.t alpha1
        
        # Store current tangent gradient for later use
        gval[0]  = g_kp1;
        gval_alpha[0] = alpha1;

        print("size derphi",np.size(derphi1))
        print("derphi",derphi1)
        
        return float(np.squeeze(derphi1));    

    # Compute the derivative w.r.t alpha at alpha=0    
    derphi0=0.    
    for index,_ in enumerate(M_0):    
        derphi0 += inner_prod(g_k[index],d_k[index],*args_IP,**kwargs_IP);    


    extra_condition2 = None

    alpha_star, phi_star, old_fval, derphi_star = scalar_search_wolfe2(
            phi, derphi, old_fval, old_old_fval, derphi0, c1, c2, amax,
            extra_condition2, maxiter=maxiter)

    if derphi_star is None:
        
        warn('The line search algorithm did not converge', LineSearchWarning)
    else:
        # derphi_star is a number (derphi) -- so use the most recently
        # calculated gradient used in computing it derphi = gfk*pk
        # this is the gradient at the next step no need to compute it
        # again in the outer loop.
        derphi_star = gval[0]

    return alpha_star, fc[0], gc[0], phi_star, old_fval, derphi_star

#------------------------------------------------------------------------------
# Pure-Python Spherical Manifold Optimisation
#------------------------------------------------------------------------------

# These three functions are specific to the update formula/Retraction

# X_k+1 = √M_0*(X_k + α_k*d_k)/|| X_k + α_k*d_k ||_2

# args, kwargs should be based as args=(a,b,c), kwargs={'a':1,'b':2,'c':3} and only unpacked when calling inner_prod(f,g,*args,**kwargs)

def transport_vector(X_k,dkm1,inner_prod,args_IP=(),kwargs_IP={}):
	"""
	Return the vector transport for an arbitrary inner product

	Inputs:
	X_k    	   - parameter vector
	dkm1 	   - previous search direction
	inner_prod - callable function: takes args_IP = (),kwargs_IP = {}

	Returns:
	T(dkm1)    - vector transported to X_k tangent plane

	"""
	L2   = np.sqrt( inner_prod(X_k,X_k ,*args_IP,**kwargs_IP) );
	return dkm1 - ( inner_prod(X_k,dkm1,*args_IP,**kwargs_IP)/(L2**2) )*X_k; 

def tangent_vector(X_k,Nab_Jk,inner_prod,args_IP=(),kwargs_IP={}):
    """
	Return the tangent vector for an arbitrary inner product

	Inputs:
	X_k    	   - parameter vector
	Nab_Jk 	   - Euclidean vector
	inner_prod - callable function: takes args_IP = (),kwargs_IP = {}

	Returns:
	gk - tangent vector
	"""
    coeff = ( inner_prod(X_k,Nab_Jk,*args_IP,**kwargs_IP)/inner_prod(X_k,X_k,*args_IP,**kwargs_IP))
    
    return Nab_Jk -coeff*X_k

def Update_vector(X_k,alpha_k,d_k,M_0,inner_prod,args_IP=(),kwargs_IP={}):

    """
    return
    Update the parameter vector in the search direction alpha_k*d_k

    Inputs:
    X_k    - vector parameter 
    alpha_k- float  step-size
    d_k    - vector search direction
    M_0    - float spherical manifold size < X_0,X_0 > = M_0
    inner_prod - callable function: takes args_IP = (),kwargs_IP = {}

    Returns:

    X_k+1  - vector new parameter vector

    Notes this is the rectraction based update see:
    Nicolas Boumal. An introduction to optimization on smooth manifolds, 2020

    """

    #Xn = np.sqrt( M_0 );
    #dn = np.sqrt( inner_prod(d_k,d_k,*args) );
    #return np.cos(alpha_k*dn)*X_k + np.sin(alpha_k*dn)*d_k*(Xn/dn); 

    f    = X_k + alpha_k*d_k;
    L2_f = inner_prod(f,f,*args_IP,**kwargs_IP);

    return f*np.sqrt(M_0/L2_f)

def Optimise_On_Multi_Sphere(X_0, M_0, f, myfprime, inner_prod, args_f = (), args_IP=(), kwargs_f = {}, kwargs_IP={}, err_tol = 1e-06, max_iters = 1000, alpha_k = 1., LS = 'LS_wolfe', CG = False, callback=None, verbose=True):
	
    """
    Function to perform the minimisation of J(X) via 
    gradient based descent  Grad_f(X) on a spherical 
    manifold <X,X> = M_0.

    Inputs:
    X_0        - list of initial parameter vector guess i.e. [x_0,x_1, .... , x_N]
    M_0        - list of spherical manifold radius of each vector i.e. [c_0,c_1, .... , c_N]
    f 	   	   - callable returns      J(X_k) takes unpacked *args_f, **kwargs_f
    Grad_f     - callable returns Grad_J(X_k) takes unpacked *args_f, **kwargs_f
    inner_prod - callable returns <F,G>		  takes unpacked *args_IP,**kwargs_IP

    Returns:

    RESIDUAL - vector of error at each iterations
    FUNCT    - vector of function of evaluations
    X_opt    - list of optimal vectors

    """
    iter=util.OptionValueD('iter',0,warn=False,**kwargs_f)

    if LS == 'LS_wolfe': 
        LS = LS_wolfe_multiple;
    elif LS == 'LS_armijo':
        LS = LS_armijo_multiple;

    error = np.ones(len(M_0)); 
    error_rescaled= np.ones(len(M_0)); 
    
    func_evals=0;
    grad_evals=0;
    alpha_max = alpha_k;

    RESIDUAL = [];
    RESIDUAL_rescaled = [];
    for val in error:
        RESIDUAL.append([]);
    for val in error_rescaled:
        RESIDUAL_rescaled.append([]);

    # Initialise the class for data handling
    R = result(len(M_0))
    file_exist=os.path.exists("optimize.txt")
    f2= open("optimize.txt", "a")
    if (not file_exist) :
        f2.write('% i  J  Res  Res_rescaled  alpha  nbr_ls \n')
    
    # Normalise X_k so that <X,X> = M_0
    J_k_old = None;
    X_k = [ x_i*np.sqrt( c_i/inner_prod(x_i,x_i,*args_IP,**kwargs_IP) ) for x_i,c_i in zip(X_0,M_0) ];
    J_k = f(X_k,*args_f,**kwargs_f); func_evals+=1;
    
    Iterations=0
    while (max(error) > err_tol) and (Iterations < max_iters):
        Iterations+=1
        
        # Reuse the gradient computed if using a strong wolfe line-search
        if (LS == LS_wolfe_multiple) and (Iterations > 1):
            g_k = derphi_star;
        else:
            Nab_Jk = myfprime(X_k,*args_f,**kwargs_f);
            g_k    = [tangent_vector(u,du,inner_prod,args_IP,kwargs_IP) for u,du in zip(X_k,Nab_Jk)];
            grad_evals +=1;

        
        # Select a search direction d_k via 
        # SD-steepest descent or CG - conjugate gradient
        if (Iterations > 1) and (CG==True):
            # Conjuagte-gradient
            beta_k_FR = 0.; 
            beta_k_PR = 0.;
            Tg_km1    = copy.deepcopy(g_k);
            Td_km1    = copy.deepcopy(g_k);

            for ii,_ in enumerate(g_k):
                
                beta_k_FR += inner_prod(g_k[ii],g_k[ii],*args_IP,**kwargs_IP)/inner_prod(g_km1[ii],g_km1[ii],*args_IP,**kwargs_IP)
                
                Tg_km1[ii] = transport_vector(X_k[ii],g_km1[ii],inner_prod,args_IP,kwargs_IP);
                beta_k_PR += ( inner_prod(g_k[ii],g_k[ii],*args_IP,**kwargs_IP) - inner_prod(g_k[ii],Tg_km1[ii],*args_IP,**kwargs_IP) )/inner_prod(g_km1[ii],g_km1[ii],*args_IP,**kwargs_IP);

                Td_km1[ii] = transport_vector(X_k[ii],d_k[ii],inner_prod,args_IP,kwargs_IP)

            # Use the Fletcher-Reeves + Polak Rib\`ere update 
            # of H. Sato Riemannian conjugate gradient methods 2021
            # to select the parameter ß_k
            beta_k = max(0.,min(beta_k_FR,beta_k_PR));    
            d_k = [-1.*g_ki + beta_k*T_ki for g_ki,T_ki in zip(g_k,Td_km1)];
        
        else:   
            # Gradient descent    
            d_k = [-1.*g_i for g_i in g_k];


        # Perform a line-search for the step-size α_k to ensure descent 
        if (R.Iterations == 1) or (LS == LS_armijo_multiple):
            alpha_k, f_evals, J_k = LS_armijo_multiple(f, inner_prod, M_0, X_k, g_k, d_k,  J_k,  args_f, args_IP,kwargs_f, kwargs_IP, alpha0 = alpha_k)
            func_evals+=f_evals;
        else:
            alpha_k, f_evals, g_evals, J_k, J_k_old, derphi_star = LS(f, myfprime, inner_prod, M_0, X_k, g_k, d_k, J_k,J_k_old, args_f, args_IP, kwargs_f, kwargs_IP,amax=alpha_max);
            grad_evals+=g_evals;
            func_evals+=f_evals;

        # Update the parameter vector - applying the norm constraints    
        for index,c_i in enumerate(M_0):
            
            if alpha_k == None:
                print("\n Couldn't find a descent direction .... Terminating \n");      
                return R;
            else:
                X_k[index]   = Update_vector(X_k[index],alpha_k,d_k[index],c_i,inner_prod,args_IP,kwargs_IP);
                error[index] = np.sqrt(inner_prod(g_k[index],g_k[index],*args_IP,**kwargs_IP));
                error_rescaled[index]=error[index]/np.sqrt(inner_prod(Nab_Jk[index],Nab_Jk[index],*args_IP,**kwargs_IP))


        # Update the optimization state        
        R.X_opt = X_k;
        R.Iterations=Iterations+iter;
        R.Function_Evals+=func_evals; 
        R.Gradient_Evals+=grad_evals; 
        
        for ii,_ in enumerate(error):
            RESIDUAL[ii].append(error[ii]);
            RESIDUAL_rescaled[ii].append(error_rescaled[ii]);

        R.Residual=RESIDUAL;
        R.Residual_rescaled=RESIDUAL_rescaled;
        R.Step_Size.append(alpha_k)
        ## R.Function_Value.append(-1.*J_k)
        R.Function_Value.append(J_k)
        
        R.counterloop_line_search.append(func_evals)
        
        f2.write(str(R.Iterations)+'  '+str(J_k)+'  '+str(error[0])+'  '+str(error_rescaled[0])+'  '+str(alpha_k)+'  '+str(func_evals)+'\n')
        f2.flush()
        
        g_km1 = copy.deepcopy(g_k);
        func_evals=0;
        grad_evals=0;

        # Save the optimisation state   
        if callback != None:
            callback(R.Iterations);
             
        # Print out the optimisation status    
        if verbose : print(R,flush=True);   
        f_txt = open("optimize_result.txt", "w")
        f_txt.write('#'+str(kwargs_f)+'\n')      
        f_txt.write(str(R))
        f_txt.flush()
        f_txt.close()
    f2.close()
    return R;
