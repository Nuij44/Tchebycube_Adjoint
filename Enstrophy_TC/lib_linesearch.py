#!/usr/bin/env python3
#
# lib_linesearch.py
# Florence Marcotte & Yannick Ponty
# version 0.4   avril 2021 (pass willis test, stratification test)
##########################################################################
import numpy as np
from math import sqrt, cos, sin
import lib_util

########################################################################## 
##########################################################################      
def line_search(fun, x, d ,direction, innerprod, super, mu0=1, **kwargs):
    """[line search function]
    
    Parameters:
    
    Args
        | fun : ([function]) [get X0 and give back J]. 
        | x : ([float vector]) [input raw vector]. 
        | d : ([float vector]) [step vector]. 
        | direction : ([int]) [1 for seek the maximun, -1 to seek the minimum]. 
        | innerprod : ([function]): [take x1, x2 give dot product]. 
        | super : ([bool]) [active super linesearch options]. 
        | mu0 : ([float])  [initial step size]. 
        
    Returns:
        [float, int] : [mu , number of line search].
    """
    mu=mu0
    J=fun(x,**kwargs)
    xnew=update_vector_rotation(x,mu,d,innerprod, **kwargs)       
    Jnew=fun(xnew,**kwargs)
    count_LS=1
    while ( direction*(J-Jnew)> 0) :
        J_mem1=Jnew
        mu_mem1=mu
        count_LS= count_LS+1
        mu=0.5*mu
        xnew=update_vector_rotation(x,mu,d, innerprod, **kwargs)
        Jnew=fun(xnew,**kwargs)
    
    # super last linesearch  
    if (count_LS > 1 and super):
        numer=J*(mu+mu_mem1)/(mu*mu_mem1) + J_mem1*mu/(mu_mem1*(mu_mem1-mu)) + Jnew*mu_mem1/(mu*(mu-mu_mem1))    
        den=2*J/(mu*mu_mem1) + 2*J_mem1/(mu_mem1*(mu_mem1-mu)) + 2*Jnew/(mu*(mu-mu_mem1))
        mutest=numer/den
        x_test=update_vector_rotation(x,mutest,d, innerprod, **kwargs)
        J_test=fun(x_test,**kwargs)
        count_LS= count_LS+1
        if (direction*(J_test-Jnew) > 0):
            mu=mutest
        
    return mu,count_LS
#####################################################################################################
def update_vector_rotation(x,mu,d,innerprod, **kwargs):
    """[update vector with rotation method]
    
    Parameters: 
    
    kwargs 
        E0 ([float]): [Energy budget (normalization constraint)]
    
    Args
        | x ([float vector]): [input raw vector].
        | mu ([type]): [step update].
        | d ([float vector]): [step vector].
        | innerprod ([function]): [take x1, x2 give dot product].

    returns
        [float vector]: [new float vector update].
    """
    E0=lib_util.OptionValue('E0',**kwargs)
    x=np.asarray(x)
    d=np.asarray(d)
    dn=d/sqrt(innerprod(d,d))*sqrt(E0)
    Xnew=x*cos(mu)+dn*sin(mu) # Rotation method Douglas 98
    return Xnew
##################################################################
#------------------------------------------------------------------------------
# Armijo line and scalar searches from  scipy/scipy/optimize/linesearch.py Commits on Nov 12, 2020
#
#  few adaptations has been made for lib_optimize.py 
#  including the rotation advance update for the vector
#------------------------------------------------------------------------------

def line_search_armijo(f, xk, pk, gfk, old_fval, direction, innerprod,  c1=1e-4, alpha0=1, **kwargs):
    """Minimize over alpha, the function ``f(xk)``.
    
    Parameters:
    
    Args
        | f  [callable] Function to be minimized.
        | xk [array_like] Current point.
        | pk [array_like] Search direction.
        | gfk [array_like] Gradient of `f` at point `xk`.
        | old_fval [float] Value of `f` at point `xk`.
        | direction [int] 
        | innerprod [callable] function of the dot product
        | c1 [float] optional Value to control stopping criterion.
        | alpha0 [scalar] optional Value of `alpha` at start of the optimization.
        
    Returns
        | alpha 
        | f_count 
        | f_val_at_alpha :
    
    Notes
    -----
    Uses the interpolation algorithm (Armijo backtracking) as suggested by
    Wright and Nocedal in 'Numerical Optimization', 1999, pp. 56-57
    """
    xk = np.atleast_1d(xk)
    fc = [0]

    def phi(alpha1, xk, pk, innerprod, **kwargs):    
        fc[0] += 1
        xnew=update_vector_rotation(xk, alpha1, pk, innerprod, **kwargs)
        return -direction*f(xnew,**kwargs)

    if old_fval is None:
        phi0 = phi(0.0, xk, pk, innerprod, **kwargs)
    else:
        phi0 = -direction*old_fval  # compute f(xk) -- done in past loop

    derphi0 = innerprod(gfk, pk)
    alpha, phi1 = scalar_search_armijo(phi, phi0, derphi0, xk, pk, innerprod, c1=c1, alpha0=alpha0, **kwargs)
    if alpha==None : 
        print('armijo failed alpha='+str(alpha)+' with loop number='+str(fc))
        alpha=1.0e-4
        print('armijo has return artificial alpha given at the value '+str(alpha))
    return alpha, fc[0]

def scalar_search_armijo(phi, phi0, derphi0,  xk, pk, innerprod, c1=1e-4, alpha0=1, amin=0, **kwargs):
    """Minimize over alpha, the function ``phi(alpha)``.
    
    Uses the interpolation algorithm (Armijo backtracking) as suggested by
    Wright and Nocedal in 'Numerical Optimization', 1999, pp. 56-57
    alpha > 0 is assumed to be a descent direction.
    Paramaters:
    
    Args
        | phi
        | phi0
        | derphi0
        | xk [array_like] Current point.
        | pk [array_like] Search direction.
        | innerprod [callable] function of the dot product
        | c1 [float] optional Value to control stopping criterion.
        | alpha0 [scalar] optional Value of `alpha` at start of the optimization.
        | amin [int] optional 
    
    Returns
        | alpha
        | phi1
    """
    phi_a0 = phi(alpha0, xk, pk, innerprod, **kwargs)
    if phi_a0 <= phi0 + c1*alpha0*derphi0:
        return alpha0, phi_a0

    # Otherwise, compute the minimizer of a quadratic interpolant:

    alpha1 = -(derphi0) * alpha0**2 / 2.0 / (phi_a0 - phi0 - derphi0 * alpha0)
    phi_a1 = phi(alpha1, xk, pk, innerprod, **kwargs)

    if (phi_a1 <= phi0 + c1*alpha1*derphi0):
        return alpha1, phi_a1
    
    # Otherwise, loop with cubic interpolation until we find an alpha which
    # satisfies the first Wolfe condition (since we are backtracking, we will
    # assume that the value of alpha is not too small and satisfies the second
    # condition.

    while alpha1 > amin:       # we are assuming alpha>0 is a descent direction
        factor = alpha0**2 * alpha1**2 * (alpha1-alpha0)
        a = alpha0**2 * (phi_a1 - phi0 - derphi0*alpha1) - \
            alpha1**2 * (phi_a0 - phi0 - derphi0*alpha0)
        a = a / factor
        b = -alpha0**3 * (phi_a1 - phi0 - derphi0*alpha1) + \
            alpha1**3 * (phi_a0 - phi0 - derphi0*alpha0)
        b = b / factor

        alpha2 = (-b + np.sqrt(abs(b**2 - 3 * a * derphi0))) / (3.0*a)
        phi_a2 = phi(alpha2, xk, pk, innerprod, **kwargs)

        if (phi_a2 <= phi0 + c1*alpha2*derphi0):    
            return alpha2, phi_a2

        if (alpha1 - alpha2) > alpha1 / 2.0 or (1 - alpha2/alpha1) < 0.96:
            alpha2 = alpha1 / 2.0

        alpha0 = alpha1
        alpha1 = alpha2
        phi_a0 = phi_a1
        phi_a1 = phi_a2

    # Failed to find a suitable step length
    return None, phi_a1
