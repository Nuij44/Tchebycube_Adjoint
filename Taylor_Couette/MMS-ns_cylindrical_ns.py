#!/usr/bin/python
from sympy import *
from sympy.utilities.codegen import codegen


def compute_linear_momentum(ur,ua,uz) :


    d_ur_dr = diff(ur, r)
    d_ua_dtheta = diff(ua, theta)
    d_ua_dr = diff(ua, r)
    d_ur_dtheta = diff(ur, theta)
    d_uz_dr = diff(uz, r)

    
    
    d2_ur_dr2 = diff(ur, r, 2)
    d2_ur_dtheta2 = diff(ur, theta, 2)
    d2_ur_dz2 = diff(ur, z, 2)

    d2_ua_dr2 = diff(ua, r, 2)
    
    d2_ua_dtheta2 = diff(ua, theta, 2)
    d2_ua_dz2 = diff(ua, z, 2)
    
    
    d2_uz_dr2 = diff(uz, r, 2)
    d2_uz_dtheta2 = diff(uz, theta, 2)
    d2_uz_dz2 = diff(uz, z, 2)

#    Lap_ur = d2_ur_dr2 + (1/r) * d_ur_dr - (ur/r**2) + (1/r**2) * d2_ur_dtheta2 + d2_ur_dz2 - (2/r**2) * d_ua_dtheta
#    Lap_ua = d2_ua_dr2 + (1/r) * d_ua_dr - (ua/r**2) + (1/r**2) * d2_ua_dtheta2 + d2_ua_dz2 + (2/r**2) * d_ur_dtheta
#    Lap_uz = d2_uz_dr2 + (1/r) * d_uz_dr + (1/r**2) * d2_uz_dtheta2 + d2_uz_dz2

    Lap_ua = (1./r)*d_ua_dr + d2_ua_dr2 + (1./r**2)*d2_ua_dtheta2 + d2_ua_dz2 - (ua/r**2) + (2./r**2)*d_ur_dtheta
    Lap_uz = (1./r)*d_uz_dr + d2_uz_dr2 + (1./r**2)*d2_uz_dtheta2 + d2_uz_dz2
    Lap_ur = (1./r)*d_ur_dr + d2_ur_dr2 + (1./r**2)*d2_ur_dtheta2 + d2_ur_dz2 - (ur/r**2) - (2./r**2)*d_ua_dtheta



    return [Lap_ur,Lap_ua,Lap_uz]

def grad_cyl(ua):

    dr_ua = diff(ua,r)
    dtheta_ua = diff(ua,theta)
    dz_ua = diff(ua,z)

    return [dtheta_ua/r,dz_ua,dr_ua]


def compute_linear_momentum_solver(ur,ua,uz) :
    
    d2_ur_dr2 = diff(ur, r, 2)
    d_ur_dr = diff(ur, r)
    d2_ur_dtheta2 = diff(ur, theta, 2)
    d2_ur_dz2 = diff(ur, z, 2)
    d_ua_dtheta = diff(ua, theta)

    d2_ua_dr2 = diff(ua, r, 2)
    d_ua_dr = diff(ua, r)
    d2_ua_dtheta2 = diff(ua, theta, 2)
    d2_ua_dz2 = diff(ua, z, 2)
    d_ur_dtheta = diff(ur, theta)
    
    d2_uz_dr2 = diff(uz, r, 2)
    d_uz_dr = diff(uz, r)
    d2_uz_dtheta2 = diff(uz, theta, 2)
    d2_uz_dz2 = diff(uz, z, 2)

    Lap_ur = d2_ur_dr2 + (1/r) * d_ur_dr + (1/r**2) * d2_ur_dtheta2 + d2_ur_dz2 - (ur/r**2) - (2/r**2) * d_ua_dtheta
    Lap_ua = d2_ua_dr2 + (1/r) * d_ua_dr + (1/r**2) * d2_ua_dtheta2 + d2_ua_dz2 - (ua/r**2) + (2/r**2) * d_ur_dtheta
    Lap_uz = d2_uz_dr2 + (1/r) * d_uz_dr + (1/r**2) * d2_uz_dtheta2 + d2_uz_dz2



    return [Lap_ur,Lap_ua,Lap_uz]

def compute_linear_scalar(fi) :
    d2S_dr2 = diff(fi, r, 2)
    dS_dr = diff(fi, r)
    d2S_dtheta2 = diff(fi, theta, 2)
    d2S_dz2 = diff(fi, z, 2)
    Lap_S = d2S_dr2 + (1/r) * dS_dr + (1/r**2) * d2S_dtheta2 + d2S_dz2
#    Lap_S = d2S_dr2 + d2S_dtheta2 + d2S_dz2
    
    return simplify(Lap_S)

def compute_nonlinear_momentum(ur,ua,uz):

    d_ur_dr     = diff(ur, r)
    d_ur_dtheta = diff(ur, theta)
    d_ur_dz     = diff(ur, z)

    d_ua_dr     = diff(ua, r)
    d_ua_dtheta = diff(ua, theta)
    d_ua_dz     = diff(ua, z)

    d_uz_dr     = diff(uz, r)
    d_uz_dtheta = diff(uz, theta)
    d_uz_dz     = diff(uz, z)

    NL_r     = ur * d_ur_dr + (ua / r) * d_ur_dtheta + uz * d_ur_dz - (ua**2 / r)
    NL_theta = ur * d_ua_dr + (ua / r) * d_ua_dtheta + uz * d_ua_dz + (ur * ua / r)
    NL_z     = ur * d_uz_dr + (ua / r) * d_uz_dtheta + uz * d_uz_dz

    NL_r     = simplify( NL_r )
    NL_theta = simplify( NL_theta )
    NL_z     = simplify( NL_z )
    
    return [NL_r,NL_theta,NL_z]



theta,z,r,t, nu, omega = symbols('theta z r t nu omega')
ri,ro,H= 1, 3 , 2*pi
d = ro-ri

# velocity field used by isa to check the accuracy in 3d
ur = cos(pi*t)*((r-1)*(r-3))**2*sin(8*theta)*sin(2.*pi*z) #+1/(2*pi  )*sin(  pi*(r-ri)/d)**2*cos(theta)*sin(2*pi*z/H)*(1+sin(2*pi*t))
ua = cos(pi*t)*sin(10.*pi*z)*((r-1)*(r-3))**2*sin(4*theta)#-1/(2*pi  )*sin(  pi*(r-ri)/d)**2*sin(theta)*sin(2*pi*z/H)*(1+sin(2*pi*t))
div_h = (1/r) * diff(r * ur, r) + (1/r) * diff(ua, theta) 
uz = ((r-1)*(r-3))**2*sin(theta)*sin(4.*pi*z)*0. + integrate(-div_h, z)
p = ((r-1)*(r-3))**2 * cos(theta) * cos(2.*pi*z)#   0.#( cos( pi*(r-ri)/d ) + cos(2*pi*z/H ) )*sin(theta)#*(1+sin(2*pi*t))
#ua = sin(theta)
#uz = sin(z)+cos(theta)+cos(r*pi)
#ur = r
phi = (sin(2.*pi*z))*sin(theta)*(r-1)*(r-3)

print("Div U = ",(1./r) * diff(r*ur,r) + (1./r) * diff(ua,theta) + diff(uz,z) )

print("values of ua on the bounds:",(ua.subs({r:ri}),ua.subs({r:ro})))
print("values of uz on the bounds:",(uz.subs({r:ri}),uz.subs({r:ro})))
print("values of ur on the bounds:",(ur.subs({r:ri}),ur.subs({r:ro})))


ST_phi = compute_linear_scalar(phi)
print('ST phi:',ST_phi)
T = sin(pi*(r-ri)/d) * cos(theta) * sin(2*pi*z/H)

div = (1/r) * diff(r * ur, r) + (1/r) * diff(ua, theta) + diff(uz, z)

grada_a,grada_z,grada_r = grad_cyl(ua)
gradz_a,gradz_z,gradz_r = grad_cyl(uz)
gradr_a,gradr_z,gradr_r = grad_cyl(ur)

NL_r, NL_a, NL_z = compute_nonlinear_momentum(ur,ua,uz)

Lap_r, Lap_a, Lap_z = compute_linear_momentum(ur,ua,uz)

Lap_T = compute_linear_scalar(T)

#
print('ua=',simplify(ua))
print('uz=',simplify(uz))
print('ur=',simplify(ur))
print(simplify(div))
print(simplify(div).is_zero)


phi_code=codegen(('udf_phi',phi), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))
ST_phi_code=codegen(('udf_st_phi',ST_phi), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))
ua_code=codegen(('udf_ua',ua), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))
uz_code=codegen(('udf_uz',uz), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))
ur_code=codegen(('udf_ur',ur), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))
p_code=codegen(('udf_p',p), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))
T_code=codegen(('udf_T',T), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))

NL_r_code=codegen(('udf_NLr',NL_r), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))
NL_a_code=codegen(('udf_NLa',NL_a), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))
NL_z_code=codegen(('udf_NLz',NL_z), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))

L_r_code=codegen(('udf_Lr',Lap_r), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))
L_a_code=codegen(('udf_La',Lap_a), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))
L_z_code=codegen(('udf_Lz',Lap_z), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))

grada_r_code=codegen(('udf_grada_r',grada_r), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))
grada_a_code=codegen(('udf_grada_a',grada_a), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))
grada_z_code=codegen(('udf_grada_z',grada_z), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))

gradz_r_code=codegen(('udf_gradz_r',gradz_r), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))
gradz_a_code=codegen(('udf_gradz_a',gradz_a), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))
gradz_z_code=codegen(('udf_gradz_z',gradz_z), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))

gradr_r_code=codegen(('udf_gradr_r',gradr_r), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))
gradr_a_code=codegen(('udf_gradr_a',gradr_a), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))
gradr_z_code=codegen(('udf_gradr_z',gradr_z), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))

Lap_T_code=codegen(('udf_Lt',Lap_T), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))

Div_code=codegen(('udf_div',div), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))

# Manufactured solutions for navier-stokes


dpdr,dpda,dpdz = diff(p,r),diff(p,theta),diff(p,z) 

d_ur_dt, d_ua_dt, d_uz_dt = diff(ur,t), diff(ua,t), diff(uz,t)

Fr = d_ur_dt +       dpdr - nu*Lap_r + NL_r
Fa = d_ua_dt + (1/r)*dpda - nu*Lap_a + NL_a
Fz = d_uz_dt +       dpdz - nu*Lap_z + NL_z

print('SA:',(Fa))
print('SZ:',(Fz))
print('SR:',(Fr))


scm_ur_code=codegen(('udf_scm_ur',Fr), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r, nu))
scm_ua_code=codegen(('udf_scm_ua',Fa), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r, nu))
scm_uz_code=codegen(('udf_scm_uz',Fz), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r, nu))

dpdr_code=codegen(('udf_dpdr',dpdr), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))
dpda_code=codegen(('udf_dpda',dpda/r), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))
dpdz_code=codegen(('udf_dpdz',dpdz), 'f95', 'my_project',header=False,argument_sequence=(t, theta, z, r))

#######








bof="""
module m_udf_mms
  implicit none
contains
"""
eof="""
end module m_udf_mms
"""

motif=""
motif=motif+bof

motif=motif+phi_code[0][1] 
motif=motif+ST_phi_code[0][1] 
motif=motif+ua_code[0][1] 
motif=motif+uz_code[0][1]
motif=motif+ur_code[0][1]
motif=motif+p_code[0][1]
motif=motif+T_code[0][1]
motif=motif+Div_code[0][1]

motif=motif+NL_r_code[0][1] 
motif=motif+NL_a_code[0][1]
motif=motif+NL_z_code[0][1]

motif=motif+L_r_code[0][1] 
motif=motif+L_a_code[0][1]
motif=motif+L_z_code[0][1]

motif=motif+grada_r_code[0][1] 
motif=motif+grada_a_code[0][1]
motif=motif+grada_z_code[0][1]

motif=motif+gradz_r_code[0][1] 
motif=motif+gradz_a_code[0][1]
motif=motif+gradz_z_code[0][1]

motif=motif+gradr_r_code[0][1] 
motif=motif+gradr_a_code[0][1]
motif=motif+gradr_z_code[0][1]

motif=motif+scm_ur_code[0][1] 
motif=motif+scm_ua_code[0][1]
motif=motif+scm_uz_code[0][1]

motif=motif+dpdr_code[0][1] 
motif=motif+dpda_code[0][1]
motif=motif+dpdz_code[0][1]


motif=motif+eof
motif = motif.replace('REAL*8 function', 'REAL*8 pure function')
motif = motif.replace('INTEGER*4 function', 'REAL*8 pure function')

f=open('m_udf_mms.f90','w')
f.write(motif)
f.close()
exit()




exit()





exit()

x,y,z,t, nu, omega = symbols('x y z t nu omega')

# velocity field used by isa to check the accuracy in 3d







lx=diff(ux,x,x)+diff(ux,y,y)+diff(ux,z,z)
ly=diff(uy,x,x)+diff(uy,y,y)+diff(uy,z,z)
lz=diff(uz,x,x)+diff(uz,y,y)+diff(uz,z,z)




#ux = sin(2*x)*sin(2*y)*sin(4*z)
#uy = sin(3*x)*cos(8*y)*cos(6*z)

# linear term
lx=diff(ux,x,x)+diff(ux,y,y)+diff(ux,z,z)
ly=diff(uy,x,x)+diff(uy,y,y)+diff(uy,z,z)
lz=diff(uz,x,x)+diff(uz,y,y)+diff(uz,z,z)

sx_cplx = - uy*omega + nu*lx
sy_cplx = + ux*omega + nu*ly


#Time derivative
dtux = diff(ux,t)
dtuy = diff(uy,t)
dtuz = diff(uz,t)

# divergence
div = diff(ux,x) + diff(uy,y) + diff(uz,z)
div=simplify(div)
print("divergence :",div)


# gradient 
gradpx = diff(p,x)
gradpy = diff(p,y)
gradpz = diff(p,z)

sfi = + (diff(p,x,x)+diff(p,y,y)+diff(p,z,z))




# non-linear term
nlx=ux*diff(ux,x)+uy*diff(ux,y)+uz*diff(ux,z)
nly=ux*diff(uy,x)+uy*diff(uy,y)+uz*diff(uy,z)
nlz=ux*diff(uz,x)+uy*diff(uz,y)+uz*diff(uz,z)



fx = dtux + nlx - nu*lx + gradpx + omega*uy
fy = dtuy + nly - nu*ly + gradpy - omega*ux
fz = dtuz + nlz - nu*lz + gradpz

fx=simplify(fx)
fy=simplify(fy)
fz=simplify(fz)


clx,cly,clz=-nlx+nu*lx+fx,-nly+nu*ly+fy,-nlz+nu*lz+fz
print(simplify(clx))
print(simplify(cly))
print(simplify(clz))

sx_cplx_code=codegen(('udf_sx_cplx',sx_cplx), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z, nu, omega))
sy_cplx_code=codegen(('udf_sy_cplx',sy_cplx), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z, nu, omega))

fx_code=codegen( ('udf_fx',fx), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z, nu, omega) )
fy_code=codegen( ('udf_fy',fy), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z, nu, omega) )
fz_code=codegen( ('udf_fz',fz), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z, nu, omega) )


ux_code=codegen(('udf_u_ex',ux), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))
uy_code=codegen(('udf_v_ex',uy), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))
uz_code=codegen(('udf_w_ex',uz), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))

dtux_code=codegen(('udf_dtu',dtux), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))
dtuy_code=codegen(('udf_dtv',dtuy), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))
dtuz_code=codegen(('udf_dtw',dtuz), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))

# check for the gradient
gradpx_code=codegen(('udf_grad_P_x',gradpx), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))
gradpy_code=codegen(('udf_grad_P_y',gradpy), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))
gradpz_code=codegen(('udf_grad_P_z',gradpz), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))

p_code=codegen(('udf_p_ex',p), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))

# check for the divergence
div_code=codegen(('udf_div', div), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))
sfi_code=codegen(('udf_sfi', sfi), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))



# check for the non-linear terms with extra linear terms
nlx_code=codegen(('udf_nlu',nlx), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))
nly_code=codegen(('udf_nlv',nly), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))
nlz_code=codegen(('udf_nlw',nlz), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))

# check for the linear terms
lx_code=codegen(('udf_lx',lx), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))
ly_code=codegen(('udf_ly',ly), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))
lz_code=codegen(('udf_lz',lz), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))






bof="""
module m_udf_mms_rotating
  implicit none
contains
"""
eof="""
end module m_udf_mms_rotating
"""

motif=""
motif=motif+bof

motif=motif+ux_code[0][1] # U_x
motif=motif+uy_code[0][1] # U_y
motif=motif+uz_code[0][1] # U_z

motif=motif+dtux_code[0][1] # dt U_x
motif=motif+dtuy_code[0][1] # dt U_y
motif=motif+dtuz_code[0][1] # dt U_z

motif=motif+gradpx_code[0][1] # dx P
motif=motif+gradpy_code[0][1] # dy P
motif=motif+gradpz_code[0][1] # dz P

motif=motif+div_code[0][1] # div U

motif=motif+lx_code[0][1] # lap(U_x)
motif=motif+ly_code[0][1] # lap(U_y)
motif=motif+lz_code[0][1] # lap(U_z)

motif=motif+nlx_code[0][1] # NL(U_x)
motif=motif+nly_code[0][1] # NL(U_y)
motif=motif+nlz_code[0][1] # NL(U_z)

motif=motif+fx_code[0][1] # 
motif=motif+fy_code[0][1] # 
motif=motif+fz_code[0][1] #

motif=motif+p_code[0][1] # 
motif=motif+sfi_code[0][1] # 

motif=motif+sx_cplx_code[0][1] #
motif=motif+sy_cplx_code[0][1] #

motif=motif+eof
motif = motif.replace('REAL*8 function', 'REAL*8 pure function')
motif = motif.replace('INTEGER*4 function', 'REAL*8 pure function')

#print(motif)
f=open('m_udf_mms_rotating.f90','w')
f.write(motif)
f.close()
exit()
