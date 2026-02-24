#!/usr/bin/python
from sympy import *
from sympy.utilities.codegen import codegen



x,y,z,t = symbols('x y z t')

# velocity field used by isa to check the accuracy in 3d
ux=-2.*(sin(pi*x))**2*sin(2.*pi*y)*sin(2.*pi*z)*(1.+(cos(4.*pi*t))**2)

uy=sin(2.*pi*x)*(sin(pi*y))**2*sin(2.*pi*z)*(1.+(cos(4.*pi*t))**2)

uz=sin(2.*pi*x)*sin(2.*pi*y)*(sin(pi*z))**2*(1.+(cos(4.*pi*t))**2)

p=(sin(pi*x))**2*(sin(pi*y))**2*(sin(pi*z))**2/pi



#Time derivative
dtux = diff(ux,t)
dtuy = diff(uy,t)
dtuz = diff(uz,t)

# divergence
div = simplify(diff(ux,x) + diff(uy,y) + diff(uz,z))
print(div)

# gradient 
gradpx = diff(p,x)
gradpy = diff(p,y)
gradpz = diff(p,z)

# linear term
lx=diff(ux,x,x)+diff(ux,y,y)+diff(ux,z,z)
ly=diff(uy,x,x)+diff(uy,y,y)+diff(uy,z,z)
lz=diff(uz,x,x)+diff(uz,y,y)+diff(uz,z,z)

# non-linear term
nlx=ux*diff(ux,x)+uy*diff(ux,y)+uz*diff(ux,z)
nly=ux*diff(uy,x)+uy*diff(uy,y)+uz*diff(uy,z)
nlz=ux*diff(uz,x)+uy*diff(uz,y)+uz*diff(uz,z)

ux_code=codegen(('udf_u_ex',ux), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))
uy_code=codegen(('udf_v_ex',uy), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))
uz_code=codegen(('udf_w_ex',uz), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))

dtux_code=codegen(('udf_dtu',dtux), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))
dtuy_code=codegen(('udf_dtv',dtuy), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))
dtuz_code=codegen(('udf_dtw',dtuz), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))

# check for the gradient
gradpx_code=codegen(('grad_P_x',gradpx), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))
gradpy_code=codegen(('grad_P_y',gradpy), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))
gradpz_code=codegen(('grad_P_z',gradpz), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))

# check for the divergence
div_code=codegen(('udf_div', div), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))

# check for the non-linear terms with extra linear terms
nlx_code=codegen(('advection_x',nlx), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))
nly_code=codegen(('advection_y',nly), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))
nlz_code=codegen(('advection_z',nlz), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))

# check for the linear terms
lx_code=codegen(('udf_lx',lx), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))
ly_code=codegen(('udf_ly',ly), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))
lz_code=codegen(('udf_lz',lz), 'f95', 'my_project',header=False,argument_sequence=(t, x, y, z))



bof="""
module m_udf_inc_ns_cart
  implicit none
contains
"""
eof="""
end module m_udf_inc_ns_cart
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
motif=motif+eof
print(motif)
f=open('m_udf_inc_ns_cart.f90','w')
f.write(motif)
f.close()
exit()
