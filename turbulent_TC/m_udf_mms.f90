
module m_udf_mms
  implicit none
contains

REAL*8 pure function udf_orr(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_orr = -0.0365873432395162d0*(6.28318530717959d0*sin( &
      6.2831853071795862d0*x) + 0.229885057471264d0*cos( &
      6.2831853071795862d0*x))*(cos(t) + 1)**2*sin(6.2831853071795862d0 &
      *x)*sin(6.2831853071795862d0*y)**2*sin(6.2831853071795862d0*z)**2

end function

REAL*8 pure function udf_lift(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_lift = -0.459770114942529d0*(cos(t) + 1)**2*sin(6.2831853071795862d0 &
      *x)**2*sin(6.2831853071795862d0*y)**2*cos(6.2831853071795862d0*z) &
      **2

end function

REAL*8 pure function udf_push(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_push = 0.229885057471264d0*(cos(t) + 1)**2*sin(6.2831853071795862d0* &
      x)**2*sin(6.2831853071795862d0*z)**2*cos(6.2831853071795862d0*y) &
      **2

end function

REAL*8 pure function udf_u_ex(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_u_ex = (1.0d0*cos(t) + 1.0d0)*sin(6.2831853071795862d0*y)*sin( &
      6.2831853071795862d0*z)*cos(6.2831853071795862d0*x)

end function

REAL*8 pure function udf_v_ex(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_v_ex = (1.0d0*cos(t) + 1.0d0)*sin(6.2831853071795862d0*x)*sin( &
      6.2831853071795862d0*z)*cos(6.2831853071795862d0*y)

end function

REAL*8 pure function udf_w_ex(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_w_ex = -2*(1.0d0*cos(t) + 1.0d0)*sin(6.2831853071795862d0*x)*sin( &
      6.2831853071795862d0*y)*cos(6.2831853071795862d0*z)

end function

REAL*8 pure function udf_dtu(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_dtu = -1.0d0*sin(t)*sin(6.2831853071795862d0*y)*sin( &
      6.2831853071795862d0*z)*cos(6.2831853071795862d0*x)

end function

REAL*8 pure function udf_dtv(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_dtv = -1.0d0*sin(t)*sin(6.2831853071795862d0*x)*sin( &
      6.2831853071795862d0*z)*cos(6.2831853071795862d0*y)

end function

REAL*8 pure function udf_dtw(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_dtw = 2.0d0*sin(t)*sin(6.2831853071795862d0*x)*sin( &
      6.2831853071795862d0*y)*cos(6.2831853071795862d0*z)

end function

REAL*8 pure function udf_grad_P_x(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_grad_P_x = cos(x)

end function

REAL*8 pure function udf_grad_P_y(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_grad_P_y = cos(y)

end function

REAL*8 pure function udf_grad_P_z(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_grad_P_z = cos(z)

end function

REAL*8 pure function udf_div(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_div = 0

end function

REAL*8 pure function udf_lx(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_lx = -118.435252813072d0*(cos(t) + 1)*sin(6.2831853071795862d0*y)* &
      sin(6.2831853071795862d0*z)*cos(6.2831853071795862d0*x)

end function

REAL*8 pure function udf_ly(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_ly = -118.435252813072d0*(cos(t) + 1)*sin(6.2831853071795862d0*x)* &
      sin(6.2831853071795862d0*z)*cos(6.2831853071795862d0*y)

end function

REAL*8 pure function udf_lz(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_lz = 236.870505626145d0*(cos(t) + 1)*sin(6.2831853071795862d0*x)*sin &
      (6.2831853071795862d0*y)*cos(6.2831853071795862d0*z)

end function

REAL*8 pure function udf_nlu(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_nlu = -6.28318530717959d0*(1.0d0*cos(t) + 1.0d0)**2*sin( &
      6.2831853071795862d0*x)*sin(6.2831853071795862d0*y)**2*sin( &
      6.2831853071795862d0*z)**2*cos(6.2831853071795862d0*x) - &
      12.5663706143592d0*(1.0d0*cos(t) + 1.0d0)**2*sin( &
      6.2831853071795862d0*x)*sin(6.2831853071795862d0*y)**2*cos( &
      6.2831853071795862d0*x)*cos(6.2831853071795862d0*z)**2 + &
      6.28318530717959d0*(1.0d0*cos(t) + 1.0d0)**2*sin( &
      6.2831853071795862d0*x)*sin(6.2831853071795862d0*z)**2*cos( &
      6.2831853071795862d0*x)*cos(6.2831853071795862d0*y)**2

end function

REAL*8 pure function udf_nlv(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_nlv = -6.28318530717959d0*(1.0d0*cos(t) + 1.0d0)**2*sin( &
      6.2831853071795862d0*x)**2*sin(6.2831853071795862d0*y)*sin( &
      6.2831853071795862d0*z)**2*cos(6.2831853071795862d0*y) - &
      12.5663706143592d0*(1.0d0*cos(t) + 1.0d0)**2*sin( &
      6.2831853071795862d0*x)**2*sin(6.2831853071795862d0*y)*cos( &
      6.2831853071795862d0*y)*cos(6.2831853071795862d0*z)**2 + &
      6.28318530717959d0*(1.0d0*cos(t) + 1.0d0)**2*sin( &
      6.2831853071795862d0*y)*sin(6.2831853071795862d0*z)**2*cos( &
      6.2831853071795862d0*x)**2*cos(6.2831853071795862d0*y)

end function

REAL*8 pure function udf_nlw(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_nlw = -25.1327412287183d0*(1.0d0*cos(t) + 1.0d0)**2*sin( &
      6.2831853071795862d0*x)**2*sin(6.2831853071795862d0*y)**2*sin( &
      6.2831853071795862d0*z)*cos(6.2831853071795862d0*z) - &
      12.5663706143592d0*(1.0d0*cos(t) + 1.0d0)**2*sin( &
      6.2831853071795862d0*x)**2*sin(6.2831853071795862d0*z)*cos( &
      6.2831853071795862d0*y)**2*cos(6.2831853071795862d0*z) - &
      12.5663706143592d0*(1.0d0*cos(t) + 1.0d0)**2*sin( &
      6.2831853071795862d0*y)**2*sin(6.2831853071795862d0*z)*cos( &
      6.2831853071795862d0*x)**2*cos(6.2831853071795862d0*z)

end function

REAL*8 pure function udf_fx(t, x, y, z, nu)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z
REAL*8, intent(in) :: nu

udf_fx = 118.435252813072d0*nu*(cos(t) + 1)*sin(6.2831853071795862d0*y)* &
      sin(6.2831853071795862d0*z)*cos(6.2831853071795862d0*x) - &
      6.28318530717959d0*(cos(t) + 1)**2*sin(6.2831853071795862d0*x)* &
      sin(6.2831853071795862d0*y)**2*sin(6.2831853071795862d0*z)**2*cos &
      (6.2831853071795862d0*x) - 12.5663706143592d0*(cos(t) + 1)**2*sin &
      (6.2831853071795862d0*x)*sin(6.2831853071795862d0*y)**2*cos( &
      6.2831853071795862d0*x)*cos(6.2831853071795862d0*z)**2 + &
      6.28318530717959d0*(cos(t) + 1)**2*sin(6.2831853071795862d0*x)* &
      sin(6.2831853071795862d0*z)**2*cos(6.2831853071795862d0*x)*cos( &
      6.2831853071795862d0*y)**2 - 1.0d0*sin(t)*sin( &
      6.2831853071795862d0*y)*sin(6.2831853071795862d0*z)*cos( &
      6.2831853071795862d0*x) + cos(x)

end function

REAL*8 pure function udf_fy(t, x, y, z, nu)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z
REAL*8, intent(in) :: nu

udf_fy = 118.435252813072d0*nu*(cos(t) + 1)*sin(6.2831853071795862d0*x)* &
      sin(6.2831853071795862d0*z)*cos(6.2831853071795862d0*y) - &
      6.28318530717959d0*(cos(t) + 1)**2*sin(6.2831853071795862d0*x)**2 &
      *sin(6.2831853071795862d0*y)*sin(6.2831853071795862d0*z)**2*cos( &
      6.2831853071795862d0*y) - 12.5663706143592d0*(cos(t) + 1)**2*sin( &
      6.2831853071795862d0*x)**2*sin(6.2831853071795862d0*y)*cos( &
      6.2831853071795862d0*y)*cos(6.2831853071795862d0*z)**2 + &
      6.28318530717959d0*(cos(t) + 1)**2*sin(6.2831853071795862d0*y)* &
      sin(6.2831853071795862d0*z)**2*cos(6.2831853071795862d0*x)**2*cos &
      (6.2831853071795862d0*y) - 1.0d0*sin(t)*sin(6.2831853071795862d0* &
      x)*sin(6.2831853071795862d0*z)*cos(6.2831853071795862d0*y) + cos( &
      y)

end function

REAL*8 pure function udf_fz(t, x, y, z, nu)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z
REAL*8, intent(in) :: nu

udf_fz = -236.870505626145d0*nu*(cos(t) + 1)*sin(6.2831853071795862d0*x) &
      *sin(6.2831853071795862d0*y)*cos(6.2831853071795862d0*z) - &
      25.1327412287183d0*(cos(t) + 1)**2*sin(6.2831853071795862d0*x)**2 &
      *sin(6.2831853071795862d0*y)**2*sin(6.2831853071795862d0*z)*cos( &
      6.2831853071795862d0*z) - 12.5663706143592d0*(cos(t) + 1)**2*sin( &
      6.2831853071795862d0*x)**2*sin(6.2831853071795862d0*z)*cos( &
      6.2831853071795862d0*y)**2*cos(6.2831853071795862d0*z) - &
      12.5663706143592d0*(cos(t) + 1)**2*sin(6.2831853071795862d0*y)**2 &
      *sin(6.2831853071795862d0*z)*cos(6.2831853071795862d0*x)**2*cos( &
      6.2831853071795862d0*z) + 2.0d0*sin(t)*sin(6.2831853071795862d0*x &
      )*sin(6.2831853071795862d0*y)*cos(6.2831853071795862d0*z) + cos(z &
      )

end function

REAL*8 pure function udf_p_ex(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_p_ex = sin(x) + sin(y) + sin(z)

end function

REAL*8 pure function udf_sfi(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_sfi = -sin(x) - sin(y) - sin(z)

end function

REAL*8 pure function udf_sx_cplx(t, x, y, z, nu, omega)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z
REAL*8, intent(in) :: nu
REAL*8, intent(in) :: omega

udf_sx_cplx = -118.435252813072d0*nu*(cos(t) + 1)*sin( &
      6.2831853071795862d0*y)*sin(6.2831853071795862d0*z)*cos( &
      6.2831853071795862d0*x) - omega*(1.0d0*cos(t) + 1.0d0)*sin( &
      6.2831853071795862d0*x)*sin(6.2831853071795862d0*z)*cos( &
      6.2831853071795862d0*y)

end function

REAL*8 pure function udf_sy_cplx(t, x, y, z, nu, omega)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z
REAL*8, intent(in) :: nu
REAL*8, intent(in) :: omega

udf_sy_cplx = -118.435252813072d0*nu*(cos(t) + 1)*sin( &
      6.2831853071795862d0*x)*sin(6.2831853071795862d0*z)*cos( &
      6.2831853071795862d0*y) + omega*(1.0d0*cos(t) + 1.0d0)*sin( &
      6.2831853071795862d0*y)*sin(6.2831853071795862d0*z)*cos( &
      6.2831853071795862d0*x)

end function

end module m_udf_mms
