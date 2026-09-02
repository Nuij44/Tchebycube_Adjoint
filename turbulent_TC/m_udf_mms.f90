
module m_udf_mms
  implicit none
contains

REAL*8 pure function udf_u_ex(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_u_ex = 1 - z**2

end function

REAL*8 pure function udf_v_ex(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_v_ex = -4*z*(1 - z**2)*sin(y)

end function

REAL*8 pure function udf_w_ex(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_w_ex = -(1 - z**2)**2*cos(y)

end function

REAL*8 pure function udf_dtu(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_dtu = 0

end function

REAL*8 pure function udf_dtv(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_dtv = 0

end function

REAL*8 pure function udf_dtw(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_dtw = 0

end function

REAL*8 pure function udf_grad_P_x(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_grad_P_x = 0

end function

REAL*8 pure function udf_grad_P_y(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_grad_P_y = 0

end function

REAL*8 pure function udf_grad_P_z(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_grad_P_z = 0

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

udf_lx = -2

end function

REAL*8 pure function udf_ly(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_ly = -4*z*(z**2 - 1)*sin(y) + 24*z*sin(y)

end function

REAL*8 pure function udf_lz(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_lz = (z**2 - 1)**2*cos(y) - 4*(3*z**2 - 1)*cos(y)

end function

REAL*8 pure function udf_nlu(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_nlu = 2*z*(1 - z**2)**2*cos(y) + (1 - z**2)**2*cos(y)

end function

REAL*8 pure function udf_nlv(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_nlv = 16*z**2*(1 - z**2)**2*sin(y)*cos(y) - (1 - z**2)**2*(8*z**2* &
      sin(y) - 4*(1 - z**2)*sin(y))*cos(y)

end function

REAL*8 pure function udf_nlw(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_nlw = -4*z*(1 - z**2)**3*sin(y)**2 - 4*z*(1 - z**2)**3*cos(y)**2

end function

REAL*8 pure function udf_fx(t, x, y, z, nu, omega)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z
REAL*8, intent(in) :: nu
REAL*8, intent(in) :: omega

udf_fx = -2*nu - omega*(z**2 - 1)**2*cos(y) + 2*z*(z**2 - 1)**2*cos(y) + &
      (z**2 - 1)**2*cos(y)

end function

REAL*8 pure function udf_fy(t, x, y, z, nu, omega)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z
REAL*8, intent(in) :: nu
REAL*8, intent(in) :: omega

udf_fy = 4*(-nu*z**3 + 7*nu*z + z**6*cos(y) - z**4*cos(y) - z**2*cos(y) &
      + cos(y))*sin(y)

end function

REAL*8 pure function udf_fz(t, x, y, z, nu, omega)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z
REAL*8, intent(in) :: nu
REAL*8, intent(in) :: omega

udf_fz = nu*z**4*cos(y) - 14*nu*z**2*cos(y) + 5*nu*cos(y) + omega*z**2 - &
      omega + 4*z**7 - 12*z**5 + 12*z**3 - 4*z

end function

REAL*8 pure function udf_p_ex(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_p_ex = 0

end function

REAL*8 pure function udf_sfi(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_sfi = 0

end function

REAL*8 pure function udf_sx_cplx(t, x, y, z, nu, omega)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z
REAL*8, intent(in) :: nu
REAL*8, intent(in) :: omega

udf_sx_cplx = -2*nu + 4*omega*z*(1 - z**2)*sin(y)

end function

REAL*8 pure function udf_sy_cplx(t, x, y, z, nu, omega)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z
REAL*8, intent(in) :: nu
REAL*8, intent(in) :: omega

udf_sy_cplx = nu*(-4*z*(z**2 - 1)*sin(y) + 24*z*sin(y)) + omega*(1 - z** &
      2)

end function

end module m_udf_mms
