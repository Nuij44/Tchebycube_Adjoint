
module m_udf_inc_ns_cart
  implicit none
contains

REAL*8 function udf_u_ex(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_u_ex = -2.0d0*(cos(31.415926535897932d0*t)**2 + 1.0d0)*sin( &
      3.1415926535897932d0*x)**2*sin(6.2831853071795865d0*y)*sin( &
      6.2831853071795865d0*z)

end function

REAL*8 function udf_v_ex(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_v_ex = (cos(31.415926535897932d0*t)**2 + 1.0d0)*sin( &
      6.2831853071795865d0*x)*sin(3.1415926535897932d0*y)**2*sin( &
      6.2831853071795865d0*z)

end function

REAL*8 function udf_w_ex(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_w_ex = (cos(31.415926535897932d0*t)**2 + 1.0d0)*sin( &
      6.2831853071795865d0*x)*sin(6.2831853071795865d0*y)*sin( &
      3.1415926535897932d0*z)**2

end function

REAL*8 function udf_dtu(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_dtu = 40.0d0*pi*sin(31.415926535897932d0*t)*sin(3.1415926535897932d0 &
      *x)**2*sin(6.2831853071795865d0*y)*sin(6.2831853071795865d0*z)* &
      cos(31.415926535897932d0*t)

end function

REAL*8 function udf_dtv(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_dtv = -20*pi*sin(31.415926535897932d0*t)*sin(6.2831853071795865d0*x) &
      *sin(3.1415926535897932d0*y)**2*sin(6.2831853071795865d0*z)*cos( &
      31.415926535897932d0*t)

end function

REAL*8 function udf_dtw(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_dtw = -20*pi*sin(31.415926535897932d0*t)*sin(6.2831853071795865d0*x) &
      *sin(6.2831853071795865d0*y)*sin(3.1415926535897932d0*z)**2*cos( &
      31.415926535897932d0*t)

end function

REAL*8 function grad_P_x(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

grad_P_x = 2*sin(3.1415926535897932d0*x)*sin(3.1415926535897932d0*y)**2* &
      sin(3.1415926535897932d0*z)**2*cos(3.1415926535897932d0*x)

end function

REAL*8 function grad_P_y(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

grad_P_y = 2*sin(3.1415926535897932d0*x)**2*sin(3.1415926535897932d0*y)* &
      sin(3.1415926535897932d0*z)**2*cos(3.1415926535897932d0*y)

end function

REAL*8 function grad_P_z(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

grad_P_z = 2*sin(3.1415926535897932d0*x)**2*sin(3.1415926535897932d0*y) &
      **2*sin(3.1415926535897932d0*z)*cos(3.1415926535897932d0*z)

end function

INTEGER*4 function udf_div(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

udf_div = 0

end function

REAL*8 function udf_lx(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_lx = 4.0d0*pi**2*(sin(3.1415926535897932d0*x)**2 - cos( &
      3.1415926535897932d0*x)**2)*(cos(31.415926535897932d0*t)**2 + &
      1.0d0)*sin(6.2831853071795865d0*y)*sin(6.2831853071795865d0*z) + &
      16.0d0*pi**2*(cos(31.415926535897932d0*t)**2 + 1.0d0)*sin( &
      3.1415926535897932d0*x)**2*sin(6.2831853071795865d0*y)*sin( &
      6.2831853071795865d0*z)

end function

REAL*8 function udf_ly(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_ly = -2*pi**2*(sin(3.1415926535897932d0*y)**2 - cos( &
      3.1415926535897932d0*y)**2)*(cos(31.415926535897932d0*t)**2 + &
      1.0d0)*sin(6.2831853071795865d0*x)*sin(6.2831853071795865d0*z) - &
      8*pi**2*(cos(31.415926535897932d0*t)**2 + 1.0d0)*sin( &
      6.2831853071795865d0*x)*sin(3.1415926535897932d0*y)**2*sin( &
      6.2831853071795865d0*z)

end function

REAL*8 function udf_lz(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_lz = -2*pi**2*(sin(3.1415926535897932d0*z)**2 - cos( &
      3.1415926535897932d0*z)**2)*(cos(31.415926535897932d0*t)**2 + &
      1.0d0)*sin(6.2831853071795865d0*x)*sin(6.2831853071795865d0*y) - &
      8*pi**2*(cos(31.415926535897932d0*t)**2 + 1.0d0)*sin( &
      6.2831853071795865d0*x)*sin(6.2831853071795865d0*y)*sin( &
      3.1415926535897932d0*z)**2

end function

REAL*8 function advection_x(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

REAL*8, parameter :: pi = 3.1415926535897932d0
advection_x = 8.0d0*pi*(cos(31.415926535897932d0*t)**2 + 1.0d0)**2*sin( &
      3.1415926535897932d0*x)**3*sin(6.2831853071795865d0*y)**2*sin( &
      6.2831853071795865d0*z)**2*cos(3.1415926535897932d0*x) - 4.0d0*pi &
      *(cos(31.415926535897932d0*t)**2 + 1.0d0)**2*sin( &
      3.1415926535897932d0*x)**2*sin(6.2831853071795865d0*x)*sin( &
      3.1415926535897932d0*y)**2*sin(6.2831853071795865d0*z)**2*cos( &
      6.2831853071795865d0*y) - 4.0d0*pi*(cos(31.415926535897932d0*t)** &
      2 + 1.0d0)**2*sin(3.1415926535897932d0*x)**2*sin( &
      6.2831853071795865d0*x)*sin(6.2831853071795865d0*y)**2*sin( &
      3.1415926535897932d0*z)**2*cos(6.2831853071795865d0*z)

end function

REAL*8 function advection_y(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

REAL*8, parameter :: pi = 3.1415926535897932d0
advection_y = -4.0d0*pi*(cos(31.415926535897932d0*t)**2 + 1.0d0)**2*sin( &
      3.1415926535897932d0*x)**2*sin(3.1415926535897932d0*y)**2*sin( &
      6.2831853071795865d0*y)*sin(6.2831853071795865d0*z)**2*cos( &
      6.2831853071795865d0*x) + 2*pi*(cos(31.415926535897932d0*t)**2 + &
      1.0d0)**2*sin(6.2831853071795865d0*x)**2*sin(3.1415926535897932d0 &
      *y)**3*sin(6.2831853071795865d0*z)**2*cos(3.1415926535897932d0*y &
      ) + 2*pi*(cos(31.415926535897932d0*t)**2 + 1.0d0)**2*sin( &
      6.2831853071795865d0*x)**2*sin(3.1415926535897932d0*y)**2*sin( &
      6.2831853071795865d0*y)*sin(3.1415926535897932d0*z)**2*cos( &
      6.2831853071795865d0*z)

end function

REAL*8 function advection_z(t, x, y, z)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: x
REAL*8, intent(in) :: y
REAL*8, intent(in) :: z

REAL*8, parameter :: pi = 3.1415926535897932d0
advection_z = -4.0d0*pi*(cos(31.415926535897932d0*t)**2 + 1.0d0)**2*sin( &
      3.1415926535897932d0*x)**2*sin(6.2831853071795865d0*y)**2*sin( &
      3.1415926535897932d0*z)**2*sin(6.2831853071795865d0*z)*cos( &
      6.2831853071795865d0*x) + 2*pi*(cos(31.415926535897932d0*t)**2 + &
      1.0d0)**2*sin(6.2831853071795865d0*x)**2*sin(3.1415926535897932d0 &
      *y)**2*sin(3.1415926535897932d0*z)**2*sin(6.2831853071795865d0*z) &
      *cos(6.2831853071795865d0*y) + 2*pi*(cos(31.415926535897932d0*t) &
      **2 + 1.0d0)**2*sin(6.2831853071795865d0*x)**2*sin( &
      6.2831853071795865d0*y)**2*sin(3.1415926535897932d0*z)**3*cos( &
      3.1415926535897932d0*z)

end function

end module m_udf_inc_ns_cart
