
module m_udf_mms
  implicit none
contains

REAL*8 pure function udf_phi(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_phi = (r - 3)*(r - 1)*sin(theta)*sin(6.2831853071795865d0*z)

end function

REAL*8 pure function udf_st_phi(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_st_phi = (r**2*(-4*pi**2*(r - 3)*(r - 1) + 2) + 2*r*(r - 2) - (r - 3 &
      )*(r - 1))*sin(theta)*sin(6.2831853071795865d0*z)/r**2

end function

REAL*8 pure function udf_ua(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_ua = 10.0d0*(r - 3)**2*(r - 1)**2*sin(theta)*sin( &
      6.2831853071795865d0*z)*cos(12.566370614359173d0*t)

end function

REAL*8 pure function udf_uz(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_uz = 5.0d0*(r - 3)**2*(r - 1)**2*cos(theta)*cos(12.566370614359173d0 &
      *t)*cos(6.2831853071795865d0*z)/(pi*r)

end function

REAL*8 pure function udf_ur(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_ur = 0

end function

REAL*8 pure function udf_p(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_p = 0.0d0

end function

REAL*8 pure function udf_T(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_T = sin(z)*sin(1.5707963267948966d0*r - 1.5707963267948966d0)*cos( &
      theta)

end function

REAL*8 pure function udf_div(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_div = 0

end function

REAL*8 pure function udf_NLr(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_NLr = -100.0d0*(r - 3)**4*(r - 1)**4*sin(theta)**2*sin( &
      6.2831853071795865d0*z)**2*cos(12.566370614359173d0*t)**2/r

end function

REAL*8 pure function udf_NLa(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_NLa = 100.0d0*(r - 3)**4*(r - 1)**4*sin(theta)*cos(theta)*cos( &
      12.566370614359173d0*t)**2/r

end function

REAL*8 pure function udf_NLz(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_NLz = -50.0d0*(r - 3)**4*(r - 1)**4*sin(6.2831853071795865d0*z)*cos( &
      12.566370614359173d0*t)**2*cos(6.2831853071795865d0*z)/(pi*r**2)

end function

REAL*8 pure function udf_Lr(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_Lr = -20.0d0*(r - 3)**2*(r - 1)**2*sin(6.2831853071795865d0*z)*cos( &
      theta)*cos(12.566370614359173d0*t)/r**2

end function

REAL*8 pure function udf_La(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_La = -40.0d0*pi**2*(r - 3)**2*(r - 1)**2*sin(theta)*sin( &
      6.2831853071795865d0*z)*cos(12.566370614359173d0*t) + (20.0d0*(r &
      - 3)**2 + 80.0d0*(r - 3)*(r - 1) + 20.0d0*(r - 1)**2)*sin(theta)* &
      sin(6.2831853071795865d0*z)*cos(12.566370614359173d0*t) + 1.0d0*( &
      10.0d0*(r - 3)**2*(2*r - 2)*sin(theta)*sin(6.2831853071795865d0*z &
      )*cos(12.566370614359173d0*t) + 10.0d0*(r - 1)**2*(2*r - 6)*sin( &
      theta)*sin(6.2831853071795865d0*z)*cos(12.566370614359173d0*t))/r &
      - 20.0d0*(r - 3)**2*(r - 1)**2*sin(theta)*sin( &
      6.2831853071795865d0*z)*cos(12.566370614359173d0*t)/r**2

end function

REAL*8 pure function udf_Lz(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_Lz = -20.0d0*pi*(r - 3)**2*(r - 1)**2*cos(theta)*cos( &
      12.566370614359173d0*t)*cos(6.2831853071795865d0*z)/r + 1.0d0*( &
      5.0d0*(r - 3)**2*(2*r - 2)*cos(theta)*cos(12.566370614359173d0*t) &
      *cos(6.2831853071795865d0*z)/(pi*r) + 5.0d0*(r - 1)**2*(2*r - 6)* &
      cos(theta)*cos(12.566370614359173d0*t)*cos(6.2831853071795865d0*z &
      )/(pi*r) - 5.0d0*(r - 3)**2*(r - 1)**2*cos(theta)*cos( &
      12.566370614359173d0*t)*cos(6.2831853071795865d0*z)/(pi*r**2))/r &
      + (10.0d0*(r - 3)**2 + 40.0d0*(r - 3)*(r - 1) + 10.0d0*(r - 1)**2 &
      - 20.0d0*(r - 3)**2*(r - 1)/r - 20.0d0*(r - 3)*(r - 1)**2/r + &
      10.0d0*(r - 3)**2*(r - 1)**2/r**2)*cos(theta)*cos( &
      12.566370614359173d0*t)*cos(6.2831853071795865d0*z)/(pi*r) - &
      5.0d0*(r - 3)**2*(r - 1)**2*cos(theta)*cos(12.566370614359173d0*t &
      )*cos(6.2831853071795865d0*z)/(pi*r**3)

end function

REAL*8 pure function udf_grada_r(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_grada_r = 10.0d0*(r - 3)**2*(2*r - 2)*sin(theta)*sin( &
      6.2831853071795865d0*z)*cos(12.566370614359173d0*t) + 10.0d0*(r - &
      1)**2*(2*r - 6)*sin(theta)*sin(6.2831853071795865d0*z)*cos( &
      12.566370614359173d0*t)

end function

REAL*8 pure function udf_grada_a(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_grada_a = 10.0d0*(r - 3)**2*(r - 1)**2*sin(6.2831853071795865d0*z)* &
      cos(theta)*cos(12.566370614359173d0*t)/r

end function

REAL*8 pure function udf_grada_z(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_grada_z = 20.0d0*pi*(r - 3)**2*(r - 1)**2*sin(theta)*cos( &
      12.566370614359173d0*t)*cos(6.2831853071795865d0*z)

end function

REAL*8 pure function udf_gradz_r(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_gradz_r = 5.0d0*(r - 3)**2*(2*r - 2)*cos(theta)*cos( &
      12.566370614359173d0*t)*cos(6.2831853071795865d0*z)/(pi*r) + &
      5.0d0*(r - 1)**2*(2*r - 6)*cos(theta)*cos(12.566370614359173d0*t) &
      *cos(6.2831853071795865d0*z)/(pi*r) - 5.0d0*(r - 3)**2*(r - 1)**2 &
      *cos(theta)*cos(12.566370614359173d0*t)*cos(6.2831853071795865d0* &
      z)/(pi*r**2)

end function

REAL*8 pure function udf_gradz_a(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_gradz_a = -5.0d0*(r - 3)**2*(r - 1)**2*sin(theta)*cos( &
      12.566370614359173d0*t)*cos(6.2831853071795865d0*z)/(pi*r**2)

end function

REAL*8 pure function udf_gradz_z(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_gradz_z = -10.0d0*(r - 3)**2*(r - 1)**2*sin(6.2831853071795865d0*z)* &
      cos(theta)*cos(12.566370614359173d0*t)/r

end function

REAL*8 pure function udf_gradr_r(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_gradr_r = 0

end function

REAL*8 pure function udf_gradr_a(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_gradr_a = 0

end function

REAL*8 pure function udf_gradr_z(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_gradr_z = 0

end function

REAL*8 pure function udf_scm_ur(t, theta, z, r, nu)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r
REAL*8, intent(in) :: nu

udf_scm_ur = 20.0d0*nu*(r - 3)**2*(r - 1)**2*sin(6.2831853071795865d0*z) &
      *cos(theta)*cos(12.566370614359173d0*t)/r**2 - 100.0d0*(r - 3)**4 &
      *(r - 1)**4*sin(theta)**2*sin(6.2831853071795865d0*z)**2*cos( &
      12.566370614359173d0*t)**2/r

end function

REAL*8 pure function udf_scm_ua(t, theta, z, r, nu)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r
REAL*8, intent(in) :: nu

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_scm_ua = -nu*(-40.0d0*pi**2*(r - 3)**2*(r - 1)**2*sin(theta)*sin( &
      6.2831853071795865d0*z)*cos(12.566370614359173d0*t) + (20.0d0*(r &
      - 3)**2 + 80.0d0*(r - 3)*(r - 1) + 20.0d0*(r - 1)**2)*sin(theta)* &
      sin(6.2831853071795865d0*z)*cos(12.566370614359173d0*t) + 1.0d0*( &
      10.0d0*(r - 3)**2*(2*r - 2)*sin(theta)*sin(6.2831853071795865d0*z &
      )*cos(12.566370614359173d0*t) + 10.0d0*(r - 1)**2*(2*r - 6)*sin( &
      theta)*sin(6.2831853071795865d0*z)*cos(12.566370614359173d0*t))/r &
      - 20.0d0*(r - 3)**2*(r - 1)**2*sin(theta)*sin( &
      6.2831853071795865d0*z)*cos(12.566370614359173d0*t)/r**2) - &
      40.0d0*pi*(r - 3)**2*(r - 1)**2*sin(theta)*sin( &
      12.566370614359173d0*t)*sin(6.2831853071795865d0*z) + 100.0d0*(r &
      - 3)**4*(r - 1)**4*sin(theta)*cos(theta)*cos(12.566370614359173d0 &
      *t)**2/r

end function

REAL*8 pure function udf_scm_uz(t, theta, z, r, nu)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r
REAL*8, intent(in) :: nu

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_scm_uz = -nu*(-20.0d0*pi*(r - 3)**2*(r - 1)**2*cos(theta)*cos( &
      12.566370614359173d0*t)*cos(6.2831853071795865d0*z)/r + 1.0d0*( &
      5.0d0*(r - 3)**2*(2*r - 2)*cos(theta)*cos(12.566370614359173d0*t) &
      *cos(6.2831853071795865d0*z)/(pi*r) + 5.0d0*(r - 1)**2*(2*r - 6)* &
      cos(theta)*cos(12.566370614359173d0*t)*cos(6.2831853071795865d0*z &
      )/(pi*r) - 5.0d0*(r - 3)**2*(r - 1)**2*cos(theta)*cos( &
      12.566370614359173d0*t)*cos(6.2831853071795865d0*z)/(pi*r**2))/r &
      + (10.0d0*(r - 3)**2 + 40.0d0*(r - 3)*(r - 1) + 10.0d0*(r - 1)**2 &
      - 20.0d0*(r - 3)**2*(r - 1)/r - 20.0d0*(r - 3)*(r - 1)**2/r + &
      10.0d0*(r - 3)**2*(r - 1)**2/r**2)*cos(theta)*cos( &
      12.566370614359173d0*t)*cos(6.2831853071795865d0*z)/(pi*r) - &
      5.0d0*(r - 3)**2*(r - 1)**2*cos(theta)*cos(12.566370614359173d0*t &
      )*cos(6.2831853071795865d0*z)/(pi*r**3)) - 20.0d0*(r - 3)**2*(r - &
      1)**2*sin(12.566370614359173d0*t)*cos(theta)*cos( &
      6.2831853071795865d0*z)/r - 50.0d0*(r - 3)**4*(r - 1)**4*sin( &
      6.2831853071795865d0*z)*cos(12.566370614359173d0*t)**2*cos( &
      6.2831853071795865d0*z)/(pi*r**2)

end function

REAL*8 pure function udf_dpdr(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_dpdr = 0

end function

REAL*8 pure function udf_dpda(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_dpda = 0

end function

REAL*8 pure function udf_dpdz(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_dpdz = 0

end function

end module m_udf_mms
