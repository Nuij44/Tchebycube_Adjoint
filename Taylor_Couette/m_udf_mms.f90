
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

udf_ua = (r - 3)**2*(r - 1)**2*sin(theta)*sin(6.2831853071795865d0*z)

end function

REAL*8 pure function udf_uz(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_uz = (1.0d0/2.0d0)*(r - 3)**2*(r - 1)**2*cos(theta)*cos( &
      6.2831853071795865d0*z)/(pi*r) - (-5.0d0/2.0d0*r**4*sin(theta)* &
      cos(6.2831853071795865d0*z)/pi + 16*r**3*sin(theta)*cos( &
      6.2831853071795865d0*z)/pi - 33*r**2*sin(theta)*cos( &
      6.2831853071795865d0*z)/pi + 24*r*sin(theta)*cos( &
      6.2831853071795865d0*z)/pi - 9.0d0/2.0d0*sin(theta)*cos( &
      6.2831853071795865d0*z)/pi)/r

end function

REAL*8 pure function udf_ur(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_ur = (r - 3)**2*(r - 1)**2*sin(theta)*sin(6.2831853071795865d0*z)

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

udf_div = (r*(r - 3)**2*(2*r - 2)*sin(theta)*sin(6.2831853071795865d0*z &
      ) + r*(r - 1)**2*(2*r - 6)*sin(theta)*sin(6.2831853071795865d0*z &
      ) + (r - 3)**2*(r - 1)**2*sin(theta)*sin(6.2831853071795865d0*z)) &
      /r - (5*r**4*sin(theta)*sin(6.2831853071795865d0*z) - 32*r**3*sin &
      (theta)*sin(6.2831853071795865d0*z) + 66*r**2*sin(theta)*sin( &
      6.2831853071795865d0*z) - 48*r*sin(theta)*sin( &
      6.2831853071795865d0*z) + 9*sin(theta)*sin(6.2831853071795865d0*z &
      ))/r

end function

REAL*8 pure function udf_NLr(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_NLr = (r - 3)**2*(r - 1)**2*(4*r*(r - 3)*(r - 2)*(r - 1)*sin(theta)* &
      sin(6.2831853071795865d0*z)**2 - (r - 3)**2*(r - 1)**2*sin(theta) &
      *sin(6.2831853071795865d0*z)**2 + (r - 3)**2*(r - 1)**2*sin( &
      6.2831853071795865d0*z)**2*cos(theta) + (5*r**4*sin(theta) - 32*r &
      **3*sin(theta) + 66*r**2*sin(theta) - 48*r*sin(theta) + (r - 3)** &
      2*(r - 1)**2*cos(theta) + 9*sin(theta))*cos(6.2831853071795865d0* &
      z)**2)*sin(theta)/r

end function

REAL*8 pure function udf_NLa(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_NLa = (r - 3)**2*(r - 1)**2*(5*r**4*sin(theta) + r**4*cos(theta) - &
      32*r**3*sin(theta) - 8*r**3*cos(theta) + 63*r**2*sin(theta) + 3* &
      sqrt(2.0d0)*r**2*sin(theta + 0.78539816339744831d0) + 19*r**2*cos &
      (theta) - 36*r*sin(theta) - 12*sqrt(2.0d0)*r*sin(theta + &
      0.78539816339744831d0) - 12*r*cos(theta) + 9*sqrt(2.0d0)*sin( &
      theta + 0.78539816339744831d0))*sin(theta)/r

end function

REAL*8 pure function udf_NLz(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_NLz = (1.0d0/2.0d0)*(r - 3)**2*(r - 1)**2*(-r**4*sin(2.0d0*theta) + &
      5*r**4*cos(2.0d0*theta) - 6*r**4 + 8*r**3*sin(2.0d0*theta) - 28*r &
      **3*cos(2.0d0*theta) + 36*r**3 - 22*r**2*sin(2.0d0*theta) + 54*r &
      **2*cos(2.0d0*theta) - 76*r**2 + 24*r*sin(2.0d0*theta) - 36*r*cos &
      (2.0d0*theta) + 60*r - 9*sin(2.0d0*theta) + 9*cos(2.0d0*theta) - &
      18)*sin(6.2831853071795865d0*z)*cos(6.2831853071795865d0*z)/(pi*r &
      **2)

end function

REAL*8 pure function udf_Lr(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_Lr = -4*pi**2*(r - 3)**2*(r - 1)**2*sin(theta)*sin( &
      6.2831853071795865d0*z) + 2*((r - 3)**2 + 4*(r - 3)*(r - 1) + (r &
      - 1)**2)*sin(theta)*sin(6.2831853071795865d0*z) + 1.0d0*((r - 3) &
      **2*(2*r - 2)*sin(theta)*sin(6.2831853071795865d0*z) + (r - 1)**2 &
      *(2*r - 6)*sin(theta)*sin(6.2831853071795865d0*z))/r - 2.0d0*(r - &
      3)**2*(r - 1)**2*sin(theta)*sin(6.2831853071795865d0*z)/r**2 - &
      2.0d0*(r - 3)**2*(r - 1)**2*sin(6.2831853071795865d0*z)*cos(theta &
      )/r**2

end function

REAL*8 pure function udf_La(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_La = -4*pi**2*(r - 3)**2*(r - 1)**2*sin(theta)*sin( &
      6.2831853071795865d0*z) + 2*((r - 3)**2 + 4*(r - 3)*(r - 1) + (r &
      - 1)**2)*sin(theta)*sin(6.2831853071795865d0*z) + 1.0d0*((r - 3) &
      **2*(2*r - 2)*sin(theta)*sin(6.2831853071795865d0*z) + (r - 1)**2 &
      *(2*r - 6)*sin(theta)*sin(6.2831853071795865d0*z))/r - 2.0d0*(r - &
      3)**2*(r - 1)**2*sin(theta)*sin(6.2831853071795865d0*z)/r**2 + &
      2.0d0*(r - 3)**2*(r - 1)**2*sin(6.2831853071795865d0*z)*cos(theta &
      )/r**2

end function

REAL*8 pure function udf_Lz(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_Lz = -2*pi*((r - 3)**2*(r - 1)**2*cos(theta) + (5*r**4 - 32*r**3 + &
      66*r**2 - 48*r + 9)*sin(theta))*cos(6.2831853071795865d0*z)/r + &
      1.0d0*((1.0d0/2.0d0)*(r - 3)**2*(2*r - 2)*cos(theta)*cos( &
      6.2831853071795865d0*z)/(pi*r) + (1.0d0/2.0d0)*(r - 1)**2*(2*r - &
      6)*cos(theta)*cos(6.2831853071795865d0*z)/(pi*r) - (-10*r**3*sin( &
      theta)*cos(6.2831853071795865d0*z)/pi + 48*r**2*sin(theta)*cos( &
      6.2831853071795865d0*z)/pi - 66*r*sin(theta)*cos( &
      6.2831853071795865d0*z)/pi + 24*sin(theta)*cos( &
      6.2831853071795865d0*z)/pi)/r - 1.0d0/2.0d0*(r - 3)**2*(r - 1)**2 &
      *cos(theta)*cos(6.2831853071795865d0*z)/(pi*r**2) + (-5.0d0/2.0d0 &
      *r**4*sin(theta)*cos(6.2831853071795865d0*z)/pi + 16*r**3*sin( &
      theta)*cos(6.2831853071795865d0*z)/pi - 33*r**2*sin(theta)*cos( &
      6.2831853071795865d0*z)/pi + 24*r*sin(theta)*cos( &
      6.2831853071795865d0*z)/pi - 9.0d0/2.0d0*sin(theta)*cos( &
      6.2831853071795865d0*z)/pi)/r**2)/r + ((r - 3)**2*cos(theta) + 4* &
      (r - 3)*(r - 1)*cos(theta) + (r - 1)**2*cos(theta) + 6*(5*r**2 - &
      16*r + 11)*sin(theta) - 2*(r - 3)**2*(r - 1)*cos(theta)/r - 2*(r &
      - 3)*(r - 1)**2*cos(theta)/r - 4*(5*r**3 - 24*r**2 + 33*r - 12)* &
      sin(theta)/r + (r - 3)**2*(r - 1)**2*cos(theta)/r**2 + (5*r**4 - &
      32*r**3 + 66*r**2 - 48*r + 9)*sin(theta)/r**2)*cos( &
      6.2831853071795865d0*z)/(pi*r) - 0.5d0*((r - 3)**2*(r - 1)**2*cos &
      (theta) + (5*r**4 - 32*r**3 + 66*r**2 - 48*r + 9)*sin(theta))*cos &
      (6.2831853071795865d0*z)/(pi*r**3)

end function

REAL*8 pure function udf_grada_r(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_grada_r = (r - 3)**2*(2*r - 2)*sin(theta)*sin(6.2831853071795865d0*z &
      ) + (r - 1)**2*(2*r - 6)*sin(theta)*sin(6.2831853071795865d0*z)

end function

REAL*8 pure function udf_grada_a(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_grada_a = (r - 3)**2*(r - 1)**2*sin(6.2831853071795865d0*z)*cos( &
      theta)/r

end function

REAL*8 pure function udf_grada_z(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_grada_z = 2*pi*(r - 3)**2*(r - 1)**2*sin(theta)*cos( &
      6.2831853071795865d0*z)

end function

REAL*8 pure function udf_gradz_r(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_gradz_r = (1.0d0/2.0d0)*(r - 3)**2*(2*r - 2)*cos(theta)*cos( &
      6.2831853071795865d0*z)/(pi*r) + (1.0d0/2.0d0)*(r - 1)**2*(2*r - &
      6)*cos(theta)*cos(6.2831853071795865d0*z)/(pi*r) - (-10*r**3*sin( &
      theta)*cos(6.2831853071795865d0*z)/pi + 48*r**2*sin(theta)*cos( &
      6.2831853071795865d0*z)/pi - 66*r*sin(theta)*cos( &
      6.2831853071795865d0*z)/pi + 24*sin(theta)*cos( &
      6.2831853071795865d0*z)/pi)/r - 1.0d0/2.0d0*(r - 3)**2*(r - 1)**2 &
      *cos(theta)*cos(6.2831853071795865d0*z)/(pi*r**2) + (-5.0d0/2.0d0 &
      *r**4*sin(theta)*cos(6.2831853071795865d0*z)/pi + 16*r**3*sin( &
      theta)*cos(6.2831853071795865d0*z)/pi - 33*r**2*sin(theta)*cos( &
      6.2831853071795865d0*z)/pi + 24*r*sin(theta)*cos( &
      6.2831853071795865d0*z)/pi - 9.0d0/2.0d0*sin(theta)*cos( &
      6.2831853071795865d0*z)/pi)/r**2

end function

REAL*8 pure function udf_gradz_a(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_gradz_a = (-1.0d0/2.0d0*(r - 3)**2*(r - 1)**2*sin(theta)*cos( &
      6.2831853071795865d0*z)/(pi*r) - (-5.0d0/2.0d0*r**4*cos(theta)* &
      cos(6.2831853071795865d0*z)/pi + 16*r**3*cos(theta)*cos( &
      6.2831853071795865d0*z)/pi - 33*r**2*cos(theta)*cos( &
      6.2831853071795865d0*z)/pi + 24*r*cos(theta)*cos( &
      6.2831853071795865d0*z)/pi - 9.0d0/2.0d0*cos(theta)*cos( &
      6.2831853071795865d0*z)/pi)/r)/r

end function

REAL*8 pure function udf_gradz_z(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_gradz_z = -(r - 3)**2*(r - 1)**2*sin(6.2831853071795865d0*z)*cos( &
      theta)/r - (5*r**4*sin(theta)*sin(6.2831853071795865d0*z) - 32*r &
      **3*sin(theta)*sin(6.2831853071795865d0*z) + 66*r**2*sin(theta)* &
      sin(6.2831853071795865d0*z) - 48*r*sin(theta)*sin( &
      6.2831853071795865d0*z) + 9*sin(theta)*sin(6.2831853071795865d0*z &
      ))/r

end function

REAL*8 pure function udf_gradr_r(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_gradr_r = (r - 3)**2*(2*r - 2)*sin(theta)*sin(6.2831853071795865d0*z &
      ) + (r - 1)**2*(2*r - 6)*sin(theta)*sin(6.2831853071795865d0*z)

end function

REAL*8 pure function udf_gradr_a(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_gradr_a = (r - 3)**2*(r - 1)**2*sin(6.2831853071795865d0*z)*cos( &
      theta)/r

end function

REAL*8 pure function udf_gradr_z(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_gradr_z = 2*pi*(r - 3)**2*(r - 1)**2*sin(theta)*cos( &
      6.2831853071795865d0*z)

end function

REAL*8 pure function udf_scm_ur(t, theta, z, r, nu)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r
REAL*8, intent(in) :: nu

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_scm_ur = -nu*(-4*pi**2*(r - 3)**2*(r - 1)**2*sin(theta)*sin( &
      6.2831853071795865d0*z) + 2*((r - 3)**2 + 4*(r - 3)*(r - 1) + (r &
      - 1)**2)*sin(theta)*sin(6.2831853071795865d0*z) + 1.0d0*((r - 3) &
      **2*(2*r - 2)*sin(theta)*sin(6.2831853071795865d0*z) + (r - 1)**2 &
      *(2*r - 6)*sin(theta)*sin(6.2831853071795865d0*z))/r - 2.0d0*(r - &
      3)**2*(r - 1)**2*sin(theta)*sin(6.2831853071795865d0*z)/r**2 - &
      2.0d0*(r - 3)**2*(r - 1)**2*sin(6.2831853071795865d0*z)*cos(theta &
      )/r**2) + (r - 3)**2*(r - 1)**2*(4*r*(r - 3)*(r - 2)*(r - 1)*sin( &
      theta)*sin(6.2831853071795865d0*z)**2 - (r - 3)**2*(r - 1)**2*sin &
      (theta)*sin(6.2831853071795865d0*z)**2 + (r - 3)**2*(r - 1)**2* &
      sin(6.2831853071795865d0*z)**2*cos(theta) + (5*r**4*sin(theta) - &
      32*r**3*sin(theta) + 66*r**2*sin(theta) - 48*r*sin(theta) + (r - &
      3)**2*(r - 1)**2*cos(theta) + 9*sin(theta))*cos( &
      6.2831853071795865d0*z)**2)*sin(theta)/r

end function

REAL*8 pure function udf_scm_ua(t, theta, z, r, nu)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r
REAL*8, intent(in) :: nu

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_scm_ua = -nu*(-4*pi**2*(r - 3)**2*(r - 1)**2*sin(theta)*sin( &
      6.2831853071795865d0*z) + 2*((r - 3)**2 + 4*(r - 3)*(r - 1) + (r &
      - 1)**2)*sin(theta)*sin(6.2831853071795865d0*z) + 1.0d0*((r - 3) &
      **2*(2*r - 2)*sin(theta)*sin(6.2831853071795865d0*z) + (r - 1)**2 &
      *(2*r - 6)*sin(theta)*sin(6.2831853071795865d0*z))/r - 2.0d0*(r - &
      3)**2*(r - 1)**2*sin(theta)*sin(6.2831853071795865d0*z)/r**2 + &
      2.0d0*(r - 3)**2*(r - 1)**2*sin(6.2831853071795865d0*z)*cos(theta &
      )/r**2) + (r - 3)**2*(r - 1)**2*(5*r**4*sin(theta) + r**4*cos( &
      theta) - 32*r**3*sin(theta) - 8*r**3*cos(theta) + 63*r**2*sin( &
      theta) + 3*sqrt(2.0d0)*r**2*sin(theta + 0.78539816339744831d0) + &
      19*r**2*cos(theta) - 36*r*sin(theta) - 12*sqrt(2.0d0)*r*sin(theta &
      + 0.78539816339744831d0) - 12*r*cos(theta) + 9*sqrt(2.0d0)*sin( &
      theta + 0.78539816339744831d0))*sin(theta)/r

end function

REAL*8 pure function udf_scm_uz(t, theta, z, r, nu)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r
REAL*8, intent(in) :: nu

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_scm_uz = -nu*(-2*pi*((r - 3)**2*(r - 1)**2*cos(theta) + (5*r**4 - 32 &
      *r**3 + 66*r**2 - 48*r + 9)*sin(theta))*cos(6.2831853071795865d0* &
      z)/r + 1.0d0*((1.0d0/2.0d0)*(r - 3)**2*(2*r - 2)*cos(theta)*cos( &
      6.2831853071795865d0*z)/(pi*r) + (1.0d0/2.0d0)*(r - 1)**2*(2*r - &
      6)*cos(theta)*cos(6.2831853071795865d0*z)/(pi*r) - (-10*r**3*sin( &
      theta)*cos(6.2831853071795865d0*z)/pi + 48*r**2*sin(theta)*cos( &
      6.2831853071795865d0*z)/pi - 66*r*sin(theta)*cos( &
      6.2831853071795865d0*z)/pi + 24*sin(theta)*cos( &
      6.2831853071795865d0*z)/pi)/r - 1.0d0/2.0d0*(r - 3)**2*(r - 1)**2 &
      *cos(theta)*cos(6.2831853071795865d0*z)/(pi*r**2) + (-5.0d0/2.0d0 &
      *r**4*sin(theta)*cos(6.2831853071795865d0*z)/pi + 16*r**3*sin( &
      theta)*cos(6.2831853071795865d0*z)/pi - 33*r**2*sin(theta)*cos( &
      6.2831853071795865d0*z)/pi + 24*r*sin(theta)*cos( &
      6.2831853071795865d0*z)/pi - 9.0d0/2.0d0*sin(theta)*cos( &
      6.2831853071795865d0*z)/pi)/r**2)/r + ((r - 3)**2*cos(theta) + 4* &
      (r - 3)*(r - 1)*cos(theta) + (r - 1)**2*cos(theta) + 6*(5*r**2 - &
      16*r + 11)*sin(theta) - 2*(r - 3)**2*(r - 1)*cos(theta)/r - 2*(r &
      - 3)*(r - 1)**2*cos(theta)/r - 4*(5*r**3 - 24*r**2 + 33*r - 12)* &
      sin(theta)/r + (r - 3)**2*(r - 1)**2*cos(theta)/r**2 + (5*r**4 - &
      32*r**3 + 66*r**2 - 48*r + 9)*sin(theta)/r**2)*cos( &
      6.2831853071795865d0*z)/(pi*r) - 0.5d0*((r - 3)**2*(r - 1)**2*cos &
      (theta) + (5*r**4 - 32*r**3 + 66*r**2 - 48*r + 9)*sin(theta))*cos &
      (6.2831853071795865d0*z)/(pi*r**3)) + (1.0d0/2.0d0)*(r - 3)**2*(r &
      - 1)**2*(-r**4*sin(2.0d0*theta) + 5*r**4*cos(2.0d0*theta) - 6*r** &
      4 + 8*r**3*sin(2.0d0*theta) - 28*r**3*cos(2.0d0*theta) + 36*r**3 &
      - 22*r**2*sin(2.0d0*theta) + 54*r**2*cos(2.0d0*theta) - 76*r**2 + &
      24*r*sin(2.0d0*theta) - 36*r*cos(2.0d0*theta) + 60*r - 9*sin( &
      2.0d0*theta) + 9*cos(2.0d0*theta) - 18)*sin(6.2831853071795865d0* &
      z)*cos(6.2831853071795865d0*z)/(pi*r**2)

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
