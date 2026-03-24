
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

udf_ua = (r - 3)**2*(r - 1)**2*sin(4.0d0*theta)*sin(31.415926535897932d0 &
      *z)*cos(3.1415926535897932d0*t)

end function

REAL*8 pure function udf_uz(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_uz = (2.0d0/5.0d0)*(r - 3)**2*(r - 1)**2*cos(4.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(31.415926535897932d0*z)/(pi*r) - ( &
      -5.0d0/2.0d0*r**4*sin(8.0d0*theta)*cos(3.1415926535897932d0*t)* &
      cos(6.2831853071795865d0*z)/pi + 16*r**3*sin(8.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(6.2831853071795865d0*z)/pi - 33*r**2* &
      sin(8.0d0*theta)*cos(3.1415926535897932d0*t)*cos( &
      6.2831853071795865d0*z)/pi + 24*r*sin(8.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(6.2831853071795865d0*z)/pi - 9.0d0/ &
      2.0d0*sin(8.0d0*theta)*cos(3.1415926535897932d0*t)*cos( &
      6.2831853071795865d0*z)/pi)/r

end function

REAL*8 pure function udf_ur(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_ur = (r - 3)**2*(r - 1)**2*sin(8.0d0*theta)*sin(6.2831853071795865d0 &
      *z)*cos(3.1415926535897932d0*t)

end function

REAL*8 pure function udf_p(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_p = (r - 3)**2*(r - 1)**2*cos(theta)*cos(6.2831853071795865d0*z)

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

udf_div = (r*(r - 3)**2*(2*r - 2)*sin(8.0d0*theta)*sin( &
      6.2831853071795865d0*z)*cos(3.1415926535897932d0*t) + r*(r - 1)** &
      2*(2*r - 6)*sin(8.0d0*theta)*sin(6.2831853071795865d0*z)*cos( &
      3.1415926535897932d0*t) + (r - 3)**2*(r - 1)**2*sin(8.0d0*theta)* &
      sin(6.2831853071795865d0*z)*cos(3.1415926535897932d0*t))/r - (5*r &
      **4*sin(8.0d0*theta)*sin(6.2831853071795865d0*z)*cos( &
      3.1415926535897932d0*t) - 32*r**3*sin(8.0d0*theta)*sin( &
      6.2831853071795865d0*z)*cos(3.1415926535897932d0*t) + 66*r**2*sin &
      (8.0d0*theta)*sin(6.2831853071795865d0*z)*cos( &
      3.1415926535897932d0*t) - 48*r*sin(8.0d0*theta)*sin( &
      6.2831853071795865d0*z)*cos(3.1415926535897932d0*t) + 9*sin(8.0d0 &
      *theta)*sin(6.2831853071795865d0*z)*cos(3.1415926535897932d0*t))/ &
      r

end function

REAL*8 pure function udf_NLr(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_NLr = (r - 3)**2*(r - 1)**2*(4*r*(r - 3)*(r - 2)*(r - 1)*sin(8.0d0* &
      theta)**2*sin(6.2831853071795865d0*z)**2 - (r - 3)**2*(r - 1)**2* &
      sin(4.0d0*theta)**2*sin(31.415926535897932d0*z)**2 + 8*(r - 3)**2 &
      *(r - 1)**2*sin(4.0d0*theta)*sin(6.2831853071795865d0*z)*sin( &
      31.415926535897932d0*z)*cos(8.0d0*theta) + (1.0d0/5.0d0)*(25*r**4 &
      *sin(8.0d0*theta)*cos(6.2831853071795865d0*z) - 160*r**3*sin( &
      8.0d0*theta)*cos(6.2831853071795865d0*z) + 330*r**2*sin(8.0d0* &
      theta)*cos(6.2831853071795865d0*z) - 240*r*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) + 4*(r - 3)**2*(r - 1)**2*cos(4.0d0*theta &
      )*cos(31.415926535897932d0*z) + 45*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z))*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z))*cos(3.1415926535897932d0*t)**2/r

end function

REAL*8 pure function udf_NLa(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_NLa = (r - 3)**2*(r - 1)**2*(4*r*(r - 3)*(r - 2)*(r - 1)*sin(8.0d0* &
      theta)*sin(6.2831853071795865d0*z)*sin(31.415926535897932d0*z) + &
      (r - 3)**2*(r - 1)**2*sin(8.0d0*theta)*sin(6.2831853071795865d0*z &
      )*sin(31.415926535897932d0*z) + 4*(r - 3)**2*(r - 1)**2*sin( &
      31.415926535897932d0*z)**2*cos(4.0d0*theta) + (25*r**4*sin(8.0d0* &
      theta)*cos(6.2831853071795865d0*z) - 160*r**3*sin(8.0d0*theta)* &
      cos(6.2831853071795865d0*z) + 330*r**2*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) - 240*r*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) + 4*(r - 3)**2*(r - 1)**2*cos(4.0d0*theta &
      )*cos(31.415926535897932d0*z) + 45*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z))*cos(31.415926535897932d0*z))*sin(4.0d0* &
      theta)*cos(3.1415926535897932d0*t)**2/r

end function

REAL*8 pure function udf_NLz(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_NLz = (1.0d0/10.0d0)*(-8*(r - 3)**2*(r - 1)**2*(-25*r**4*cos(8.0d0* &
      theta)*cos(6.2831853071795865d0*z) + 160*r**3*cos(8.0d0*theta)* &
      cos(6.2831853071795865d0*z) - 330*r**2*cos(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) + 240*r*cos(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) + 2*(r - 3)**2*(r - 1)**2*sin(4.0d0*theta &
      )*cos(31.415926535897932d0*z) - 45*cos(8.0d0*theta)*cos( &
      6.2831853071795865d0*z))*sin(4.0d0*theta)*sin( &
      31.415926535897932d0*z) + (r - 3)**2*(r - 1)**2*(-25*r**4*sin( &
      8.0d0*theta)*cos(6.2831853071795865d0*z) + 160*r**3*sin(8.0d0* &
      theta)*cos(6.2831853071795865d0*z) - 330*r**2*sin(8.0d0*theta)* &
      cos(6.2831853071795865d0*z) + 4*r*(25*r**3*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) - 120*r**2*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) + 165*r*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) + 2*(r - 3)**2*(r - 1)*cos(4.0d0*theta)* &
      cos(31.415926535897932d0*z) + 2*(r - 3)*(r - 1)**2*cos(4.0d0* &
      theta)*cos(31.415926535897932d0*z) - 60*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z)) + 240*r*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) - 4*(r - 3)**2*(r - 1)**2*cos(4.0d0*theta &
      )*cos(31.415926535897932d0*z) - 45*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z))*sin(8.0d0*theta)*sin( &
      6.2831853071795865d0*z) - (5*r**4*sin(8.0d0*theta)*sin( &
      6.2831853071795865d0*z) - 32*r**3*sin(8.0d0*theta)*sin( &
      6.2831853071795865d0*z) + 66*r**2*sin(8.0d0*theta)*sin( &
      6.2831853071795865d0*z) - 48*r*sin(8.0d0*theta)*sin( &
      6.2831853071795865d0*z) + 4*(r - 3)**2*(r - 1)**2*sin( &
      31.415926535897932d0*z)*cos(4.0d0*theta) + 9*sin(8.0d0*theta)*sin &
      (6.2831853071795865d0*z))*(25*r**4*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) - 160*r**3*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) + 330*r**2*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) - 240*r*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) + 4*(r - 3)**2*(r - 1)**2*cos(4.0d0*theta &
      )*cos(31.415926535897932d0*z) + 45*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z)))*cos(3.1415926535897932d0*t)**2/(pi*r**2 &
      )

end function

REAL*8 pure function udf_Lr(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_Lr = -4*pi**2*(r - 3)**2*(r - 1)**2*sin(8.0d0*theta)*sin( &
      6.2831853071795865d0*z)*cos(3.1415926535897932d0*t) + 2*((r - 3) &
      **2 + 4*(r - 3)*(r - 1) + (r - 1)**2)*sin(8.0d0*theta)*sin( &
      6.2831853071795865d0*z)*cos(3.1415926535897932d0*t) + 1.0d0*((r - &
      3)**2*(2*r - 2)*sin(8.0d0*theta)*sin(6.2831853071795865d0*z)*cos( &
      3.1415926535897932d0*t) + (r - 1)**2*(2*r - 6)*sin(8.0d0*theta)* &
      sin(6.2831853071795865d0*z)*cos(3.1415926535897932d0*t))/r - &
      65.0d0*(r - 3)**2*(r - 1)**2*sin(8.0d0*theta)*sin( &
      6.2831853071795865d0*z)*cos(3.1415926535897932d0*t)/r**2 - 8.0d0* &
      (r - 3)**2*(r - 1)**2*sin(31.415926535897932d0*z)*cos(4.0d0*theta &
      )*cos(3.1415926535897932d0*t)/r**2

end function

REAL*8 pure function udf_La(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_La = -100*pi**2*(r - 3)**2*(r - 1)**2*sin(4.0d0*theta)*sin( &
      31.415926535897932d0*z)*cos(3.1415926535897932d0*t) + 2*((r - 3) &
      **2 + 4*(r - 3)*(r - 1) + (r - 1)**2)*sin(4.0d0*theta)*sin( &
      31.415926535897932d0*z)*cos(3.1415926535897932d0*t) + 1.0d0*((r - &
      3)**2*(2*r - 2)*sin(4.0d0*theta)*sin(31.415926535897932d0*z)*cos( &
      3.1415926535897932d0*t) + (r - 1)**2*(2*r - 6)*sin(4.0d0*theta)* &
      sin(31.415926535897932d0*z)*cos(3.1415926535897932d0*t))/r - &
      17.0d0*(r - 3)**2*(r - 1)**2*sin(4.0d0*theta)*sin( &
      31.415926535897932d0*z)*cos(3.1415926535897932d0*t)/r**2 + 16.0d0 &
      *(r - 3)**2*(r - 1)**2*sin(6.2831853071795865d0*z)*cos(8.0d0* &
      theta)*cos(3.1415926535897932d0*t)/r**2

end function

REAL*8 pure function udf_Lz(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_Lz = -2*pi*(20*(r - 3)**2*(r - 1)**2*cos(4.0d0*theta)*cos( &
      31.415926535897932d0*z) + (5*r**4 - 32*r**3 + 66*r**2 - 48*r + 9) &
      *sin(8.0d0*theta)*cos(6.2831853071795865d0*z))*cos( &
      3.1415926535897932d0*t)/r + 1.0d0*((2.0d0/5.0d0)*(r - 3)**2*(2*r &
      - 2)*cos(4.0d0*theta)*cos(3.1415926535897932d0*t)*cos( &
      31.415926535897932d0*z)/(pi*r) + (2.0d0/5.0d0)*(r - 1)**2*(2*r - &
      6)*cos(4.0d0*theta)*cos(3.1415926535897932d0*t)*cos( &
      31.415926535897932d0*z)/(pi*r) - (-10*r**3*sin(8.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(6.2831853071795865d0*z)/pi + 48*r**2* &
      sin(8.0d0*theta)*cos(3.1415926535897932d0*t)*cos( &
      6.2831853071795865d0*z)/pi - 66*r*sin(8.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(6.2831853071795865d0*z)/pi + 24*sin( &
      8.0d0*theta)*cos(3.1415926535897932d0*t)*cos(6.2831853071795865d0 &
      *z)/pi)/r - 2.0d0/5.0d0*(r - 3)**2*(r - 1)**2*cos(4.0d0*theta)* &
      cos(3.1415926535897932d0*t)*cos(31.415926535897932d0*z)/(pi*r**2 &
      ) + (-5.0d0/2.0d0*r**4*sin(8.0d0*theta)*cos(3.1415926535897932d0* &
      t)*cos(6.2831853071795865d0*z)/pi + 16*r**3*sin(8.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(6.2831853071795865d0*z)/pi - 33*r**2* &
      sin(8.0d0*theta)*cos(3.1415926535897932d0*t)*cos( &
      6.2831853071795865d0*z)/pi + 24*r*sin(8.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(6.2831853071795865d0*z)/pi - 9.0d0/ &
      2.0d0*sin(8.0d0*theta)*cos(3.1415926535897932d0*t)*cos( &
      6.2831853071795865d0*z)/pi)/r**2)/r + 2*((2.0d0/5.0d0)*(r - 3)**2 &
      *cos(4.0d0*theta)*cos(31.415926535897932d0*z) + (8.0d0/5.0d0)*(r &
      - 3)*(r - 1)*cos(4.0d0*theta)*cos(31.415926535897932d0*z) + ( &
      2.0d0/5.0d0)*(r - 1)**2*cos(4.0d0*theta)*cos(31.415926535897932d0 &
      *z) + 3*(5*r**2 - 16*r + 11)*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) - 4.0d0/5.0d0*(r - 3)**2*(r - 1)*cos( &
      4.0d0*theta)*cos(31.415926535897932d0*z)/r - 4.0d0/5.0d0*(r - 3)* &
      (r - 1)**2*cos(4.0d0*theta)*cos(31.415926535897932d0*z)/r - 2*(5* &
      r**3 - 24*r**2 + 33*r - 12)*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z)/r + (2.0d0/5.0d0)*(r - 3)**2*(r - 1)**2* &
      cos(4.0d0*theta)*cos(31.415926535897932d0*z)/r**2 + (1.0d0/2.0d0) &
      *(5*r**4 - 32*r**3 + 66*r**2 - 48*r + 9)*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z)/r**2)*cos(3.1415926535897932d0*t)/(pi*r) &
      - 32.0d0*((1.0d0/5.0d0)*(r - 3)**2*(r - 1)**2*cos(4.0d0*theta)* &
      cos(31.415926535897932d0*z) + (5*r**4 - 32*r**3 + 66*r**2 - 48*r &
      + 9)*sin(8.0d0*theta)*cos(6.2831853071795865d0*z))*cos( &
      3.1415926535897932d0*t)/(pi*r**3)

end function

REAL*8 pure function udf_grada_r(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_grada_r = (r - 3)**2*(2*r - 2)*sin(4.0d0*theta)*sin( &
      31.415926535897932d0*z)*cos(3.1415926535897932d0*t) + (r - 1)**2* &
      (2*r - 6)*sin(4.0d0*theta)*sin(31.415926535897932d0*z)*cos( &
      3.1415926535897932d0*t)

end function

REAL*8 pure function udf_grada_a(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_grada_a = 4*(r - 3)**2*(r - 1)**2*sin(31.415926535897932d0*z)*cos( &
      4.0d0*theta)*cos(3.1415926535897932d0*t)/r

end function

REAL*8 pure function udf_grada_z(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_grada_z = 10*pi*(r - 3)**2*(r - 1)**2*sin(4.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(31.415926535897932d0*z)

end function

REAL*8 pure function udf_gradz_r(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_gradz_r = (2.0d0/5.0d0)*(r - 3)**2*(2*r - 2)*cos(4.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(31.415926535897932d0*z)/(pi*r) + ( &
      2.0d0/5.0d0)*(r - 1)**2*(2*r - 6)*cos(4.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(31.415926535897932d0*z)/(pi*r) - (-10 &
      *r**3*sin(8.0d0*theta)*cos(3.1415926535897932d0*t)*cos( &
      6.2831853071795865d0*z)/pi + 48*r**2*sin(8.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(6.2831853071795865d0*z)/pi - 66*r*sin &
      (8.0d0*theta)*cos(3.1415926535897932d0*t)*cos( &
      6.2831853071795865d0*z)/pi + 24*sin(8.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(6.2831853071795865d0*z)/pi)/r - 2.0d0 &
      /5.0d0*(r - 3)**2*(r - 1)**2*cos(4.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(31.415926535897932d0*z)/(pi*r**2) + ( &
      -5.0d0/2.0d0*r**4*sin(8.0d0*theta)*cos(3.1415926535897932d0*t)* &
      cos(6.2831853071795865d0*z)/pi + 16*r**3*sin(8.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(6.2831853071795865d0*z)/pi - 33*r**2* &
      sin(8.0d0*theta)*cos(3.1415926535897932d0*t)*cos( &
      6.2831853071795865d0*z)/pi + 24*r*sin(8.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(6.2831853071795865d0*z)/pi - 9.0d0/ &
      2.0d0*sin(8.0d0*theta)*cos(3.1415926535897932d0*t)*cos( &
      6.2831853071795865d0*z)/pi)/r**2

end function

REAL*8 pure function udf_gradz_a(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_gradz_a = (-8.0d0/5.0d0*(r - 3)**2*(r - 1)**2*sin(4.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(31.415926535897932d0*z)/(pi*r) - (-20 &
      *r**4*cos(8.0d0*theta)*cos(3.1415926535897932d0*t)*cos( &
      6.2831853071795865d0*z)/pi + 128*r**3*cos(8.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(6.2831853071795865d0*z)/pi - 264*r**2 &
      *cos(8.0d0*theta)*cos(3.1415926535897932d0*t)*cos( &
      6.2831853071795865d0*z)/pi + 192*r*cos(8.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(6.2831853071795865d0*z)/pi - 36*cos( &
      8.0d0*theta)*cos(3.1415926535897932d0*t)*cos(6.2831853071795865d0 &
      *z)/pi)/r)/r

end function

REAL*8 pure function udf_gradz_z(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_gradz_z = -4*(r - 3)**2*(r - 1)**2*sin(31.415926535897932d0*z)*cos( &
      4.0d0*theta)*cos(3.1415926535897932d0*t)/r - (5*r**4*sin(8.0d0* &
      theta)*sin(6.2831853071795865d0*z)*cos(3.1415926535897932d0*t) - &
      32*r**3*sin(8.0d0*theta)*sin(6.2831853071795865d0*z)*cos( &
      3.1415926535897932d0*t) + 66*r**2*sin(8.0d0*theta)*sin( &
      6.2831853071795865d0*z)*cos(3.1415926535897932d0*t) - 48*r*sin( &
      8.0d0*theta)*sin(6.2831853071795865d0*z)*cos(3.1415926535897932d0 &
      *t) + 9*sin(8.0d0*theta)*sin(6.2831853071795865d0*z)*cos( &
      3.1415926535897932d0*t))/r

end function

REAL*8 pure function udf_gradr_r(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_gradr_r = (r - 3)**2*(2*r - 2)*sin(8.0d0*theta)*sin( &
      6.2831853071795865d0*z)*cos(3.1415926535897932d0*t) + (r - 1)**2* &
      (2*r - 6)*sin(8.0d0*theta)*sin(6.2831853071795865d0*z)*cos( &
      3.1415926535897932d0*t)

end function

REAL*8 pure function udf_gradr_a(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_gradr_a = 8*(r - 3)**2*(r - 1)**2*sin(6.2831853071795865d0*z)*cos( &
      8.0d0*theta)*cos(3.1415926535897932d0*t)/r

end function

REAL*8 pure function udf_gradr_z(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_gradr_z = 2*pi*(r - 3)**2*(r - 1)**2*sin(8.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(6.2831853071795865d0*z)

end function

REAL*8 pure function udf_scm_ur(t, theta, z, r, nu)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r
REAL*8, intent(in) :: nu

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_scm_ur = -nu*(-4*pi**2*(r - 3)**2*(r - 1)**2*sin(8.0d0*theta)*sin( &
      6.2831853071795865d0*z)*cos(3.1415926535897932d0*t) + 2*((r - 3) &
      **2 + 4*(r - 3)*(r - 1) + (r - 1)**2)*sin(8.0d0*theta)*sin( &
      6.2831853071795865d0*z)*cos(3.1415926535897932d0*t) + 1.0d0*((r - &
      3)**2*(2*r - 2)*sin(8.0d0*theta)*sin(6.2831853071795865d0*z)*cos( &
      3.1415926535897932d0*t) + (r - 1)**2*(2*r - 6)*sin(8.0d0*theta)* &
      sin(6.2831853071795865d0*z)*cos(3.1415926535897932d0*t))/r - &
      65.0d0*(r - 3)**2*(r - 1)**2*sin(8.0d0*theta)*sin( &
      6.2831853071795865d0*z)*cos(3.1415926535897932d0*t)/r**2 - 8.0d0* &
      (r - 3)**2*(r - 1)**2*sin(31.415926535897932d0*z)*cos(4.0d0*theta &
      )*cos(3.1415926535897932d0*t)/r**2) - pi*(r - 3)**2*(r - 1)**2* &
      sin(8.0d0*theta)*sin(3.1415926535897932d0*t)*sin( &
      6.2831853071795865d0*z) + (r - 3)**2*(2*r - 2)*cos(theta)*cos( &
      6.2831853071795865d0*z) + (r - 1)**2*(2*r - 6)*cos(theta)*cos( &
      6.2831853071795865d0*z) + (r - 3)**2*(r - 1)**2*(4*r*(r - 3)*(r - &
      2)*(r - 1)*sin(8.0d0*theta)**2*sin(6.2831853071795865d0*z)**2 - ( &
      r - 3)**2*(r - 1)**2*sin(4.0d0*theta)**2*sin(31.415926535897932d0 &
      *z)**2 + 8*(r - 3)**2*(r - 1)**2*sin(4.0d0*theta)*sin( &
      6.2831853071795865d0*z)*sin(31.415926535897932d0*z)*cos(8.0d0* &
      theta) + (1.0d0/5.0d0)*(25*r**4*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) - 160*r**3*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) + 330*r**2*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) - 240*r*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) + 4*(r - 3)**2*(r - 1)**2*cos(4.0d0*theta &
      )*cos(31.415926535897932d0*z) + 45*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z))*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z))*cos(3.1415926535897932d0*t)**2/r

end function

REAL*8 pure function udf_scm_ua(t, theta, z, r, nu)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r
REAL*8, intent(in) :: nu

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_scm_ua = -nu*(-100*pi**2*(r - 3)**2*(r - 1)**2*sin(4.0d0*theta)*sin( &
      31.415926535897932d0*z)*cos(3.1415926535897932d0*t) + 2*((r - 3) &
      **2 + 4*(r - 3)*(r - 1) + (r - 1)**2)*sin(4.0d0*theta)*sin( &
      31.415926535897932d0*z)*cos(3.1415926535897932d0*t) + 1.0d0*((r - &
      3)**2*(2*r - 2)*sin(4.0d0*theta)*sin(31.415926535897932d0*z)*cos( &
      3.1415926535897932d0*t) + (r - 1)**2*(2*r - 6)*sin(4.0d0*theta)* &
      sin(31.415926535897932d0*z)*cos(3.1415926535897932d0*t))/r - &
      17.0d0*(r - 3)**2*(r - 1)**2*sin(4.0d0*theta)*sin( &
      31.415926535897932d0*z)*cos(3.1415926535897932d0*t)/r**2 + 16.0d0 &
      *(r - 3)**2*(r - 1)**2*sin(6.2831853071795865d0*z)*cos(8.0d0* &
      theta)*cos(3.1415926535897932d0*t)/r**2) - pi*(r - 3)**2*(r - 1) &
      **2*sin(4.0d0*theta)*sin(3.1415926535897932d0*t)*sin( &
      31.415926535897932d0*z) + (r - 3)**2*(r - 1)**2*(4*r*(r - 3)*(r - &
      2)*(r - 1)*sin(8.0d0*theta)*sin(6.2831853071795865d0*z)*sin( &
      31.415926535897932d0*z) + (r - 3)**2*(r - 1)**2*sin(8.0d0*theta)* &
      sin(6.2831853071795865d0*z)*sin(31.415926535897932d0*z) + 4*(r - &
      3)**2*(r - 1)**2*sin(31.415926535897932d0*z)**2*cos(4.0d0*theta) &
      + (25*r**4*sin(8.0d0*theta)*cos(6.2831853071795865d0*z) - 160*r** &
      3*sin(8.0d0*theta)*cos(6.2831853071795865d0*z) + 330*r**2*sin( &
      8.0d0*theta)*cos(6.2831853071795865d0*z) - 240*r*sin(8.0d0*theta) &
      *cos(6.2831853071795865d0*z) + 4*(r - 3)**2*(r - 1)**2*cos(4.0d0* &
      theta)*cos(31.415926535897932d0*z) + 45*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z))*cos(31.415926535897932d0*z))*sin(4.0d0* &
      theta)*cos(3.1415926535897932d0*t)**2/r - (r - 3)**2*(r - 1)**2* &
      sin(theta)*cos(6.2831853071795865d0*z)/r

end function

REAL*8 pure function udf_scm_uz(t, theta, z, r, nu)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r
REAL*8, intent(in) :: nu

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_scm_uz = -nu*(-2*pi*(20*(r - 3)**2*(r - 1)**2*cos(4.0d0*theta)*cos( &
      31.415926535897932d0*z) + (5*r**4 - 32*r**3 + 66*r**2 - 48*r + 9) &
      *sin(8.0d0*theta)*cos(6.2831853071795865d0*z))*cos( &
      3.1415926535897932d0*t)/r + 1.0d0*((2.0d0/5.0d0)*(r - 3)**2*(2*r &
      - 2)*cos(4.0d0*theta)*cos(3.1415926535897932d0*t)*cos( &
      31.415926535897932d0*z)/(pi*r) + (2.0d0/5.0d0)*(r - 1)**2*(2*r - &
      6)*cos(4.0d0*theta)*cos(3.1415926535897932d0*t)*cos( &
      31.415926535897932d0*z)/(pi*r) - (-10*r**3*sin(8.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(6.2831853071795865d0*z)/pi + 48*r**2* &
      sin(8.0d0*theta)*cos(3.1415926535897932d0*t)*cos( &
      6.2831853071795865d0*z)/pi - 66*r*sin(8.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(6.2831853071795865d0*z)/pi + 24*sin( &
      8.0d0*theta)*cos(3.1415926535897932d0*t)*cos(6.2831853071795865d0 &
      *z)/pi)/r - 2.0d0/5.0d0*(r - 3)**2*(r - 1)**2*cos(4.0d0*theta)* &
      cos(3.1415926535897932d0*t)*cos(31.415926535897932d0*z)/(pi*r**2 &
      ) + (-5.0d0/2.0d0*r**4*sin(8.0d0*theta)*cos(3.1415926535897932d0* &
      t)*cos(6.2831853071795865d0*z)/pi + 16*r**3*sin(8.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(6.2831853071795865d0*z)/pi - 33*r**2* &
      sin(8.0d0*theta)*cos(3.1415926535897932d0*t)*cos( &
      6.2831853071795865d0*z)/pi + 24*r*sin(8.0d0*theta)*cos( &
      3.1415926535897932d0*t)*cos(6.2831853071795865d0*z)/pi - 9.0d0/ &
      2.0d0*sin(8.0d0*theta)*cos(3.1415926535897932d0*t)*cos( &
      6.2831853071795865d0*z)/pi)/r**2)/r + 2*((2.0d0/5.0d0)*(r - 3)**2 &
      *cos(4.0d0*theta)*cos(31.415926535897932d0*z) + (8.0d0/5.0d0)*(r &
      - 3)*(r - 1)*cos(4.0d0*theta)*cos(31.415926535897932d0*z) + ( &
      2.0d0/5.0d0)*(r - 1)**2*cos(4.0d0*theta)*cos(31.415926535897932d0 &
      *z) + 3*(5*r**2 - 16*r + 11)*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) - 4.0d0/5.0d0*(r - 3)**2*(r - 1)*cos( &
      4.0d0*theta)*cos(31.415926535897932d0*z)/r - 4.0d0/5.0d0*(r - 3)* &
      (r - 1)**2*cos(4.0d0*theta)*cos(31.415926535897932d0*z)/r - 2*(5* &
      r**3 - 24*r**2 + 33*r - 12)*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z)/r + (2.0d0/5.0d0)*(r - 3)**2*(r - 1)**2* &
      cos(4.0d0*theta)*cos(31.415926535897932d0*z)/r**2 + (1.0d0/2.0d0) &
      *(5*r**4 - 32*r**3 + 66*r**2 - 48*r + 9)*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z)/r**2)*cos(3.1415926535897932d0*t)/(pi*r) &
      - 32.0d0*((1.0d0/5.0d0)*(r - 3)**2*(r - 1)**2*cos(4.0d0*theta)* &
      cos(31.415926535897932d0*z) + (5*r**4 - 32*r**3 + 66*r**2 - 48*r &
      + 9)*sin(8.0d0*theta)*cos(6.2831853071795865d0*z))*cos( &
      3.1415926535897932d0*t)/(pi*r**3)) - 2*pi*(r - 3)**2*(r - 1)**2* &
      sin(6.2831853071795865d0*z)*cos(theta) - 2.0d0/5.0d0*(r - 3)**2*( &
      r - 1)**2*sin(3.1415926535897932d0*t)*cos(4.0d0*theta)*cos( &
      31.415926535897932d0*z)/r - ((5.0d0/2.0d0)*r**4*sin(8.0d0*theta)* &
      sin(3.1415926535897932d0*t)*cos(6.2831853071795865d0*z) - 16*r**3 &
      *sin(8.0d0*theta)*sin(3.1415926535897932d0*t)*cos( &
      6.2831853071795865d0*z) + 33*r**2*sin(8.0d0*theta)*sin( &
      3.1415926535897932d0*t)*cos(6.2831853071795865d0*z) - 24*r*sin( &
      8.0d0*theta)*sin(3.1415926535897932d0*t)*cos(6.2831853071795865d0 &
      *z) + (9.0d0/2.0d0)*sin(8.0d0*theta)*sin(3.1415926535897932d0*t)* &
      cos(6.2831853071795865d0*z))/r + (1.0d0/10.0d0)*(-8*(r - 3)**2*(r &
      - 1)**2*(-25*r**4*cos(8.0d0*theta)*cos(6.2831853071795865d0*z) + &
      160*r**3*cos(8.0d0*theta)*cos(6.2831853071795865d0*z) - 330*r**2* &
      cos(8.0d0*theta)*cos(6.2831853071795865d0*z) + 240*r*cos(8.0d0* &
      theta)*cos(6.2831853071795865d0*z) + 2*(r - 3)**2*(r - 1)**2*sin( &
      4.0d0*theta)*cos(31.415926535897932d0*z) - 45*cos(8.0d0*theta)* &
      cos(6.2831853071795865d0*z))*sin(4.0d0*theta)*sin( &
      31.415926535897932d0*z) + (r - 3)**2*(r - 1)**2*(-25*r**4*sin( &
      8.0d0*theta)*cos(6.2831853071795865d0*z) + 160*r**3*sin(8.0d0* &
      theta)*cos(6.2831853071795865d0*z) - 330*r**2*sin(8.0d0*theta)* &
      cos(6.2831853071795865d0*z) + 4*r*(25*r**3*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) - 120*r**2*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) + 165*r*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) + 2*(r - 3)**2*(r - 1)*cos(4.0d0*theta)* &
      cos(31.415926535897932d0*z) + 2*(r - 3)*(r - 1)**2*cos(4.0d0* &
      theta)*cos(31.415926535897932d0*z) - 60*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z)) + 240*r*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) - 4*(r - 3)**2*(r - 1)**2*cos(4.0d0*theta &
      )*cos(31.415926535897932d0*z) - 45*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z))*sin(8.0d0*theta)*sin( &
      6.2831853071795865d0*z) - (5*r**4*sin(8.0d0*theta)*sin( &
      6.2831853071795865d0*z) - 32*r**3*sin(8.0d0*theta)*sin( &
      6.2831853071795865d0*z) + 66*r**2*sin(8.0d0*theta)*sin( &
      6.2831853071795865d0*z) - 48*r*sin(8.0d0*theta)*sin( &
      6.2831853071795865d0*z) + 4*(r - 3)**2*(r - 1)**2*sin( &
      31.415926535897932d0*z)*cos(4.0d0*theta) + 9*sin(8.0d0*theta)*sin &
      (6.2831853071795865d0*z))*(25*r**4*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) - 160*r**3*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) + 330*r**2*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) - 240*r*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z) + 4*(r - 3)**2*(r - 1)**2*cos(4.0d0*theta &
      )*cos(31.415926535897932d0*z) + 45*sin(8.0d0*theta)*cos( &
      6.2831853071795865d0*z)))*cos(3.1415926535897932d0*t)**2/(pi*r**2 &
      )

end function

REAL*8 pure function udf_dpdr(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_dpdr = (r - 3)**2*(2*r - 2)*cos(theta)*cos(6.2831853071795865d0*z) + &
      (r - 1)**2*(2*r - 6)*cos(theta)*cos(6.2831853071795865d0*z)

end function

REAL*8 pure function udf_dpda(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_dpda = -(r - 3)**2*(r - 1)**2*sin(theta)*cos(6.2831853071795865d0*z) &
      /r

end function

REAL*8 pure function udf_dpdz(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_dpdz = -2*pi*(r - 3)**2*(r - 1)**2*sin(6.2831853071795865d0*z)*cos( &
      theta)

end function

end module m_udf_mms
