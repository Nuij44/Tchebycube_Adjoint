
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

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_ua = -1.0d0/2.0d0*(sin(6.2831853071795865d0*t) + 1)*sin(theta)*sin( &
      6.2831853071795865d0*r)**2*sin(6.2831853071795865d0*z)/pi

end function

REAL*8 pure function udf_uz(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_uz = -1.0d0/4.0d0*(sin(6.2831853071795865d0*t) + 1)*sin( &
      6.2831853071795865d0*r)**2*cos(theta)*cos(6.2831853071795865d0*z) &
      /(pi**2*r) - (-r*sin(6.2831853071795865d0*r)*sin( &
      6.2831853071795865d0*t)*cos(theta)*cos(6.2831853071795865d0*r)* &
      cos(6.2831853071795865d0*z)/pi - r*sin(6.2831853071795865d0*r)* &
      cos(theta)*cos(6.2831853071795865d0*r)*cos(6.2831853071795865d0*z &
      )/pi - 1.0d0/4.0d0*sin(6.2831853071795865d0*r)**2*sin( &
      6.2831853071795865d0*t)*cos(theta)*cos(6.2831853071795865d0*z)/pi &
      **2 - 1.0d0/4.0d0*sin(6.2831853071795865d0*r)**2*cos(theta)*cos( &
      6.2831853071795865d0*z)/pi**2)/r

end function

REAL*8 pure function udf_ur(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_ur = (1.0d0/2.0d0)*(sin(6.2831853071795865d0*t) + 1)*sin( &
      6.2831853071795865d0*r)**2*sin(6.2831853071795865d0*z)*cos(theta) &
      /pi

end function

REAL*8 pure function udf_p(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_p = cos(6.2831853071795865d0*r)

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

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_div = (2*r*(sin(6.2831853071795865d0*t) + 1)*sin( &
      6.2831853071795865d0*r)*sin(6.2831853071795865d0*z)*cos(theta)* &
      cos(6.2831853071795865d0*r) + (1.0d0/2.0d0)*(sin( &
      6.2831853071795865d0*t) + 1)*sin(6.2831853071795865d0*r)**2*sin( &
      6.2831853071795865d0*z)*cos(theta)/pi)/r - (2*r*sin( &
      6.2831853071795865d0*r)*sin(6.2831853071795865d0*t)*sin( &
      6.2831853071795865d0*z)*cos(theta)*cos(6.2831853071795865d0*r) + &
      2*r*sin(6.2831853071795865d0*r)*sin(6.2831853071795865d0*z)*cos( &
      theta)*cos(6.2831853071795865d0*r) + (1.0d0/2.0d0)*sin( &
      6.2831853071795865d0*r)**2*sin(6.2831853071795865d0*t)*sin( &
      6.2831853071795865d0*z)*cos(theta)/pi + (1.0d0/2.0d0)*sin( &
      6.2831853071795865d0*r)**2*sin(6.2831853071795865d0*z)*cos(theta) &
      /pi)/r

end function

REAL*8 pure function udf_NLr(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_NLr = (sin(6.2831853071795865d0*t) + 1)**2*sin(6.2831853071795865d0* &
      r)**3*cos(theta)**2*cos(6.2831853071795865d0*r)/pi

end function

REAL*8 pure function udf_NLa(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_NLa = -(sin(6.2831853071795865d0*t) + 1)**2*sin(theta)*sin( &
      6.2831853071795865d0*r)**3*cos(theta)*cos(6.2831853071795865d0*r) &
      /pi

end function

REAL*8 pure function udf_NLz(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_NLz = (1.0d0/2.0d0)*(sin(6.2831853071795865d0*t) + 1)**2*(-pi*r*cos( &
      2.0d0*theta) - pi*r + (1.0d0/4.0d0)*sin(12.566370614359173d0*r) - &
      1.0d0/8.0d0*sin(12.566370614359173d0*r - 2.0d0*theta) - 1.0d0/ &
      8.0d0*sin(12.566370614359173d0*r + 2.0d0*theta))*sin( &
      6.2831853071795865d0*r)**2*sin(6.2831853071795865d0*z)*cos( &
      6.2831853071795865d0*z)/(pi**2*r)

end function

REAL*8 pure function udf_Lr(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_Lr = -4*pi*(sin(6.2831853071795865d0*r)**2 - cos( &
      6.2831853071795865d0*r)**2)*(sin(6.2831853071795865d0*t) + 1)*sin &
      (6.2831853071795865d0*z)*cos(theta) - 2*pi*(sin( &
      6.2831853071795865d0*t) + 1)*sin(6.2831853071795865d0*r)**2*sin( &
      6.2831853071795865d0*z)*cos(theta) + 2.0d0*(sin( &
      6.2831853071795865d0*t) + 1)*sin(6.2831853071795865d0*r)*sin( &
      6.2831853071795865d0*z)*cos(theta)*cos(6.2831853071795865d0*r)/r

end function

REAL*8 pure function udf_La(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_La = 4*pi*(sin(6.2831853071795865d0*r)**2 - cos(6.2831853071795865d0 &
      *r)**2)*(sin(6.2831853071795865d0*t) + 1)*sin(theta)*sin( &
      6.2831853071795865d0*z) + 2*pi*(sin(6.2831853071795865d0*t) + 1)* &
      sin(theta)*sin(6.2831853071795865d0*r)**2*sin( &
      6.2831853071795865d0*z) - 2.0d0*(sin(6.2831853071795865d0*t) + 1) &
      *sin(theta)*sin(6.2831853071795865d0*r)*sin(6.2831853071795865d0* &
      z)*cos(6.2831853071795865d0*r)/r

end function

REAL*8 pure function udf_Lz(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_Lz = 1.0d0*(-(sin(6.2831853071795865d0*t) + 1)*sin( &
      6.2831853071795865d0*r)*cos(theta)*cos(6.2831853071795865d0*r)* &
      cos(6.2831853071795865d0*z)/(pi*r) - (2*r*sin( &
      6.2831853071795865d0*r)**2*sin(6.2831853071795865d0*t)*cos(theta) &
      *cos(6.2831853071795865d0*z) + 2*r*sin(6.2831853071795865d0*r)**2 &
      *cos(theta)*cos(6.2831853071795865d0*z) - 2*r*sin( &
      6.2831853071795865d0*t)*cos(theta)*cos(6.2831853071795865d0*r)**2 &
      *cos(6.2831853071795865d0*z) - 2*r*cos(theta)*cos( &
      6.2831853071795865d0*r)**2*cos(6.2831853071795865d0*z) - 2*sin( &
      6.2831853071795865d0*r)*sin(6.2831853071795865d0*t)*cos(theta)* &
      cos(6.2831853071795865d0*r)*cos(6.2831853071795865d0*z)/pi - 2* &
      sin(6.2831853071795865d0*r)*cos(theta)*cos(6.2831853071795865d0*r &
      )*cos(6.2831853071795865d0*z)/pi)/r + (1.0d0/4.0d0)*(sin( &
      6.2831853071795865d0*t) + 1)*sin(6.2831853071795865d0*r)**2*cos( &
      theta)*cos(6.2831853071795865d0*z)/(pi**2*r**2) + (-r*sin( &
      6.2831853071795865d0*r)*sin(6.2831853071795865d0*t)*cos(theta)* &
      cos(6.2831853071795865d0*r)*cos(6.2831853071795865d0*z)/pi - r* &
      sin(6.2831853071795865d0*r)*cos(theta)*cos(6.2831853071795865d0*r &
      )*cos(6.2831853071795865d0*z)/pi - 1.0d0/4.0d0*sin( &
      6.2831853071795865d0*r)**2*sin(6.2831853071795865d0*t)*cos(theta) &
      *cos(6.2831853071795865d0*z)/pi**2 - 1.0d0/4.0d0*sin( &
      6.2831853071795865d0*r)**2*cos(theta)*cos(6.2831853071795865d0*z) &
      /pi**2)/r**2)/r + (-4*pi*r*sin(6.2831853071795865d0*t)*cos( &
      6.2831853071795865d0*r) - 4*pi*r*cos(6.2831853071795865d0*r) + ( &
      sin(6.2831853071795865d0*t) + 1)*sin(6.2831853071795865d0*r) - &
      sin(6.2831853071795865d0*r)*sin(6.2831853071795865d0*t) - sin( &
      6.2831853071795865d0*r))*sin(6.2831853071795865d0*r)*cos(theta)* &
      cos(6.2831853071795865d0*z)/r + (-16*pi*r*sin( &
      6.2831853071795865d0*r)*sin(6.2831853071795865d0*t)*cos( &
      6.2831853071795865d0*r) - 16*pi*r*sin(6.2831853071795865d0*r)*cos &
      (6.2831853071795865d0*r) + 2*(sin(6.2831853071795865d0*t) + 1)* &
      sin(6.2831853071795865d0*r)**2 - 2*(sin(6.2831853071795865d0*t) + &
      1)*cos(6.2831853071795865d0*r)**2 - 6*sin(6.2831853071795865d0*r) &
      **2*sin(6.2831853071795865d0*t) - 6*sin(6.2831853071795865d0*r)** &
      2 + 6*sin(6.2831853071795865d0*t)*cos(6.2831853071795865d0*r)**2 &
      + 6*cos(6.2831853071795865d0*r)**2 + 2*(sin(6.2831853071795865d0* &
      t) + 1)*sin(6.2831853071795865d0*r)*cos(6.2831853071795865d0*r)/( &
      pi*r) - 4*(-r*sin(6.2831853071795865d0*r)**2*sin( &
      6.2831853071795865d0*t) - r*sin(6.2831853071795865d0*r)**2 + r* &
      sin(6.2831853071795865d0*t)*cos(6.2831853071795865d0*r)**2 + r* &
      cos(6.2831853071795865d0*r)**2 + sin(6.2831853071795865d0*r)*sin( &
      6.2831853071795865d0*t)*cos(6.2831853071795865d0*r)/pi + sin( &
      6.2831853071795865d0*r)*cos(6.2831853071795865d0*r)/pi)/r - 1.0d0 &
      /2.0d0*(sin(6.2831853071795865d0*t) + 1)*sin(6.2831853071795865d0 &
      *r)**2/(pi**2*r**2) + (1.0d0/2.0d0)*(4*r*sin(6.2831853071795865d0 &
      *t)*cos(6.2831853071795865d0*r) + 4*r*cos(6.2831853071795865d0*r &
      ) + sin(6.2831853071795865d0*r)*sin(6.2831853071795865d0*t)/pi + &
      sin(6.2831853071795865d0*r)/pi)*sin(6.2831853071795865d0*r)/(pi*r &
      **2))*cos(theta)*cos(6.2831853071795865d0*z)/r + 1.0d0*(-r*sin( &
      6.2831853071795865d0*t)*cos(6.2831853071795865d0*r) - r*cos( &
      6.2831853071795865d0*r) + (1.0d0/4.0d0)*(sin(6.2831853071795865d0 &
      *t) + 1)*sin(6.2831853071795865d0*r)/pi - 1.0d0/4.0d0*sin( &
      6.2831853071795865d0*r)*sin(6.2831853071795865d0*t)/pi - 1.0d0/ &
      4.0d0*sin(6.2831853071795865d0*r)/pi)*sin(6.2831853071795865d0*r) &
      *cos(theta)*cos(6.2831853071795865d0*z)/(pi*r**3)

end function

REAL*8 pure function udf_grada_r(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_grada_r = -2*(sin(6.2831853071795865d0*t) + 1)*sin(theta)*sin( &
      6.2831853071795865d0*r)*sin(6.2831853071795865d0*z)*cos( &
      6.2831853071795865d0*r)

end function

REAL*8 pure function udf_grada_a(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_grada_a = -1.0d0/2.0d0*(sin(6.2831853071795865d0*t) + 1)*sin( &
      6.2831853071795865d0*r)**2*sin(6.2831853071795865d0*z)*cos(theta) &
      /(pi*r)

end function

REAL*8 pure function udf_grada_z(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_grada_z = -(sin(6.2831853071795865d0*t) + 1)*sin(theta)*sin( &
      6.2831853071795865d0*r)**2*cos(6.2831853071795865d0*z)

end function

REAL*8 pure function udf_gradz_r(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_gradz_r = -(sin(6.2831853071795865d0*t) + 1)*sin( &
      6.2831853071795865d0*r)*cos(theta)*cos(6.2831853071795865d0*r)* &
      cos(6.2831853071795865d0*z)/(pi*r) - (2*r*sin( &
      6.2831853071795865d0*r)**2*sin(6.2831853071795865d0*t)*cos(theta) &
      *cos(6.2831853071795865d0*z) + 2*r*sin(6.2831853071795865d0*r)**2 &
      *cos(theta)*cos(6.2831853071795865d0*z) - 2*r*sin( &
      6.2831853071795865d0*t)*cos(theta)*cos(6.2831853071795865d0*r)**2 &
      *cos(6.2831853071795865d0*z) - 2*r*cos(theta)*cos( &
      6.2831853071795865d0*r)**2*cos(6.2831853071795865d0*z) - 2*sin( &
      6.2831853071795865d0*r)*sin(6.2831853071795865d0*t)*cos(theta)* &
      cos(6.2831853071795865d0*r)*cos(6.2831853071795865d0*z)/pi - 2* &
      sin(6.2831853071795865d0*r)*cos(theta)*cos(6.2831853071795865d0*r &
      )*cos(6.2831853071795865d0*z)/pi)/r + (1.0d0/4.0d0)*(sin( &
      6.2831853071795865d0*t) + 1)*sin(6.2831853071795865d0*r)**2*cos( &
      theta)*cos(6.2831853071795865d0*z)/(pi**2*r**2) + (-r*sin( &
      6.2831853071795865d0*r)*sin(6.2831853071795865d0*t)*cos(theta)* &
      cos(6.2831853071795865d0*r)*cos(6.2831853071795865d0*z)/pi - r* &
      sin(6.2831853071795865d0*r)*cos(theta)*cos(6.2831853071795865d0*r &
      )*cos(6.2831853071795865d0*z)/pi - 1.0d0/4.0d0*sin( &
      6.2831853071795865d0*r)**2*sin(6.2831853071795865d0*t)*cos(theta) &
      *cos(6.2831853071795865d0*z)/pi**2 - 1.0d0/4.0d0*sin( &
      6.2831853071795865d0*r)**2*cos(theta)*cos(6.2831853071795865d0*z) &
      /pi**2)/r**2

end function

REAL*8 pure function udf_gradz_a(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_gradz_a = ((1.0d0/4.0d0)*(sin(6.2831853071795865d0*t) + 1)*sin(theta &
      )*sin(6.2831853071795865d0*r)**2*cos(6.2831853071795865d0*z)/(pi &
      **2*r) - (r*sin(theta)*sin(6.2831853071795865d0*r)*sin( &
      6.2831853071795865d0*t)*cos(6.2831853071795865d0*r)*cos( &
      6.2831853071795865d0*z)/pi + r*sin(theta)*sin( &
      6.2831853071795865d0*r)*cos(6.2831853071795865d0*r)*cos( &
      6.2831853071795865d0*z)/pi + (1.0d0/4.0d0)*sin(theta)*sin( &
      6.2831853071795865d0*r)**2*sin(6.2831853071795865d0*t)*cos( &
      6.2831853071795865d0*z)/pi**2 + (1.0d0/4.0d0)*sin(theta)*sin( &
      6.2831853071795865d0*r)**2*cos(6.2831853071795865d0*z)/pi**2)/r)/ &
      r

end function

REAL*8 pure function udf_gradz_z(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_gradz_z = (1.0d0/2.0d0)*(sin(6.2831853071795865d0*t) + 1)*sin( &
      6.2831853071795865d0*r)**2*sin(6.2831853071795865d0*z)*cos(theta) &
      /(pi*r) - (2*r*sin(6.2831853071795865d0*r)*sin( &
      6.2831853071795865d0*t)*sin(6.2831853071795865d0*z)*cos(theta)* &
      cos(6.2831853071795865d0*r) + 2*r*sin(6.2831853071795865d0*r)*sin &
      (6.2831853071795865d0*z)*cos(theta)*cos(6.2831853071795865d0*r) + &
      (1.0d0/2.0d0)*sin(6.2831853071795865d0*r)**2*sin( &
      6.2831853071795865d0*t)*sin(6.2831853071795865d0*z)*cos(theta)/pi &
      + (1.0d0/2.0d0)*sin(6.2831853071795865d0*r)**2*sin( &
      6.2831853071795865d0*z)*cos(theta)/pi)/r

end function

REAL*8 pure function udf_gradr_r(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_gradr_r = 2*(sin(6.2831853071795865d0*t) + 1)*sin( &
      6.2831853071795865d0*r)*sin(6.2831853071795865d0*z)*cos(theta)* &
      cos(6.2831853071795865d0*r)

end function

REAL*8 pure function udf_gradr_a(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_gradr_a = -1.0d0/2.0d0*(sin(6.2831853071795865d0*t) + 1)*sin(theta)* &
      sin(6.2831853071795865d0*r)**2*sin(6.2831853071795865d0*z)/(pi*r)

end function

REAL*8 pure function udf_gradr_z(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

udf_gradr_z = (sin(6.2831853071795865d0*t) + 1)*sin(6.2831853071795865d0 &
      *r)**2*cos(theta)*cos(6.2831853071795865d0*z)

end function

REAL*8 pure function udf_scm_ur(t, theta, z, r, nu)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r
REAL*8, intent(in) :: nu

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_scm_ur = -nu*(-4*pi*(sin(6.2831853071795865d0*r)**2 - cos( &
      6.2831853071795865d0*r)**2)*(sin(6.2831853071795865d0*t) + 1)*sin &
      (6.2831853071795865d0*z)*cos(theta) - 2*pi*(sin( &
      6.2831853071795865d0*t) + 1)*sin(6.2831853071795865d0*r)**2*sin( &
      6.2831853071795865d0*z)*cos(theta) + 2.0d0*(sin( &
      6.2831853071795865d0*t) + 1)*sin(6.2831853071795865d0*r)*sin( &
      6.2831853071795865d0*z)*cos(theta)*cos(6.2831853071795865d0*r)/r &
      ) + (sin(6.2831853071795865d0*t) + 1)**2*sin(6.2831853071795865d0 &
      *r)**3*cos(theta)**2*cos(6.2831853071795865d0*r)/pi + sin( &
      6.2831853071795865d0*r)**2*sin(6.2831853071795865d0*z)*cos(theta) &
      *cos(6.2831853071795865d0*t) - 2*pi*sin(6.2831853071795865d0*r)

end function

REAL*8 pure function udf_scm_ua(t, theta, z, r, nu)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r
REAL*8, intent(in) :: nu

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_scm_ua = -nu*(4*pi*(sin(6.2831853071795865d0*r)**2 - cos( &
      6.2831853071795865d0*r)**2)*(sin(6.2831853071795865d0*t) + 1)*sin &
      (theta)*sin(6.2831853071795865d0*z) + 2*pi*(sin( &
      6.2831853071795865d0*t) + 1)*sin(theta)*sin(6.2831853071795865d0* &
      r)**2*sin(6.2831853071795865d0*z) - 2.0d0*(sin( &
      6.2831853071795865d0*t) + 1)*sin(theta)*sin(6.2831853071795865d0* &
      r)*sin(6.2831853071795865d0*z)*cos(6.2831853071795865d0*r)/r) - ( &
      sin(6.2831853071795865d0*t) + 1)**2*sin(theta)*sin( &
      6.2831853071795865d0*r)**3*cos(theta)*cos(6.2831853071795865d0*r) &
      /pi - sin(theta)*sin(6.2831853071795865d0*r)**2*sin( &
      6.2831853071795865d0*z)*cos(6.2831853071795865d0*t)

end function

REAL*8 pure function udf_scm_uz(t, theta, z, r, nu)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r
REAL*8, intent(in) :: nu

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_scm_uz = -nu*(1.0d0*(-(sin(6.2831853071795865d0*t) + 1)*sin( &
      6.2831853071795865d0*r)*cos(theta)*cos(6.2831853071795865d0*r)* &
      cos(6.2831853071795865d0*z)/(pi*r) - (2*r*sin( &
      6.2831853071795865d0*r)**2*sin(6.2831853071795865d0*t)*cos(theta) &
      *cos(6.2831853071795865d0*z) + 2*r*sin(6.2831853071795865d0*r)**2 &
      *cos(theta)*cos(6.2831853071795865d0*z) - 2*r*sin( &
      6.2831853071795865d0*t)*cos(theta)*cos(6.2831853071795865d0*r)**2 &
      *cos(6.2831853071795865d0*z) - 2*r*cos(theta)*cos( &
      6.2831853071795865d0*r)**2*cos(6.2831853071795865d0*z) - 2*sin( &
      6.2831853071795865d0*r)*sin(6.2831853071795865d0*t)*cos(theta)* &
      cos(6.2831853071795865d0*r)*cos(6.2831853071795865d0*z)/pi - 2* &
      sin(6.2831853071795865d0*r)*cos(theta)*cos(6.2831853071795865d0*r &
      )*cos(6.2831853071795865d0*z)/pi)/r + (1.0d0/4.0d0)*(sin( &
      6.2831853071795865d0*t) + 1)*sin(6.2831853071795865d0*r)**2*cos( &
      theta)*cos(6.2831853071795865d0*z)/(pi**2*r**2) + (-r*sin( &
      6.2831853071795865d0*r)*sin(6.2831853071795865d0*t)*cos(theta)* &
      cos(6.2831853071795865d0*r)*cos(6.2831853071795865d0*z)/pi - r* &
      sin(6.2831853071795865d0*r)*cos(theta)*cos(6.2831853071795865d0*r &
      )*cos(6.2831853071795865d0*z)/pi - 1.0d0/4.0d0*sin( &
      6.2831853071795865d0*r)**2*sin(6.2831853071795865d0*t)*cos(theta) &
      *cos(6.2831853071795865d0*z)/pi**2 - 1.0d0/4.0d0*sin( &
      6.2831853071795865d0*r)**2*cos(theta)*cos(6.2831853071795865d0*z) &
      /pi**2)/r**2)/r + (-4*pi*r*sin(6.2831853071795865d0*t)*cos( &
      6.2831853071795865d0*r) - 4*pi*r*cos(6.2831853071795865d0*r) + ( &
      sin(6.2831853071795865d0*t) + 1)*sin(6.2831853071795865d0*r) - &
      sin(6.2831853071795865d0*r)*sin(6.2831853071795865d0*t) - sin( &
      6.2831853071795865d0*r))*sin(6.2831853071795865d0*r)*cos(theta)* &
      cos(6.2831853071795865d0*z)/r + (-16*pi*r*sin( &
      6.2831853071795865d0*r)*sin(6.2831853071795865d0*t)*cos( &
      6.2831853071795865d0*r) - 16*pi*r*sin(6.2831853071795865d0*r)*cos &
      (6.2831853071795865d0*r) + 2*(sin(6.2831853071795865d0*t) + 1)* &
      sin(6.2831853071795865d0*r)**2 - 2*(sin(6.2831853071795865d0*t) + &
      1)*cos(6.2831853071795865d0*r)**2 - 6*sin(6.2831853071795865d0*r) &
      **2*sin(6.2831853071795865d0*t) - 6*sin(6.2831853071795865d0*r)** &
      2 + 6*sin(6.2831853071795865d0*t)*cos(6.2831853071795865d0*r)**2 &
      + 6*cos(6.2831853071795865d0*r)**2 + 2*(sin(6.2831853071795865d0* &
      t) + 1)*sin(6.2831853071795865d0*r)*cos(6.2831853071795865d0*r)/( &
      pi*r) - 4*(-r*sin(6.2831853071795865d0*r)**2*sin( &
      6.2831853071795865d0*t) - r*sin(6.2831853071795865d0*r)**2 + r* &
      sin(6.2831853071795865d0*t)*cos(6.2831853071795865d0*r)**2 + r* &
      cos(6.2831853071795865d0*r)**2 + sin(6.2831853071795865d0*r)*sin( &
      6.2831853071795865d0*t)*cos(6.2831853071795865d0*r)/pi + sin( &
      6.2831853071795865d0*r)*cos(6.2831853071795865d0*r)/pi)/r - 1.0d0 &
      /2.0d0*(sin(6.2831853071795865d0*t) + 1)*sin(6.2831853071795865d0 &
      *r)**2/(pi**2*r**2) + (1.0d0/2.0d0)*(4*r*sin(6.2831853071795865d0 &
      *t)*cos(6.2831853071795865d0*r) + 4*r*cos(6.2831853071795865d0*r &
      ) + sin(6.2831853071795865d0*r)*sin(6.2831853071795865d0*t)/pi + &
      sin(6.2831853071795865d0*r)/pi)*sin(6.2831853071795865d0*r)/(pi*r &
      **2))*cos(theta)*cos(6.2831853071795865d0*z)/r + 1.0d0*(-r*sin( &
      6.2831853071795865d0*t)*cos(6.2831853071795865d0*r) - r*cos( &
      6.2831853071795865d0*r) + (1.0d0/4.0d0)*(sin(6.2831853071795865d0 &
      *t) + 1)*sin(6.2831853071795865d0*r)/pi - 1.0d0/4.0d0*sin( &
      6.2831853071795865d0*r)*sin(6.2831853071795865d0*t)/pi - 1.0d0/ &
      4.0d0*sin(6.2831853071795865d0*r)/pi)*sin(6.2831853071795865d0*r) &
      *cos(theta)*cos(6.2831853071795865d0*z)/(pi*r**3)) - (-2*r*sin( &
      6.2831853071795865d0*r)*cos(theta)*cos(6.2831853071795865d0*r)* &
      cos(6.2831853071795865d0*t)*cos(6.2831853071795865d0*z) - 1.0d0/ &
      2.0d0*sin(6.2831853071795865d0*r)**2*cos(theta)*cos( &
      6.2831853071795865d0*t)*cos(6.2831853071795865d0*z)/pi)/r + ( &
      1.0d0/2.0d0)*(sin(6.2831853071795865d0*t) + 1)**2*(-pi*r*cos( &
      2.0d0*theta) - pi*r + (1.0d0/4.0d0)*sin(12.566370614359173d0*r) - &
      1.0d0/8.0d0*sin(12.566370614359173d0*r - 2.0d0*theta) - 1.0d0/ &
      8.0d0*sin(12.566370614359173d0*r + 2.0d0*theta))*sin( &
      6.2831853071795865d0*r)**2*sin(6.2831853071795865d0*z)*cos( &
      6.2831853071795865d0*z)/(pi**2*r) - 1.0d0/2.0d0*sin( &
      6.2831853071795865d0*r)**2*cos(theta)*cos(6.2831853071795865d0*t) &
      *cos(6.2831853071795865d0*z)/(pi*r)

end function

REAL*8 pure function udf_dpdr(t, theta, z, r)
implicit none
REAL*8, intent(in) :: t
REAL*8, intent(in) :: theta
REAL*8, intent(in) :: z
REAL*8, intent(in) :: r

REAL*8, parameter :: pi = 3.1415926535897932d0
udf_dpdr = -2*pi*sin(6.2831853071795865d0*r)

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
