module m_navier_stokes_cart
  use m_operator_base
  use m_mesh_base
  implicit none

  INTERFACE DIV
     MODULE PROCEDURE NS_INC_COLOCATED_CARTESIAN_DIVERGENCE
  END INTERFACE DIV

  INTERFACE GRAD
     MODULE PROCEDURE NS_INC_COLOCATED_CARTESIAN_GRADIENT
  END INTERFACE GRAD

  INTERFACE CROSS
     MODULE PROCEDURE NS_INC_COLOCATED_CARTESIAN_CROSS_PRODUCT
  END INTERFACE CROSS

  INTERFACE CURL
     MODULE PROCEDURE NS_INC_COLOCATED_CARTESIAN_CURL
  END INTERFACE CURL

  INTERFACE VORTICITY
     MODULE PROCEDURE NS_INC_COLOCATED_CARTESIAN_VORTICITY
  END INTERFACE VORTICITY

  INTERFACE COMPUTE_NON_LINEAR_TERMS
     MODULE PROCEDURE M_NAVIER_STOKES_CART_NLINEAR_MOMENTUM_conv
     MODULE PROCEDURE M_NAVIER_STOKES_CART_NLINEAR_SCALAR_SKEW
  END INTERFACE COMPUTE_NON_LINEAR_TERMS
  
  INTERFACE COMPUTE_LINEAR_TERMS
     MODULE PROCEDURE M_NAVIER_STOKES_CART_LINEAR_MOMENTUM
     MODULE PROCEDURE M_NAVIER_STOKES_CART_LINEAR_SCALAR
  END INTERFACE COMPUTE_LINEAR_TERMS
  
contains
  
  subroutine ns_inc_colocated_cartesian_divergence( op_x1, op_x2, op_x3, u1, u2, u3, div, dg1, dg2 , dg3 )
    implicit none
    class(t_operator_base) :: op_x1, op_x2, op_x3
    real(kind=8),dimension(:,:,:),allocatable :: u1, u2, u3, div, dg1, dg2,dg3
    
    call op_x1%d1(u1,dg1)
    call op_x2%d1(u2,dg2)
    call op_x3%d1(u3,dg3)
    div = dg1 + dg2 + dg3
  end subroutine ns_inc_colocated_cartesian_divergence
  
  subroutine ns_inc_colocated_cartesian_gradient(&
       op_x1, op_x2, op_x3, fi, fi_dx1, fi_dx2, fi_dx3)
    implicit none
    class(t_operator_base) :: op_x1, op_x2, op_x3
    real(kind=8),dimension(:,:,:),allocatable :: fi, fi_dx1, fi_dx2, fi_dx3
    
    call op_x1%d1(fi,fi_dx1)
    call op_x2%d1(fi,fi_dx2)
    call op_x3%d1(fi,fi_dx3)
    
  end subroutine ns_inc_colocated_cartesian_gradient  

  subroutine ns_inc_colocated_cartesian_cross_product(&
       op_x1, op_x2, op_x3, ax, ay, az, bx, by, bz, cx, cy, cz)
    implicit none
    class(t_operator_base) :: op_x1, op_x2, op_x3
    real(kind=8),dimension(:,:,:),allocatable :: ax,ay,az
    real(kind=8),dimension(:,:,:),allocatable :: bx,by,bz
    real(kind=8),dimension(:,:,:),allocatable :: cx,cy,cz
    
    cx = ay*bz - az*by
    cy = az*bx - ax*bz
    cz = ax*by - ay*bx
    
  end subroutine ns_inc_colocated_cartesian_cross_product

  subroutine ns_inc_colocated_cartesian_curl(&
       op_x1, op_x2, op_x3, ax, ay, az, cx, cy, cz,dg1,dg2,dg3)
    implicit none
    class(t_operator_base) :: op_x1, op_x2, op_x3
    real(kind=8),dimension(:,:,:),allocatable :: ax,ay,az
    real(kind=8),dimension(:,:,:),allocatable :: cx,cy,cz
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: DG1,DG2,DG3

    CALL OP_X2%D1(az,dg1)
    CALL OP_X3%D1(ay,dg2)
    cx = dg1 - dg2

    CALL OP_X3%D1(ax,dg1)
    CALL OP_X1%D1(az,dg2)
    cy = dg1 - dg2

    CALL OP_X1%D1(ay,dg1)
    CALL OP_X2%D1(ax,dg2)
    cz = dg1 - dg2

  end subroutine ns_inc_colocated_cartesian_curl

  subroutine m_navier_stokes_cart_nlinear_momentum_skew(&
       op_x1,op_x2,op_x3, u, v, w, hu, hv, hw,&
       dg1, dg2, dg3, dg4, dg5, dg6, dg7, dg8, dg9)
    implicit none
    
    class(t_operator_base)     :: op_x1, op_x2, op_x3
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: U,HU
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: V,HV
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: W,HW
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: DG1,DG2,DG3
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: DG4,DG5,DG6
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: DG7,DG8,DG9

    
    
    CALL OP_X1%D1(U,DG1)
    CALL OP_X2%D1(U,DG2)
    CALL OP_X3%D1(U,DG3)
    
    DG4 = U*DG1 
    DG5 = V*DG2 
    DG6 = W*DG3 
    
    HU = (DG4 + DG5 + DG6)*0.5
    
    DG4 = U*U 
    DG5 = V*U 
    DG6 = W*U 
    
    CALL OP_X1%D1(DG4,DG7)
    CALL OP_X2%D1(DG5,DG8)
    CALL OP_X3%D1(DG6,DG9)
    HU = HU + (dg7 + dg8 + dg9)*0.5

    CALL OP_X1%D1(V,DG1)
    CALL OP_X2%D1(V,DG2)
    CALL OP_X3%D1(V,DG3)
    
    DG4 = U*DG1 
    DG5 = V*DG2 
    DG6 = W*DG3 
    
    HV = (DG4 + DG5 + DG6)*0.5
    
    DG4 = U*V 
    DG5 = V*V 
    DG6 = W*V 
    
    CALL OP_X1%D1(DG4,DG7)
    CALL OP_X2%D1(DG5,DG8)
    CALL OP_X3%D1(DG6,DG9)
    HV = HV + (DG7 + DG8 + DG9)*0.5

    !>
    CALL OP_X1%D1(W,DG1)
    CALL OP_X2%D1(W,DG2)
    CALL OP_X3%D1(W,DG3)
    
    DG4 = U*DG1 
    DG5 = V*DG2 
    DG6 = W*DG3 
    
    HW = (DG4 + DG5 + DG6)*0.5
    
    DG4 = U*W 
    DG5 = V*W 
    DG6 = W*W 
    
    CALL OP_X1%D1(DG4,DG7)
    CALL OP_X2%D1(DG5,DG8)
    CALL OP_X3%D1(DG6,DG9)
    HW = HW + (DG7 + DG8 + DG9)*0.5
    
  END subroutine m_navier_stokes_cart_nlinear_momentum_skew



    subroutine m_navier_stokes_cart_nlinear_momentum_conv(&
       op_x1,op_x2,op_x3, u, v, w, hu, hv, hw,&
       dg1, dg2, dg3, dg4, dg5, dg6, dg7, dg8, dg9)
    implicit none
    
    class(t_operator_base)     :: op_x1, op_x2, op_x3
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: U,HU
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: V,HV
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: W,HW
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: DG1,DG2,DG3
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: DG4,DG5,DG6
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: DG7,DG8,DG9

    
    
    CALL OP_X1%D1(U,DG1)
    CALL OP_X2%D1(U,DG2)
    CALL OP_X3%D1(U,DG3)
    
    DG4 = U*DG1 
    DG5 = V*DG2 
    DG6 = W*DG3 
    
    HU = (DG4 + DG5 + DG6)
    

    CALL OP_X1%D1(V,DG1)
    CALL OP_X2%D1(V,DG2)
    CALL OP_X3%D1(V,DG3)
    
    DG4 = U*DG1 
    DG5 = V*DG2 
    DG6 = W*DG3 
    
    HV = (DG4 + DG5 + DG6)
    

    !>
    CALL OP_X1%D1(W,DG1)
    CALL OP_X2%D1(W,DG2)
    CALL OP_X3%D1(W,DG3)
    
    DG4 = U*DG1 
    DG5 = V*DG2 
    DG6 = W*DG3 
    
    HW = (DG4 + DG5 + DG6)
    

    
  END subroutine m_navier_stokes_cart_nlinear_momentum_conv



  subroutine m_navier_stokes_cart_nlinear_scalar_skew(&
       op_x1,op_x2,op_x3, u, v, w, T , ht, &
       dg1, dg2, dg3, dg4, dg5, dg6, dg7, dg8, dg9)
    implicit none
    
    class(t_operator_base)     :: op_x1, op_x2, op_x3
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: U,v,w
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: t,ht
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: DG1,DG2,DG3
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: DG4,DG5,DG6
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: DG7,DG8,DG9
    
    CALL OP_X1%D1(T,DG1)
    CALL OP_X2%D1(T,DG2)
    CALL OP_X3%D1(T,DG3)
    
    DG4 = U*DG1 
    DG5 = V*DG2 
    DG6 = W*DG3 
    
    HT = (DG4 + DG5 + DG6)*0.5
   
    DG4 = U*T 
    DG5 = V*T 
    DG6 = W*T 
    
    CALL OP_X1%D1(DG4,DG7)
    CALL OP_X2%D1(DG5,DG8)
    CALL OP_X3%D1(DG6,DG9)
    HT = HT + (dg7 + dg8 + dg9)*0.5
!!$

    
    

  END subroutine m_navier_stokes_cart_nlinear_scalar_skew

  subroutine m_navier_stokes_cart_linear_momentum(&
       op_x1,op_x2,op_x3, nu,u, v, w, lu, lv, lw,&
       dg1, dg2, dg3 )
    implicit none
    
    class(t_operator_base)     :: op_x1, op_x2, op_x3
    REAL(KIND=8)               :: nu
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: U,lU
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: V,lV
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: W,lW
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: DG1,DG2,DG3
    
    CALL OP_X1%D2(U,DG1)
    CALL OP_X2%D2(U,DG2)
    CALL OP_X3%D2(U,DG3)    
    LU = NU*(DG1+DG2+DG3)

    CALL OP_X1%D2(V,DG1)
    CALL OP_X2%D2(V,DG2)
    CALL OP_X3%D2(V,DG3)    
    LV = NU*(DG1+DG2+DG3)

    CALL OP_X1%D2(W,DG1)
    CALL OP_X2%D2(W,DG2)
    CALL OP_X3%D2(W,DG3)    
    LW = NU*(DG1+DG2+DG3)
    
  END subroutine m_navier_stokes_cart_linear_momentum

  
  subroutine m_navier_stokes_cart_linear_scalar(&
       op_x1,op_x2,op_x3,nu, t, lt,&
       dg1, dg2, dg3 )
    implicit none
    
    class(t_operator_base)     :: op_x1, op_x2, op_x3
    REAL(KIND=8)                 :: nu
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: T,lT
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: DG1,DG2,DG3
    
    CALL OP_X1%D2(T,DG1)
    CALL OP_X2%D2(T,DG2)
    CALL OP_X3%D2(T,DG3)    
    LT = NU*(DG1+DG2+DG3)

    
  END subroutine m_navier_stokes_cart_linear_scalar

  
  SUBROUTINE NS_INC_COLOCATED_CARTESIAN_VORTICITY(OPX,OPY,OPZ,U,V,W,VORT_U,VORT_V,VORT_W,DG1,DG2,DG3,DG4,DG5,DG6)
    implicit none
    class(t_operator_base)  ::  opx,opy,opz
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE,INTENT(IN) :: U,V,W
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE,INTENT(OUT) :: VORT_U,VORT_V,VORT_W
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: DG1,DG2,DG3,DG4,DG5,DG6

    CALL OPY%D1(W,DG1)
    CALL OPZ%D1(V,DG2)

    VORT_U = DG1 - DG2

    CALL OPZ%D1(U,DG3)
    CALL OPX%D1(W,DG4)

    VORT_V = DG3 - DG4

    CALL OPX%D1(V,DG5)
    CALL OPY%D1(W,DG6)

    VORT_W = DG5 - DG6
  END SUBROUTINE NS_INC_COLOCATED_CARTESIAN_VORTICITY

  
end module m_navier_stokes_cart
