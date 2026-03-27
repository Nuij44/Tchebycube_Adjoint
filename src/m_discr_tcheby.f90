module m_discr_tcheby
  use m_mesh_base
  use m_fourier_transform
  use m_numerics
  
  implicit none

  type t_discr_tcheby
     integer :: N ! nb
     real(dp),allocatable :: A(:,:), TAU(:,:)

     real(dp),allocatable :: d(:,:)
     real(dp),allocatable :: d2(:,:)
     real(dp),allocatable :: x(:) ! $ X $
     real(dp),allocatable :: g(:) ! $x^{-1}_X(X)$
     real(dp),allocatable :: h(:) ! $ to do $
     
   contains

     procedure :: init_discr_tcheby , check_discr_tcheby , free_discr_tcheby
     procedure :: eval_discr_tcheby_x,eval_discr_tcheby_y,eval_discr_tcheby_z
     procedure :: eval_discr_tcheby_x_clean,eval_discr_tcheby_y_clean,eval_discr_tcheby_z_clean

  end type t_discr_tcheby

  public ::  init_discr_tcheby , check_discr_tcheby , free_discr_tcheby
  private :: tensor_dpoduct_i,tensor_dpoduct_j,tensor_dpoduct_k

  real(kind=8),dimension(:,:,:),allocatable,private :: mem_dg1y,mem_dg2y
  real(kind=8),dimension(:,:,:),allocatable,private :: mem_dg1z,mem_dg2z
      

  
contains
  
  
  subroutine init_discr_tcheby(this,msh,opt)
    use decomp_2d
    use decomp_2d_mpi   
    use m_mesh_base

    implicit none
    !.. args
    class(t_discr_tcheby ),target :: this
    type(t_mesh_base) :: msh
    integer opt
    !.. variables
    integer :: i,j,n
    real(dp) :: c_i,c_j
    real(dp),allocatable :: nu(:)
    real(dp), allocatable:: d(:,:),d_tilde(:,:)
    ! penser à décaler les 0:n -> is:ie



    this%n = msh%n
    n = msh%n

    
    allocate( this%x(1:n+1), source = msh%x )
    
    allocate( d (1:n+1,1:n+1) )
    
    allocate( nu(0:n) )
    
    
    do i=0,N
       nu(i) = cos(i*pi/n)
    end do
    
    
    
    
    do i=0,n
       do j=0,n
   
          c_j = 1._dp
          if (j==0) c_j = 2._dp
          if (j==N) c_j = 2._dp
          
          c_i = 1._dp
          if (i==0) c_i = 2._dp
          if (i==N) c_i = 2._dp

          if (i==j) then
             
             if (i==0) then
                d(i+1,j+1) = +(2*n**2+1._dp)/6._dp
             else if (i==n) then
                d(i+1,j+1) = -(2*n**2+1._dp)/6._dp
             else
                d(i+1,j+1) = 0.5*nu(i)/(nu(i)**2-1)
             end if
             
          else

             d(i+1,j+1) = (c_i/c_j)*(-1)**(i+j)/(nu(i)-nu(j))
             
          end if
          
       end do
    end do


    

    
    d = -d*2./(msh%x(n+1)-msh%x(1))
    
    if (opt==100) then
       d=0
       do i=1,n+1
          d(i,i) = 1
       end do
       allocate( this%A(1:n+1,1:n+1) , source=d ) ! 1000
       allocate( this%TAU(1:n+1,1:n+1) , source=d ) 
    else if (opt==101) then
       allocate( this%A(1:n+1,1:n+1) , source=d )
    else if (opt==102) then
       allocate( this%A(1:n+1,1:n+1) , source=matmul(d,d))
    end if
    
    deallocate(nu)
    deallocate(d)


    if (.not.allocated(mem_dg1y) ) call alloc_y( mem_dg1y , opt_global=.true. )
    if (.not.allocated(mem_dg2y) ) call alloc_y( mem_dg2y , opt_global=.true. )
    if (.not.allocated(mem_dg1z) ) call alloc_z( mem_dg1z , opt_global=.true. )
    if (.not.allocated(mem_dg2z) ) call alloc_z( mem_dg2z , opt_global=.true. )
    
  end subroutine init_discr_tcheby

  subroutine check_discr_tcheby(this)
    implicit none
    !.. 
    class(t_discr_tcheby ) :: this
    !.. variables
      REAL(dp), allocatable  :: z(:),x_scale(:)
    real(kind=dp),allocatable :: Tk(:),Tk_der(:),dg(:),errD1(:)
    REAL(dp), allocatable  :: Tk_der2(:), dg2(:), errD2(:)
    integer :: k,n
    REAL(dp) :: xmin,xmax
    




    return

    
    n = this%n
    xmin = this%x(1)
    xmax = this%x(n+1)
    
    allocate( x_scale(1:N+1), source=0._dp)
    allocate( z(1:N+1), source=0._dp)
    x_scale = -2/(xmax-xmin)*this%x + (xmax+xmin)/(xmax-xmin)
    z = acos( x_scale )
    
    allocate( dg(1:n+1) ,tk(1:n+1),tk_der(1:n+1),errD1(0:n), source = 0._dp)
    !.. Checking for first order derivation on Tchebyshev polynomials
    do k=0,N
       
       tk = cos( k * z )
       tk_der = matmul( this%d , tk)
       
       dg(1)   =             k**2
       dg(n+1) = (-1)**(k+1)*k**2
       dg(2:n) = k * sin(k*z(2:n)) / sin(z(2:n))

       !dg = 2*this%x
       errD1(k) = maxval(abs(dg-tk_der))
       print*,k,errD1(k)
       !       do i=0,N
       !          print*,z(i),this%x(i),tk_der(i),dg(i),4*this%x(i)
       !       end do
    end do
    print'("check err D1 -> ",e15.8 )',maxval(errD1)
    
    allocate(dg2(1:N+1), Tk_der2(1:N+1), errD2(0:N))
    
    !.. Checking second order derivation on Tchebyshev polynomials
    
    do k=0,N
       tk = cos( k * z )
       Tk_der2 = matmul( this%d2 , tk)
       
       dg2(1) =             ( k**2 ) * ( k**2 - 1._dp ) / 3._dp
       dg2(n+1) = ((-1)**(k))*( k**2 ) * ( k**2 - 1._dp ) / 3._dp
       dg2(2:n ) = k * ( sin(k*z(2:N)) * cos(z(2:N)) - k * cos(k*z(2:N)) * sin(z(2:N))) / ( sin(z(2:N))**3 )
       errD2(k) = maxval(abs(dg2-tk_der2))

    end do

    print'("check err D2 -> ",e15.8 )',maxval(errD2)
    
    deallocate(dg,tk,tk_der,errD1)
    deallocate(Tk_der2,dg2,errD2)
    deallocate(x_scale,z)
 
  end subroutine check_discr_tcheby
  
  
  subroutine free_discr_tcheby(this)
    implicit none
    class(t_discr_tcheby ) :: this
    deallocate(this%A)
  end subroutine free_discr_tcheby
  
  SUBROUTINE eval_discr_tcheby_x(this,fi,dfi)
    use decomp_2d
    implicit none
    CLASS(t_discr_tcheby) :: this
    real(kind=8),dimension(:,:,:),allocatable :: fi,dfi
    TYPE(DECOMP_INFO)  :: dd
    
    CALL GET_DECOMP_INFO(DD)
    call tensor_dpoduct_i(THIS%A,FI,DFI,DD%XST,DD%XEN)
    
  END SUBROUTINE EVAL_DISCR_TCHEBY_X


  SUBROUTINE eval_discr_tcheby_y(this,fi,dfi,DG1Y,DG2Y)
    use decomp_2d
    use decomp_2d_fft
    use decomp_2d_constants
   use decomp_2d_mpi
   implicit none
   CLASS(t_discr_tcheby) :: this
    real(mytype),dimension(:,:,:),allocatable :: fi,dfi
    real(mytype),dimension(:,:,:),allocatable,optional :: DG1Y,DG2Y
    TYPE(DECOMP_INFO)  :: dd
    
    CALL GET_DECOMP_INFO(DD)
    
    CALL TRANSPOSE_X_TO_Y(     FI, MEM_DG1Y )
    CALL TENSOR_dpODUCT_J( THIS%A, MEM_DG1Y, MEM_DG2Y, DD%YST, DD%YEN )
    CALL TRANSPOSE_Y_TO_X(   MEM_DG2Y,  DFI )
    
  END SUBROUTINE EVAL_DISCR_TCHEBY_Y
  
  SUBROUTINE eval_discr_tcheby_z(this,fi,dfi,DG1Y,DG2Y,DG1Z,DG2Z)
    use decomp_2d
    implicit none
    CLASS(t_discr_tcheby) :: this
    real(kind=8),dimension(:,:,:),allocatable :: fi,dfi
    real(kind=8),dimension(:,:,:),allocatable,optional :: DG1Y,DG2Y
    real(kind=8),dimension(:,:,:),allocatable,optional :: DG1Z,DG2Z
    TYPE(DECOMP_INFO)  :: dd
    
    CALL GET_DECOMP_INFO(DD)
    
    CALL TRANSPOSE_X_TO_Y(   FI, MEM_DG1Y )
    CALL TRANSPOSE_Y_TO_Z( MEM_DG1Y, MEM_DG1Z )
    CALL TENSOR_dpODUCT_K(THIS%A, MEM_DG1Z, MEM_DG2Z, DD%ZST, DD%ZEN )
    CALL TRANSPOSE_Z_TO_Y( MEM_DG2Z , MEM_DG2Y  )
    CALL TRANSPOSE_Y_TO_X( MEM_DG2Y ,  DFI  )
    
  END SUBROUTINE EVAL_DISCR_TCHEBY_Z

  SUBROUTINE eval_discr_tcheby_x_clean(this,fi,dfi)
    use decomp_2d
    implicit none
    CLASS(t_discr_tcheby) :: this
    real(kind=8),dimension(:,:,:),allocatable :: fi,dfi
    TYPE(DECOMP_INFO)  :: dd
    
    CALL GET_DECOMP_INFO(DD)
    call tensor_dpoduct_i(THIS%A,FI,DFI,DD%XST,DD%XEN)
    
  END SUBROUTINE EVAL_DISCR_TCHEBY_X_clean


  SUBROUTINE eval_discr_tcheby_y_clean(this,fi,dfi,DG1Y,DG2Y)
    use decomp_2d
    use decomp_2d_fft
    use decomp_2d_constants
   use decomp_2d_mpi
   implicit none
   CLASS(t_discr_tcheby) :: this
    real(mytype),dimension(:,:,:),allocatable :: fi,dfi
    real(mytype),dimension(:,:,:),allocatable,optional :: DG1Y,DG2Y
    TYPE(DECOMP_INFO)  :: dd
    
    CALL GET_DECOMP_INFO(DD)
    
    CALL TRANSPOSE_X_TO_Y(     FI, MEM_DG1Y )
    CALL TENSOR_dpODUCT_J( THIS%A, MEM_DG1Y, MEM_DG2Y, DD%YST, DD%YEN )
    CALL TRANSPOSE_Y_TO_X(   MEM_DG2Y,  DFI )
    
  END SUBROUTINE EVAL_DISCR_TCHEBY_Y_clean
  
  SUBROUTINE eval_discr_tcheby_z_clean(this,fi,dfi,DG1Y,DG2Y,DG1Z,DG2Z)
    use decomp_2d
    implicit none
    CLASS(t_discr_tcheby) :: this
    real(kind=8),dimension(:,:,:),allocatable :: fi,dfi
    real(kind=8),dimension(:,:,:),allocatable,optional :: DG1Y,DG2Y
    real(kind=8),dimension(:,:,:),allocatable,optional :: DG1Z,DG2Z
    TYPE(DECOMP_INFO)  :: dd
    INTEGER :: N,K_max,k
    
    CALL GET_DECOMP_INFO(DD)
    
    CALL TRANSPOSE_X_TO_Y(   FI, MEM_DG1Y )
    CALL TRANSPOSE_Y_TO_Z( MEM_DG1Y, MEM_DG1Z )


    N = DD%ZEN(3)-DD%ZST(3)+1
    
    k_max = floor(2._DP/3._DP*N)
    print*,k_max
    CALL dct_1m_z(MEM_DG1Z,MEM_DG2Z,dct_plan_fwd_z,N)
    CALL TENSOR_dpODUCT_K(THIS%TAU, MEM_DG2Z, MEM_DG1Z, DD%ZST, DD%ZEN )

    do k = k_max, N
       MEM_DG1Z(:,:,K) = 0._DP
    END do
    
    print*,'k_max:',k_max

    do k = 1, n
       print*,MEM_DG1Z(:,:,k)
    end do

    CALL dct_1M_z(MEM_DG2Z,MEM_DG1Z,dct_plan_bck_z,N)
    MEM_DG2Z=MEM_DG1Z

    CALL TRANSPOSE_Z_TO_Y( MEM_DG2Z , MEM_DG2Y  )
    CALL TRANSPOSE_Y_TO_X( MEM_DG2Y ,  DFI  )
    
  END SUBROUTINE EVAL_DISCR_TCHEBY_Z_clean

  
  subroutine tensor_dpoduct_i(DX,FI,DFI,XST,XEN)
    implicit none
    INTEGER :: XST(3),XEN(3)
    real(dp) ::  DX(XST(1):XEN(1),XST(1):XEN(1))
    real(dp) ::  FI(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
    real(dp) :: DFI(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
    
    real(dp) ::  ALPHA,BETA
    INTEGER :: N(3),NI,NJ,NK
    
    N = XEN-XST + 1
    NI = N(1)
    NJ = N(2)
    NK = N(3)
    ALPHA = 1._dp
    BETA = 0.0_dp
    
    CALL DGEMM( 'N', 'N', NI, NJ*NK, NI, ALPHA, DX, NI, FI, NI, BETA, DFI, NI )
    
  end subroutine tensor_dpoduct_i
  
  SUBROUTINE TENSOR_dpODUCT_J(DY,FI,DFI,YST,YEN)
    IMPLICIT NONE
        INTEGER :: YST(3),YEN(3)
        real(dp) :: DY(YST(2):YEN(2),YST(2):YEN(2))
        real(dp) ::  FI(YST(1):YEN(1),YST(2):YEN(2),YST(3):YEN(3))
        real(dp) :: DFI(YST(1):YEN(1),YST(2):YEN(2),YST(3):YEN(3))
        real(dp) ::  ALPHA,BETA
        INTEGER :: K,N(3),NI,NJ,NK
        
        ALPHA = 1._dp
        BETA = 0.0_dp
        N = YEN - YST + 1
        NI = N(1)
        NJ = N(2)
        NK = N(3)
        DO K=YST(3),YEN(3)
           CALL DGEMM('N','T',NI,NJ,NJ,ALPHA,FI(YST(1),YST(2),K),NI,DY,NJ,BETA,DFI(YST(1),YST(2),K),NI)
        END DO
        
      END SUBROUTINE TENSOR_dpODUCT_J


      
      
      SUBROUTINE TENSOR_dpODUCT_K(DZ,FI,DFI,ZST,ZEN)
        IMPLICIT NONE
        INTEGER :: ZST(3),ZEN(3)
        real(dp) :: DZ(ZST(3):ZEN(3),ZST(3):ZEN(3)) 
        real(dp) ::  FI(ZST(1):ZEN(1),ZST(2):ZEN(2),ZST(3):ZEN(3))
        real(dp) :: DFI(ZST(1):ZEN(1),ZST(2):ZEN(2),ZST(3):ZEN(3))

        INTEGER :: N(3),NI,NJ,NK
        real(dp) ::  ALPHA,BETA
        
        ALPHA = 1._dp
        BETA = 0._dp
        N = ZEN - ZST + 1
        NI = N(1)
        NJ = N(2)
        NK = N(3)
        CALL DGEMM('N','T',NI*NJ,NK,NK, ALPHA, FI, NI*NJ , DZ, NK , BETA , DFI, NI*NJ )
        
      END SUBROUTINE TENSOR_dpODUCT_K



  
end module m_discr_tcheby
