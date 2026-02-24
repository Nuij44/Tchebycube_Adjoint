module m_fourier_transform
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'

  complex(kind=8), allocatable,public,save :: dg_x(:,:,:),zg_x(:,:,:)
  complex(kind=8), allocatable,public,save :: zg1_y(:,:,:),zg2_y(:,:,:)  
  complex(kind=8), allocatable,public,save :: zg1_z(:,:,:),zg2_z(:,:,:)  
  real(kind=8),allocatable,public,save :: dg_z(:,:,:)

  ! for computing the 
  type(C_PTR),public :: dct_fwd_x, dct_bck_x
  type(C_PTR),public :: dct_fwd_y, dct_bck_y
  type(C_PTR),public :: dct_fwd_z, dct_bck_z
  
  type(C_PTR),public :: plan_fwd_x, plan_bck_x
  type(C_PTR),public :: plan_fwd_y, plan_bck_y
  type(C_PTR),public :: plan_fwd_z, plan_bck_z

  type(C_PTR),public :: plan_r2c_x, plan_c2r_x
  type(C_PTR),public :: plan_r2c_y, plan_c2r_y
  type(C_PTR),public :: plan_r2c_c, plan_c2r_z

  type(C_PTR),public,save :: plan_r2hc_x, plan_hc2r_x
  type(C_PTR),public,save :: plan_r2hc_y, plan_hc2r_y

  type(C_PTR),public,save :: dct_plan_fwd_z,dct_plan_bck_z
  
  real(kind=8),public :: normalization_factor 
  interface
     subroutine fftw_execute(plan) bind(C, name="fftw_execute")
       use, intrinsic :: iso_c_binding
       type(C_PTR), value :: plan
     end subroutine fftw_execute
  end interface
  
  integer, save :: plan_type = FFTW_ESTIMATE !FFTW_MEASURE
  real(kind=8),public :: norm_factor(3),norm_factor_dct(3)
  
contains


  subroutine test_dct(&
       xst,xen,nx,&
       yst,yen,ny,&
       zst,zen,nz )
    use decomp_2d
    use decomp_2d_mpi
    use m_tensor_product
    implicit none

    integer :: xst(3),xen(3),nx
    integer :: yst(3),yen(3),ny
    integer :: zst(3),zen(3),nz

    integer :: i,j,k
    real(kind=8) :: hx,hy,hz
    integer :: ierr
    real(kind=8),allocatable :: wrk(:,:,:), wrk1(:,:,:)
    real(kind=8),allocatable :: X(:,:,:), Y(:,:,:),Z(:,:,:)
    real(kind=8),allocatable :: quad_z(:)
    real(kind=8) :: quad
    real(kind=8) :: pi=acos(-1.)
    

    
    call dct_1m_x_plan(dct_fwd_x,xst,xen, FFTW_FORWARD )
    call dct_1m_y_plan(dct_fwd_y,yst,yen, FFTW_FORWARD )
    call dct_1m_z_plan(dct_fwd_z,zst,zen, FFTW_FORWARD )
    
    ALLOCATE( X(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3)) )
    ALLOCATE( Y(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3)) )
    ALLOCATE( Z(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3)) )

    FORALL(I=XST(1):XEN(1),J=XST(2):XEN(2),K=XST(3):XEN(3))
       X(I,J,K) = COS( DBLE(I-1)/DBLE(NX-1) * PI  )
       Y(I,J,K) = COS( DBLE(J-1)/DBLE(NY-1) * PI  )
       Z(I,J,K) = COS( DBLE(K-1)/DBLE(NZ-1) * PI  )
    END FORALL
    

    ALLOCATE( WRK (XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3)) )
    ALLOCATE( WRK1(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3)) )
    WRK = 0
    WRK1 = 0
    
    
    FORALL(I=XST(1):XEN(1),J=XST(2):XEN(2),K=XST(3):XEN(3))
       WRK(I,J,K) = 1 !+COS( 3*ACOS(Z(I,J,K)) )
    END FORALL
    !CALL FFTW_EXECUTE_R2R(DCT_FWD_X, WRK,WRK1)
    
    
!!$    CALL TRANSPOSE_X_TO_Y(WRK,DG1_Y)
!!$    CALL TRANSPOSE_Y_TO_Z(DG1_Y,DG1_Z)
!!$    call fftw_execute_r2r(dct_fwd_z,DG1_Z,DG2_Z)
!!$    CALL TRANSPOSE_Z_TO_Y(DG2_Z,DG1_Y)
!!$    CALL TRANSPOSE_Y_TO_X(DG1_Y,WRK1)
!!$    wrk1=wrk1/(nz-1)
!!$    wrk1(:,:,xst(3))= wrk1(:,:,xst(3))*.5
!!$
!!$    allocate( quad_z( zst(3):zen(3) ) )
!!$    quad_z = 0
!!$    do k=zst(3),zen(3)
!!$       if (mod(k-1,2)==0) then
!!$          quad_z(k) = 2./((k-1)**2-1)
!!$       end if
!!$    end do
!!$    
!!$    do i=xst(3),xen(3)
!!$       do j=yst(3),yen(3)
!!$          quad=0
!!$          do k=zst(3),zen(3)
!!$             quad = quad + wrk1(i,j,k)*quad_z(k) 
!!$          end do
!!$          do k=zst(3),zen(3)
!!$             wrk1(i,j,k) = quad
!!$          end do
!!$       end do
!!$    end do
    
!!$    if (nrank==0) then
!!$       do k=xst(3),xen(3)
!!$          print*,z(3,3,k),wrk1(3,3,k)
!!$       end do
!!$    end if
  end subroutine test_dct

  
  subroutine start_fourier_transforms(&
       xst,xen,nx,&
       yst,yen,ny,&
       zst,zen,nz )
    use decomp_2d
    use decomp_2d_mpi
    implicit none
    

    integer :: xst(3),xen(3),nx
    integer :: yst(3),yen(3),ny
    integer :: zst(3),zen(3),nz

    integer :: i,j,k
    real(kind=8) :: hx,hy,hz
    integer :: ierr
    real(kind=8),allocatable :: wrk(:,:,:), wrk1(:,:,:)
    real(kind=8) :: diff

    norm_factor(1) = 1./sqrt( dble(nx) )
    call c2c_1m_x_plan(plan_fwd_x,xst,xen, FFTW_FORWARD )
    call c2c_1m_x_plan(plan_bck_x,xst,xen, FFTW_BACKWARD)
    allocate(dg_x(xst(1):xen(1),xst(2):xen(2),xst(3):xen(3)))
    allocate(zg_x(xst(1):xen(1),xst(2):xen(2),xst(3):xen(3)))
    
    norm_factor(2) = 1./sqrt( dble(ny) )

    call c2c_1m_y_plan(plan_fwd_y,yst,yen, FFTW_FORWARD )
    call c2c_1m_y_plan(plan_bck_y,yst,yen, FFTW_BACKWARD)
    allocate(zg1_y(yst(1):yen(1),yst(2):yen(2),yst(3):yen(3)))
    allocate(zg2_y(yst(1):yen(1),yst(2):yen(2),yst(3):yen(3)))

    
    norm_factor(3) = 1./sqrt( dble(nz) )
    norm_factor_dct(3) = 1./( dble(nz-1) )

    call c2c_1m_z_plan(plan_fwd_z,zst,zen, FFTW_FORWARD )
    call c2c_1m_z_plan(plan_bck_z,zst,zen, FFTW_BACKWARD)
    CALL dct_1m_z_plan(dct_plan_fwd_z,zst,zen, FFTW_FORWARD )
    CALL dct_1m_z_plan(dct_plan_bck_z,zst,zen, FFTW_BACKWARD)
    allocate(zg1_z(zst(1):zen(1),zst(2):zen(2),zst(3):zen(3)))
    allocate(zg2_z(zst(1):zen(1),zst(2):zen(2),zst(3):zen(3)))
    allocate(dg_z(zst(1):zen(1),zst(2):zen(2),zst(3):zen(3)))

    
    hx = 2*acos(-1d0)/nx
    hy = 2*acos(-1d0)/ny
    hz = 2*acos(-1d0)/nz
    
    !> init in x-pencil cplx
    forall(i=xst(1):xen(1),j=xst(2):xen(2),k=xst(3):xen(3))
       dg_x(i,j,k) = cos( 2*(j-1)*hy ) + sin( 4*(k-1)*hz )
    end forall

    allocate(wrk(xst(1):xen(1),xst(2):xen(2),xst(3):xen(3)))
    allocate(wrk1(xst(1):xen(1),xst(2):xen(2),xst(3):xen(3)))
    call random_number(wrk)
    dg_x = wrk
    
    call c2c_1m_x(dg_x,plan_fwd_x)
    call transpose_x_to_y(dg_x,zg1_y)
    call c2c_1m_y(zg1_y,plan_fwd_y)
    call transpose_y_to_z(zg1_y,zg1_z)
    call c2c_1m_z(zg1_z,plan_fwd_z)
    
    call c2c_1m_z(zg1_z,plan_bck_z)
    call transpose_z_to_y(zg1_z,zg1_y)
    call c2c_1m_y(zg1_y,plan_bck_y)
    call transpose_y_to_x(zg1_y,dg_x)
    call c2c_1m_x(dg_x,plan_bck_x)

    diff=maxval(abs(wrk - dg_x))
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,diff,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,IERR)
    if (nrank==0) then
       print'("error on fft  base change",e15.8)',diff
    end if
    
    
    call r2r_1m_x_plan(plan_r2hc_x, xst,xen,+1)
    call r2r_1m_x_plan(plan_hc2r_x, xst,xen,-1)
    
    forall(i=xst(1):xen(1),j=xst(2):xen(2),k=xst(3):xen(3))
       wrk(i,j,k) = 1+cos(2*(i-xst(1))*hx)+sin(2*(i-xst(1))*hx)
    end forall
    call fftw_execute_r2r(plan_r2hc_x, wrk,wrk1)
    call fftw_execute_r2r(plan_hc2r_x, wrk1,wrk)
    
    deallocate(wrk )
    deallocate(wrk1)

    
    allocate(wrk(yst(1):yen(1),yst(2):yen(2),yst(3):yen(3)))
    allocate(wrk1(yst(1):yen(1),yst(2):yen(2),yst(3):yen(3)))

    
    forall(i=yst(1):yen(1),j=yst(2):yen(2),k=yst(3):yen(3))
       wrk(i,j,k) = 1+cos(2*(j-1)*hy)+sin(2*(j-1)*hy)
    end forall
    call r2r_1m_y_plan(plan_r2hc_y, yst,yen,+1)
    call r2r_1m_y_plan(plan_hc2r_y, yst,yen,-1)

    call r2r_1m_y(wrk,wrk1,plan_r2hc_y)
    
    return
!!$    do j=yst(1),yen(2)
!!$       print*,wrk1(1,j,1)
!!$    end do

    
    
  end subroutine start_fourier_transforms

  ! for chebychev transform 

   subroutine dct_1m_z(in, out, plan1,nz)
     
      implicit none
      
      real(kind=8), dimension(:, :, :), intent(IN) :: in
      integer, intent(in) :: nz
      real(kind=8), dimension(:, :, :), intent(out) :: out
      type(C_PTR), intent(IN) :: plan1
      
      DG_Z = IN
      DG_Z(:,:,1 ) = DG_Z(:,:,1 )*sqrt(dble(2))
      DG_Z(:,:,nz) = DG_Z(:,:,nz)*sqrt(dble(2))
      call dfftw_execute_r2r(plan1, DG_Z, out)
!      out = out*norm_factor_DCT(3)
      out = out*dble(0.5)/sqrt(dble(nz-1))
      out(:,:,2:nz-1) = out(:,:,2:nz-1)*sqrt(dble(2))
      return
    end subroutine dct_1m_z
    

  subroutine dct_1m_x_plan(plan1,xst,xen, isign)
    
    implicit none

    type(C_PTR), intent(OUT) :: plan1
    integer, intent(IN) :: xst(3),xen(3)
    integer, intent(IN) :: isign
    
    real(C_DOUBLE), pointer :: a1(:, :, :)
    real(C_DOUBLE), pointer :: a1o(:, :, :)

    type(C_PTR) :: a1_p
    integer(C_SIZE_T) :: sz
    integer :: xsz(3)

    !complex(kind=8), allocatable, dimension(:, :, :) :: a1

    
    xsz = xen-xst+1
    sz = xsz(1)*xsz(2)*xsz(3)
    a1_p = fftw_alloc_real(sz)

    call c_f_pointer(a1_p, a1 , [xsz(1),xsz(2),xsz(3)] )
    call c_f_pointer(a1_p, a1o, [xsz(1),xsz(2),xsz(3)] )
    
    plan1 = fftw_plan_many_r2r(&
         1, [xsz(1)], xsz(2) * xsz(3), &
         a1 , xsz(1), 1, xsz(1), &
         a1o, xsz(1), 1, xsz(1), &
         [ FFTW_REDFT00 ], FFTW_ESTIMATE)

    call fftw_free(a1_p)
    
    return
  end subroutine dct_1m_x_plan



  subroutine dct_1m_y_plan(plan1, yst,yen, isign)
    use, intrinsic :: iso_c_binding
    implicit none
    
    type(C_PTR) :: plan1
    integer :: yst(3),yen(3)
    integer, intent(IN) :: isign
    
    real(C_DOUBLE), pointer :: a1(:, :)
    real(C_DOUBLE), pointer :: a1o(:, :)
    
    type(C_PTR) :: a1_p
    integer(C_SIZE_T) :: sz
    integer :: ysz(3)

    ysz = yen-yst+1
    sz = ysz(1) * ysz(2)
    a1_p = fftw_alloc_real(sz)

    call c_f_pointer(a1_p, a1 , [ysz(1),ysz(2)])
    call c_f_pointer(a1_p, a1o, [ysz(1),ysz(2)])

    plan1 = fftw_plan_many_r2r(1, ysz(2), ysz(1), &
         a1 , ysz(2), ysz(1), 1, &
         a1o, ysz(2), ysz(1), 1, &
         [ FFTW_REDFT00 ], FFTW_ESTIMATE)
    
    call fftw_free(a1_p)
    
    return
  end subroutine dct_1m_y_plan
  

  subroutine dct_1m_z_plan(plan1,zst,zen, isign)

    implicit none
    
    type(C_PTR), intent(OUT) :: plan1
    integer, intent(IN) :: isign
    integer, intent(IN) :: zst(3),zen(3)
    
    REAL(C_DOUBLE), allocatable, dimension(:, :, :) :: a1,a1o
    integer :: zsz(3)

    zsz = zen-zst+1
    allocate (a1(zsz(1), zsz(2), zsz(3)))
    allocate (a1o(zsz(1), zsz(2), zsz(3)))
    
    plan1 = fftw_plan_many_r2r( 1, zsz(3), zsz(1) * zsz(2),&
         a1,  zsz(3), zsz(1) * zsz(2), 1, &
         a1o, zsz(3), zsz(1) * zsz(2), 1, &
         [ FFTW_REDFT00 ] , plan_type)
    
    deallocate (a1)
    deallocate (a1o)
    
    return
  end subroutine dct_1m_z_plan  

  subroutine r2r_1m_x(in,out, plan1)
    
    implicit none
    
    real(kind=8), dimension(:, :, :),allocatable  :: in
    real(kind=8), dimension(:, :, :),allocatable  :: out
    type(C_PTR), intent(IN) :: plan1
    
    call fftw_execute_r2r(plan1, in,out)
    return
  end subroutine r2r_1m_x


  subroutine r2r_1m_y(in,out, plan1)
    
    implicit none
    
    real(kind=8), dimension(:, :, :),allocatable  :: in
    real(kind=8), dimension(:, :, :),allocatable  :: out
    type(C_PTR), intent(IN) :: plan1
    
    integer :: k,ks,ke,is,ie,js,je

    is = lbound(in, 1)
    ie = ubound(in, 1) 

    js = lbound(in, 2)
    je = ubound(in, 2) 

    
    ks = lbound(in, 3)
    ke = ubound(in, 3) 
    
    ! transform on one Z-plane at a time
    do k = ks, ke
       call fftw_execute_r2r(plan1, in(:, :, k), out(:,:, k))
    end do
    
    return
  end subroutine r2r_1m_y
   

  
  
  subroutine c2c_1m_x(inout, plan1)

    implicit none
    
    complex(kind=8), dimension(:, :, :), intent(INOUT) :: inout
    type(C_PTR) :: plan1
    
    call fftw_execute_dft(plan1, inout, inout)
    inout = inout*norm_factor(1)
    return
  end subroutine c2c_1m_x
  
  subroutine c2c_1m_y(inout, plan1)
    
      implicit none

      complex(kind=8), dimension(:, :, :), intent(INOUT) :: inout
      type(C_PTR), intent(IN) :: plan1
      
      integer :: k,ks,ke

      ks = lbound(inout, 3)
      ke = ubound(inout, 3) 
      
      ! transform on one Z-plane at a time
      do k = ks, ke
         call dfftw_execute_dft(plan1, inout(:, :, k), inout(:, :, k))
      end do
      inout = inout*norm_factor(2)
      return
   end subroutine c2c_1m_y
   
   subroutine c2c_1m_z(inout, plan1)
     
      implicit none
      
      complex(kind=8), dimension(:, :, :), intent(INOUT) :: inout
      type(C_PTR), intent(IN) :: plan1
      
      call dfftw_execute_dft(plan1, inout, inout)
      inout = inout*norm_factor(3)
      return
   end subroutine c2c_1m_z
   
   
  subroutine c2c_1m_x_plan(plan1,xst,xen, isign)
    
    implicit none

    type(C_PTR) :: plan1
    integer :: xst(3),xen(3)
    integer, intent(IN) :: isign
    
    complex(C_DOUBLE_COMPLEX), pointer :: a1(:, :, :)
    complex(C_DOUBLE_COMPLEX), pointer :: a1o(:, :, :)

    type(C_PTR) :: a1_p
    integer(C_SIZE_T) :: sz
    integer :: xsz(3)

    xsz = xen-xst+1
    sz = xsz(1)*xsz(2)*xsz(3)
    a1_p = fftw_alloc_complex(sz)
    call c_f_pointer(a1_p, a1 , [xsz(1),xsz(2),xsz(3)] )
    call c_f_pointer(a1_p, a1o, [xsz(1),xsz(2),xsz(3)] )
    
    plan1 = fftw_plan_many_dft(&
         1, xsz(1), xsz(2) * xsz(3), &
         a1 , xsz(1), 1, xsz(1), &
         a1o, xsz(1), 1, xsz(1), &
         isign, plan_type)
    call fftw_free(a1_p)
    
    return
  end subroutine c2c_1m_x_plan
  
  subroutine c2c_1m_y_plan(plan1, yst,yen, isign)
    use, intrinsic :: iso_c_binding
    implicit none
    
    type(C_PTR) :: plan1
    integer :: yst(3),yen(3)
    integer, intent(IN) :: isign
    
    complex(C_DOUBLE_COMPLEX), pointer :: a1(:, :)
    complex(C_DOUBLE_COMPLEX), pointer :: a1o(:, :)
    
    type(C_PTR) :: a1_p
    integer(C_SIZE_T) :: sz
    integer :: ysz(3)

    ysz = yen-yst+1
    sz = ysz(1) * ysz(2)
    a1_p = fftw_alloc_complex(sz)

    call c_f_pointer(a1_p, a1 , [ysz(1),ysz(2)])
    call c_f_pointer(a1_p, a1o, [ysz(1),ysz(2)])

    plan1 = fftw_plan_many_dft(1, ysz(2), ysz(1), &
         a1 , ysz(2), ysz(1), 1, &
         a1o, ysz(2), ysz(1), 1, &
         isign, plan_type)

    call fftw_free(a1_p)
    
    return
  end subroutine c2c_1m_y_plan
  
  
  
  subroutine c2c_1m_z_plan(plan1,zst,zen, isign)

    implicit none
    
    type(C_PTR), intent(OUT) :: plan1
    integer, intent(IN) :: isign
    integer, intent(IN) :: zst(3),zen(3)
    
    complex(kind=8), allocatable, dimension(:, :, :) :: a1
    integer :: zsz(3)

    zsz = zen-zst+1
    allocate (a1(zsz(1), zsz(2), zsz(3)))
    
    call dfftw_plan_many_dft(plan1, 1, zsz(3), zsz(1) * zsz(2),&
         a1, zsz(3), zsz(1) * zsz(2), 1, &
         a1, zsz(3), zsz(1) * zsz(2), 1, &
         isign, plan_type)
    
    deallocate (a1)
    
    return
  end subroutine c2c_1m_z_plan

  integer function half_complex_wp(i,n)
    implicit none
    integer :: i,n
      
    if (i<=n/2+1) then
       half_complex_wp = mod(i-1,n)
    else
       half_complex_wp = n-i+1
    end if
  end function half_complex_wp



  subroutine r2r_1m_x_plan(plan1, xst,xen,sign )
    
    implicit none
    type(C_PTR), intent(OUT) :: plan1
    integer, intent(IN) :: xst(3),xen(3)
    integer :: xsz(3)
    integer :: sign
    integer(C_FFTW_R2R_KIND), dimension(1) :: kind
    
    real(kind=8), allocatable, dimension(:, :, :) :: a1
    real(kind=8), allocatable, dimension(:, :, :) :: a2

    xsz = xen-xst+1
    allocate (a1( xsz(1), xsz(2), xsz(3) ) )
    allocate (a2( xsz(1), xsz(2), xsz(3) ) )

    if (sign==1) then
       kind(1) = FFTW_R2HC
    else
       kind(1) = FFTW_HC2R
    end if

    plan1 = fftw_plan_many_r2r(1,[xsz(1)],xsz(2) * xsz(3) ,&
         a1,  xsz(1), 1,  xsz(1), &
         a2,  xsz(1), 1,  xsz(1), kind, plan_type)
    deallocate (a1, a2)
    
    return
  end subroutine r2r_1m_x_plan


  subroutine r2r_1m_y_plan(plan1, yst,yen,sign )
    
    implicit none
    type(C_PTR), intent(OUT) :: plan1
    integer, intent(IN) :: yst(3),yen(3)
    integer :: ysz(3)
    integer :: sign
    integer(C_FFTW_R2R_KIND), dimension(1) :: kind

    REAL(C_DOUBLE_COMPLEX), pointer :: a1(:, :)
    REAL(C_DOUBLE_COMPLEX), pointer :: a1o(:, :)
    
    type(C_PTR) :: a1_p
    integer(C_SIZE_T) :: sz
    

    ysz = yen-yst+1
    sz = ysz(1) * ysz(2)
    a1_p = fftw_alloc_real(sz)
    call c_f_pointer(a1_p, a1 , [ysz(1),ysz(2)])
    call c_f_pointer(a1_p, a1o, [ysz(1),ysz(2)])
    
    if (sign==1) then
       kind(1) = FFTW_R2HC
    else
       kind(1) = FFTW_HC2R
    end if
    plan1 = fftw_plan_many_r2r(1,[ysz(2)],ysz(1) ,&
         a1,  ysz(2), ysz(1), 1, &
         a1o,  ysz(2), ysz(1), 1, kind, plan_type)
    call fftw_free(a1_p)
    
    return
  end subroutine r2r_1m_y_plan



  
  
end module M_FOURIER_TRANSFORM
