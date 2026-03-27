module m_quadrature
  use decomp_2d
  use m_numerics
  use mpi
  use m_fourier_transform
  implicit none
  type t_quadrature
     real(kind=8)::xmin(3),xmax(3)
     real(kind=8),allocatable::quad_x(:),quad_y(:),quad_z(:),mean(:),cheb_coeffs(:),clenshaw_curtis_w(:,:,:)
     integer :: rank
     REAL(DP),ALLOCATABLE,DIMENSION(:,:,:) :: DGY,DGZ
     real,ALLOCATABLE,dimension(:,:) :: wkU
  end type t_quadrature
  
contains

  SUBROUTINE INIT_quadrature_hhi(this,xmin,xmax,xst,xen,nx,yst,yen,ny,zst,zen,nz)
    use m_fourier_transform
    use m_tensor_product
    implicit none
    !.. args
    class(t_quadrature),target :: this
    real(kind=8) :: xmin(3),xmax(3)

    integer :: xst(3),xen(3),nx
    integer :: yst(3),yen(3),ny
    integer :: zst(3),zen(3),nz

    integer :: i,j,k
    REAL(KIND=8),ALLOCATABLE,DIMENSION(:) :: mean,cheb_coeffs
    REAL(DP) :: wk,dn

    call mpi_comm_rank(mpi_comm_world,this%rank,i)

    call alloc_y(this%dgy, OPT_GLOBAL=.TRUE.)
    call alloc_z(this%dgz, OPT_GLOBAL=.TRUE.)

    allocate(this%wkU(ZST(1):ZEN(1),ZST(2):ZEN(2)))
    
    ALLOCATE(this%mean(1:nz+1),this%cheb_coeffs(1:Nz+1))


    ALLOCATE(THIS%QUAD_X(XST(1):XEN(1)))
    ALLOCATE(THIS%QUAD_Y(YST(2):YEN(2)))
    ALLOCATE(THIS%QUAD_Z(ZST(3):ZEN(3)))

    ! natural quadrature rule on [+1,-1]
    do k=zst(3),zen(3)
       if (mod(k-1,2)==0) then
          THIS%quad_z(k) = 2./((k-1)**2-1)
       else
          THIS%quad_z(k) = 0
       end if
    end do
    
    ! scaling form fftw computation of  
    THIS%QUAD_Z(ZST(3)) = THIS%QUAD_Z(ZST(3))*0.5
    DO K=ZST(3),ZEN(3)
       THIS%QUAD_Z(K) = THIS%QUAD_Z(K)/NZ
    END DO

    ! scaling from [+1,-1] -> [z_min,z_max]
    ! Z(z) = 1-2*(z-z_min)/(z_max-z_min)
    ! dZ = -2*dz/(z_max-z_min)
    ! dz = -(z_max-z_min)/2*dZ
    DO K=ZST(3),ZEN(3)
       THIS%QUAD_Z(K) = THIS%QUAD_Z(K)*(-(XMAX(3)-XMIN(3))*0.5)
    END DO

    DO I=XST(1),XEN(1)
       THIS%QUAD_X(I) = (XMAX(1)-XMIN(1))/(NX+1)
    END DO

    DO J=YST(2),YEN(2)
       THIS%QUAD_Y(J) = (XMAX(2)-XMIN(2))/(NY+1)
    END DO
    

!    THIS%QUAD_X = THIS%QUAD_X/(XMAX(1)-XMIN(1))
!    THIS%QUAD_Y = THIS%QUAD_Y/(XMAX(2)-XMIN(2))
!    THIS%QUAD_Z = THIS%QUAD_Z/(XMAX(3)-XMIN(3))


    ! CLENSHAW-CURTS WEIGHT
    call alloc_z(this%clenshaw_curtis_w, OPT_GLOBAL=.TRUE.)

    DO k = ZST(3),ZEN(3)
       wk = 0._DP
       do i = ZST(3),ZEN(3),2
          if (i == zst(3) .OR. i == zen(3)) then
             dn = 1._DP
          else
             dn = 2._DP
          end if
          wk = wk + dn/(1._DP - (i - 1)**2) * cos((i-1)*(k-1)*pi/(ZEN(3) - 1))
       end do
       if (k == ZST(3) .OR. k == ZEN(3)) then
          wk = wk / (ZEN(3)-1)
       else
          wk = 2._DP * wk / (ZEN(3)-1)
       end if

       this%clenshaw_curtis_w(:,:,k) = wk
       
    end DO

    
  end subroutine INIT_quadrature_hhi
  
  subroutine get_quadrature_hhi(this,fi,quad,xst,xen,nx,yst,yen,ny,zst,zen,nz)
    use m_fourier_transform
    use m_tensor_product
    use decomp_2d
    implicit none
    type(t_quadrature) :: this
    real(kind=8),allocatable ::fi(:,:,:), quad(:,:,:)
    integer :: xst(3),xen(3),nx
    integer :: yst(3),yen(3),ny
    integer :: zst(3),zen(3),nz
    real(kind=8) :: tmp
    integer ::i,j,k 

    
    forall(j=xst(2):xen(2),k=xst(3):xen(3))
       DG1_X(1,J,K) = SUM(FI(1:NX+1,j,k)*THIS%QUAD_X(1:NX+1))
    end forall
    forall(i=xst(1):xen(1),j=xst(2):xen(2),k=xst(3):xen(3))
       DG2_X(i,J,K) = DG1_X(1,J,K) 
    end forall
    CALL TRANSPOSE_X_TO_Y(DG2_X,DG2_Y)

    forall(i=yst(1):yen(1),k=yst(3):yen(3))
       DG1_Y(I,1,K) = SUM(DG2_Y(I,1:NY+1,k)*THIS%QUAD_Y(1:NY+1))
    end forall
    forall(i=yst(1):yen(1),j=yst(2):yen(2),k=yst(3):yen(3))
       DG2_Y(i,J,K) = DG1_Y(I,1,K) 
    end forall
    
    CALL TRANSPOSE_Y_TO_Z(DG2_Y,DG1_Z)
    
    CALL FFTW_EXECUTE_R2R(DCT_FWD_Z,DG1_Z,DG2_Z)
    FORALL(I=ZST(1):ZEN(1),J=ZST(2):ZEN(2))
       DG1_Z(I,J,1) = SUM( DG2_Z(I,J,1:NZ+1) * THIS%QUAD_Z(1:NZ+1) )
    END FORALL
    FORALL(I=ZST(1):ZEN(1),J=ZST(2):ZEN(2),K=ZST(3):ZEN(3))
       DG2_Z(I,J,K) = DG1_Z(I,J,1) 
    END FORALL
    CALL TRANSPOSE_Z_TO_Y(DG2_Z,DG1_Y)
    CALL TRANSPOSE_Y_TO_X(DG1_Y, QUAD)
        
  end subroutine get_quadrature_hhi
  
  subroutine normalize(this,U,U_tilde,PH,N)
    type(t_quadrature) :: this
    TYPE(DECOMP_INFO) :: ph
    real(kind=8),allocatable ::U(:,:,:), U_TILDE(:,:,:)
    INTEGER, DIMENSION(3) :: N

    REAL(DP),allocatable :: DG1(:,:,:)
    REAL(DP) :: norm

    INTEGER :: i,j,k

    U_tilde = U**2

    call alloc_x(DG1 , OPT_GLOBAL=.TRUE.) ; DG1 = 0

    CALL get_quadrature_hhi(this,u_tilde,DG1,ph%xst,ph%xen,n(1),ph%yst,ph%yen,n(2),ph%zst,ph%zen,n(3))

    norm = DG1(PH%XST(1),PH%XST(2),PH%XST(3))

    U_tilde = norm*U
  END subroutine normalize
  

  subroutine integrate_spec(this,UA,VOL,PH,NA,NZ,NR,xmax,xmin)
    type(t_quadrature) :: this
    real(kind=8),allocatable,intent(in) ::UA(:,:,:)
    real(kind=8) :: vol
    integer :: na,nz,nr
    REAL(kind=8),dimension(3) :: xmax,xmin
    integer :: i,j,k,n,rank_mode0
    type(decomp_info) :: ph
    REAL(kind=8) :: weight,integral_z,LX,LY,LZ
    integer :: ierr
    
    DG_x = UA
    call c2c_1m_x(dg_x,plan_fwd_x)
    zg_X = 0._DP
    zg_X(1,:,:) = DG_X(1,:,:)
    call c2c_1m_x(zg_x,plan_bck_x)
    
    CALL TRANSPOSE_X_TO_Y(zg_X,zg1_Y)
    call c2c_1m_y(zg1_y,plan_fwd_y)
    zg2_Y = 0._DP
    zg2_Y(:,1,:) = zg1_Y(:,1,:)
    call c2c_1m_y(zg2_y,plan_bck_y)

    CALL TRANSPOSE_Y_TO_Z(zg2_Y,zg1_Z)

    CALL integrate_spec_r(this,ZG1_Z,ZG2_Z,xmin,xmax,na,nz,nr,ph)


    integral_z = zg2_z(PH%zst(1),PH%zst(2),PH%zst(3))

    ! ----------------------------------------------------------------
    ! Etape 4 : Assemblage final
    !           I = Lx * Ly * Lz * ∫ mean_xy(ξ) dξ
    ! ----------------------------------------------------------------
    Lx = xmax(1) - xmin(1)
    LY = xmax(2) - xmin(2)
  
    vol = Lx * Ly * integral_z
  end subroutine integrate_spec
  
  
  subroutine integrate_spec_r(this,U,vol_az,xmin,xmax,na,nz,nr,ph)
    type(t_quadrature) :: this
    COMPLEX(kind=8), allocatable :: U(:,:,:)
    COMPLEX(kind=8), allocatable :: vol_az(:,:,:)
    integer :: na, nz, nr
    REAL(kind=8), dimension(3) :: xmax, xmin
    integer :: i,j,k
    type(decomp_info) :: ph

!    call transpose_x_to_y(U,this%DGY)
!    call transpose_y_to_z(this%DGY,this%DGZ)
    this%wkU = 0._DP
    do k = 1, nr+1
       this%wkU = this%wkU + U(:,:,K)*this%clenshaw_curtis_w(:,:,K)
    end do

    do i = PH%ZST(1),PH%ZEN(1)
       do j = PH%ZST(2),PH%ZEN(2)
          VOL_AZ(i,j,:) = this%wkU(i,j) * 0.5_DP * (Xmax(3) - Xmin(3))
       end do
    end do

    
!    call transpose_z_to_y(this%DGZ,this%DGY)
!    call transpose_y_to_x(this%DGY,vol_az)


    
  end subroutine integrate_spec_r
  
end module m_quadrature
