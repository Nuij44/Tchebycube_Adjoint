module m_quadrature
  use decomp_2d
  use m_numerics
  use mpi
  use m_fourier_transform
  implicit none
  type t_quadrature
     real(kind=8)::xmin(3),xmax(3)
     real(kind=8),allocatable::quad_x(:),quad_y(:),quad_z(:),mean(:),cheb_coeffs(:)
     integer :: rank
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

    call mpi_comm_rank(mpi_comm_world,this%rank,i)
    
    
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
    
    CALL TRANSPOSE_X_TO_Y(zg_X,zg1_Y)
    call c2c_1m_y(zg1_y,plan_fwd_y)
    zg2_Y = 0._DP
    zg2_Y(:,1,:) = zg1_Y(:,1,:)

    CALL TRANSPOSE_Y_TO_Z(zg2_Y,zg1_Z)
    
    !Recherche du rang du processus possedant zg1_Z(1,1,:)

    rank_mode0 = 0

    if (PH%ZST(1)==1 .AND. PH%ZST(2)==1) then
       rank_mode0 = this%rank
       this%mean(1:nr+1) = zg1_Z(1,1,1:NR+1)/sqrt(dble((na+1)*(nz+1)))
    end if

    CALL MPI_ALLREDUCE(MPI_IN_PLACE,rank_mode0,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERR)

    call MPI_Bcast(this%mean, nr+1, MPI_DOUBLE_PRECISION,rank_mode0,MPI_COMM_WORLD,ierr)

    ! ----------------------------------------------------
    ! Etape 2 : Coefficients de Tchebychev de mean(r)
    !           via DCT-II (formule directe)
    ! ----------------------------------------------------

    do n = 0, NR
       this%cheb_coeffs(n+1) = 0.0d0
       do k = 1, NR+1
          this%cheb_coeffs(n+1) = this%cheb_coeffs(n+1) + this%mean(k) * cos(PI * dble(n) * dble(k-1) / dble(NR))  ! ← dble(NR) correct
       end do
       this%cheb_coeffs(n+1) = this%cheb_coeffs(n+1) * 2.0d0 / dble(NR)  ! ← dble(NR) correct
    end do

    ! Correction du mode 0
    this%cheb_coeffs(1) = this%cheb_coeffs(1) * 0.5d0
    this%cheb_coeffs(NR+1) = this%cheb_coeffs(NR+1) * 0.5d0

  ! ----------------------------------------------------------------
  ! Etape 3 : Intégration Tchebychev
  !           ∫_{-1}^{1} T_n(ξ) dξ = 2/(1-n²)  si n pair, 0 sinon
  ! ----------------------------------------------------------------
  LZ = 0.5d0 * (xmax(3) - xmin(3))    ! jacobien dξ → dz

  integral_z = 0.0d0
  do n = 0, Nr
    if (mod(n,2) == 0) then     ! modes pairs uniquement
      if (n == 0) then
        weight = 2.0d0          ! ∫T_0 dξ = 2
      else
        weight = 2.0d0 / (1.0d0 - dble(n*n))
      end if
      integral_z = integral_z + this%cheb_coeffs(n+1) * weight
    end if
  end do

  ! ----------------------------------------------------------------
  ! Etape 4 : Assemblage final
  !           I = Lx * Ly * Lz * ∫ mean_xy(ξ) dξ
  ! ----------------------------------------------------------------
  Lx = xmax(1) - xmin(1)
  LY = xmax(2) - xmin(2)
  
  vol = Lx * Ly * Lz * integral_z
    
  end subroutine integrate_spec
  
end module m_quadrature
