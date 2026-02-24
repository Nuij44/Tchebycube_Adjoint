module m_solver_diag_cart_hhh
  use M_FOURIER_TRANSFORM
    
  implicit none
  ! nu
  type t_solver_diag_cart_hhh
     
     integer :: n(3)
     integer :: xst(3),xen(3)  ! x-distributed array size 
     integer :: yst(3),yen(3)  ! y-distributed array size 
     integer :: zst(3),zen(3)  ! z-distributed array size 
     
     REAL(KIND=8) :: nu(3), sigma ! diffusion coefficient
     
     
     COMPLEX(KIND=8), dimension(:)  , allocatable :: LBDX ! 
     
     COMPLEX(KIND=8), dimension(:)  , allocatable :: LBDY !
     
     COMPLEX(KIND=8), dimension(:)  , allocatable :: LBDZ
     
   contains

     PROCEDURE :: SOLVE => IFCE_SOLVE
     PROCEDURE :: SOLVE_POISSON => IFCE_SOLVE_POISSON
     PROCEDURE :: INITIALISE => IFCE_INITIALISE
     PROCEDURE :: SET_PARAMS => SOLVER_DIAG_CART_HHH_SET_PARAMS

     
  end type t_solver_diag_cart_hhh
  
  REAL(KIND=8), dimension(:,:,:), allocatable :: LBD
  
contains

  ! common interface for 


  subroutine ifce_initialise(this,grid,opx,opy,opz,ph)
    use m_mesh_base
    use m_operator_tcheby
    use m_operator_fourier_dft
    use decomp_2d
    implicit none
    class(t_solver_diag_cart_hhh):: this
    type(t_mesh_base)            :: grid(3)
    type(t_operator_fourier_dft) :: opx,opy,opz    
    real(kind=8) :: sigma
    real(kind=8) :: nu(3),bcl(3,2),bcr(3,2)
    integer :: n(3)
    type(decomp_info) :: ph

    ! get the global size. ni := nx+1 
    n = [ph%xen(1)-ph%xst(1),ph%yen(2)-ph%yst(2),ph%zen(3)-ph%zst(3)]
    BCL = 0
    BCR = 0
    THIS%N = N
    
    CALL SOLVER_DIAG_HHH_BUILD (THIS, &
         GRID(1)%X, THIS%NU(1), THIS%N(1), PH%XST,  PH%XEN, &
         GRID(2)%X, THIS%NU(2), THIS%N(2), PH%YST,  PH%YEN, &
         GRID(3)%X, THIS%NU(3), THIS%N(3), PH%ZST,  PH%ZEN)
    
  end subroutine ifce_initialise

  
  
  subroutine ifce_solve(this,fi,sfi,sigma,ph)
    use decomp_2d
    implicit none
    class(t_solver_diag_cart_hhh) :: this
    real(kind=8),dimension(:,:,:),allocatable :: fi,sfi
    real(kind=8) :: sigma
    integer :: n(3)
    type(decomp_info) :: ph


    call solver_diag_cart_hhh_solve(fi,sfi,sigma,&
         this%lbdx, this%n(1), ph%xst, ph%xen, &
         this%lbdy, this%n(2), ph%yst, ph%yen, &
         this%lbdz, this%n(3), ph%zst, ph%zen)
      
  end subroutine ifce_solve


  
  SUBROUTINE SOLVER_DIAG_CART_HHH_SET_PARAMS(THIS,NU,SIGMA )
    IMPLICIT NONE
    CLASS(T_SOLVER_DIAG_CART_HHH) :: THIS
    REAL(KIND=8) :: NU(3)
    REAL(KIND=8) :: SIGMA

    THIS%NU = NU
    THIS%SIGMA = SIGMA
    
  END SUBROUTINE SOLVER_DIAG_CART_HHH_SET_PARAMS


  subroutine solver_diag_hhh_build(this,&
       x, nux, nx, xst, xen, &
       y, nuy, ny, yst, yen, &
       z, nuz, nz, zst, zen)
    use decomp_2d_mpi
    implicit none
    class(t_solver_diag_cart_hhh) :: this
    integer :: nx,xst(3),xen(3)
    REAL(KIND=8) :: x(nx+1),nux
    
    integer :: ny,yst(3),yen(3)
    REAL(KIND=8) :: y(ny+1),nuy
      
    integer :: nz,zst(3),zen(3)
    REAL(KIND=8) :: z(nz+1),nuz

      
    real(kind=8) :: pi=acos(-1.)
    real(kind=8) :: xmin,xmax,hx,lx,wx
    real(kind=8) :: ymin,ymax,hy,ly,wy
    real(kind=8) :: zmin,zmax,hz,lz,wz

    integer :: nxp,nyp,nzp
    
    integer :: i,im,jm,j,k,km

    this%n(1) = nx ; this%xst = xst ; this%xen = xen ;
    this%n(2) = ny ; this%yst = yst ; this%yen = yen ;
    this%n(3) = nz ; this%zst = zst ; this%zen = zen ;

    ! wavenumber in the direction-z := \varphi
    nxp = nx+1
    hx = x(2)-x(1)
    xmin = x(1)
    xmax = x(nx+1)+hx
    Lx = xmax-xmin 
    wx = (2*pi)/Lx

    nyp = ny+1
    hy = y(2)-y(1)
    ymin = y(1)
    ymax = y(ny+1)+hy
    Ly = ymax-ymin 
    wy = (2*pi)/Ly
    
    nzp = nz + 1
    hz = z(2)- z(1)
    zmin = z(1)
    zmax = z(nz+1)+hz
    Lz = zmax-zmin 
    wz = (2*pi)/Lz

    this%nu=[nux,nuy,nuz]

    allocate( this%lbdx( 1:nx+1 ) ); this%lbdx = 0
    do i = 1 , nx+1
       im = half_complex_wp(i,nx+1)
       this%lbdx(i) = - nux*(im*wx)**2
    end do    

    allocate( this%lbdy( 1:ny+1 ) ); this%lbdy = 0
    do j = 1 , ny+1
       jm = half_complex_wp(j,ny+1)
       this%lbdy(j) = - nuy*(jm*wy)**2
    end do

    allocate( this%lbdz( 1:nz+1 ) ); this%lbdz = 0
    do k = 1 , nz+1
       km = half_complex_wp(k,nz+1)
       this%lbdz(k) = - nuz*(km*wz)**2
    end do

    if (.not. allocated(LBD)) then
       allocate( LBD(zst(1):zen(1),zst(2):zen(2),zst(3):zen(3) ) )
    end if

  end subroutine solver_diag_hhh_build
  
  subroutine solver_diag_cart_hhh_solve(fi,sfi,sigma,&
       lbdx, nx, xst, xen, &
       lbdy, ny, yst, yen, &
       lbdz, nz, zst, zen)
    use m_fourier_transform
    use decomp_2d
    implicit none
    real(kind=8) :: sigma
    integer :: nx,xst(3),xen(3)
    integer :: ny,yst(3),yen(3)
    integer :: nz,zst(3),zen(3)
    real(kind=8), dimension(xst(1):xen(1),xst(2):xen(2),xst(3):xen(3) ) :: fi,sfi
    complex(kind=8), dimension( 1:nx+1 ) :: lbdx 
    complex(kind=8), dimension( 1:ny+1 ) :: lbdy
    complex(kind=8), dimension( 1:nz+1 ) :: lbdz
    
    integer :: i,j,k

    
    forall(i=xst(1):xen(1),j=xst(2):xen(2),k=xst(3):xen(3))
       dg_x(i,j,k) = sfi(i,j,k)
    end forall

    forall(i=zst(1):zen(1),j=zst(2):zen(2),k=zst(3):zen(3))
       LBD(I,J,K) = ( lbdx(i) + lbdy(j) +  lbdz(k) + sigma )**(-1)
    end forall


    
    call c2c_1m_x(dg_x,plan_fwd_x)
    call transpose_x_to_y(dg_x, zg1_y)
    call c2c_1m_y(zg1_y,plan_fwd_y)
    call transpose_y_to_z(zg1_y,zg1_z)
    call c2c_1m_z(zg1_z,plan_fwd_z)

    
    forall(i=zst(1):zen(1),j=zst(2):zen(2),k=zst(3):zen(3))
       zg1_z(i,j,k) = zg2_z(i,j,k) * LBD(i,j,k)
    end forall
    
    
    call c2c_1m_z(zg1_z,plan_bck_z)
    call transpose_z_to_y(zg2_z,zg1_y)
    call c2c_1m_y(zg1_y,plan_bck_y)
    call transpose_y_to_x(zg1_y,dg_x)
    call c2c_1m_x(dg_x,plan_bck_x)
    
    forall(i=xst(1):xen(1),j=xst(2):xen(2),k=xst(3):xen(3))
       fi(i,j,k) = dg_x(i,j,k)
    end forall
  end subroutine solver_diag_cart_hhh_solve

  subroutine ifce_solve_poisson(this,fi,sfi,nu,ph)
    use decomp_2d
    implicit none
    class(t_solver_diag_cart_hhh) :: this
    real(kind=8),dimension(:,:,:),allocatable :: fi,sfi
    real(kind=8) :: nu(3)
    integer :: n(3)
    type(decomp_info) :: ph


    call solver_diag_cart_hhh_solve_poisson(fi,sfi,&
         nu(1),this%lbdx, this%n(1), ph%xst, ph%xen, &
         nu(2),this%lbdy, this%n(2), ph%yst, ph%yen, &
         nu(3),this%lbdz, this%n(3), ph%zst, ph%zen)
      
  end subroutine ifce_solve_poisson


  subroutine solver_diag_cart_hhh_solve_poisson(fi,sfi,&
       nux,lbdx, nx, xst, xen, &
       nuy,lbdy, ny, yst, yen, &
       nuz,lbdz, nz, zst, zen)
    use m_fourier_transform
    use decomp_2d
    implicit none
    real(kind=8) :: nux,nuy,nuz
    integer :: nx,xst(3),xen(3)
    integer :: ny,yst(3),yen(3)
    integer :: nz,zst(3),zen(3)
    real(kind=8), dimension(xst(1):xen(1),xst(2):xen(2),xst(3):xen(3) ) :: fi,sfi
    complex(kind=8), dimension( 1:nx+1 ) :: lbdx 
    complex(kind=8), dimension( 1:ny+1 ) :: lbdy
    complex(kind=8), dimension( 1:nz+1 ) :: lbdz
    
    integer :: i,j,k

    
    dg_x = sfi 
    
    DO K = ZST(3),ZEN(3)
       DO J = ZST(2),ZEN(2) 
          DO I=ZST(1),ZEN(1)
             LBD(I,J,K) = ( nux*lbdx(i) + nuy*lbdy(j) +  nuz*lbdz(k)  )**(-1)
          END DO
       END DO
    END DO
    
    
    call c2c_1m_x(dg_x,plan_fwd_x)
    call transpose_x_to_y(dg_x, zg1_y)
    call c2c_1m_y(zg1_y,plan_fwd_y)
    call transpose_y_to_z(zg1_y,zg1_z)
    call c2c_1m_z(zg1_z,plan_fwd_z)

    forall(i=zst(1):zen(1),j=zst(2):zen(2),k=zst(3):zen(3))
       zg1_z(i,j,k) = zg2_z(i,j,k) * LBD(i,j,k)
    end forall
    
    call c2c_1m_z(zg1_z,plan_bck_z)
    call transpose_z_to_y(zg2_z,zg1_y)
    call c2c_1m_y(zg1_y,plan_bck_y)
    call transpose_y_to_x(zg1_y,dg_x)
    call c2c_1m_x(dg_x,plan_bck_x)
    
    forall(i=xst(1):xen(1),j=xst(2):xen(2),k=xst(3):xen(3))
       fi(i,j,k) = dg_x(i,j,k)
    end forall
  end subroutine solver_diag_cart_hhh_solve_poisson



end module m_solver_diag_cart_hhh
