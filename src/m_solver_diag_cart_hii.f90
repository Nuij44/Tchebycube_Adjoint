module m_solver_diag_cart_hii
  use M_FOURIER_TRANSFORM
    
  implicit none
  
  type t_solver_diag_cart_hii
     
     integer :: n(3)
     integer :: xst(3),xen(3)  ! x-distributed array size 
     integer :: yst(3),yen(3)  ! y-distributed array size 
     integer :: zst(3),zen(3)  ! z-distributed array size 
     
     REAL(KIND=8) :: nu(3) ! diffusion coefficient
     
     REAL(KIND=8) :: alf_l_y, bet_l_y 
     REAL(KIND=8) :: alf_r_y, bet_r_y

     REAL(KIND=8) :: alf_l_z, bet_l_z 
     REAL(KIND=8) :: alf_r_z, bet_r_z  

     
     REAL(KIND=8), dimension(:,:)  , allocatable :: bc_val_l_z, bc_val_r_z
     REAL(KIND=8), dimension(:,:)  , allocatable :: bc_val_l_y, bc_val_r_y
     
     COMPLEX(KIND=8), dimension(:)  , allocatable :: LBDX ! 

     COMPLEX(KIND=8), dimension(:)  , allocatable :: LBDY !
     COMPLEX(KIND=8), dimension(:,:), allocatable :: PY,PYM1
     REAL(KIND=8)   , dimension(:)  , allocatable :: FLY,FRY
     COMPLEX(KIND=8), dimension(:,:), allocatable :: CCY
     REAL(KIND=8) :: cly(2,2)

     
     COMPLEX(KIND=8), dimension(:)  , allocatable :: LBDZ
     COMPLEX(KIND=8), dimension(:,:), allocatable :: PZ,PZM1
     REAL(KIND=8)   , dimension(:)  , allocatable :: FLZ,FRZ
     COMPLEX(KIND=8), dimension(:,:)  , allocatable :: CCZ

     REAL(KIND=8) :: clz(2,2)
     
  end type t_solver_diag_cart_hii
  
  COMPLEX(KIND=8), dimension(:,:,:), allocatable :: dg_cl_z
  COMPLEX(KIND=8), dimension(:,:,:), allocatable :: dg_cl_y
  REAL(KIND=8), dimension(:,:,:), allocatable :: dg_LBD
  
contains

  ! common interface for 
  
  
  subroutine build(this,grid,opx,opy,opz,sigma,nu,bcl,bcr,ph,n)
    use m_mesh_base
    use m_operator_tcheby
    use m_operator_fourier
    use decomp_2d
    implicit none
    class(t_solver_diag_cart_hii) :: this
    type(t_mesh_base)        :: grid(3)
    type(t_operator_fourier) :: opx
    type(t_operator_tcheby)  :: opy,opz
    
    real(kind=8) :: sigma
    real(kind=8) :: nu(3),bcl(3,2),bcr(3,2)
    integer :: n(3)
    type(decomp_info) :: ph

    
    call solver_diag_hhi_build (this, &
         grid(1)%x, nu(1), n(1), ph%xst,  ph%xen, &
         grid(2)%x, nu(2), n(2), ph%yst,  ph%yen, opy%get_d2(), opy%get_d1(), bcl(2,:), bcr(2,:) , &
         grid(3)%x, nu(3), n(3), ph%zst,  ph%zen, opz%get_d2(), opz%get_d1(), bcl(3,:), bcr(3,:) )
    
  end subroutine build
  

  subroutine set_bval(this,grid,udf_bc_left_y,udf_bc_right_y,udf_bc_left_z,udf_bc_right_z,ph)
    use m_mesh_base
    use decomp_2d

    implicit none
    class(t_solver_diag_cart_hii) :: this
    type(t_mesh_base)        :: grid(3)
    PROCEDURE(UDF_TIMESPACE)      :: UDF_BC_LEFT_Y,UDF_BC_RIGHT_Y
    PROCEDURE(UDF_TIMESPACE)      :: UDF_BC_LEFT_Z,UDF_BC_RIGHT_Z
    integer :: n(3)
    type(decomp_info) :: ph
    
    call solver_diag_hhi_set_boundary_field(&
         grid(1)%x, grid(2)%x, grid(3)%x, ph%xst, ph%xen, this%n, &
         udf_bc_left_y, udf_bc_right_y, this%bc_val_l_y, this%bc_val_r_y,&
         udf_bc_left_z, udf_bc_right_z, this%bc_val_l_z, this%bc_val_r_z )
    
  end subroutine set_bval


  subroutine solve(this,fi,sfi,sigma,ph)
    use decomp_2d
    implicit none
    class(t_solver_diag_cart_hii) :: this
    real(kind=8),dimension(:,:,:),allocatable :: fi,sfi
    real(kind=8) :: sigma
    integer :: n(3)
    type(decomp_info) :: ph
    
    call solver_diag_cart_hii_solve(fi,sfi,sigma,&
         this%lbdx, this%n(1), ph%xst, ph%xen, &
         ! xxxx
         this%lbdy, this%n(2), ph%yst, ph%yen, this%py,  this%pym1, &
         this%cly, this%ccy,  this%fly,  this%fry, this%bc_val_l_y, this%bc_val_r_y,&
         ! xxxx
         this%lbdz, this%n(3), ph%zst, ph%zen, this%pz,  this%pzm1,  &
         this%clz, this%ccz,  this%flz,  this%frz, this%bc_val_l_z, this%bc_val_r_z)
      
  end subroutine solve
  !> basic routine 
  
  subroutine reduction_operator(D2,BL,BR,P,PM1,LBD,CL,CC,FL,FR,N)
    use m_numerics
    implicit none
    !... args in
    integer :: n
    real(kind=8),intent(in) :: d2(1:n+1,1:n+1)
    real(kind=8),intent(in) :: bl(1:n+1,1:n+1)
      real(kind=8),intent(in) :: br(1:n+1,1:n+1)
      !... args out
      real(kind=8),intent(out) :: P(1:N+1,1:N+1),PM1(1:N+1,1:N+1)
      real(kind=8),intent(out) :: LBD(1:N+1)
      real(kind=8),intent(out) :: FL(1:N+1),FR(1:N+1)
      real(kind=8),intent(out) :: CL(1:2,1:2),CC(1:2,1:N+1)
      !... variables
      real(kind=8) :: D2_(2:N,2:N)
      real(kind=8) :: P_(2:N,2:N),PM1_(2:N,2:N)
      real(kind=8) :: LBD_(2:N)
      
      integer :: i,j
      
      ! inverse of cl attention au det
      ! CL(1,1:2) = BL(0,0) , BL(0,N)
      ! CL(2,1:2) = BR(N,0) , BR(N,N)
      
      cl(1,1:2) = [   br(n+1,n+1) , - bl(1,n+1) ]
      cl(2,1:2) = [ - br(n+1,  1) ,   bl(1,  1) ] 
      cl = cl / ( bl(1,1)*br(n+1,n+1) - br(n+1,1)*bl(1,n+1)  )
      
      forall(i=2:n,j=2:n)
         d2_(i,j) = d2(i,j) &
              - d2(i,1  )*( cl(1,1)*bl(1,j) + cl(1,2)*br(n+1,j)  ) &
              - d2(i,n+1)*( cl(2,1)*bl(1,j) + cl(2,2)*br(n+1,j)  ) 
      end forall
      
      call diagonalise_real( d2_(2:n,2:n), p_(2:n,2:n), pm1_(2:n,2:n), lbd_(2:n), n-1 )
      
      
      p=0
      pm1=0
      !> 
      forall(i=1:n+1) p(i,i) = 1d0
      forall(i=1:n+1) pm1(i,i) = 1d0
      
      forall(i=2:n,j=2:n)
         p(i,j) = p_(i,j)
         pm1(i,j) = pm1_(i,j)
      end forall
      
      CC(1,  1) = 1
      CC(1,2:n) = -bl(1,2:n)
      CC(1,n+1) = 0
      
      CC(2,  1) = 0
      CC(2,2:n) = -br(n+1,2:n)
      CC(2,n+1) = 1
      
      
      forall(i=1:n+1) lbd(i) = 1e16 
      forall(i=2:n  ) lbd(i) = lbd_(i)
      
      forall(i=2:n)
         fl(i) = d2(i,1)*cl(1,1) + d2(i,n+1)*cl(2,1)
         fr(i) = d2(i,1)*cl(1,2) + d2(i,n+1)*cl(2,2)
      end forall
      
    end subroutine reduction_operator
    

    subroutine solver_diag_hhi_build(this,&
         x, nux, nx, xst, xen, &
         y, nuy, ny, yst, yen, d2y, dy, cl_left_y, cl_right_y , &
         z, nuz, nz, zst, zen, d2z, dz, cl_left_z, cl_right_z )
      use decomp_2d_mpi
      implicit none
      class(t_solver_diag_cart_hii) :: this
      integer :: nx,xst(3),xen(3)
      REAL(KIND=8) :: x(nx+1),nux
      
      integer :: ny,yst(3),yen(3)
      REAL(KIND=8) :: y(ny+1),nuy,d2y(1:ny+1,1:ny+1),dy(1:ny+1,1:ny+1)
      REAL(KIND=8) :: cl_left_y(2), cl_right_y(2)
      
      integer :: nz,zst(3),zen(3)
      REAL(KIND=8) :: z(nz+1),nuz,d2z(1:nz+1,1:nz+1),dz(1:nz+1,1:nz+1)
      REAL(KIND=8) :: cl_left_z(2), cl_right_z(2)
      
      real(kind=8) :: pi=acos(-1.)
      real(kind=8) :: xmin,xmax,hx,lx,wx
    real(kind=8) :: ymin,ymax,hy,ly,wy

    integer :: nxp,nyp,nzp
    
    integer :: i,im,jm,j,k

    real(kind=8),allocatable :: d2(:,:),bl(:,:),br(:,:),id(:,:)
    
    REAL(KIND=8), dimension(:,:), allocatable :: P,PM1
    REAL(KIND=8), dimension(:)  , allocatable :: LBD
    REAL(KIND=8), dimension(:)  , allocatable :: FL,FR
    REAL(KIND=8), dimension(:,:)  , allocatable :: CC
    REAL(KIND=8) :: cl(2,2)

    

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

    this%nu=[nux,nuy,nuz]
    
    !
    nyp = ny + 1
    
    allocate( p  ( 1:ny+1 , 1:ny+1 ) ) ; p   = 0
    allocate( pm1( 1:ny+1 , 1:ny+1 ) ) ; pm1 = 0
    allocate( lbd( 1:ny+1          ) ) ; lbd = 0 

    allocate( fl ( 1:ny+1         ) ) ; fl = 0
    allocate( fr ( 1:ny+1         ) ) ; fr = 0
    allocate( cc ( 1:2   , 1:ny+1 ) ) ; cc = 0
    
    allocate( id ( 1:ny+1 , 1:ny+1 ) ) ; id = 0
    allocate( d2 ( 1:ny+1 , 1:ny+1 ) ) ; d2 = 0
    allocate( bl ( 1:ny+1 , 1:ny+1 ) ) ; bl = 0
    allocate( br ( 1:ny+1 , 1:ny+1 ) ) ; br = 0
    
    id=0
    do i=1,ny+1
       id(i,i)=1
    end do
    d2 = 0

    bl = cl_left_y(1) * id + cl_left_y(2) * dy
    br = cl_right_y(1) * id + cl_right_y(2) * dy
    !
    
    d2 = nuy*d2y 
    call reduction_operator( &
         d2, bl, br, p(1:nyp,1:nyp), pm1(1:nyp,1:nyp), lbd(1:nyp), & 
         cl, cc, fl(1:nyp),fr(1:nyp), ny )
    
    allocate( this%py  ( 1:ny+1 , 1:ny+1 ) ) ; this%py = p
    allocate( this%pym1( 1:ny+1 , 1:ny+1 ) ) ; this%pym1 = pm1
    allocate( this%lbdy( 1:ny+1          ) ) ; this%lbdy = lbd

    allocate( this%fly ( 1:ny+1          ) ) ; this%fly = fl
    allocate( this%fry ( 1:ny+1          ) ) ; this%fry = fr
    allocate( this%ccy ( 1:2    , 1:ny+1 ) ) ; this%ccy = cc


    
    this%cly = cl

    deallocate( d2,bl,br,id )
    deallocate( p,pm1,lbd,fl,fr,cc )
    
    nzp = nz + 1
    
    allocate( p  ( 1:nz+1 , 1:nz+1 ) ) ; p = 0
    allocate( pm1( 1:nz+1 , 1:nz+1 ) ) ; pm1 = 0
    allocate( lbd( 1:nz+1          ) ) ; lbd = 0 

    allocate( fl ( 1:nz+1         ) ) ; fl = 0
    allocate( fr ( 1:nz+1         ) ) ; fr = 0
    allocate( cc ( 1:2   , 1:nz+1 ) ) ; cc = 0


    allocate( id ( 1:nz+1 , 1:nz+1 ) ) ; id = 0
    allocate( d2 ( 1:nz+1 , 1:nz+1 ) ) ; d2 = 0
    allocate( bl ( 1:nz+1 , 1:nz+1 ) ) ; bl = 0
    allocate( br ( 1:nz+1 , 1:nz+1 ) ) ; br = 0
    
    id=0
    do i=1,nz+1
       id(i,i)=1
    end do
    d2 = 0

    bl = cl_left_z(1) * id + cl_left_z(2) * dz
    br = cl_right_z(1) * id + cl_right_z(2) * dz
    !
    
    d2 = nuz*d2z 
    call reduction_operator( &
         d2, bl, br, p(1:nzp,1:nzp), pm1(1:nzp,1:nzp), lbd(1:nzp), & 
         cl, cc, fl(1:nzp),fr(1:nzp), nz )
    
    allocate( this%pz  ( 1:nz+1 , 1:nz+1 ) ) ; this%pz = p
    allocate( this%pzm1( 1:nz+1 , 1:nz+1 ) ) ; this%pzm1 = pm1
    allocate( this%lbdz( 1:nz+1          ) ) ; this%lbdz = lbd 

    allocate( this%flz ( 1:nz+1          ) ) ; this%flz = fl
    allocate( this%frz ( 1:nz+1          ) ) ; this%frz = fr
    allocate( this%ccz ( 1:2    , 1:nz+1 ) ) ; this%ccz = cc

    this%clz = cl

    deallocate( d2,bl,br,id )
    deallocate( p,pm1,lbd,fl,fr,cc )


    

    allocate( this%lbdx( 1:nx+1 ) ); this%lbdx = 0
    do i = 1 , nx+1
       im = half_complex_wp(i,nx+1)
       this%lbdx(i) = - nux*(im*wx)**2
    end do
    
    allocate( this%bc_val_l_y ( xst(1):xen(1) , xst(3):xen(3) ) ) ; this%bc_val_l_y = 0
    allocate( this%bc_val_r_y ( xst(1):xen(1) , xst(3):xen(3) ) ) ; this%bc_val_r_y = 0
    allocate( dg_cl_y(yst(1):yen(1),1:2,yst(3):yen(3))) 
    
    allocate( this%bc_val_l_z ( xst(1):xen(1) , xst(2):xen(2) ) ) ; this%bc_val_l_z = 0
    allocate( this%bc_val_r_z ( xst(1):xen(1) , xst(2):xen(2) ) ) ; this%bc_val_r_z = 0
    allocate( dg_cl_z(zst(1):zen(1),zst(2):zen(2),1:2 ))

    
    allocate( DG_LBD(zst(1):zen(1),zst(2):zen(2),zst(3):zen(3) ) )
!!$    if (nrank==0) then
!!$       print*,"---"
!!$       do i = 1 , nx+1
!!$          print*, this%lbdx(i)
!!$       end do
!!$       print*,"---"
!!$       do j = 1 , ny+1
!!$          print*, this%lbdy(j)
!!$       end do
!!$        print*,"---"
!!$       do k = 1 , nz+1
!!$          print*, this%lbdz(k)
!!$       end do
!!$
!!$    end if

    
  end subroutine solver_diag_hhi_build

  subroutine solver_diag_hhi_set_boundary_field(x,y,z,xst,xen,n,&
       udf_bc_left_y, udf_bc_right_y, bvyl, bvyr ,&
       udf_bc_left_z, udf_bc_right_z, bvzl, bvzr  )
    use m_numerics
    implicit none
    integer :: xst(3),xen(3),n(3)
    real(kind=8), dimension(:),allocatable :: x,y,z
    real(kind=8), dimension(xst(1):xen(1),xst(3):xen(3)) :: bvyl,bvyr
    real(kind=8), dimension(xst(1):xen(1),xst(2):xen(2)) :: bvzl,bvzr
    PROCEDURE(UDF_TIMESPACE) :: UDF_BC_LEFT_Y , UDF_BC_RIGHT_Y
    PROCEDURE(UDF_TIMESPACE) :: UDF_BC_LEFT_Z , UDF_BC_RIGHT_Z
    
    integer :: i,j,k
    integer :: nx,ny,nz
    real(kind=8) :: ymin,ymax
    real(kind=8) :: zmin,zmax

    ny = n(2)
    ymin = y(1)
    ymax = y(ny+1)

    nz = n(3)
    zmin = z(1)
    zmax = z(nz+1)

    
    forall(i=xst(1):xen(1),j=xst(2):xen(2))
       bvzl(i,j) = udf_bc_left_z(0d0,X(I),Y(J),ZMIN)
    end forall
    
    forall(i=xst(1):xen(1),j=xst(2):xen(2))
       bvzr(i,j) = udf_bc_right_z(0d0,X(I),Y(J),ZMAX)
    end forall

    
    forall(i=xst(1):xen(1),k=xst(3):xen(3))
       bvyl(i,k) = udf_bc_left_y(0d0,X(I),YMIN,Z(K))
    end forall
    
    forall(i=xst(1):xen(1),k=xst(3):xen(3))
       bvyr(i,k) = udf_bc_right_y(0d0,X(I),YMAX,Z(K))
    end forall

    
    
    
  end subroutine solver_diag_hhi_set_boundary_field
  
  subroutine solver_diag_cart_hii_solve(fi,sfi,sigma,&
       lbdx, nx, xst, xen, &
       lbdy, ny, yst, yen, py, pym1, cly, ccy, fyl, fyr, bvyl, bvyr, &
       lbdz, nz, zst, zen, pz, pzm1, clz, ccz, fzl, fzr, bvzl, bvzr)
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
    complex(kind=8), dimension( 1:ny+1,1:ny+1 ) :: py,pym1
    complex(kind=8), dimension(1:ny+1,1:2) :: ccy
    real(kind=8)::cly(2,2)
    real(kind=8), dimension(1:ny+1) :: fyl,fyr
    real(kind=8), dimension(xst(1):xen(1),xst(3):xen(3)) :: bvyl
    real(kind=8), dimension(xst(1):xen(1),xst(3):xen(3)) :: bvyr

    complex(kind=8), dimension( 1:nz+1 ) :: lbdz
    complex(kind=8), dimension( 1:nz+1,1:nz+1 ) :: pz,pzm1
    complex(kind=8), dimension(1:nz+1,1:2) :: ccz
    real(kind=8)::clz(2,2)
    real(kind=8), dimension(1:nz+1) :: fzl,fzr
    real(kind=8), dimension(xst(1):xen(1),xst(2):xen(2)) :: bvzl
    real(kind=8), dimension(xst(1):xen(1),xst(2):xen(2)) :: bvzr

    integer :: is(3),ie(3)
    integer :: i,j,k


    
    FORALL(I=XST(1):XEN(1),J=XST(2):XEN(2),K=XST(3):XEN(3))
       DG_X(I,J,K) = SFI(I,J,K) &
            - (FYL(J)*BVYL(I,K) + FYR(J)*BVYR(I,K)) &
            - (FZL(K)*BVZL(I,J) + FZR(K)*BVZR(I,J))
    END FORALL
    
    
    DO K = ZST(3),ZEN(3)
       DO J = ZST(2),ZEN(2) 
          DO I=ZST(1),ZEN(1)
             !
             DG_LBD(I,J,K) = ( LBDX(I) + LBDY(J) +  LBDZ(K) + SIGMA )**(-1)
             !
             IF (J==   1) DG_LBD(I,J,K) = 1
             IF (J==NY+1) DG_LBD(I,J,K) = 1
             
             IF (K==   1) DG_LBD(I,J,K) = 1
             IF (K==NZ+1) DG_LBD(I,J,K) = 1
          END DO
       END DO
    END DO
    



    IS = GET_IS_X([0,0,0])
    IE = GET_IE_X([0,0,0])
    if (xst(2)==1) then
       forall(i=IS(1):IE(1),k=IS(3):IE(3))
          dg_x(i,1  ,k) = bvyl(i,k)
       end forall
    end if
    
    if (xen(2)==ny+1) then
       forall(i=is(1):ie(1),k=is(3):ie(3))
          dg_x(i,ny+1,k) = bvyr(i,k)
       end forall
    end if
    
    IS = GET_IS_X([0,1,0])
    IE = GET_IE_X([0,1,0])
    if (xst(3)==1) then
       forall(i=is(1):ie(1),j=is(2):ie(2))
          dg_x(i,j,   1) = bvzl(i,j)
       end forall
    end if
    
    if (xen(3)==nz+1) then
       forall(i=is(1):ie(1),j=is(2):ie(2))
          dg_x(i,j,nz+1) = bvzr(i,j)
       end forall
    end if


    call c2c_1m_x(dg_x,plan_fwd_x)
    
    call transpose_x_to_y(dg_x, zg1_y)
    call tensor_product_j(pym1,zg1_y,zg2_y,yst,yen)
    
    call transpose_y_to_z(zg2_y,zg1_z)
    call tensor_product_k(pzm1,zg1_z,zg2_z,zst,zen)
    
    forall(i=zst(1):zen(1),j=zst(2):zen(2),k=zst(3):zen(3))
       zg1_z(i,j,k) = zg2_z(i,j,k) * DG_LBD(i,j,k)
    end forall
    
    call tensor_product_k(pz,zg1_z,zg2_z,zst,zen)
    call tensor_product_k_cl(ccz, zg2_z, dg_cl_z, zst, zen , nz )


    IS = GET_IS_Z([0,1,0])
    IE = GET_IE_Z([0,1,0])
    forall(i=is(1):ie(1),j=is(2):ie(2))
       zg2_z(i,j,   1) = dg_cl_z(i,j,1)*clz(1,1) + dg_cl_z(i,j,2)*clz(1,2)
       zg2_z(i,j,nz+1) = dg_cl_z(i,j,1)*clz(2,1) + dg_cl_z(i,j,2)*clz(2,2)
    end forall
    
    call transpose_z_to_y(zg2_z,zg1_y)
    call tensor_product_j(py,zg1_y,zg2_y,yst,yen)


    call tensor_product_j_cl(ccy, zg2_y, dg_cl_y, yst, yen , ny )
    IS = GET_IS_Y([0,0,0])
    IE = GET_IE_Y([0,0,0])    
    forall(i=is(1):ie(1),k=is(3):ie(3))
       zg2_y(i,1   ,k) = (dg_cl_y(i,1,k)*cly(1,1) + dg_cl_y(i,2,k)*cly(1,2))
       zg2_y(i,ny+1,k) = (dg_cl_y(i,1,k)*cly(2,1) + dg_cl_y(i,2,k)*cly(2,2))
    end forall
    
    call transpose_y_to_x(zg2_y,dg_x)
    call c2c_1m_x(dg_x,plan_bck_x)
    
    forall(i=xst(1):xen(1),j=xst(2):xen(2),k=xst(3):xen(3))
       fi(i,j,k) = dg_x(i,j,k)
    end forall

  contains

    function get_is_x(dec) result(res)
      implicit none
      integer res(3),dec(3)

      res = xst 
      if (xst(1) == 1) res(1) = res(1)+dec(1)
      if (xst(2) == 1) res(2) = res(2)+dec(2)
      if (xst(3) == 1) res(3) = res(3)+dec(3)
      
    end function get_is_x

    function get_ie_x(dec) result(res)
      implicit none
      integer res(3),dec(3)
      res = xen
      if (xen(1) == nx+1) res(1) = res(1)-dec(1)
      if (xen(2) == ny+1) res(2) = res(2)-dec(2)
      if (xen(3) == nz+1) res(3) = res(3)-dec(3)
    end function get_ie_x
    
    
    function get_is_y(dec) result(res)
      implicit none
      integer res(3),dec(3)

      res = yst 
      if (yst(1) == 1) res(1) = res(1)+dec(1)
      if (yst(2) == 1) res(2) = res(2)+dec(2)
      if (yst(3) == 1) res(3) = res(3)+dec(3)
      
    end function get_is_y
    
    function get_ie_y(dec) result(res)
      implicit none
      integer res(3),dec(3)
      res = yen
      if (yen(1) == nx+1) res(1) = res(1)-dec(1)
      if (yen(2) == ny+1) res(2) = res(2)-dec(2)
      if (yen(3) == nz+1) res(3) = res(3)-dec(3)
    end function get_ie_y


    function get_is_z(dec) result(res)
      implicit none
      integer res(3),dec(3)
      res = zst 
      if (zst(1) == 1) res(1) = res(1)+dec(1)
      if (zst(2) == 1) res(2) = res(2)+dec(2)
      if (zst(3) == 1) res(3) = res(3)+dec(3)
    end function get_is_z
    
    function get_ie_z(dec) result(res)
      implicit none
      integer res(3),dec(3)
      res = zen
      if (zen(1) == nx+1) res(1) = res(1)-dec(1)
      if (zen(2) == ny+1) res(2) = res(2)-dec(2)
      if (zen(3) == nz+1) res(3) = res(3)-dec(3)
    end function get_ie_z

    
    
  end subroutine solver_diag_cart_hii_solve

  SUBROUTINE TENSOR_PRODUCT_J(DY,FI,DFI,YST,YEN)
    IMPLICIT NONE
    INTEGER :: YST(3),YEN(3)
    complex(kind=8) :: DY(YST(2):YEN(2),YST(2):YEN(2))
    complex(kind=8) ::  FI(YST(1):YEN(1),YST(2):YEN(2),YST(3):YEN(3))
    complex(kind=8) :: DFI(YST(1):YEN(1),YST(2):YEN(2),YST(3):YEN(3))
    complex(kind=8) ::  ALPHA,BETA
    INTEGER :: K,N(3),NI,NJ,NK,is,js,ks,ke
        
    ALPHA = dble(1.)
    BETA = dble(0.0)
    N = YEN - YST + 1
    NI = N(1)
    NJ = N(2)
    NK = N(3)
    
    is=lbound(dfi,1)
    js=lbound(dfi,2)
    ks=lbound(dfi,3)
    ke=ubound(dfi,3)
    
    DO K=KS,KE
       CALL ZGEMM('N','T',NI,NJ,NJ, ALPHA ,FI(IS,JS,K), NI , DY, NJ , BETA , DFI(IS,JS,K), NI )
    END DO
    
  END SUBROUTINE TENSOR_PRODUCT_J
  
  
  
  SUBROUTINE TENSOR_PRODUCT_K(DZ, FI, DFI, ZST, ZEN)
    IMPLICIT NONE
    INTEGER         :: ZST(3), ZEN(3)
    complex(KIND=8) :: DZ(ZST(3):ZEN(3), ZST(3):ZEN(3), ZST(1):ZEN(1))
    complex(KIND=8) :: FI(ZST(1):ZEN(1), ZST(2):ZEN(2), ZST(3):ZEN(3))
    complex(KIND=8) :: DFI(ZST(1):ZEN(1), ZST(2):ZEN(2), ZST(3):ZEN(3))
    complex(KIND=8) :: ALPHA, BETA
    INTEGER :: I, N(3), NI, NJ, NK,info
    
    ! Paramètres DGEMM
    ALPHA = 1D0
    BETA = 0D0
    
    ! Dimensions locales
    N = ZEN - ZST + 1
    NI = N(1)
    NJ = N(2)
    NK = N(3)
    
    
    CALL ZGEMM('N','T',NI*NJ,NK,NK, ALPHA, FI, NI*NJ , DZ, NK , BETA , DFI, NI*NJ )
    
  END SUBROUTINE TENSOR_PRODUCT_K
  
  SUBROUTINE TENSOR_PRODUCT_K_CL(CC,FI,DFI,ZST,ZEN,NZ)
    IMPLICIT NONE
    INTEGER :: ZST(3),ZEN(3),NZ
    complex(kind=8) :: CC(1:2,1:NZ+1)
    complex(kind=8) ::  FI(ZST(1):ZEN(1),ZST(2):ZEN(2),ZST(3):ZEN(3))
    complex(kind=8) :: DFI(ZST(1):ZEN(1),ZST(2):ZEN(2),1:2)
    
    INTEGER :: N(3),NI,NJ,NK
    complex(kind=8) ::  ALPHA,BETA
    
    ALPHA = 1d0
    BETA = 0d0
    N = ZEN - ZST + 1
    NI = N(1)
    NJ = N(2)
    NK = N(3) ! NZ + 1 car Z-pencil
    CALL ZGEMM('N','T',NI*NJ,2,NK, ALPHA, FI, NI*NJ , CC, 2 , BETA , DFI, NI*NJ )
    
  END SUBROUTINE TENSOR_PRODUCT_K_CL
  
  SUBROUTINE TENSOR_PRODUCT_J_CL(CC,FI,DFI,YST,YEN,NY)
    IMPLICIT NONE
    INTEGER  :: YST(3),YEN(3),NY
    complex(kind=8) :: CC(1:2,1:NY+1)
    complex(kind=8) ::  FI(YST(1):YEN(1),YST(2):YEN(2),YST(3):YEN(3))
    complex(kind=8) :: DFI(YST(1):YEN(1),1:2,YST(3):YEN(3))
    complex(kind=8) ::  ALPHA,BETA
    INTEGER  :: K,N(3),NI,NJ,NK
    
    ALPHA = cmplx(1.0,0.0,kind=8)
    BETA =  cmplx(0.0,0.0,kind=8)
    N = YEN - YST + 1
    NI = N(1)
    NJ = N(2)
    DO K=YST(3),YEN(3)
       CALL ZGEMM('N','T',NI,2,NJ,ALPHA,FI(YST(1),YST(2),K),NI,CC,2,BETA,DFI(YST(1),1,K),NI)
    END DO
    
  END SUBROUTINE TENSOR_PRODUCT_J_CL
   
  
end module m_solver_diag_cart_hii
