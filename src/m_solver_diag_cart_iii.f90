module m_solver_diag_cart_iii
  use M_FOURIER_TRANSFORM
    
  implicit none
  
  type t_solver_diag_cart_iii
     
     integer :: n(3)
     integer :: xst(3),xen(3)  ! x-distributed array size 
     integer :: yst(3),yen(3)  ! y-distributed array size 
     integer :: zst(3),zen(3)  ! z-distributed array size 
     
     REAL(KIND=8) :: nu(3),sigma ! diffusion coefficient


     REAL(KIND=8) :: alf_l_x, bet_l_x
     REAL(KIND=8) :: alf_r_x, bet_r_x
     
     REAL(KIND=8) :: alf_l_y, bet_l_y 
     REAL(KIND=8) :: alf_r_y, bet_r_y

     REAL(KIND=8) :: alf_l_z, bet_l_z 
     REAL(KIND=8) :: alf_r_z, bet_r_z  


     REAL(KIND=8), dimension(:,:)  , allocatable :: bc_val_l_x, bc_val_r_x
     REAL(KIND=8), dimension(:,:)  , allocatable :: bc_val_l_y, bc_val_r_y
     REAL(KIND=8), dimension(:,:)  , allocatable :: bc_val_l_z, bc_val_r_z

     
     COMPLEX(KIND=8), dimension(:)  , allocatable :: LBDX ! 
     COMPLEX(KIND=8), dimension(:,:), allocatable :: PX,PXM1
     REAL(KIND=8)   , dimension(:)  , allocatable :: FLX,FRX
     COMPLEX(KIND=8), dimension(:,:), allocatable :: CCX
     REAL(KIND=8) :: clX(2,2)
     
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

     logical :: scale
   contains

     

     PROCEDURE :: SOLVE => IFCE_SOLVE
     PROCEDURE :: SOLVE_POISSON =>IFCE_SOLVE_PRESSURE
     PROCEDURE :: INITIALISE => IFCE_INITIALISE
     PROCEDURE :: SET_BCS => SOLVER_DIAG_CART_III_SET_BCS
     PROCEDURE :: SET_BVS => SOLVER_DIAG_CART_III_SET_BVS
     PROCEDURE :: SET_PARAMS => SOLVER_DIAG_CART_III_SET_PARAMS
     
     
  end type t_solver_diag_cart_iii
  
  COMPLEX(KIND=8), dimension(:,:,:), allocatable :: dg_cl_x
  COMPLEX(KIND=8), dimension(:,:,:), allocatable :: dg_cl_z
  COMPLEX(KIND=8), dimension(:,:,:), allocatable :: dg_cl_y
  complex(KIND=8), dimension(:,:,:), allocatable :: dg_LBD
  
contains

  ! common interface for 
  subroutine ifce_initialise(this,grid,opx,opy,opz,ph)
    use m_mesh_base
    use m_operator_tcheby
    use m_operator_fourier_dft
    use decomp_2d
    implicit none
    class(t_solver_diag_cart_iii):: this
    type(t_mesh_base)            :: grid(3)
    type(t_operator_tcheby)      :: opx,opy,opz
    
    real(kind=8) :: sigma
    real(kind=8) :: nu(3),bcl(3,2),bcr(3,2)
    integer :: n(3)
    type(decomp_info) :: ph

    ! get the global size. ni := nx+1 
    n = [ph%xen(1)-ph%xst(1),ph%yen(2)-ph%yst(2),ph%zen(3)-ph%zst(3)]
    BCL = 0
    BCR = 0
    THIS%N = N

    
    BCL(1,1:2) = [ THIS%ALF_L_X , THIS%BET_L_X ]
    BCR(1,1:2) = [ THIS%ALF_R_X , THIS%BET_R_X ]
    
    BCL(2,1:2) = [ THIS%ALF_L_Y , THIS%BET_L_Y ]
    BCR(2,1:2) = [ THIS%ALF_R_Y , THIS%BET_R_Y ]

    BCL(3,1:2) = [ THIS%ALF_L_Z , THIS%BET_L_Z ]
    BCR(3,1:2) = [ THIS%ALF_R_Z , THIS%BET_R_Z ]
    
    CALL SOLVER_DIAG_III_BUILD (THIS, &
         GRID(1)%X, THIS%NU(1), THIS%N(1), PH%XST,  PH%XEN, OPX%GET_D2(), OPX%GET_D1(), BCL(1,:), BCR(1,:), &
         GRID(2)%X, THIS%NU(2), THIS%N(2), PH%YST,  PH%YEN, OPY%GET_D2(), OPY%GET_D1(), BCL(2,:), BCR(2,:), &
         GRID(3)%X, THIS%NU(3), THIS%N(3), PH%ZST,  PH%ZEN, OPZ%GET_D2(), OPZ%GET_D1(), BCL(3,:), BCR(3,:)  )
    
  end subroutine ifce_initialise

  
  SUBROUTINE SOLVER_DIAG_CART_III_SET_BVS(THIS,MESH,AXIS,UDF_MINUS,UDF_PLUS)
    use decomp_2d
    use m_mesh_base
    use mpi
    IMPLICIT NONE
    CLASS(T_SOLVER_DIAG_CART_III) :: THIS
    TYPE(T_MESH_BASE)             :: MESH(3)
    INTEGER                       :: AXIS
    PROCEDURE(UDF_TIMESPACE)      :: UDF_PLUS
    PROCEDURE(UDF_TIMESPACE)      :: UDF_MINUS
    
    TYPE(DECOMP_INFO) :: PH
    REAL(KIND=8) :: TIME
    integer :: ierr
    
    CALL GET_DECOMP_INFO(PH)
    
    
    if (axis==1) then
       CALL SOLVER_DIAG_III_SET_BOUNDARY_FIELD_X(&
            MESH(1)%X, MESH(2)%X, MESH(3)%X,PH%XST,PH%XEN,THIS%N ,&
            UDF_MINUS,UDF_PLUS,THIS%BC_VAL_L_X,THIS%BC_VAL_R_X)
       
    else if (axis==2) then
       CALL SOLVER_DIAG_III_SET_BOUNDARY_FIELD_Y(&
            MESH(1)%X, MESH(2)%X, MESH(3)%X,PH%XST,PH%XEN,THIS%N, &
            UDF_MINUS,UDF_PLUS,THIS%BC_VAL_L_Y,THIS%BC_VAL_R_Y)
       
    else if (axis==3) then
       CALL SOLVER_DIAG_III_SET_BOUNDARY_FIELD_Z(&
            MESH(1)%X, MESH(2)%X, MESH(3)%X,PH%XST,PH%XEN,THIS%N, &
            UDF_MINUS,UDF_PLUS,THIS%BC_VAL_L_Z,THIS%BC_VAL_R_Z)
       
    else
       
       call MPI_Abort(mpi_comm_world,0,ierr)
    end if
  end subroutine SOLVER_DIAG_CART_III_SET_BVS



  SUBROUTINE SOLVER_DIAG_CART_III_SET_BCS(THIS,AXIS,BCS_MINUS,BCS_PLUS)
    use decomp_2d
    use mpi
    IMPLICIT NONE
    CLASS(T_SOLVER_DIAG_CART_III) :: THIS
    INTEGER                       :: AXIS
    real(kind=8)                  :: BCS_MINUS(2)
    real(kind=8)                  :: BCS_PLUS(2)
    integer :: ierr
    TYPE(DECOMP_INFO) :: PH
    
    CALL GET_DECOMP_INFO(PH)
    IF (AXIS==1) THEN

       THIS%ALF_L_X = BCS_MINUS(1)
       THIS%BET_L_X = BCS_MINUS(2)

       THIS%ALF_R_X = BCS_PLUS(1)
       THIS%BET_R_X = BCS_PLUS(2)
       
    ELSE IF (AXIS==2) THEN

       THIS%ALF_L_Y = BCS_MINUS(1)
       THIS%BET_L_Y = BCS_MINUS(2)

       THIS%ALF_R_Y = BCS_PLUS(1)
       THIS%BET_R_Y = BCS_PLUS(2)
       
    ELSE IF (AXIS==3) THEN
       
       THIS%ALF_L_Z = BCS_MINUS(1)
       THIS%BET_L_Z = BCS_MINUS(2)

       THIS%ALF_R_Z = BCS_PLUS(1)
       THIS%BET_R_Z = BCS_PLUS(2)

    END IF

    
  end subroutine SOLVER_DIAG_CART_III_SET_BCS

  
  
  subroutine ifce_solve(this,fi,sfi,sigma,ph)
    use decomp_2d
    implicit none
    class(t_solver_diag_cart_iii) :: this
    real(kind=8),dimension(:,:,:),allocatable :: fi,sfi
    real(kind=8) :: sigma
    integer :: n(3)
    type(decomp_info) :: ph
    integer :: i


    
    call solver_diag_cart_iii_solve(fi,sfi,sigma,&
         ! xxxx
         this%lbdx, this%n(1), ph%xst, ph%xen, this%px,  this%pxm1, &
         this%clx, this%ccx,  this%flx,  this%frx, this%bc_val_l_x, this%bc_val_r_x,&
         ! xxxx
         this%lbdy, this%n(2), ph%yst, ph%yen, this%py,  this%pym1, &
         this%cly, this%ccy,  this%fly,  this%fry, this%bc_val_l_y, this%bc_val_r_y,&
         ! xxxx
         this%lbdz, this%n(3), ph%zst, ph%zen, this%pz,  this%pzm1,  &
         this%clz, this%ccz,  this%flz,  this%frz, this%bc_val_l_z, this%bc_val_r_z)
    
  end subroutine ifce_solve


  subroutine ifce_solve_pressure(this,fi,sfi,nu,ph)
    use decomp_2d
    implicit none
    class(t_solver_diag_cart_iii) :: this
    real(kind=8),dimension(:,:,:),allocatable :: fi,sfi
    real(kind=8) :: nu(3)
    integer :: n(3)
    type(decomp_info) :: ph
    integer :: i


    
    call solver_diag_cart_iii_solve_scaled(fi,sfi,0d0,&
         ! xxxx
         nu(1),this%lbdx, this%n(1), ph%xst, ph%xen, this%px,  this%pxm1, &
         this%clx, this%ccx,  this%flx,  this%frx, this%bc_val_l_x, this%bc_val_r_x,&
         ! xxxx
         nu(2),this%lbdy, this%n(2), ph%yst, ph%yen, this%py,  this%pym1, &
         this%cly, this%ccy,  this%fly,  this%fry, this%bc_val_l_y, this%bc_val_r_y,&
         ! xxxx
         nu(3),this%lbdz, this%n(3), ph%zst, ph%zen, this%pz,  this%pzm1,  &
         this%clz, this%ccz,  this%flz,  this%frz, this%bc_val_l_z, this%bc_val_r_z)
    
  end subroutine ifce_solve_pressure

  
  !> basic routine 

  SUBROUTINE SOLVER_DIAG_CART_III_SET_PARAMS(THIS,NU,SIGMA )
    IMPLICIT NONE
    CLASS(T_SOLVER_DIAG_CART_III) :: THIS
    REAL(KIND=8) :: NU(3)
    REAL(KIND=8) :: SIGMA

    THIS%NU = NU
    THIS%SIGMA = SIGMA
    
  END SUBROUTINE SOLVER_DIAG_CART_III_SET_PARAMS

  
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
    

    subroutine solver_diag_iii_build(this,&
         x, nux, nx, xst, xen, d2x, dx, cl_left_x, cl_right_x ,&
         y, nuy, ny, yst, yen, d2y, dy, cl_left_y, cl_right_y , &
         z, nuz, nz, zst, zen, d2z, dz, cl_left_z, cl_right_z )
      use decomp_2d_mpi
      implicit none
      class(t_solver_diag_cart_iii) :: this
      integer :: nx,xst(3),xen(3)
      REAL(KIND=8) :: x(nx+1),nux,d2x(1:nx+1,1:nx+1),dx(1:nx+1,1:nx+1)
      REAL(KIND=8) :: cl_left_x(2), cl_right_x(2)
      
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
    
    this%nu=[nux,nuy,nuz]
    

    ! OX ------------------------- !
    nxp = nx + 1
    
    allocate( p  ( 1:nx+1 , 1:nx+1 ) ) ; p   = 0
    allocate( pm1( 1:nx+1 , 1:nx+1 ) ) ; pm1 = 0
    allocate( lbd( 1:nx+1          ) ) ; lbd = 0 

    allocate( fl ( 1:nx+1         ) ) ; fl = 0
    allocate( fr ( 1:nx+1         ) ) ; fr = 0
    allocate( cc ( 1:2   , 1:nx+1 ) ) ; cc = 0
    
    allocate( id ( 1:nx+1 , 1:nx+1 ) ) ; id = 0
    allocate( d2 ( 1:nx+1 , 1:nx+1 ) ) ; d2 = 0
    allocate( bl ( 1:nx+1 , 1:nx+1 ) ) ; bl = 0
    allocate( br ( 1:nx+1 , 1:nx+1 ) ) ; br = 0
    
    id=0
    do i=1,nx+1
       id(i,i)=1
    end do
    d2 = 0

    bl = cl_left_x(1) * id + cl_left_x(2) * dx
    br = cl_right_x(1) * id + cl_right_x(2) * dx
    !
    
    d2 = nux*d2x 

    call reduction_operator( &
         d2, bl, br, p(1:nxp,1:nxp), pm1(1:nxp,1:nxp), lbd(1:nxp), & 
         cl, cc, fl(1:nxp),fr(1:nxp), nx )

    allocate( this%px  ( 1:nx+1 , 1:nx+1 ) ) ; this%px = p
    allocate( this%pxm1( 1:nx+1 , 1:nx+1 ) ) ; this%pxm1 = pm1
    allocate( this%lbdx( 1:nx+1          ) ) ; this%lbdx = lbd

    allocate( this%flx ( 1:nx+1          ) ) ; this%flx = fl
    allocate( this%frx ( 1:nx+1          ) ) ; this%frx = fr
    allocate( this%ccx ( 1:2    , 1:nx+1 ) ) ; this%ccx = cc
    
    this%clx = cl

    deallocate( d2,bl,br,id )
    deallocate( p,pm1,lbd,fl,fr,cc )
    
    ! OY ------------------------- !
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
    

    allocate( this%bc_val_l_x ( xst(2):xen(2) , xst(3):xen(3) ) ) ; this%bc_val_l_x = 0
    allocate( this%bc_val_r_x ( xst(2):xen(2) , xst(3):xen(3) ) ) ; this%bc_val_r_x = 0
    
    
    allocate( this%bc_val_l_y ( xst(1):xen(1) , xst(3):xen(3) ) ) ; this%bc_val_l_y = 0
    allocate( this%bc_val_r_y ( xst(1):xen(1) , xst(3):xen(3) ) ) ; this%bc_val_r_y = 0
    
    
    allocate( this%bc_val_l_z ( xst(1):xen(1) , xst(2):xen(2) ) ) ; this%bc_val_l_z = 0
    allocate( this%bc_val_r_z ( xst(1):xen(1) , xst(2):xen(2) ) ) ; this%bc_val_r_z = 0


    if (.not.allocated(dg_cl_x)) then
       
       allocate( dg_cl_x(1:2,xst(2):xen(2),xst(3):xen(3))) 
       allocate( dg_cl_y(yst(1):yen(1),1:2,yst(3):yen(3))) 
       allocate( dg_cl_z(zst(1):zen(1),zst(2):zen(2),1:2 ))
       allocate( DG_LBD(zst(1):zen(1),zst(2):zen(2),zst(3):zen(3) ) )

       dg_cl_x = 0
       dg_cl_y = 0
       dg_cl_z = 0
       DG_LBD = 0
       
    end if
    
  end subroutine solver_diag_iii_build

  SUBROUTINE SOLVER_DIAG_III_SET_BOUNDARY_FIELD_X(X,Y,Z,XST,XEN,N,&
       UDF_BC_LEFT_X, UDF_BC_RIGHT_X, BVXL, BVXR )
    USE M_NUMERICS
    IMPLICIT NONE
    INTEGER :: XST(3),XEN(3),N(3)
    REAL(KIND=8), DIMENSION(:),ALLOCATABLE :: X,Y,Z
    REAL(KIND=8), DIMENSION(XST(2):XEN(2),XST(3):XEN(3)) :: BVXL,BVXR
    PROCEDURE(UDF_TIMESPACE) :: UDF_BC_LEFT_X , UDF_BC_RIGHT_X

    INTEGER :: I,J,K
    INTEGER :: NX,NY,NZ
    REAL(KIND=8) :: XMIN, XMAX
    
    NX = N(1)
    XMIN = X(1)
    XMAX = X(NX+1)
    
    FORALL(J=XST(2):XEN(2),K=XST(3):XEN(3))
       BVXL(J,K) = UDF_BC_LEFT_X(0D0,XMIN,Y(J),Z(K))
    END FORALL
    
    FORALL(J=XST(2):XEN(2),K=XST(3):XEN(3))
       BVXR(J,K) = UDF_BC_RIGHT_X(0D0,XMAX,Y(J),Z(K))
    END FORALL
    
  END SUBROUTINE SOLVER_DIAG_III_SET_BOUNDARY_FIELD_X

  SUBROUTINE SOLVER_DIAG_III_SET_BOUNDARY_FIELD_Y(X,Y,Z,XST,XEN,N,&
       UDF_BC_LEFT_Y, UDF_BC_RIGHT_Y, BVYL, BVYR )
    USE M_NUMERICS
    IMPLICIT NONE
    INTEGER :: XST(3),XEN(3),N(3)
    REAL(KIND=8), DIMENSION(:),ALLOCATABLE :: X,Y,Z
    REAL(KIND=8), DIMENSION(XST(1):XEN(1),XST(3):XEN(3)) :: BVYL,BVYR
    PROCEDURE(UDF_TIMESPACE) :: UDF_BC_LEFT_Y , UDF_BC_RIGHT_Y
    
    INTEGER :: I,J,K
    INTEGER :: NX,NY,NZ
    REAL(KIND=8) :: YMIN, YMAX
    
    NY = N(2)
    YMIN = Y(1)
    YMAX = Y(NY+1)
    
    forall(i=xst(1):xen(1),k=xst(3):xen(3))
       bvyl(i,k) = udf_bc_left_y(0d0,X(I),YMIN,Z(K))
    end forall
    
    forall(i=xst(1):xen(1),k=xst(3):xen(3))
       bvyr(i,k) = udf_bc_right_y(0d0,X(I),YMAX,Z(K))
    end forall
    
    
  END SUBROUTINE SOLVER_DIAG_III_SET_BOUNDARY_FIELD_Y

  SUBROUTINE SOLVER_DIAG_III_SET_BOUNDARY_FIELD_Z(X,Y,Z,XST,XEN,N,&
       UDF_BC_LEFT_Z, UDF_BC_RIGHT_Z, BVZL, BVZR )
    USE M_NUMERICS
    IMPLICIT NONE
    INTEGER :: XST(3),XEN(3),N(3)
    REAL(KIND=8), DIMENSION(:),ALLOCATABLE :: X,Y,Z
    REAL(KIND=8), DIMENSION(XST(1):XEN(1),XST(2):XEN(2)) :: BVZL,BVZR
    PROCEDURE(UDF_TIMESPACE) :: UDF_BC_LEFT_Z , UDF_BC_RIGHT_Z
    
    INTEGER :: I,J,K
    INTEGER :: NX,NY,NZ
    REAL(KIND=8) :: ZMIN, ZMAX
    
    NZ = N(3)
    ZMIN = Z(1)
    ZMAX = Z(NZ+1)
    
    forall(i=xst(1):xen(1),j=xst(2):xen(2))
       bvzl(i,j) = udf_bc_left_z(0d0,X(I),Y(J),ZMIN)
    end forall
    
    forall(i=xst(1):xen(1),j=xst(2):xen(2))
       bvzr(i,j) = udf_bc_right_z(0d0,X(I),Y(J),ZMAX)
    end forall
    
    
  END SUBROUTINE SOLVER_DIAG_III_SET_BOUNDARY_FIELD_Z

  

  subroutine solver_diag_iii_set_boundary_field(x,y,z,xst,xen,n,&
       udf_bc_left_x, udf_bc_right_x, bvxl, bvxr ,&
       udf_bc_left_y, udf_bc_right_y, bvyl, bvyr ,&
       udf_bc_left_z, udf_bc_right_z, bvzl, bvzr  )
    use m_numerics
    implicit none
    integer :: xst(3),xen(3),n(3)
    real(kind=8), dimension(:),allocatable :: x,y,z
    real(kind=8), dimension(xst(2):xen(2),xst(3):xen(3)) :: bvxl,bvxr
    real(kind=8), dimension(xst(1):xen(1),xst(3):xen(3)) :: bvyl,bvyr
    real(kind=8), dimension(xst(1):xen(1),xst(2):xen(2)) :: bvzl,bvzr
    PROCEDURE(UDF_TIMESPACE) :: UDF_BC_LEFT_X , UDF_BC_RIGHT_X
    PROCEDURE(UDF_TIMESPACE) :: UDF_BC_LEFT_Y , UDF_BC_RIGHT_Y
    PROCEDURE(UDF_TIMESPACE) :: UDF_BC_LEFT_Z , UDF_BC_RIGHT_Z
    
    integer :: i,j,k
    integer :: nx,ny,nz
    real(kind=8) :: xmin, xmax
    real(kind=8) :: ymin, ymax
    real(kind=8) :: zmin, zmax


    nx = n(1)
    xmin = x(1)
    xmax = x(nx+1)

    
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

    
    forall(j=xst(2):xen(2),k=xst(3):xen(3))
       bvxl(j,k) = udf_bc_left_x(0d0,XMIN,Y(J),Z(K))
    end forall
    
    forall(j=xst(2):xen(2),k=xst(3):xen(3))
       bvxr(j,k) = udf_bc_right_x(0d0,XMAX,Y(J),Z(K))
    end forall

    
    
  end subroutine solver_diag_iii_set_boundary_field
  
  subroutine solver_diag_cart_iii_solve(fi,sfi,sigma,&
       lbdx, nx, xst, xen, px, pxm1, clx, ccx, fxl, fxr, bvxl, bvxr, &
       lbdy, ny, yst, yen, py, pym1, cly, ccy, fyl, fyr, bvyl, bvyr, &
       lbdz, nz, zst, zen, pz, pzm1, clz, ccz, fzl, fzr, bvzl, bvzr)
    use m_fourier_transform
    use decomp_2d
    use mpi
    implicit none
    real(kind=8) :: sigma
    integer :: nx,xst(3),xen(3)
    integer :: ny,yst(3),yen(3)
    integer :: nz,zst(3),zen(3)
    real(kind=8), dimension(xst(1):xen(1),xst(2):xen(2),xst(3):xen(3) ) :: fi,sfi

    complex(kind=8), dimension( 1:nx+1 ) :: lbdx
    complex(kind=8), dimension( 1:nx+1,1:nx+1 ) :: px,pxm1
    complex(kind=8), dimension(1:nx+1,1:2) :: ccx
    real(kind=8) :: clx(2,2)
    real(kind=8), dimension(1:nx+1) :: fxl,fxr
    real(kind=8), dimension(xst(2):xen(2),xst(3):xen(3)) :: bvxl
    real(kind=8), dimension(xst(2):xen(2),xst(3):xen(3)) :: bvxr
    
    complex(kind=8), dimension( 1:ny+1 ) :: lbdy
    complex(kind=8), dimension( 1:ny+1,1:ny+1 ) :: py,pym1
    complex(kind=8), dimension(1:ny+1,1:2) :: ccy
    real(kind=8) :: cly(2,2)
    real(kind=8), dimension(1:ny+1) :: fyl,fyr
    real(kind=8), dimension(xst(1):xen(1),xst(3):xen(3)) :: bvyl
    real(kind=8), dimension(xst(1):xen(1),xst(3):xen(3)) :: bvyr

    complex(kind=8), dimension( 1:nz+1 ) :: lbdz
    complex(kind=8), dimension( 1:nz+1,1:nz+1 ) :: pz,pzm1
    complex(kind=8), dimension(1:nz+1,1:2) :: ccz
    real(kind=8) :: clz(2,2)
    real(kind=8), dimension(1:nz+1) :: fzl,fzr
    real(kind=8), dimension(xst(1):xen(1),xst(2):xen(2)) :: bvzl
    real(kind=8), dimension(xst(1):xen(1),xst(2):xen(2)) :: bvzr

    integer :: is(3),ie(3)
    integer :: i,j,k
    integer :: rank,ierr


    FORALL(I=XST(1):XEN(1),J=XST(2):XEN(2),K=XST(3):XEN(3))
       DG_X(I,J,K) = SFI(I,J,K) &
            - (FXL(I)*BVXL(J,K) + FXR(I)*BVXR(J,K)) &
            - (FYL(J)*BVYL(I,K) + FYR(J)*BVYR(I,K)) &
            - (FZL(K)*BVZL(I,J) + FZR(K)*BVZR(I,J))
    END FORALL
    
    
    DO K = ZST(3),ZEN(3)
       DO J = ZST(2),ZEN(2) 
          DO I=ZST(1),ZEN(1)
             !
             DG_LBD(I,J,K) = ( LBDX(I) + LBDY(J) +  LBDZ(K) + SIGMA )**(-1)

             IF (I==   1) DG_LBD(I,J,K) = 1
             IF (I==NX+1) DG_LBD(I,J,K) = 1
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
    !if (xst(1)==1) then
       forall(J=IS(2):IE(2),k=IS(3):IE(3))
          dg_x(1   ,j,k) = bvxl(j,k)
       end forall
    !end if
    
    !if (xen(1)==nx+1) then
       forall(j=is(2):ie(2),k=is(3):ie(3))
          dg_x(nx+1,j,k) = bvxr(j,k)
       end forall

    !end if
    
    IS = GET_IS_X([1,1,0])
    IE = GET_IE_X([1,1,0])
    
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
    
    IS = GET_IS_X([1,1,1])
    IE = GET_IE_X([1,1,1])
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


    
    call tensor_product_i(pxm1,dg_x,zg_x,xst,xen)
    
    call transpose_x_to_y(zg_x, zg1_y)
    call tensor_product_j(pym1,zg1_y,zg2_y,yst,yen)
    
    call transpose_y_to_z(zg2_y,zg1_z)
    call tensor_product_k(pzm1,zg1_z,zg2_z,zst,zen)
    
    forall(i=zst(1):zen(1),j=zst(2):zen(2),k=zst(3):zen(3))
       zg1_z(i,j,k) = zg2_z(i,j,k) * DG_LBD(i,j,k)
    end forall
    
    call tensor_product_k(pz,zg1_z,zg2_z,zst,zen)
    call tensor_product_k_cl(ccz, zg2_z, dg_cl_z, zst, zen , nz )


    IS = GET_IS_Z([1,1,1])
    IE = GET_IE_Z([1,1,1])
    forall(i=is(1):ie(1),j=is(2):ie(2))
       zg2_z(i,j,   1) = dg_cl_z(i,j,1)*clz(1,1) + dg_cl_z(i,j,2)*clz(1,2)
       zg2_z(i,j,nz+1) = dg_cl_z(i,j,1)*clz(2,1) + dg_cl_z(i,j,2)*clz(2,2)
    end forall
    
    CALL TRANSPOSE_Z_TO_Y(ZG2_Z,ZG1_Y)
    CALL TENSOR_PRODUCT_J(PY,ZG1_Y,ZG2_Y,YST,YEN)
    CALL TENSOR_PRODUCT_J_CL(CCY, ZG2_Y, DG_CL_Y, YST, YEN , NY )
    IS = GET_IS_Y([1,1,0])
    IE = GET_IE_Y([1,1,0])    
    FORALL(I=IS(1):IE(1),K=IS(3):IE(3))
       ZG2_Y(I,1   ,K) = (DG_CL_Y(I,1,K)*CLY(1,1) + DG_CL_Y(I,2,K)*CLY(1,2))
       ZG2_Y(I,NY+1,K) = (DG_CL_Y(I,1,K)*CLY(2,1) + DG_CL_Y(I,2,K)*CLY(2,2))
    END FORALL

    
    call transpose_y_to_x(zg2_y,zg_x)
    call tensor_product_i(px,zg_x,dg_x,xst,xen)
    call tensor_product_i_cl(ccx, dg_x, dg_cl_x, xst, xen )
    
    
    IS = GET_IS_X([0,0,0])
    IE = GET_IE_X([0,0,0])
    forall(j=is(2):ie(2),k=is(3):ie(3))
       dg_x(1   ,j,k) = dg_cl_x(1,j,k)*clx(1,1) + dg_cl_x(2,j,k)*clx(1,2)
       dg_x(nx+1,j,k) = dg_cl_x(1,j,k)*clx(2,1) + dg_cl_x(2,j,k)*clx(2,2)
    end forall
    
    call mpi_comm_rank(mpi_comm_world,rank,ierr)
    
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

    
    
  end subroutine solver_diag_cart_iii_solve


  subroutine solver_diag_cart_iii_solve_scaled(fi,sfi,sigma,&
         nu_x,lbdx, nx, xst, xen, px, pxm1, clx, ccx, fxl, fxr, bvxl, bvxr, &
         nu_y,lbdy, ny, yst, yen, py, pym1, cly, ccy, fyl, fyr, bvyl, bvyr, &
         nu_z,lbdz, nz, zst, zen, pz, pzm1, clz, ccz, fzl, fzr, bvzl, bvzr)
      use m_fourier_transform
      use decomp_2d
      use mpi
      implicit none
      real(kind=8) :: sigma,nu_x,nu_y,nu_z
      integer :: nx,xst(3),xen(3)
      integer :: ny,yst(3),yen(3)
      integer :: nz,zst(3),zen(3)
      real(kind=8), dimension(xst(1):xen(1),xst(2):xen(2),xst(3):xen(3) ) :: fi,sfi

    complex(kind=8), dimension( 1:nx+1 ) :: lbdx
    complex(kind=8), dimension( 1:nx+1,1:nx+1 ) :: px,pxm1
    complex(kind=8), dimension(1:nx+1,1:2) :: ccx
    real(kind=8) :: clx(2,2)
    real(kind=8), dimension(1:nx+1) :: fxl,fxr
    real(kind=8), dimension(xst(2):xen(2),xst(3):xen(3)) :: bvxl
    real(kind=8), dimension(xst(2):xen(2),xst(3):xen(3)) :: bvxr
    
    complex(kind=8), dimension( 1:ny+1 ) :: lbdy
    complex(kind=8), dimension( 1:ny+1,1:ny+1 ) :: py,pym1
    complex(kind=8), dimension(1:ny+1,1:2) :: ccy
    real(kind=8) :: cly(2,2)
    real(kind=8), dimension(1:ny+1) :: fyl,fyr
    real(kind=8), dimension(xst(1):xen(1),xst(3):xen(3)) :: bvyl
    real(kind=8), dimension(xst(1):xen(1),xst(3):xen(3)) :: bvyr

    complex(kind=8), dimension( 1:nz+1 ) :: lbdz
    complex(kind=8), dimension( 1:nz+1,1:nz+1 ) :: pz,pzm1
    complex(kind=8), dimension(1:nz+1,1:2) :: ccz
    real(kind=8) :: clz(2,2)
    real(kind=8), dimension(1:nz+1) :: fzl,fzr
    real(kind=8), dimension(xst(1):xen(1),xst(2):xen(2)) :: bvzl
    real(kind=8), dimension(xst(1):xen(1),xst(2):xen(2)) :: bvzr

    integer :: is(3),ie(3)
    integer :: i,j,k
    integer :: ierr

    FORALL(I=XST(1):XEN(1),J=XST(2):XEN(2),K=XST(3):XEN(3))
       DG_X(I,J,K) = SFI(I,J,K) &
            - (FXL(I)*BVXL(J,K) + FXR(I)*BVXR(J,K)) &
            - (FYL(J)*BVYL(I,K) + FYR(J)*BVYR(I,K)) &
            - (FZL(K)*BVZL(I,J) + FZR(K)*BVZR(I,J))
    END FORALL
    
    
    DO K = ZST(3),ZEN(3)
       DO J = ZST(2),ZEN(2) 
          DO I=ZST(1),ZEN(1)
             !
             DG_LBD(I,J,K) = ( nu_x*LBDX(I) + nu_y*LBDY(J) +  nu_z*LBDZ(K) + SIGMA )**(-1)

             IF (I==   1) DG_LBD(I,J,K) = 1
             IF (I==NX+1) DG_LBD(I,J,K) = 1
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
    !if (xst(1)==1) then
       forall(J=IS(2):IE(2),k=IS(3):IE(3))
          dg_x(1   ,j,k) = bvxl(j,k)
       end forall
    !end if
    
    !if (xen(1)==nx+1) then
       forall(j=is(2):ie(2),k=is(3):ie(3))
          dg_x(nx+1,j,k) = bvxr(j,k)
       end forall

    !end if
    
    IS = GET_IS_X([1,1,0])
    IE = GET_IE_X([1,1,0])
    
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
    
    IS = GET_IS_X([1,1,1])
    IE = GET_IE_X([1,1,1])
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


    
    call tensor_product_i(pxm1,dg_x,zg_x,xst,xen)
    
    call transpose_x_to_y(zg_x, zg1_y)
    call tensor_product_j(pym1,zg1_y,zg2_y,yst,yen)
    
    call transpose_y_to_z(zg2_y,zg1_z)
    call tensor_product_k(pzm1,zg1_z,zg2_z,zst,zen)
    
    forall(i=zst(1):zen(1),j=zst(2):zen(2),k=zst(3):zen(3))
       zg1_z(i,j,k) = zg2_z(i,j,k) * DG_LBD(i,j,k)
    end forall
    
    call tensor_product_k(pz,zg1_z,zg2_z,zst,zen)
    call tensor_product_k_cl(ccz, zg2_z, dg_cl_z, zst, zen , nz )


    IS = GET_IS_Z([1,1,1])
    IE = GET_IE_Z([1,1,1])
    forall(i=is(1):ie(1),j=is(2):ie(2))
       zg2_z(i,j,   1) = dg_cl_z(i,j,1)*clz(1,1) + dg_cl_z(i,j,2)*clz(1,2)
       zg2_z(i,j,nz+1) = dg_cl_z(i,j,1)*clz(2,1) + dg_cl_z(i,j,2)*clz(2,2)
    end forall
    
    CALL TRANSPOSE_Z_TO_Y(ZG2_Z,ZG1_Y)
    CALL TENSOR_PRODUCT_J(PY,ZG1_Y,ZG2_Y,YST,YEN)
    CALL TENSOR_PRODUCT_J_CL(CCY, ZG2_Y, DG_CL_Y, YST, YEN , NY )
    IS = GET_IS_Y([1,1,0])
    IE = GET_IE_Y([1,1,0])    
    FORALL(I=IS(1):IE(1),K=IS(3):IE(3))
       ZG2_Y(I,1   ,K) = (DG_CL_Y(I,1,K)*CLY(1,1) + DG_CL_Y(I,2,K)*CLY(1,2))
       ZG2_Y(I,NY+1,K) = (DG_CL_Y(I,1,K)*CLY(2,1) + DG_CL_Y(I,2,K)*CLY(2,2))
    END FORALL

    
    call transpose_y_to_x(zg2_y,zg_x)
    call tensor_product_i(px,zg_x,dg_x,xst,xen)
    call tensor_product_i_cl(ccx, dg_x, dg_cl_x, xst, xen )
    
    
    IS = GET_IS_X([0,0,0])
    IE = GET_IE_X([0,0,0])
    forall(j=is(2):ie(2),k=is(3):ie(3))
       dg_x(1   ,j,k) = dg_cl_x(1,j,k)*clx(1,1) + dg_cl_x(2,j,k)*clx(1,2)
       dg_x(nx+1,j,k) = dg_cl_x(1,j,k)*clx(2,1) + dg_cl_x(2,j,k)*clx(2,2)
    end forall
    
!    call mpi_comm_rank(mpi_comm_world,rank,ierr)
    
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

    
    
  end subroutine solver_diag_cart_iii_solve_scaled


  
  subroutine tensor_product_i(DX,FI,DFI,XST,XEN)
    implicit none
    INTEGER :: XST(3),XEN(3)
    complex(kind=8) ::  DX(XST(1):XEN(1),XST(1):XEN(1))
    complex(kind=8) ::  FI(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
    complex(kind=8) :: DFI(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
        
    complex(kind=8) ::  ALPHA,BETA
    INTEGER :: N(3),NI,NJ,NK

    N = XEN-XST + 1
    NI = N(1)
    NJ = N(2)
    NK = N(3)
    ALPHA = 1.
    BETA = 0.
    
    CALL ZGEMM( 'N', 'N', NI, NJ*NK, NI, ALPHA, DX, NI, FI, NI, BETA, DFI, NI )
    
  end subroutine tensor_product_i
  
    SUBROUTINE TENSOR_PRODUCT_I_CL(CC,FI,DFI,XST,XEN)
    IMPLICIT NONE
    INTEGER      :: XST(3),XEN(3)
    complex(KIND=8) :: CC(1:2,XST(1):XEN(1))
    complex(KIND=8) :: FI(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
    complex(KIND=8) :: DFI(1:2,XST(2):XEN(2),XST(3):XEN(3))
    
    complex(KIND=8) ::  ALPHA,BETA
    INTEGER      :: N(3),NI,NJ,NK

    N = XEN-XST + 1
    NI = N(1)
    NJ = N(2)
    NK = N(3)
    ALPHA = 1.
    BETA = 0.
    
    CALL ZGEMM( 'N', 'N', 2, NJ*NK, NI, ALPHA, CC, 2, FI, NI, BETA, DFI, 2 )
    
  END SUBROUTINE TENSOR_PRODUCT_I_CL

  
  
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
   
  
end module m_solver_diag_cart_iii
