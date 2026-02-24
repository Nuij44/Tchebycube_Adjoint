module m_solver_diag_cart_ttt
  use decomp_2d
  use m_solver_diag_base
  use m_lap1d_tcheby
  use m_boundary_conditions
use m_operator_tcheby
  use decomp_2d_fft
  use decomp_2d_mpi
  implicit none

  type , extends(t_solver_diag_base) :: t_solver_diag_cart_ttt
     type(t_lap1d_tcheby) :: lap_x,lap_y,lap_z
     real(kind=8),dimension(:,:,:),allocatable :: lbd
     real(kind=8),dimension(:,:,:),allocatable :: cnt
   contains 
     procedure :: INITIALISE => INIT_SOLVER_DIAG_CART_TTT
     procedure :: SET_BCS => SOLVER_DIAG_CART_TTT_SET_BCS
     procedure :: SET_BVS => SOLVER_DIAG_CART_TTT_SET_BVS_VALUES_TIME
     procedure :: SET_BVS_TIME => SOLVER_DIAG_CART_TTT_SET_BVS_VALUES_UNSTEADY
     procedure :: SET_UBVS => SOLVER_DIAG_CART_TTT_SET_BVS_VALUES_TIME_UNSTEADY
     procedure :: UPDATE_BCS => SOLVER_DIAG_CART_TTT_UPDATE_BCS_CONTRIB
     procedure :: SOLVE => SOLVER_DIAG_CART_TTT_SOLVE
     procedure :: SET_PARAMS => SOLVER_DIAG_CART_TTT_SET_PARAMS
     procedure :: SET_PARAMS_VEC => SOLVER_DIAG_CART_TTT_SET_PARAMS_NU
  end type t_solver_diag_cart_ttt

  real(dp),DIMENSION(:,:,:),ALLOCATABLE :: DG1X, DG2X, DG1Y, DG2Y, DG1Z, DG2Z ,DG_CL_X,DG_CL_Y,DG_CL_Z
  PRIVATE :: DG1X, DG2X, DG1Y, DG2Y, DG1Z, DG2Z 

 
  
contains 
  SUBROUTINE SOLVER_DIAG_CART_TTT_SET_PARAMS(THIS,DIM,NU,sigma)
    IMPLICIT NONE
    CLASS(T_SOLVER_DIAG_CART_TTT) :: THIS
    integer                       :: dim(3)
    real(dp)                      :: nu
    real(dp)                      :: sigma
    
    THIS%N = DIM
    THIS%NU = nu
    THIS%sigma = sigma
  end SUBROUTINE SOLVER_DIAG_CART_TTT_SET_PARAMS

  SUBROUTINE SOLVER_DIAG_CART_TTT_SET_PARAMS_nu(THIS,DIM,NU,sigma)
    IMPLICIT NONE
    CLASS(T_SOLVER_DIAG_CART_TTT) :: THIS
    integer                       :: dim(3)
    real(dp)                      :: nu(3)
    real(dp)                      :: sigma
    integer :: ierr
    
    THIS%N = DIM
    THIS%NU = nu
    THIS%sigma = sigma

    
  end SUBROUTINE SOLVER_DIAG_CART_TTT_SET_PARAMS_NU


  

    
  SUBROUTINE SOLVER_DIAG_CART_TTT_SET_BCS(THIS,AXIS,BCS_MINUS,BCS_PLUS)
    use m_boundary_conditions
      use m_lap1d_tcheby
      IMPLICIT NONE
      CLASS(T_SOLVER_DIAG_CART_TTT) :: THIS
      INTEGER                       :: AXIS
      real(dp)                      :: BCS_MINUS(2)
      real(dp)                      :: BCS_PLUS(2)
      
      TYPE(DECOMP_INFO) :: PH

      
      CALL GET_DECOMP_INFO(PH)

      IF (AXIS==1) THEN
         CALL THIS%BCS_X%INITIALIZE(AXIS,BCS_MINUS(1),BCS_MINUS(2),BCS_PLUS(1),BCS_PLUS(2),PH%XST,PH%XEN,THIS%N) ! X
      ELSE IF (AXIS==2) THEN
         CALL THIS%BCS_Y%INITIALIZE(AXIS,BCS_MINUS(1),BCS_MINUS(2),BCS_PLUS(1),BCS_PLUS(2),PH%XST,PH%XEN,THIS%N) ! X
      ELSE IF (AXIS==3) THEN
         CALL THIS%BCS_Z%INITIALIZE(AXIS,BCS_MINUS(1),BCS_MINUS(2),BCS_PLUS(1),BCS_PLUS(2),PH%XST,PH%XEN,THIS%N) ! X
      END IF
      
    END SUBROUTINE SOLVER_DIAG_CART_TTT_SET_BCS
    
  
  SUBROUTINE INIT_SOLVER_DIAG_CART_TTT(THIS,OPX,OPY,OPZ)
    use m_boundary_conditions
    use m_lap1d_tcheby
    IMPLICIT NONE
    CLASS(T_SOLVER_DIAG_CART_TTT) :: THIS
    TYPE(T_OPERATOR_TCHEBY) ::  OPX,OPY,OPZ
    INTEGER N(3)
    TYPE(DECOMP_INFO) :: PH
    real(kind=8), dimension(3) :: nu
    
    real(kind=8), dimension(3) :: alfa_l, beta_l
    real(kind=8), dimension(3) :: alfa_r, beta_r
    integer :: i,j,k
    
    CALL GET_DECOMP_INFO(PH)
    
    ALFA_L = [ THIS%BCS_X%ALPHA_L, THIS%BCS_Y%ALPHA_L, THIS%BCS_Z%ALPHA_L ]
    BETA_L = [ THIS%BCS_X%BETA_L , THIS%BCS_Y%BETA_L , THIS%BCS_Z%BETA_L ]

    ALFA_R = [ THIS%BCS_X%ALPHA_R, THIS%BCS_Y%ALPHA_R, THIS%BCS_Z%ALPHA_R ]
    BETA_R = [ THIS%BCS_X%BETA_R , THIS%BCS_Y%BETA_R , THIS%BCS_Z%BETA_R ]

    nu = this%nu

    ! INITIALISATION DE CHAQUE DIRECTION
    CALL THIS%LAP_X%INITIALIZE(OPX,NU(1),ALFA_L(1),BETA_L(1),ALFA_R(1),BETA_R(1),THIS%N(1)) 
    CALL THIS%LAP_Y%INITIALIZE(OPY,NU(2),ALFA_L(2),BETA_L(2),ALFA_R(2),BETA_R(2),THIS%N(2))
    CALL THIS%LAP_Z%INITIALIZE(OPZ,NU(3),ALFA_L(3),BETA_L(3),ALFA_R(3),BETA_R(3),THIS%N(3))
    
    CALL ALLOC_Z( THIS%LBD , OPT_GLOBAL=.TRUE. )
    N = THIS%N
    DO K = PH%ZST(3),PH%ZEN(3)
       DO J = PH%ZST(2),PH%ZEN(2) 
          DO I=PH%ZST(1),PH%ZEN(1)
             THIS%LBD(I,J,K) =  ( this%sigma + (THIS%LAP_X%LBD(I) + THIS%LAP_Y%LBD(J) +  THIS%LAP_Z%LBD(K) ) )**(-1)
             IF (I==     1) THIS%LBD(I,J,K)=1
             IF (I==N(1)+1) THIS%LBD(I,J,K)=1
             IF (J==     1) THIS%LBD(I,J,K)=1
             IF (J==N(2)+1) THIS%LBD(I,J,K)=1
             IF (K==     1) THIS%LBD(I,J,K)=1
             IF (K==N(3)+1) THIS%LBD(I,J,K)=1
          END DO
       END DO
    END DO
    

    
    CALL ALLOC_X( THIS%CNT , OPT_GLOBAL=.TRUE. )

    if (.not.allocated(DG1X)) CALL ALLOC_X( DG1X , OPT_GLOBAL=.TRUE. )
    if (.not.allocated(DG2X)) CALL ALLOC_X( DG2X , OPT_GLOBAL=.TRUE. )

    if (.not.allocated(DG1Y)) CALL ALLOC_Y( DG1Y , OPT_GLOBAL=.TRUE. )
    if (.not.allocated(DG2Y)) CALL ALLOC_Y( DG2Y , OPT_GLOBAL=.TRUE. )

    if (.not.allocated(DG1Z)) CALL ALLOC_Z( DG1Z , OPT_GLOBAL=.TRUE. )
    if (.not.allocated(DG2Z)) CALL ALLOC_Z( DG2Z , OPT_GLOBAL=.TRUE. )

    if (.not.allocated(DG_CL_X)) allocate( DG_CL_X(1:2,PH%XST(2):PH%XEN(2),PH%XST(3):PH%XEN(3)) )
    if (.not.allocated(DG_CL_Y)) allocate( DG_CL_Y(PH%YST(1):PH%YEN(1),1:2,PH%YST(3):PH%YEN(3)) )
    if (.not.allocated(DG_CL_Z)) allocate( DG_CL_Z(PH%ZST(1):PH%ZEN(1),PH%ZST(2):PH%ZEN(2),1:2) )
    
  END SUBROUTINE INIT_SOLVER_DIAG_CART_TTT

  SUBROUTINE SOLVER_DIAG_CART_TTT_SOLVE(THIS,FI,SFI)
    use decomp_2d
    implicit none
    CLASS(T_SOLVER_DIAG_CART_TTT) :: THIS
    real(dp),DIMENSION(:,:,:),ALLOCATABLE :: FI
    real(dp),DIMENSION(:,:,:),ALLOCATABLE :: SFI

    TYPE(DECOMP_INFO) :: PH

    CALL GET_DECOMP_INFO(PH)

    call SOLVER_DIAG_CART_TTT_ADD_CONTRIB( &
         this%lap_x%fl,  this%lap_x%fr, this%bcs_x%bv_l,  this%bcs_x%bv_r, ph%xst, ph%xen, this%n(1), &
         this%lap_y%fl,  this%lap_y%fr, this%bcs_y%bv_l,  this%bcs_y%bv_r, ph%xst, ph%xen, this%n(2), &
         this%lap_z%fl,  this%lap_z%fr, this%bcs_z%bv_l,  this%bcs_z%bv_r, ph%xst, ph%xen, this%n(3), &
         THIS%cnt )
    
    CALL KERNEL_SOLVE_TTT(&
         THIS%LAP_X%P, THIS%LAP_X%PM1, THIS%LAP_X%LBD, THIS%LAP_X%CL, THIS%BCS_X%BV_L,  THIS%BCS_X%BV_R, PH%XST, PH%XEN, THIS%N(1), &
         THIS%LAP_y%P, THIS%LAP_y%PM1, THIS%LAP_y%LBD, THIS%LAP_y%CL, THIS%BCS_y%BV_L,  THIS%BCS_y%BV_R, PH%YST, PH%YEN, THIS%N(2), &
         THIS%LAP_z%P, THIS%LAP_z%PM1, THIS%LAP_z%LBD, THIS%LAP_z%CL, THIS%BCS_z%BV_L,  THIS%BCS_z%BV_R, PH%ZST, PH%ZEN, THIS%N(3), &
         THIS%SIGMA,THIS%LBD, THIS%CNT, THIS%LAP_X%CC, THIS%LAP_Y%CC, THIS%LAP_Z%CC, FI, SFI, DG1X, DG2X, DG1Y, DG2Y, DG1Z, DG2Z, DG_CL_X, DG_CL_Y,DG_CL_Z)
    
    
  END SUBROUTINE SOLVER_DIAG_CART_TTT_SOLVE
  
  
  subroutine SOLVER_DIAG_CART_TTT_ADD_CONTRIB( &
       fxl,fxr,bvxl,bvxr,xst,xen,nx, &
       fyl,fyr,bvyl,bvyr,yst,yen,ny, &
       fzl,fzr,bvzl,bvzr,zst,zen,nz, &
       fxyz )
    implicit none
      integer :: xst(3),xen(3),nx
      integer :: yst(3),yen(3),ny
      integer :: zst(3),zen(3),nz
      real(dp) :: fxl(1:nx+1), fxr(1:nx+1), bvxl(xst(2):xen(2),xst(3):xen(3)), bvxr(xst(2):xen(2),xst(3):xen(3))
      real(dp) :: fyl(1:ny+1), fyr(1:ny+1), bvyl(xst(1):xen(1),xst(3):xen(3)), bvyr(xst(1):xen(1),xst(3):xen(3))
      real(dp) :: fzl(1:nz+1), fzr(1:nz+1), bvzl(xst(1):xen(1),xst(2):xen(2)), bvzr(xst(1):xen(1),xst(2):xen(2))
      real(dp) :: fxyz(xst(1):xen(1),xst(2):xen(2),xst(3):xen(3))
      
      integer :: i,j,k

      forall(i=xst(1):xen(1),j=xst(2):xen(2),k=xst(3):xen(3))
         fxyz(i,j,k) = &
              + fxl(i)*bvxl(j,k) + fxr(i)*bvxr(j,k) &
              + fyl(j)*bvyl(i,k) + fyr(j)*bvyr(i,k) &
              + fzl(k)*bvzl(i,j) + fzr(k)*bvzr(i,j)
      end forall
      
    end subroutine SOLVER_DIAG_CART_TTT_ADD_CONTRIB
    


    SUBROUTINE SOLVER_DIAG_CART_TTT_SET_BVS_VALUES_UNSTEADY(THIS,time,MSH,AXIS,UDF_MINUS,UDF_PLUS)
      use m_boundary_conditions
      use m_lap1d_tcheby
      
      IMPLICIT NONE
      CLASS(T_SOLVER_DIAG_CART_TTT) :: THIS
      TYPE(T_MESH_BASE)             :: MSH(3)
      INTEGER                       :: AXIS
      PROCEDURE(UDF_TIMESPACE)      :: UDF_plus
      PROCEDURE(UDF_TIMESPACE)      :: UDF_minus
      
      TYPE(DECOMP_INFO) :: PH
      REAL(KIND=8) :: TIME
      
      CALL GET_DECOMP_INFO(PH)



    if (axis==1) then
       
       CALL SOLVER_DIAG_TTT_SET_UNSTEADY_BOUNDARY_FIELD(&
            MESH(1)%X, MESH(2)%X, MESH(3)%X, TC, &
            UDF_MINUS,UDF_PLUS,THIS%BCS_X%BV_L,THIS%BCS_X%BV_R,&
            PH%XST,PH%XEN,PH%YST,PH%YEN,PH%ZST,PH%ZEN )
    else
       
       call MPI_Abort(mpi_comm_world,0,ierr)
    end if

    if (axis==2) then
       
       CALL SOLVER_DIAG_TTT_SET_UNSTEADY_BOUNDARY_FIELD(&
            MESH(1)%X, MESH(2)%X, MESH(3)%X, TC, &
            UDF_MINUS,UDF_PLUS,THIS%BCS_Y%BV_L,THIS%BCS_Y%BV_R,&
            PH%XST,PH%XEN,PH%YST,PH%YEN,PH%ZST,PH%ZEN )
    else
       
       call MPI_Abort(mpi_comm_world,0,ierr)
    end if

    if (axis==3) then
       
       CALL SOLVER_DIAG_TTT_SET_UNSTEADY_BOUNDARY_FIELD(&
            MESH(1)%X, MESH(2)%X, MESH(3)%X, TC, &
            UDF_MINUS,UDF_PLUS,THIS%BCS_Z%BV_L,THIS%BCS_Z%BV_R,&
            PH%XST,PH%XEN,PH%YST,PH%YEN,PH%ZST,PH%ZEN )
    else
       
       call MPI_Abort(mpi_comm_world,0,ierr)
    end if
    
  end subroutine SOLVER_DIAG_CART_TTT_SET_BVS_VALUES_UNSTEADY
    
  subroutine solver_diag_ttt_set_unsteady_boundary_field(x,y,z,tc,udf_bc_left,udf_bc_right,bvzl,bvzr,xst,xen,yst,yen,zst,zen)
    use m_numerics
    implicit none
    integer :: xst(3),xen(3)
    integer :: yst(3),yen(3)
    integer :: zst(3),zen(3)
    real(kind=8) :: tc
    real(kind=8), dimension(:),allocatable :: x,y,z
    real(kind=8), dimension(xst(1):xen(1),xst(2):xen(2)) :: bvzl
    real(kind=8), dimension(xst(1):xen(1),xst(2):xen(2)) :: bvzr
    PROCEDURE(UDF_TIMESPACE)      :: UDF_BC_LEFT
    PROCEDURE(UDF_TIMESPACE)      :: UDF_BC_RIGHT
    integer :: i,j
    integer :: nx,ny,nz,n(3)
    real(kind=8) :: zmin,zmax

    n = xen - xst
    nz = zen(3) - zst(3) + 1
    
    zmin = z(1)
    zmax = z(nz)
    
    
    forall(i=xst(1):xen(1),j=xst(2):xen(2))
       bvzl(i,j) = udf_bc_left(tc,X(I),Y(J),ZMIN)
    end forall

    
    forall(i=xst(1):xen(1),j=xst(2):xen(2))
       bvzr(i,j) = udf_bc_right(tc,X(I),Y(J),ZMAX)
    end forall
    
  end subroutine solver_diag_hhi_set_unsteady_boundary_field
    


    SUBROUTINE SOLVER_DIAG_CART_TTT_SET_BVS_VALUES_TIME(THIS,MSH,AXIS,UDF_MINUS,UDF_PLUS)
      use m_boundary_conditions
      use m_lap1d_tcheby
      
      IMPLICIT NONE
      CLASS(T_SOLVER_DIAG_CART_TTT) :: THIS
      TYPE(T_MESH_BASE)             :: MSH(3)
      INTEGER                       :: AXIS
      PROCEDURE(UDF_TIMESPACE)      :: UDF_plus
      PROCEDURE(UDF_TIMESPACE)      :: UDF_minus
      
      TYPE(DECOMP_INFO) :: PH
      REAL(KIND=8) :: TIME
      
      CALL GET_DECOMP_INFO(PH)



      
      if (axis==1) then
         CALL THIS%BCS_X%SET_BCS(AXIS,-1,TIME,MSH,UDF_MINUS,PH%XST,PH%XEN,THIS%N)
         CALL THIS%BCS_X%SET_BCS(AXIS,+1,TIME,MSH,UDF_PLUS,PH%XST,PH%XEN,THIS%N)
      else if (axis==2) then
         CALL THIS%BCS_Y%SET_BCS(AXIS,-1,TIME,MSH,UDF_MINUS,PH%XST,PH%XEN,THIS%N)
         CALL THIS%BCS_Y%SET_BCS(AXIS,+1,TIME,MSH,UDF_PLUS,PH%XST,PH%XEN,THIS%N)
      else if (axis==3) then
         CALL THIS%BCS_Z%SET_BCS(AXIS,-1,TIME,MSH,UDF_MINUS,PH%XST,PH%XEN,THIS%N)
         CALL THIS%BCS_Z%SET_BCS(AXIS,+1,TIME,MSH,UDF_PLUS,PH%XST,PH%XEN,THIS%N)
      end if

      
      
    end subroutine SOLVER_DIAG_CART_TTT_SET_BVS_VALUES_TIME
    
    SUBROUTINE SOLVER_DIAG_CART_TTT_UPDATE_BCS_CONTRIB(this)
      IMPLICIT NONE
      CLASS(T_SOLVER_DIAG_CART_TTT) :: THIS
      TYPE(DECOMP_INFO) :: PH
      CALL GET_DECOMP_INFO(PH)
      call SOLVER_DIAG_CART_TTT_ADD_CONTRIB( &
           this%lap_x%fl,  this%lap_x%fr, this%bcs_x%bv_l,  this%bcs_x%bv_r, ph%xst, ph%xen, this%n(1), &
           this%lap_y%fl,  this%lap_y%fr, this%bcs_y%bv_l,  this%bcs_y%bv_r, ph%xst, ph%xen, this%n(2), &
           this%lap_z%fl,  this%lap_z%fr, this%bcs_z%bv_l,  this%bcs_z%bv_r, ph%xst, ph%xen, this%n(3), &
           THIS%cnt )

    END SUBROUTINE SOLVER_DIAG_CART_TTT_UPDATE_BCS_CONTRIB
    
    !>
    SUBROUTINE KERNEL_SOLVE_TTT(&
         PX, PXM1, LBD_X, CL_X, BC_VAL_L_X, BC_VAL_R_X, XST, XEN, NX, &
         PY, PYM1, LBD_Y, CL_Y, BC_VAL_L_Y, BC_VAL_R_Y, YST, YEN, NY, &
         PZ, PZM1, LBD_Z, CL_Z, BC_VAL_L_Z, BC_VAL_R_Z, ZST, ZEN, NZ, &
         sigma, LBD_XYZ, F_XYZ, CC_X, CC_y, CC_z, FI, SFI, DG1X, DG2X, DG1Y, DG2Y, DG1Z, DG2Z ,DG_CL_X,DG_CL_Y,DG_CL_Z)
      use mpi
      use decomp_2d
      use decomp_2d_mpi
      implicit none
      INTEGER :: XST(3) , XEN(3) , NX
      INTEGER :: YST(3) , YEN(3) , NY
      INTEGER :: ZST(3) , ZEN(3) , NZ
      
      real(dp) :: PX(1:NX+1,1:NX+1) , PXM1(1:NX+1,1:NX+1) , LBD_X(1:NX+1)
      real(dp) :: PY(1:NY+1,1:NY+1) , PYM1(1:NY+1,1:NY+1) , LBD_Y(1:NY+1)
      real(dp) :: PZ(1:NZ+1,1:NZ+1) , PZM1(1:NZ+1,1:NZ+1) , LBD_Z(1:NZ+1)
      real(dp) :: CL_X(2,2)
      real(dp) :: CL_Y(2,2)
      real(dp) :: CL_Z(2,2)
      
      real(dp) :: BC_VAL_L_X(XST(2):XEN(2),XST(3):XEN(3)),BC_VAL_R_X(XST(2):XEN(2),XST(3):XEN(3))
      real(dp) :: BC_VAL_L_Y(XST(1):XEN(1),XST(3):XEN(3)),BC_VAL_R_Y(XST(1):XEN(1),XST(3):XEN(3))
      real(dp) :: BC_VAL_L_Z(XST(1):XEN(1),XST(2):XEN(2)),BC_VAL_R_Z(XST(1):XEN(1),XST(2):XEN(2))
      
      real(dp) ::  FI(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
      real(dp) :: SFI(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))

      real(dp) :: sigma
      real(dp) :: lbd_xyz(ZST(1):ZEN(1),ZST(2):ZEN(2),ZST(3):ZEN(3))
      real(dp) ::   f_xyz(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
      
      REAL(kind=8) :: DG1X(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
      REAL(kind=8) :: DG2X(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
      
      REAL(kind=8) :: DG1Y(YST(1):YEN(1),YST(2):YEN(2),YST(3):YEN(3))
      REAL(kind=8) :: DG2Y(YST(1):YEN(1),YST(2):YEN(2),YST(3):YEN(3))

      REAL(kind=8) :: DG1Z(ZST(1):ZEN(1),ZST(2):ZEN(2),ZST(3):ZEN(3))
      REAL(kind=8) :: DG2Z(ZST(1):ZEN(1),ZST(2):ZEN(2),ZST(3):ZEN(3))

      real(dp) :: CC_X(1:2,1:NX+1)
      real(dp) :: CC_Y(1:2,1:NY+1)
      real(dp) :: CC_Z(1:2,1:NZ+1)
      
      real(dp) :: DG_CL_X(1:2,XST(2):XEN(2),XST(3):XEN(3))
      real(dp) :: DG_CL_Y(YST(1):YEN(1),1:2,YST(3):YEN(3))
      real(dp) :: DG_CL_Z(ZST(1):ZEN(1),ZST(2):ZEN(2),1:2)

      integer :: i,j,k,is(3),ie(3)


      DO K = ZST(3),ZEN(3)
         DO J = ZST(2),ZEN(2) 
            DO I=ZST(1),ZEN(1)
               lbd_xyz(I,J,K) = (sigma + LBD_x(I) + LBD_y(J) +  LBD_z(K))**(-1)
               IF (I==   1) LBD_XYZ(I,J,K) = 1
               IF (I==NX+1) LBD_XYZ(I,J,K) = 1
               IF (J==   1) LBD_XYZ(I,J,K) = 1
               IF (J==NY+1) LBD_XYZ(I,J,K) = 1
               IF (K==   1) LBD_XYZ(I,J,K) = 1
               IF (K==NZ+1) LBD_XYZ(I,J,K) = 1
            END DO
         END DO
      END DO
      
      DG1X = SFI - F_XYZ

      !> 
      IS = GET_IS_X([0,0,0])
      IE = GET_IE_X([0,0,0])
      
      FORALL(J=IS(2):IE(2),K=IS(3):IE(3))
         DG1X(   1,J,K) = BC_VAL_L_X(J,K)
      END FORALL
      
      FORALL(J=IS(2):IE(2),K=IS(3):IE(3))
         DG1X(NX+1,J,K) = BC_VAL_R_X(J,K)
      END FORALL
      
      
      IS = GET_IS_X([1,1,0])
      IE = GET_IE_X([1,1,0])
      
      IF (XST(2)==1) THEN
         FORALL(I=IS(1):IE(1),K=IS(3):IE(3))
          DG1X(I,   1,K) = BC_VAL_L_Y(I,K)
       END FORALL
       
    END IF
    
    IF (XEN(2)==NY+1) THEN
       
       FORALL(I=IS(1):IE(1),K=IS(3):IE(3))
          DG1X(I,NY+1,K) = BC_VAL_R_Y(I,K)
       END FORALL
    END IF


    IS = GET_IS_X([1,1,1])
    IE = GET_IE_X([1,1,1])
    IF (XST(3)==1) THEN
       
       FORALL(I=IS(1):IE(1),J=IS(2):IE(2))
          DG1X(I,J,   1) = BC_VAL_L_Z(I,J)
       END FORALL
       
    END IF
    
    IF (XEN(3)==NZ+1) THEN
       
       FORALL(I=IS(1):IE(1),J=IS(2):IE(2))
          DG1X(I,J,NZ+1) = BC_VAL_R_Z(I,J)
       END FORALL
       
    END IF
    
    
    CALL TENSOR_dpODUCT_I( PXM1, DG1X, DG2X, XST, XEN, NX )
    CALL TRANSPOSE_X_TO_Y( DG2X, DG1Y )
    
    CALL TENSOR_dpODUCT_J( PYM1, DG1Y, DG2Y, YST, YEN, NY )
    
    CALL TRANSPOSE_Y_TO_Z( DG2Y, DG1Z )
    CALL TENSOR_dpODUCT_K( PZM1, DG1Z, DG2Z, ZST, ZEN, NZ )
    
    
    
    FORALL(I=ZST(1):ZEN(1),J=ZST(2):ZEN(2),K=ZST(3):ZEN(3))
       DG1Z(I,J,K) = LBD_XYZ(I,J,K)*DG2Z(I,J,K)
    END FORALL
    
    
    CALL TENSOR_dpODUCT_K( PZ   , DG1Z, DG2Z, ZST, ZEN, NZ )
    call TENSOR_dpODUCT_CL_K(CC_Z, DG2Z, DG_CL_Z, ZST, ZEN , NZ )

    IS = GET_IS_Z([1,1,1])
    IE = GET_IE_Z([1,1,1])
    !> Z-PENCIL
    FORALL(I=IS(1):IE(1),J=IS(2):IE(2))
       DG2Z(I,J,   1) = DG_CL_Z(I,J,1)*CL_Z(1,1) + DG_CL_Z(I,J,2)*CL_Z(1,2)
       DG2Z(I,J,NZ+1) = DG_CL_Z(I,J,1)*CL_Z(2,1) + DG_CL_Z(I,J,2)*CL_Z(2,2)
    END FORALL
    
    CALL TRANSPOSE_Z_TO_Y( DG2Z , DG1Y )
    
    
    !> Y-PENCIL
    CALL TENSOR_dpODUCT_J( PY  , DG1Y, DG2Y, YST, YEN, NY )
    CALL TENSOR_dpODUCT_CL_J(CC_Y, DG2Y , DG_CL_Y, YST,YEN,NY)
    
    IS = GET_IS_Y([1,1,0])
    IE = GET_IE_Y([1,1,0])
    FORALL(I=IS(1):IE(1),K=IS(3):IE(3))
       DG2Y(I,   1,K) = DG_CL_Y(I,1,K)*CL_Y(1,1) + DG_CL_Y(I,2,K)*CL_Y(1,2)
       DG2Y(I,NY+1,K) = DG_CL_Y(I,1,K)*CL_Y(2,1) + DG_CL_Y(I,2,K)*CL_Y(2,2)
    END FORALL
    
    
    CALL TRANSPOSE_Y_TO_X( DG2Y, DG1X )
    
    !> X-PENCIL
    
    CALL TENSOR_dpODUCT_I( PX  , DG1X , FI, XST, XEN, NX )
    CALL TENSOR_dpODUCT_CL_I(CC_X,FI,DG_CL_X,XST,XEN,NX)
    IS = GET_IS_X([0,0,0])
    IE = GET_IE_X([0,0,0])
    FORALL(J=IS(2):IE(2),K=IS(3):IE(3))
       FI(   1,J,K) = DG_CL_X(1,J,K)*CL_X(1,1) + DG_CL_X(2,J,K)*CL_X(1,2)
       FI(NX+1,J,K) = DG_CL_X(1,J,K)*CL_X(2,1) + DG_CL_X(2,J,K)*CL_X(2,2)
    END FORALL
    
    
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
    
    
    
    SUBROUTINE TENSOR_dpODUCT_CL_I(CC,FI,DFI,XST,XEN,NX)
      IMPLICIT NONE
        INTEGER :: XST(3),XEN(3),NX
        real(dp) ::  CC(1:2,1:NX+1)
        real(dp) ::  FI(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
        real(dp) :: DFI(1:2,XST(2):XEN(2),XST(3):XEN(3))
        
        real(dp) ::  ALPHA,BETA
        INTEGER :: N(3),NI,NJ,NK
        
        N = XEN-XST + 1
        NI = N(1)
        NJ = N(2)
        NK = N(3)
        ALPHA = 1._dp
        BETA = 0.0_dp

        CALL DGEMM( 'N', 'N', 2, NJ*NK, NI, ALPHA, CC, 2, FI, NI, BETA, DFI, 2 )
        
      END SUBROUTINE TENSOR_dpODUCT_CL_I
      
      SUBROUTINE TENSOR_dpODUCT_CL_J(CC,FI,DFI,YST,YEN,NY)
        IMPLICIT NONE
        INTEGER  :: YST(3),YEN(3),NY
        real(dp) :: CC(1:2,1:NY+1)
        real(dp) ::  FI(YST(1):YEN(1),YST(2):YEN(2),YST(3):YEN(3))
        real(dp) :: DFI(YST(1):YEN(1),1:2,YST(3):YEN(3))
        real(dp) ::  ALPHA,BETA
        INTEGER  :: K,N(3),NI,NJ,NK
        
        ALPHA = 1._dp
        BETA = 0.0_dp
        N = YEN - YST + 1
        NI = N(1)
        NJ = N(2)
        DO K=YST(3),YEN(3)
           CALL DGEMM('N','T',NI,2,NJ,ALPHA,FI(YST(1),YST(2),K),NI,CC,2,BETA,DFI(YST(1),1,K),NI)
        END DO
        
      END SUBROUTINE TENSOR_dpODUCT_CL_J


      SUBROUTINE TENSOR_dpODUCT_CL_K(CC,FI,DFI,ZST,ZEN,NZ)
        IMPLICIT NONE
        INTEGER :: ZST(3),ZEN(3),NZ
        real(dp) :: CC(1:2,1:NZ+1)
        real(dp) ::  FI(ZST(1):ZEN(1),ZST(2):ZEN(2),ZST(3):ZEN(3))
        real(dp) :: DFI(ZST(1):ZEN(1),ZST(2):ZEN(2),1:2)
        
        INTEGER :: N(3),NI,NJ,NK
        real(dp) ::  ALPHA,BETA
        
        ALPHA = 1._dp
        BETA = 0._dp
        N = ZEN - ZST + 1
        NI = N(1)
        NJ = N(2)
        NK = N(3) ! NZ + 1 car Z-pencil
        CALL DGEMM('N','T',NI*NJ,2,NK, ALPHA, FI, NI*NJ , CC, 2 , BETA , DFI, NI*NJ )
        
      END SUBROUTINE TENSOR_dpODUCT_CL_K
      
      
      subroutine tensor_dpoduct_i(DX,FI,DFI,XST,XEN,NX)
        implicit none
        INTEGER :: XST(3),XEN(3),NX
        real(dp) ::  DX(1:NX+1,1:NX+1)
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
      
      SUBROUTINE TENSOR_dpODUCT_J(DY,FI,DFI,YST,YEN,NY)
        IMPLICIT NONE
        INTEGER :: YST(3),YEN(3),NY
        real(dp) :: DY(1:NY+1,1:NY+1)
        real(dp) ::  FI(YST(1):YEN(1),YST(2):YEN(2),YST(3):YEN(3))
        real(dp) :: DFI(YST(1):YEN(1),YST(2):YEN(2),YST(3):YEN(3))
        real(dp) ::  ALPHA,BETA
        INTEGER :: K,N(3),NI,NJ,NK,is,js,ks,ke,info
        
        ALPHA = 1._dp
        BETA = 0.0_dp
        N = YEN - YST + 1
        NI = N(1)
        NJ = N(2)
        NK = N(3)

        ni = ubound(dfi,1)-lbound(dfi,1)+1
        nj = ubound(dfi,2)-lbound(dfi,2)+1
        nk = ubound(dfi,3)-lbound(dfi,3)+1

        is=lbound(dfi,1)
        js=lbound(dfi,2)
        ks=lbound(dfi,3)
        ke=ubound(dfi,3)
        
        do k=ks,ke
           CALL DGEMM('N','T',NI,NJ,NJ, ALPHA ,FI(is,js,k), NI , DY, NJ , beta , DFI(is,js,k), NI )
        end do
        
!        DO K=YST(3),YEN(3)
!           CALL DGEMM('N','T',NI,NJ,NJ,ALPHA,FI(YST(1),YST(2),K),NI,DY,NJ,BETA,DFI(YST(1),YST(2),K),NI)
!        END DO
        
      END SUBROUTINE TENSOR_dpODUCT_J


      
      
      SUBROUTINE TENSOR_dpODUCT_K(DZ,FI,DFI,ZST,ZEN,NZ)
        IMPLICIT NONE
        INTEGER :: ZST(3),ZEN(3),NZ
        real(dp) :: DZ(1:NZ+1,1:NZ+1)
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
      
      
    end subroutine KERNEL_SOLVE_TTT

    
    
  end module m_solver_diag_cart_ttt
