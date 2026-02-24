module m_solver_diag_cart_fft_cplx
  use decomp_2d
  use m_operator_tcheby
  use m_solver_diag_base
  use m_lap1d_tcheby
  use m_lap1d_fourier
  use m_boundary_conditions
  use decomp_2d_fft
  use decomp_2d_mpi
  implicit none
  
  type , extends(t_solver_diag_base) :: t_solver_diag_cart_fft_cplx

     type(t_lap1d_fourier) :: lap_x
     type(t_lap1d_fourier) :: lap_y
     type(t_lap1d_tcheby)  :: lap_z
     
     complex(kind=8),dimension(:,:,:),allocatable :: lbd
     real(kind=8),dimension(:,:,:),allocatable :: cnt_u
     real(kind=8),dimension(:,:,:),allocatable :: cnt_v
     
     real(kind=8) :: omega
     
     !> fork for the complex
     COMPLEX(DP),DIMENSION(:,:),ALLOCATABLE :: PX,PXM1
     COMPLEX(DP),DIMENSION(:,:),ALLOCATABLE :: PY,PYM1
     COMPLEX(DP),DIMENSION(:,:),ALLOCATABLE :: PZ,PZM1,CC_Z

     ! we declare new boundary conditions
     type(t_boundary_conditions) :: bcs_u_z
     type(t_boundary_conditions) :: bcs_v_z
     
     !>
   contains 
     procedure :: INITIALISE => INIT_SOLVER_DIAG_CART_FFT_CPLX
     procedure :: SET_BCS => SOLVER_DIAG_CART_FFT_CPLX_SET_BCS
     procedure :: SET_BVS => SOLVER_DIAG_CART_FFT_CPLX_SET_BVS_VALUES_TIME
     procedure :: SET_BVS_U => SOLVER_DIAG_CART_FFT_CPLX_SET_BVS_VALUES_TIME_REAL
     procedure :: SET_BVS_V => SOLVER_DIAG_CART_FFT_CPLX_SET_BVS_VALUES_TIME_IMAG

     procedure :: UPDATE_BCS => SOLVER_DIAG_CART_FFT_CPLX_UPDATE_BCS_CONTRIB
     procedure :: SOLVE => SOLVER_DIAG_CART_FFT_CPLX_SOLVE
     procedure :: SET_PARAMS => SOLVER_DIAG_CART_FFT_CPLX_SET_PARAMS
  end type t_solver_diag_cart_fft_cplx

  
  COMPLEX(DP),DIMENSION(:,:,:),ALLOCATABLE,PRIVATE :: DG1X, DG2X, DG1Y, DG2Y, DG1Z, DG2Z ,DG_CL_X,DG_CL_Y,DG_CL_Z
  COMPLEX(DP),DIMENSION(:,:,:),ALLOCATABLE,PRIVATE :: FI,SFI
  
contains 

  SUBROUTINE SOLVER_DIAG_CART_FFT_CPLX_SET_PARAMS(THIS,DIM,NU,SIGMA,OMEGA)
    IMPLICIT NONE
    CLASS(T_SOLVER_DIAG_CART_FFT_CPLX) :: THIS
    INTEGER                            :: DIM(3)
    REAL(DP)                           :: NU
    REAL(DP)                           :: SIGMA
    REAL(DP)                           :: OMEGA
    
    THIS%N = DIM
    THIS%NU = NU
    THIS%OMEGA = OMEGA
    THIS%SIGMA = SIGMA
    
  END SUBROUTINE SOLVER_DIAG_CART_FFT_CPLX_SET_PARAMS
  
  SUBROUTINE SOLVER_DIAG_CART_FFT_CPLX_SET_PARAMS_PHYSICAL(THIS,SIGMA,OMEGA)
    IMPLICIT NONE
    CLASS(T_SOLVER_DIAG_CART_FFT_CPLX) :: THIS
    INTEGER                            :: DIM(3)
    REAL(DP)                           :: SIGMA
    REAL(DP)                           :: OMEGA

    THIS%OMEGA = OMEGA
    THIS%SIGMA = SIGMA
    
  END SUBROUTINE SOLVER_DIAG_CART_FFT_CPLX_SET_PARAMS_PHYSICAL
  
  SUBROUTINE SOLVER_DIAG_CART_FFT_CPLX_SET_BCS(THIS,AXIS,BCS_MINUS,BCS_PLUS)
    USE M_BOUNDARY_CONDITIONS
      USE M_LAP1D_TCHEBY
      IMPLICIT NONE
      CLASS(T_SOLVER_DIAG_CART_FFT_CPLX) :: THIS
      INTEGER                            :: AXIS
      REAL(DP)                           :: BCS_MINUS(2)
      REAL(DP)                           :: BCS_PLUS(2)
      INTEGER :: IERR
      TYPE(DECOMP_INFO) :: PH
      
      CALL GET_DECOMP_INFO(PH)
      IF ((AXIS==1).OR.(AXIS==2)) THEN
         CALL MPI_ABORT(MPI_COMM_WORLD,0,IERR)
      ELSE IF (AXIS==3) THEN
         CALL THIS%BCS_U_Z%INITIALIZE(AXIS,BCS_MINUS(1),BCS_MINUS(2),BCS_PLUS(1),BCS_PLUS(2),PH%XST,PH%XEN,THIS%N) ! Z
         CALL THIS%BCS_V_Z%INITIALIZE(AXIS,BCS_MINUS(1),BCS_MINUS(2),BCS_PLUS(1),BCS_PLUS(2),PH%XST,PH%XEN,THIS%N) ! Z
      END IF

    END SUBROUTINE SOLVER_DIAG_CART_FFT_CPLX_SET_BCS

  
  SUBROUTINE INIT_SOLVER_DIAG_CART_FFT_CPLX(THIS,OPX,OPY,OPZ)
    IMPLICIT NONE
    CLASS(T_SOLVER_DIAG_CART_FFT_CPLX) :: THIS
    TYPE(T_OPERATOR_FOURIER) ::  OPX
    TYPE(T_OPERATOR_FOURIER) ::  OPY
    TYPE(T_OPERATOR_TCHEBY)  ::  OPZ
    INTEGER N(3)
    TYPE(DECOMP_INFO) :: PH
    REAL(KIND=8), DIMENSION(3) :: NU
    
    REAL(KIND=8) :: ALFA_Z_L, BETA_Z_L
    REAL(KIND=8) :: ALFA_Z_R, BETA_Z_R
    
    INTEGER :: I,J,K
    
    CALL GET_DECOMP_INFO(PH)
    
    ALFA_Z_L = THIS%BCS_U_Z%ALPHA_L
    ALFA_Z_R = THIS%BCS_U_Z%ALPHA_R
    
    BETA_Z_L = THIS%BCS_U_Z%BETA_L
    BETA_Z_R = THIS%BCS_U_Z%BETA_R
    
    NU = THIS%NU
    
    CALL THIS%LAP_X%INITIALIZE(OPX,NU(1),THIS%N(1)) 
    CALL THIS%LAP_Y%INITIALIZE(OPY,NU(2),THIS%N(2))
    CALL THIS%LAP_Z%INITIALIZE(OPZ,NU(3),ALFA_Z_L,BETA_Z_L,ALFA_Z_R,BETA_z_R,THIS%N(3))

    
    CALL ALLOC_Z( THIS%LBD , OPT_GLOBAL=.TRUE. )
    N = THIS%N
!!$    DO K = PH%ZST(3),PH%ZEN(3)
!!$       DO J = PH%ZST(2),PH%ZEN(2) 
!!$          DO I=PH%ZST(1),PH%ZEN(1)
!!$             THIS%LBD(I,J,K) =  ( this%sigma + (THIS%LAP_Z%LBD(I) + THIS%LAP_Y%LBD(J) +  THIS%LAP_Z%LBD(K) ) )**(-1)
!!$             IF (K==     1) THIS%LBD(I,J,K)=1
!!$             IF (K==N(3)+1) THIS%LBD(I,J,K)=1
!!$          END DO
!!$       END DO
!!$    END DO
    

    
    CALL ALLOC_X( THIS%CNT_U , OPT_GLOBAL=.TRUE. )
    CALL ALLOC_X( THIS%CNT_V , OPT_GLOBAL=.TRUE. )


    if (.not.allocated(  FI )) CALL ALLOC_X(  FI , OPT_GLOBAL=.TRUE. )
    if (.not.allocated( SFI )) CALL ALLOC_X( SFI , OPT_GLOBAL=.TRUE. )

    
    if (.not.allocated(DG1X)) CALL ALLOC_X( DG1X , OPT_GLOBAL=.TRUE. )
    if (.not.allocated(DG2X)) CALL ALLOC_X( DG2X , OPT_GLOBAL=.TRUE. )

    if (.not.allocated(DG1Y)) CALL ALLOC_Y( DG1Y , OPT_GLOBAL=.TRUE. )
    if (.not.allocated(DG2Y)) CALL ALLOC_Y( DG2Y , OPT_GLOBAL=.TRUE. )

    if (.not.allocated(DG1Z)) CALL ALLOC_Z( DG1Z , OPT_GLOBAL=.TRUE. )
    if (.not.allocated(DG2Z)) CALL ALLOC_Z( DG2Z , OPT_GLOBAL=.TRUE. )

    if (.not.allocated(DG_CL_Z)) allocate( DG_CL_Z(PH%ZST(1):PH%ZEN(1),PH%ZST(2):PH%ZEN(2),1:2) )

    ALLOCATE(THIS%PX  ,SOURCE=DCMPLX(THIS%LAP_X%P  ))
    ALLOCATE(THIS%PXM1,SOURCE=DCMPLX(THIS%LAP_X%PM1))
    ALLOCATE(THIS%PY  ,SOURCE=DCMPLX(THIS%LAP_Y%P  ))
    ALLOCATE(THIS%PYM1,SOURCE=DCMPLX(THIS%LAP_Y%PM1))
    ALLOCATE(THIS%PZ  ,SOURCE=DCMPLX(THIS%LAP_Z%P  ))
    ALLOCATE(THIS%PZM1,SOURCE=DCMPLX(THIS%LAP_Z%PM1))

    ALLOCATE(THIS%CC_Z  ,SOURCE=DCMPLX(THIS%LAP_Z%CC  ))
    
  END SUBROUTINE INIT_SOLVER_DIAG_CART_FFT_CPLX

  SUBROUTINE SOLVER_DIAG_CART_FFT_CPLX_SOLVE(THIS,U,SU,V,SV)
    use decomp_2d
    implicit none
    CLASS(T_SOLVER_DIAG_CART_fft_cplx) :: THIS
    real(dp),DIMENSION(:,:,:),ALLOCATABLE :: U,SU
    real(dp),DIMENSION(:,:,:),ALLOCATABLE :: V,SV
    INTEGER :: I,J,K
    TYPE(DECOMP_INFO) :: PH
    real(kind=8) :: OMEGA,SIGMA

    CALL GET_DECOMP_INFO(PH)
    
    FORALL( I=PH%XST(1):PH%XEN(1) , J=PH%XST(2):PH%XEN(2) , K=PH%XST(3):PH%XEN(3) )
       SFI(I,J,K) = CMPLX( SU(I,J,K) , SV(I,J,K) )
    END FORALL

    SIGMA = THIS%SIGMA
    OMEGA = THIS%OMEGA
    
    call SOLVER_DIAG_CART_FFT_CPLX_ADD_CONTRIB( &
         THIS%CNT_U, ph%xst, ph%xen, this%n(1), this%n(2), this%n(3), this%lap_z%fl,  this%lap_z%fr, this%bcs_u_z%bv_l, this%bcs_u_z%bv_r )
    
    call SOLVER_DIAG_CART_FFT_CPLX_ADD_CONTRIB( &
         THIS%CNT_V, ph%xst, ph%xen, this%n(1), this%n(2), this%n(3), this%lap_z%fl,  this%lap_z%fr, this%bcs_v_z%bv_l, this%bcs_v_z%bv_r )
    
    CALL KERNEL_SOLVE_FFT_CPLX(&
         THIS%PX, THIS%PXM1, THIS%LAP_X%LBD, PH%XST, PH%XEN, THIS%N(1), &
         THIS%PY, THIS%PYM1, THIS%LAP_Y%LBD, PH%YST, PH%YEN, THIS%N(2), &
         THIS%PZ, THIS%PZM1, THIS%LAP_Z%LBD, THIS%LAP_Z%CL, &
         THIS%BCS_U_Z%BV_L,  THIS%BCS_U_Z%BV_R,  THIS%BCS_V_Z%BV_L,  THIS%BCS_V_Z%BV_R, PH%ZST, PH%ZEN, THIS%N(3), &
         THIS%SIGMA, THIS%OMEGA, THIS%LBD, THIS%CNT_U, THIS%CNT_V, THIS%CC_Z, FI, SFI, &
         DG1X, DG2X, DG1Y, DG2Y, DG1Z, DG2Z, DG_CL_Z)
    
    FORALL( I=PH%XST(1):PH%XEN(1) , J=PH%XST(2):PH%XEN(2) , K=PH%XST(3):PH%XEN(3) )
       U(I,J,K) = REAL( FI(I,J,K) )
       V(I,J,K) = AIMAG( FI(I,J,K) )
    END FORALL

  END SUBROUTINE SOLVER_DIAG_CART_FFT_CPLX_SOLVE
  
  
  subroutine SOLVER_DIAG_CART_FFT_CPLX_ADD_CONTRIB( &
       fxyz,xst,xen,nx,ny,nz,fzl,fzr,bvzl,bvzr)
    implicit none
    integer :: xst(3),xen(3),nx,ny,nz
    real(dp) :: fxyz(xst(1):xen(1),xst(2):xen(2),xst(3):xen(3))
    real(dp) :: fzl(1:nz+1), fzr(1:nz+1), bvzl(xst(1):xen(1),xst(2):xen(2)), bvzr(xst(1):xen(1),xst(2):xen(2))
    integer :: i,j,k
    
    forall(i=xst(1):xen(1),j=xst(2):xen(2),k=xst(3):xen(3))
       fxyz(i,j,k) = + fzl(k)*bvzl(i,j) + fzr(k)*bvzr(i,j)
    end forall
    
  end subroutine SOLVER_DIAG_CART_FFT_CPLX_ADD_CONTRIB



    SUBROUTINE SOLVER_DIAG_CART_FFT_CPLX_SET_BVS_VALUES_TIME_REAL(THIS,TC,MSH,AXIS,UDF_PLUS,UDF_MINUS)
      use m_boundary_conditions
      use m_lap1d_tcheby
      IMPLICIT NONE
      CLASS(T_SOLVER_DIAG_CART_fft_cplx) :: THIS
      TYPE(T_MESH_BASE)                  :: MSH(3)
      INTEGER                            :: AXIS
      PROCEDURE(UDF_TIMESPACE)           :: UDF_plus
      PROCEDURE(UDF_TIMESPACE)           :: UDF_minus
      
      TYPE(DECOMP_INFO) :: PH
      REAL(KIND=8) :: TIME,TC
      integer :: ierr
      
      CALL GET_DECOMP_INFO(PH)
      TIME = TC

      if (AXIS==1) THEN
         CALL MPI_ABORT(MPI_COMM_WORLD,0,IERR)
      else if (axis==2) then
         CALL MPI_ABORT(MPI_COMM_WORLD,0,IERR)
      else if (axis==3) then
         CALL THIS%BCS_U_Z%SET_BCS(AXIS,-1,TIME,MSH,UDF_MINUS,PH%XST,PH%XEN,THIS%N)
         CALL THIS%BCS_U_Z%SET_BCS(AXIS,+1,TIME,MSH,UDF_PLUS,PH%XST,PH%XEN,THIS%N)
      END if
      
    end subroutine SOLVER_DIAG_CART_FFT_CPLX_SET_BVS_VALUES_TIME_REAL


    SUBROUTINE SOLVER_DIAG_CART_FFT_CPLX_SET_BVS_VALUES_TIME_IMAG(THIS,TC,MSH,AXIS,UDF_PLUS,UDF_MINUS)
      use m_boundary_conditions
      use m_lap1d_tcheby
      IMPLICIT NONE
      CLASS(T_SOLVER_DIAG_CART_fft_cplx) :: THIS
      TYPE(T_MESH_BASE)                  :: MSH(3)
      INTEGER                            :: AXIS
      PROCEDURE(UDF_TIMESPACE)           :: UDF_plus
      PROCEDURE(UDF_TIMESPACE)           :: UDF_minus
      
      TYPE(DECOMP_INFO) :: PH
      REAL(KIND=8) :: TIME,TC
      integer :: ierr
      
      CALL GET_DECOMP_INFO(PH)
      TIME = TC

      if (AXIS==1) THEN
         CALL MPI_ABORT(MPI_COMM_WORLD,0,IERR)
      else if (axis==2) then
         CALL MPI_ABORT(MPI_COMM_WORLD,0,IERR)
      else if (axis==3) then
         CALL THIS%BCS_V_Z%SET_BCS(AXIS,-1,TIME,MSH,UDF_MINUS,PH%XST,PH%XEN,THIS%N)
         CALL THIS%BCS_V_Z%SET_BCS(AXIS,+1,TIME,MSH,UDF_PLUS,PH%XST,PH%XEN,THIS%N)
      END if
      
    end subroutine SOLVER_DIAG_CART_FFT_CPLX_SET_BVS_VALUES_TIME_IMAG
    

    SUBROUTINE SOLVER_DIAG_CART_FFT_CPLX_SET_BVS_VALUES_TIME(THIS,MSH,AXIS,UDF_PLUS,UDF_MINUS)
      use m_boundary_conditions
      use m_lap1d_tcheby
      IMPLICIT NONE
      CLASS(T_SOLVER_DIAG_CART_fft_cplx) :: THIS
      TYPE(T_MESH_BASE)             :: MSH(3)
      INTEGER                       :: AXIS
      PROCEDURE(UDF_TIMESPACE)      :: UDF_plus
      PROCEDURE(UDF_TIMESPACE)      :: UDF_minus
      
      TYPE(DECOMP_INFO) :: PH
      REAL(KIND=8) :: TIME
      integer :: ierr
      
      CALL GET_DECOMP_INFO(PH)
      TIME = 0

      if (AXIS==1) THEN
         CALL MPI_ABORT(MPI_COMM_WORLD,0,IERR)
      else if (axis==2) then
         CALL MPI_ABORT(MPI_COMM_WORLD,0,IERR)
      else if (axis==3) then
         CALL MPI_ABORT(MPI_COMM_WORLD,0,IERR)
      END if
      
    end subroutine SOLVER_DIAG_CART_FFT_CPLX_SET_BVS_VALUES_TIME
    
    SUBROUTINE SOLVER_DIAG_CART_FFT_CPLX_UPDATE_BCS_CONTRIB(this)
      IMPLICIT NONE
      CLASS(T_SOLVER_DIAG_CART_fft_cplx) :: THIS
      TYPE(DECOMP_INFO) :: PH
      CALL GET_DECOMP_INFO(PH)
!!$      call SOLVER_DIAG_CART_FFT_CPLX_ADD_CONTRIB( &
!!$           this%lap_x%fl,  this%lap_x%fr, this%bcs_x%bv_l,  this%bcs_x%bv_r, ph%xst, ph%xen, this%n(1), &
!!$           this%lap_y%fl,  this%lap_y%fr, this%bcs_y%bv_l,  this%bcs_y%bv_r, ph%xst, ph%xen, this%n(2), &
!!$           this%lap_z%fl,  this%lap_z%fr, this%bcs_z%bv_l,  this%bcs_z%bv_r, ph%xst, ph%xen, this%n(3), &
!!$           THIS%cnt )
    END SUBROUTINE SOLVER_DIAG_CART_FFT_CPLX_UPDATE_BCS_CONTRIB
    
    !>
    SUBROUTINE KERNEL_SOLVE_FFT_CPLX(&
         PX, PXM1, LBD_X, XST, XEN, NX, &
         PY, PYM1, LBD_Y, YST, YEN, NY, &
         PZ, PZM1, LBD_Z, CL_Z, BC_VAL_L_Z_REAL, BC_VAL_R_Z_REAL, BC_VAL_L_Z_IMAG, BC_VAL_R_Z_IMAG, ZST, ZEN, NZ, &
         SIGMA, OMEGA, LBD_XYZ, FU_XYZ, FV_XYZ, CC_Z, FI, SFI, DG1X, DG2X, DG1Y, DG2Y, DG1Z, DG2Z ,DG_CL_Z)
      use mpi
      use decomp_2d
      use decomp_2d_mpi
      implicit none
      INTEGER :: XST(3) , XEN(3) , NX
      INTEGER :: YST(3) , YEN(3) , NY
      INTEGER :: ZST(3) , ZEN(3) , NZ
      
      COMPLEX(KIND=8) ::  PX(1:NX+1,1:NX+1) , PXM1(1:NX+1,1:NX+1)
      real(dp) :: LBD_X(1:NX+1)

      COMPLEX(KIND=8) :: PY(1:NY+1,1:NY+1) , PYM1(1:NY+1,1:NY+1)
      real(dp) :: LBD_Y(1:NY+1)

      COMPLEX(KIND=8) :: PZ(1:NZ+1,1:NZ+1) , PZM1(1:NZ+1,1:NZ+1)
      real(dp) :: LBD_Z(1:NZ+1)
      
      real(dp) :: CL_X(2,2)
      real(dp) :: CL_Y(2,2)
      real(dp) :: CL_Z(2,2)
      
      real(dp) :: BC_VAL_L_X_REAL(XST(2):XEN(2),XST(3):XEN(3)),BC_VAL_R_X_REAL(XST(2):XEN(2),XST(3):XEN(3))
      real(dp) :: BC_VAL_L_Y_REAL(XST(1):XEN(1),XST(3):XEN(3)),BC_VAL_R_Y_REAL(XST(1):XEN(1),XST(3):XEN(3))
      real(dp) :: BC_VAL_L_Z_REAL(XST(1):XEN(1),XST(2):XEN(2)),BC_VAL_R_Z_REAL(XST(1):XEN(1),XST(2):XEN(2))

      real(dp) :: BC_VAL_L_X_IMAG(XST(2):XEN(2),XST(3):XEN(3)),BC_VAL_R_X_IMAG(XST(2):XEN(2),XST(3):XEN(3))
      real(dp) :: BC_VAL_L_Y_IMAG(XST(1):XEN(1),XST(3):XEN(3)),BC_VAL_R_Y_IMAG(XST(1):XEN(1),XST(3):XEN(3))
      real(dp) :: BC_VAL_L_Z_IMAG(XST(1):XEN(1),XST(2):XEN(2)),BC_VAL_R_Z_IMAG(XST(1):XEN(1),XST(2):XEN(2))


      
      complex(dp) :: fi(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
      complex(dp) ::sfi(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
      
      real(dp) :: sigma,OMEGA
      complex(dp) :: lbd_xyz(ZST(1):ZEN(1),ZST(2):ZEN(2),ZST(3):ZEN(3))
      real(dp) ::   fu_xyz(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
      real(dp) ::   fv_xyz(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
      
      complex(kind=8) :: DG1X(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
      complex(kind=8) :: DG2X(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
      
      complex(kind=8) :: DG1Y(YST(1):YEN(1),YST(2):YEN(2),YST(3):YEN(3))
      complex(kind=8) :: DG2Y(YST(1):YEN(1),YST(2):YEN(2),YST(3):YEN(3))

      complex(kind=8) :: DG1Z(ZST(1):ZEN(1),ZST(2):ZEN(2),ZST(3):ZEN(3))
      complex(kind=8) :: DG2Z(ZST(1):ZEN(1),ZST(2):ZEN(2),ZST(3):ZEN(3))
      
      complex(dp) :: CC_Z(1:2,1:NZ+1)
      
      complex(dp) :: DG_CL_Z(ZST(1):ZEN(1),ZST(2):ZEN(2),1:2)

      integer :: i,j,k,is(3),ie(3)
      complex(dp) :: aa


      
      DO K = ZST(3),ZEN(3)
         DO J = ZST(2),ZEN(2) 
            DO I=ZST(1),ZEN(1)

               AA = cmplx(0d0,1d0)*OMEGA + SIGMA + LBD_X(I) + LBD_Y(J) +  LBD_Z(K) 

               LBD_XYZ(I,J,K) = AA**(-1) 
               
               IF (K==   1) LBD_XYZ(I,J,K) = CMPLX(1.0,0.0)
               IF (K==NZ+1) LBD_XYZ(I,J,K) = CMPLX(1.0,0.0)
               
            END DO
         END DO
      END DO
      
      DG1X = SFI - CMPLX(FU_XYZ,FV_XYZ )
      
      !> 
      IS = GET_IS_X([0,0,0])
      IE = GET_IE_X([0,0,0])
      
      IF (XST(3)==1) THEN
         FORALL(I=IS(1):IE(1),J=IS(2):IE(2))
            DG1X(I,J,   1) = CMPLX(  BC_VAL_L_Z_REAL(I,J) ,  BC_VAL_L_Z_IMAG(I,J) )
         END FORALL
      END IF
      
      IF (XEN(3)==NZ+1) THEN
         FORALL(I=IS(1):IE(1),J=IS(2):IE(2))
            DG1X(I,J,NZ+1) = CMPLX(  BC_VAL_R_Z_REAL(I,J) ,  BC_VAL_R_Z_IMAG(I,J) )
         END FORALL
      END IF
      
      
    CALL TENSOR_PRODUCT_I( PXM1, DG1X, DG2X, XST, XEN, NX )
    CALL TRANSPOSE_X_TO_Y( DG2X, DG1Y )
    
    CALL TENSOR_PRODUCT_J( PYM1, DG1Y, DG2Y, YST, YEN, NY )
    
    CALL TRANSPOSE_Y_TO_Z( DG2Y, DG1Z )
    CALL TENSOR_PRODUCT_K( PZM1, DG1Z, DG2Z, ZST, ZEN, NZ )
    
    
    
    FORALL(I=ZST(1):ZEN(1),J=ZST(2):ZEN(2),K=ZST(3):ZEN(3))
       DG1Z(I,J,K) = LBD_XYZ(I,J,K)*DG2Z(I,J,K)
    END FORALL
    
    
    CALL TENSOR_PRODUCT_K( PZ   , DG1Z, DG2Z, ZST, ZEN, NZ )
    call TENSOR_PRODUCT_CL_K(CC_Z, DG2Z, DG_CL_Z, ZST, ZEN , NZ )

    IS = GET_IS_Z([0,0,0])
    IE = GET_IE_Z([0,0,0])
    !> Z-PENCIL
    FORALL(I=IS(1):IE(1),J=IS(2):IE(2))
       DG2Z(I,J,   1) = DG_CL_Z(I,J,1)*CL_Z(1,1) + DG_CL_Z(I,J,2)*CL_Z(1,2)
       DG2Z(I,J,NZ+1) = DG_CL_Z(I,J,1)*CL_Z(2,1) + DG_CL_Z(I,J,2)*CL_Z(2,2)
    END FORALL
    
    CALL TRANSPOSE_Z_TO_Y( DG2Z , DG1Y )
    !> Y-PENCIL
    CALL TENSOR_PRODUCT_J( PY  , DG1Y, DG2Y, YST, YEN, NY )
    CALL TRANSPOSE_Y_TO_X( DG2Y, DG1X )
    !> X-PENCIL
    CALL TENSOR_PRODUCT_I( PX  , DG1X , FI, XST, XEN, NX )
    
    
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
    
    SUBROUTINE TENSOR_PRODUCT_CL_K(CC,FI,DFI,ZST,ZEN,NZ)
      IMPLICIT NONE
      INTEGER :: ZST(3),ZEN(3),NZ
      COMPLEX(dp) :: CC(1:2,1:NZ+1)
      COMPLEX(dp) ::  FI(ZST(1):ZEN(1),ZST(2):ZEN(2),ZST(3):ZEN(3))
      COMPLEX(dp) :: DFI(ZST(1):ZEN(1),ZST(2):ZEN(2),1:2)
        
      INTEGER :: N(3),NI,NJ,NK
      COMPLEX(dp) ::  ALPHA,BETA
        
      ALPHA = 1._dp
      BETA = 0._dp
      N = ZEN - ZST + 1
      NI = N(1)
      NJ = N(2)
      NK = N(3) ! NZ + 1 car Z-pencil
      CALL ZGEMM('N','T',NI*NJ,2,NK, ALPHA, FI, NI*NJ , CC, 2 , BETA , DFI, NI*NJ )
      
    END SUBROUTINE TENSOR_PRODUCT_CL_K
    
    
    subroutine tensor_product_i(DX,FI,DFI,XST,XEN,NX)
        implicit none
        INTEGER :: XST(3),XEN(3),NX
        COMPLEX(dp) ::  DX(1:NX+1,1:NX+1)
        COMPLEX(dp) ::  FI(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
        COMPLEX(dp) :: DFI(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
        
        COMPLEX(dp) ::  ALPHA,BETA
        INTEGER :: N(3),NI,NJ,NK
        
        N = XEN-XST + 1
        NI = N(1)
        NJ = N(2)
        NK = N(3)
        ALPHA = 1._dp
        BETA = 0.0_dp

        CALL ZGEMM( 'N', 'N', NI, NJ*NK, NI, ALPHA, DX, NI, FI, NI, BETA, DFI, NI )
        
      end subroutine tensor_product_i
      
      SUBROUTINE TENSOR_PRODUCT_J(DY,FI,DFI,YST,YEN,NY)
        IMPLICIT NONE
        INTEGER :: YST(3),YEN(3),NY
        COMPLEX(dp) :: DY(1:NY+1,1:NY+1)
        COMPLEX(dp) ::  FI(YST(1):YEN(1),YST(2):YEN(2),YST(3):YEN(3))
        COMPLEX(dp) :: DFI(YST(1):YEN(1),YST(2):YEN(2),YST(3):YEN(3))
        COMPLEX(dp) ::  ALPHA,BETA
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
           CALL ZGEMM('N','T',NI,NJ,NJ, ALPHA ,FI(is,js,k), NI , DY, NJ , beta , DFI(is,js,k), NI )
        end do
        
!        DO K=YST(3),YEN(3)
!           CALL DGEMM('N','T',NI,NJ,NJ,ALPHA,FI(YST(1),YST(2),K),NI,DY,NJ,BETA,DFI(YST(1),YST(2),K),NI)
!        END DO
        
      END SUBROUTINE TENSOR_PRODUCT_J


      
      
      SUBROUTINE TENSOR_PRODUCT_K(DZ,FI,DFI,ZST,ZEN,NZ)
        IMPLICIT NONE
        INTEGER :: ZST(3),ZEN(3),NZ
        COMPLEX(dp) :: DZ(1:NZ+1,1:NZ+1)
        COMPLEX(dp) ::  FI(ZST(1):ZEN(1),ZST(2):ZEN(2),ZST(3):ZEN(3))
        COMPLEX(dp) :: DFI(ZST(1):ZEN(1),ZST(2):ZEN(2),ZST(3):ZEN(3))

        INTEGER :: N(3),NI,NJ,NK
        COMPLEX(dp) ::  ALPHA,BETA
        
        ALPHA = 1._dp
        BETA = 0._dp
        N = ZEN - ZST + 1
        NI = N(1)
        NJ = N(2)
        NK = N(3)
        CALL ZGEMM('N','T',NI*NJ,NK,NK, ALPHA, FI, NI*NJ , DZ, NK , BETA , DFI, NI*NJ )
        
      END SUBROUTINE TENSOR_PRODUCT_K
      
      
    end subroutine KERNEL_SOLVE_FFT_CPLX
    
    
    
  end module m_solver_diag_cart_fft_cplx
