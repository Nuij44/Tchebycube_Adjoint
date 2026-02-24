module m_solver_diag_cart_fft_dft
  use FFTW3
  use m_solver_diag_base
  use m_lap1d_tcheby
  use m_lap1d_fourier_dft
  use m_boundary_conditions

  use decomp_2d
  use decomp_2d_fft
  use decomp_2d_mpi
  use mpi
  implicit none
  
  type , extends(t_solver_diag_base) :: t_solver_diag_cart_fft_DFT
     type(t_lap1d_fourier_DFT) :: lap_x
     type(t_lap1d_fourier_DFT) :: lap_y
     type(t_lap1d_tcheby)  :: lap_z

     integer*8 :: X_FW,X_BW,Y_FW,Y_BW
     
     real(dp),dimension(:,:,:),allocatable :: lbd
     real(dp),dimension(:,:,:),allocatable :: cnt
     
   contains 
     procedure :: INITIALISE => INIT_SOLVER_DIAG_CART_FFT
     procedure :: SET_BCS => SOLVER_DIAG_CART_FFT_SET_BCS
     procedure :: SET_BVS => SOLVER_DIAG_CART_FFT_SET_BVS_VALUES_TIME
     procedure :: UPDATE_BCS => SOLVER_DIAG_CART_FFT_UPDATE_BCS_CONTRIB
     procedure :: SOLVE => SOLVER_DIAG_CART_FFT_SOLVE
     procedure :: SET_PARAMS => SOLVER_DIAG_CART_FFT_SET_PARAMS
  end type t_solver_diag_cart_fft_DFT

  COMPLEX(DP),DIMENSION(:,:,:),ALLOCATABLE :: DG1X, DG2X, DG1Y, DG2Y, DG1Z, DG2Z 
  COMPLEX(DP),DIMENSION(:,:,:),ALLOCATABLE    :: DG_CL_X,DG_CL_Y,DG_CL_Z
  PRIVATE :: DG1X, DG2X, DG1Y, DG2Y, DG1Z, DG2Z 

contains 
  SUBROUTINE SOLVER_DIAG_CART_FFT_SET_PARAMS(THIS,DIM,NU,sigma)
    IMPLICIT NONE
    CLASS(T_SOLVER_DIAG_CART_FFT_DFT) :: THIS
    integer                       :: dim(3)
    real(dp)                      :: nu
    real(dp)                      :: sigma
    
    THIS%N = DIM
    THIS%NU = nu
    THIS%sigma = sigma
    
  end SUBROUTINE SOLVER_DIAG_CART_FFT_SET_PARAMS
  

    
  SUBROUTINE SOLVER_DIAG_CART_FFT_SET_BCS(THIS,AXIS,BCS_MINUS,BCS_PLUS)
    use m_boundary_conditions
    IMPLICIT NONE
    CLASS(T_SOLVER_DIAG_CART_FFT_DFT) :: THIS
    INTEGER                       :: AXIS
    real(dp)                      :: BCS_MINUS(2)
    real(dp)                      :: BCS_PLUS(2)
    integer :: ierr
    TYPE(DECOMP_INFO) :: PH
    
    CALL GET_DECOMP_INFO(PH)
    if ((AXIS==1).or.(axis==2)) THEN
       call MPI_Abort(mpi_comm_world,0,ierr)
    END if

    CALL THIS%BCS_Z%INITIALIZE(AXIS,BCS_MINUS(1),BCS_MINUS(2),BCS_PLUS(1),BCS_PLUS(2),PH%XST,PH%XEN,THIS%N) ! X
    
  end subroutine SOLVER_DIAG_CART_FFT_SET_BCS
  
  
  SUBROUTINE INIT_SOLVER_DIAG_CART_FFT(THIS,OPX,OPY,OPZ)
    use m_boundary_conditions
    IMPLICIT NONE
    CLASS(T_SOLVER_DIAG_CART_FFT_DFT) :: THIS
    TYPE(T_OPERATOR_FOURIER_DFT) ::  OPX
    TYPE(T_OPERATOR_FOURIER_DFT) ::  OPY
    TYPE(T_OPERATOR_TCHEBY)  ::  OPZ
    INTEGER N(3)
    TYPE(DECOMP_INFO) :: PH
    real(dp), dimension(3) :: nu
    
    real(dp) :: alfa_z_l, beta_z_l
    real(dp) :: alfa_z_r, beta_z_r
    integer :: i,j,k
    complex(dp), allocatable:: inx(:,:,:),iny(:,:,:)
    complex(dp), allocatable:: outx(:,:,:),outy(:,:,:)
    integer :: X(3),Y(3)
    
    CALL GET_DECOMP_INFO(PH)

    allocate( inX(PH%XST(1):PH%XEN(1),PH%XST(2):PH%XEN(2),PH%XST(3):PH%XEN(3)) )
    allocate( inY(PH%YST(1):PH%YEN(1),PH%YST(2):PH%YEN(2),PH%YST(3):PH%YEN(3)) )

    allocate( outX(PH%XST(1):PH%XEN(1),PH%XST(2):PH%XEN(2),PH%XST(3):PH%XEN(3)) )
    allocate( outY(PH%YST(1):PH%YEN(1),PH%YST(2):PH%YEN(2),PH%YST(3):PH%YEN(3)) )

    
    ALFA_Z_L = THIS%BCS_Z%ALPHA_L
    ALFA_Z_R = THIS%BCS_Z%ALPHA_R
    
    BETA_Z_L = THIS%BCS_Z%BETA_L
    BETA_Z_R = THIS%BCS_Z%BETA_R

    nu = this%nu
    N = THIS%N
    
    ! INITIALISATION DE CHAQUE DIRECTION
    X = PH%XEN-PH%XST + 1
    Y  = PH%YEN-PH%YST + 1


    call dfftw_plan_many_dft(THIS%X_FW, 1, n(1)+1, X(2)*X(3), inX , n(1)+1    , 1, n(1)+1    , outX , n(1)+1, 1, n(1)+1 , FFTW_FORWARD, FFTW_MEASURE )
    call dfftw_plan_many_dft(THIS%X_BW, 1, n(1)+1, X(2)*X(3), outX, n(1)+1    , 1, n(1)+1    , inX  , n(1)+1, 1, n(1)+1 , FFTW_BACKWARD, FFTW_MEASURE)

    call dfftw_plan_many_dft(THIS%Y_FW, 1, n(2)+1, Y(1), inY(:,:,1) , n(2)+1    , Y(1), 1    , outY(:,:,1) , n(2)+1, Y(1), 1 , FFTW_FORWARD, FFTW_MEASURE )
    call dfftw_plan_many_dft(THIS%Y_BW, 1, n(2)+1, Y(1), outY(:,:,1), n(2)+1    , Y(1), 1    , inY(:,:,1)  , n(2)+1, Y(1), 1 , FFTW_BACKWARD, FFTW_MEASURE )



    
    CALL THIS%LAP_X%INITIALIZE(OPX,NU(1),THIS%N(1)) 

    CALL THIS%LAP_Y%INITIALIZE(OPY,NU(2),THIS%N(2))

    CALL THIS%LAP_Z%INITIALIZE(OPZ,NU(3),ALFA_Z_L,BETA_Z_L,ALFA_Z_R,BETA_Z_R,THIS%N(3))


     
    CALL ALLOC_Z( THIS%LBD , OPT_GLOBAL=.TRUE. )
    
    DO K = PH%ZST(3),PH%ZEN(3)
       DO J = PH%ZST(2),PH%ZEN(2) 
          DO I=PH%ZST(1),PH%ZEN(1)
             THIS%LBD(I,J,K) =  ( this%sigma + (  THIS%LAP_Z%LBD(K) ) )**(-1)
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
    
  END SUBROUTINE INIT_SOLVER_DIAG_CART_FFT

  SUBROUTINE SOLVER_DIAG_CART_FFT_SOLVE(THIS,FI,SFI)
    use decomp_2d
    implicit none
    CLASS(T_SOLVER_DIAG_CART_FFT_DFT) :: THIS
    real(dp),DIMENSION(:,:,:),ALLOCATABLE :: FI
    real(dp),DIMENSION(:,:,:),ALLOCATABLE :: SFI

    TYPE(DECOMP_INFO) :: PH

    CALL GET_DECOMP_INFO(PH)

      call SOLVER_DIAG_CART_FFT_ADD_CONTRIB( &
           THIS%cnt, ph%xst, ph%xen, this%n(1),  this%n(2), this%n(3), &
           this%lap_z%fl,  this%lap_z%fr, this%bcs_z%bv_l,  this%bcs_z%bv_r )    
    
 
    CALL KERNEL_SOLVE_FFT(&
         THIS%X_FW, THIS%X_BW, THIS%Y_FW, THIS%Y_BW, &
         THIS%lap_X%wave,PH%XST, PH%XEN, THIS%N(1),  &
         THIS%lap_Y%wave,PH%YST, PH%YEN, THIS%N(2),  &
         THIS%LAP_Z%P, THIS%LAP_Z%PM1, THIS%LAP_Z%LBD, PH%ZST, PH%ZEN, THIS%N(3), &
         THIS%LAP_Z%CL, THIS%BCS_Z%BV_L,  THIS%BCS_Z%BV_R, &
         THIS%SIGMA, THIS%NU, THIS%LBD, THIS%CNT, THIS%LAP_Z%CC, FI, SFI, DG1X, DG2X, DG1Y, DG2Y, DG1Z, DG2Z, DG_CL_X, DG_CL_Y,DG_CL_Z)
    
  END SUBROUTINE SOLVER_DIAG_CART_FFT_SOLVE
  
  
  subroutine SOLVER_DIAG_CART_FFT_ADD_CONTRIB( &
       fxyz,xst,xen,nx,ny,nz, &
       fzl,fzr,bvzl,bvzr )
    implicit none

    INTEGER :: XST(3),XEN(3),NX,NY,NZ
    real(dp) :: FXYZ(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
    real(dp) :: FZL(1:NZ+1), FZR(1:NZ+1), BVZL(XST(1):XEN(1),XST(2):XEN(2)), BVZR(XST(1):XEN(1),XST(2):XEN(2))
    INTEGER :: I,J,K
    
    forall(i=xst(1):xen(1),j=xst(2):xen(2),k=xst(3):xen(3))
       fxyz(i,j,k) = + fzl(k)*bvzl(i,j) + fzr(k)*bvzr(i,j)
    end forall

  end subroutine SOLVER_DIAG_CART_FFT_ADD_CONTRIB
  
  
  SUBROUTINE SOLVER_DIAG_CART_FFT_SET_BVS_VALUES_TIME(THIS,MSH,AXIS,UDF_PLUS,UDF_MINUS)
      use m_boundary_conditions
      use m_lap1d_tcheby
      IMPLICIT NONE
      CLASS(T_SOLVER_DIAG_CART_FFT_DFT) :: THIS
      TYPE(T_MESH_BASE)             :: MSH(3)
      INTEGER                       :: AXIS
      PROCEDURE(UDF_TIMESPACE)      :: UDF_plus
      PROCEDURE(UDF_TIMESPACE)      :: UDF_minus
      
      TYPE(DECOMP_INFO) :: PH
      real(dp) :: TIME
      integer :: ierr
      
      CALL GET_DECOMP_INFO(PH)
      TIME = 0

      if ((AXIS==1).or.(axis==2)) THEN
         call MPI_Abort(mpi_comm_world,0,ierr)
      END if
      CALL THIS%BCS_Z%SET_BCS(AXIS,-1,TIME,MSH,UDF_MINUS,PH%XST,PH%XEN,THIS%N)
      CALL THIS%BCS_Z%SET_BCS(AXIS,+1,TIME,MSH,UDF_PLUS ,PH%XST,PH%XEN,THIS%N)
      
    end subroutine SOLVER_DIAG_CART_FFT_SET_BVS_VALUES_TIME
    
    SUBROUTINE SOLVER_DIAG_CART_FFT_UPDATE_BCS_CONTRIB(this)
      IMPLICIT NONE
      CLASS(T_SOLVER_DIAG_CART_FFT_DFT) :: THIS
      TYPE(DECOMP_INFO) :: PH
      CALL GET_DECOMP_INFO(PH)

      call SOLVER_DIAG_CART_FFT_ADD_CONTRIB( &
           THIS%cnt, ph%xst, ph%xen, this%n(1),  this%n(2), this%n(3), &
           this%lap_z%fl,  this%lap_z%fr, this%bcs_z%bv_l,  this%bcs_z%bv_r )
    END SUBROUTINE SOLVER_DIAG_CART_FFT_UPDATE_BCS_CONTRIB
    
    !>
    SUBROUTINE KERNEL_SOLVE_FFT(&
         X_FW, X_BW, Y_FW, Y_BW, &
         wave_X,XST,XEN,NX, &
         wave_Y,YST,YEN,NY, &
         PZ, PZM1, LBD_Z, ZST,ZEN,NZ, CL_Z, BC_VAL_L_Z, BC_VAL_R_Z, &
         SIGMA, NU, LBD_XYZ, F_XYZ, CC_Z, FI, SFI, &
         DG1X, DG2X, DG1Y, DG2Y, DG1Z, DG2Z, DG_CL_X, DG_CL_Y, DG_CL_Z )
      use mpi
      use decomp_2d
      use decomp_2d_mpi
      implicit none

      INTEGER*8 :: X_FW,X_BW,Y_FW,Y_BW

      INTEGER :: XST(3) , XEN(3) , NX
      INTEGER :: YST(3) , YEN(3) , NY
      INTEGER :: ZST(3) , ZEN(3) , NZ
      
      REAL(dp) :: LX, LY

      REAL(dp) :: wave_X(1:NX),wave_Y(1:NY)

      real(dp) :: PZ(1:NZ+1,1:NZ+1) , PZM1(1:NZ+1,1:NZ+1) , LBD_Z(1:NZ+1)
      real(dp) :: CL_Z(2,2)
      
      real(dp) :: BC_VAL_L_Z(XST(1):XEN(1),XST(2):XEN(2)),BC_VAL_R_Z(XST(1):XEN(1),XST(2):XEN(2))
      
      real(dp) ::  FI(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
      real(dp) :: SFI(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))

      real(dp) :: sigma,nu(3)
      real(dp) :: lbd_xyz(ZST(1):ZEN(1),ZST(2):ZEN(2),ZST(3):ZEN(3))
      real(dp) ::   f_xyz(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
      
      COMPLEX(dp) :: DG1X(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
      COMPLEX(dp) :: DG2X(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
      
      COMPLEX(dp) :: DG1Y(YST(1):YEN(1),YST(2):YEN(2),YST(3):YEN(3))
      COMPLEX(dp) :: DG2Y(YST(1):YEN(1),YST(2):YEN(2),YST(3):YEN(3))

      COMPLEX(dp) :: DG1Z(ZST(1):ZEN(1),ZST(2):ZEN(2),ZST(3):ZEN(3))
      COMPLEX(dp) :: DG2Z(ZST(1):ZEN(1),ZST(2):ZEN(2),ZST(3):ZEN(3))

      real(dp) :: CC_Z(1:2,1:NZ+1)
      
      COMPLEX(DP) :: DG_CL_X(1:2,XST(2):XEN(2),XST(3):XEN(3))
      COMPLEX(DP) :: DG_CL_Y(YST(1):YEN(1),1:2,YST(3):YEN(3))
      COMPLEX(DP) :: DG_CL_Z(ZST(1):ZEN(1),ZST(2):ZEN(2),1:2)

      COMPLEX(DP) :: stock

      integer :: i,j,k,is(3),ie(3)


      DO K = ZST(3),ZEN(3)
         DO J = ZST(2),ZEN(2) 
            DO I=ZST(1),ZEN(1)
               lbd_xyz(I,J,K) = (sigma - nu(1)*wave_X(I)**2 - nu(2)*wave_Y(J)**2 +  LBD_z(K))**(-1)
               IF (K==   1) LBD_XYZ(I,J,K) = 1
               IF (K==NZ+1) LBD_XYZ(I,J,K) = 1
             END DO
         END DO
      END DO


      
      DG1X = SFI - F_XYZ
      
       IS = GET_IS_X([0,0,1])
       IE = GET_IE_X([0,0,1])
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
      
      CALL FFT_I(DG1X,DG2X,XST,XEN,X_BW)

      CALL TRANSPOSE_X_TO_Y( DG2X, DG1Y )
      

      CALL FFT_J(DG1Y,DG2Y,YST,YEN,Y_BW)
      
      CALL TRANSPOSE_Y_TO_Z( DG2Y, DG1Z )

      CALL TENSOR_dpODUCT_K( PZM1, DG1Z, DG2Z, ZST, ZEN, NZ )
      
      
      FORALL(I=ZST(1):ZEN(1),J=ZST(2):ZEN(2),K=ZST(3):ZEN(3))
         DG1Z(I,J,K) = LBD_XYZ(I,J,K)*DG2Z(I,J,K)
      END FORALL
      
      CALL TENSOR_dpODUCT_K( PZ , DG1Z, DG2Z, ZST, ZEN, NZ )
      

      call TENSOR_dpODUCT_CL_K(CC_Z, DG2Z, DG_CL_Z, ZST, ZEN , NZ )
      
      !> Z-PENCIL
      IS = GET_IS_Z([0,0,1])
      IE = GET_IE_Z([0,0,1])

      FORALL(I=IS(1):IE(1),J=IS(2):IE(2))
         DG2Z(I,J,   1) = DG_CL_Z(I,J,1)*CL_Z(1,1) + DG_CL_Z(I,J,2)*CL_Z(1,2)
         DG2Z(I,J,NZ+1) = DG_CL_Z(I,J,1)*CL_Z(2,1) + DG_CL_Z(I,J,2)*CL_Z(2,2) 
      END FORALL

      CALL TRANSPOSE_Z_TO_Y( DG2Z , DG2Y )
      !> Y-PENCIL
      CALL FFT_J(DG2Y,DG1Y,YST,YEN,Y_FW)


      CALL TRANSPOSE_Y_TO_X( DG1Y, DG2X )

      !> X-PENCIL
      CALL FFT_I (DG2X,DG1X,XST,XEN,X_FW)

      !.. Normalisation
      DG1X = DG1X/((NX+1)*(NY+1))
      FI = REAL(DG1X)
      
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
    
    SUBROUTINE TENSOR_dpODUCT_CL_K(CC,FI,DFI,ZST,ZEN,NZ)
      IMPLICIT NONE
      INTEGER :: ZST(3),ZEN(3),NZ
      real(dp) :: CC(1:2,1:NZ+1)
      COMPLEX(DP) ::  FI(ZST(1):ZEN(1),ZST(2):ZEN(2),ZST(3):ZEN(3))
      COMPLEX(DP) :: DFI(ZST(1):ZEN(1),ZST(2):ZEN(2),1:2)
        
      INTEGER :: N(3),NI,NJ,NK
      real(dp) ::  ALPHA,BETA
        
      ALPHA = 1._dp
      BETA = 0._dp
      N = ZEN - ZST + 1
      NI = N(1)
      NJ = N(2)
      NK = N(3) ! NZ + 1 car Z-pencil

      CALL ZGEMM('N','T',NI*NJ,2,NK, ALPHA, FI, NI*NJ , CMPLX(CC,0._dp), 2 , BETA , DFI, NI*NJ )

    END SUBROUTINE TENSOR_dpODUCT_CL_K

    subroutine FFT_I(FI,DFI,XST,XEN,plan)
      INTEGER :: XST(3),XEN(3)
      COMPLEX(DP) ::  FI(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
      COMPLEX(dp) :: DFI(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
      integer*8 :: plan
      integer :: is(3), ie(3)
      
!      IS = GET_IS_X([0,0,1])
!      IE = GET_IE_X([0,0,1])
       
!      call dfftw_execute_dft(plan, fi(:,:,is(3):ie(3)), dfi(:,:,is(3):ie(3)))
      call dfftw_execute_dft(plan, fi, dfi)
      
    end subroutine FFT_I

    subroutine FFT_J(fi,dfi,YST,YEN,plan)
      INTEGER :: YST(3),YEN(3)
      COMPLEX(DP) ::  FI(YST(1):YEN(1),YST(2):YEN(2),YST(3):YEN(3))
      COMPLEX(dp) :: DFI(YST(1):YEN(1),YST(2):YEN(2),YST(3):YEN(3))
      integer*8 :: plan
      integer :: k,is(3),ie(3)
      
!      IS = GET_IS_Y([0,0,1])
!      IE = GET_IE_Y([0,0,1])
      
      do k=YST(3),YEN(3)!is(3),ie(3)
         call dfftw_execute_dft(plan, FI(:,:,k), DFI(:,:,k))
      end do
      
    end subroutine FFT_J


    SUBROUTINE TENSOR_dpODUCT_K(DZ,FI,DFI,ZST,ZEN,NZ)
      IMPLICIT NONE
      INTEGER :: ZST(3),ZEN(3),NZ
      real(dp) :: DZ(1:NZ+1,1:NZ+1)
      COMPLEX(DP) ::  FI(ZST(1):ZEN(1),ZST(2):ZEN(2),ZST(3):ZEN(3))
      COMPLEX(DP) :: DFI(ZST(1):ZEN(1),ZST(2):ZEN(2),ZST(3):ZEN(3))
      
      INTEGER :: N(3),NI,NJ,NK
      COMPLEX(DP) ::  ALPHA,BETA
      
      ALPHA = 1._dp
      BETA = 0._dp
      N = ZEN - ZST + 1
      NI = N(1)
      NJ = N(2)
      NK = N(3)
      CALL ZGEMM('N','T',NI*NJ,NK,NK, ALPHA, FI, NI*NJ , CMPLX(DZ,0._dp), NK , BETA , DFI, NI*NJ )
      
    END SUBROUTINE TENSOR_dpODUCT_K
  

  end subroutine KERNEL_SOLVE_FFT
    
end module m_solver_diag_cart_fft_DFT
