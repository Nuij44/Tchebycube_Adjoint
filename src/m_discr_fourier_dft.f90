module m_discr_fourier_DFT
  use m_mesh_base
  use m_fourier_transform
  use m_tensor_product
  use decomp_2d  
  use decomp_2d_mpi
  implicit none
  

  type t_discr_fourier_DFT
     integer :: N(3)
     integer :: order_of_derivative
     real(kind=8) , allocatable :: w(:)
   contains
     procedure :: init_discr_fourier_dft 
     procedure :: eval_discr_fourier_dft_x,eval_discr_fourier_dft_y,eval_discr_fourier_dft_z
     procedure :: eval_discr_fourier_dft_x_clean,eval_discr_fourier_dft_y_clean,eval_discr_fourier_dft_z_clean
  end type t_discr_fourier_DFT
  
contains

  
  SUBROUTINE INIT_DISCR_FOURIER_DFT(this,xmin,xmax,n,OPT)
    implicit none
    !.. args
    class(t_discr_fourier_DFT ),target :: this
    real(kind=8) :: xmin,xmax,wm
    integer :: n,opt
    real(kind=8),parameter :: pi = acos(-1d0)

    !.. variables
    integer :: i,im
    
    if (opt==100) then
       this%order_of_derivative = 0
    else if (opt==101) then
       this%order_of_derivative = 1
    else if (opt==102) then
       this%order_of_derivative = 2
    end if
    ! the mesh size is n+1
    
    ALLOCATE( THIS%W(1:N+1) ) ; THIS%W = 0.0
    
    IF (OPT==100) THEN

       DO I=1,N+1
          THIS%W(I) = 1._DP
       END DO
       
    ELSE IF (OPT==101) THEN
       
       DO I=1,N+1
          IM = HALF_COMPLEX_WP(I,N+1)
          THIS%W(I) = IM*2*PI/(XMAX-XMIN) / (n+1)
       END DO
       
    ELSE IF (OPT==102) THEN
       
       DO I=1,N+1
          IM = HALF_COMPLEX_WP(I,N+1)
          THIS%W(I) = -(IM*2*PI/(XMAX-XMIN))**2 / (n+1)
       END DO

    END IF

    THIS%W(n+1)=0
    
  end subroutine init_discr_fourier_DFT
  
  SUBROUTINE eval_discr_fourier_DFT_x(this,fi,dfi)
    implicit none
    CLASS(T_DISCR_FOURIER_DFT) :: THIS
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: FI,DFI
    TYPE(DECOMP_INFO) :: PH
    INTEGER ::I,J,K,n,im

    
    CALL GET_DECOMP_INFO(PH)
    !CALL FFTW_EXECUTE_R2R(PLAN_R2HC_X,FI,DG1_X)
    CALL R2R_1M_X(FI,DG1_X,PLAN_R2HC_X)
    
    n = PH%XEN(1)-PH%XST(1)+1
    
    if (THIS%order_of_derivative==1) then
       
       FORALL(I=1:1,J=PH%XST(2):PH%XEN(2),K=PH%XST(3):PH%XEN(3))
          DG2_X(I,J,K) = 0d0
       END FORALL
       FORALL(I=2:N/2+1,J=PH%XST(2):PH%XEN(2),K=PH%XST(3):PH%XEN(3))
          DG2_X(    I,J,K) = -DG1_X(N+2-I,J,K)*THIS%W(I)
          DG2_X(N+2-I,J,K) =  DG1_X(    I,J,K)*THIS%W(I)
       END FORALL
       
    elseif (THIS%order_of_derivative==2) then

       FORALL(I=1:N,J=PH%XST(2):PH%XEN(2),K=PH%XST(3):PH%XEN(3))
          DG2_X( I,J,K) = DG1_X(I,J,K)*THIS%W(I)
       END FORALL
       
    end if
    
!    CALL FFTW_EXECUTE_R2R(PLAN_HC2R_X,DG2_X,DFI)
    CALL R2R_1M_X(DG2_x,DFI,PLAN_HC2R_X)

  END SUBROUTINE EVAL_DISCR_FOURIER_DFT_X


  SUBROUTINE eval_discr_fourier_DFT_y(this,fi,dfi)
    implicit none
    CLASS(t_discr_fourier_DFT) :: this
    real(kind=8),dimension(:,:,:),allocatable :: fi,dfi
    integer :: i,j,k,nj,jm
    TYPE(DECOMP_INFO) :: PH

    CALL GET_DECOMP_INFO(PH)
    CALL TRANSPOSE_X_TO_Y(FI,DG2_Y)
    
    CALL R2R_1M_Y(DG2_Y,DG1_Y,PLAN_R2HC_Y)

    NJ = PH%YEN(2)-PH%YST(2)+1
    
    if (THIS%order_of_derivative==1) then

       FORALL(J=1:1,I=PH%YST(1):PH%YEN(1),K=PH%YST(3):PH%YEN(3))
          DG2_Y(I,J,K) = 0d0
       END FORALL
       FORALL(J=2:NJ/2+1,I=PH%YST(1):PH%YEN(1),K=PH%YST(3):PH%YEN(3))
          DG2_Y(I,     J,K) = -DG1_Y(I,NJ+2-J,K)*THIS%W(J)
          DG2_Y(I,NJ+2-J,K) =  DG1_Y(I,     J,K)*THIS%W(J)
       END FORALL

    else if (THIS%order_of_derivative==2) then

       FORALL(J=1:NJ,I=PH%YST(1):PH%YEN(1),K=PH%YST(3):PH%YEN(3))
          DG2_Y(I,     J,K) = DG1_Y(I,J,K)*THIS%W(J)
       END FORALL
    end if
       
    CALL R2R_1M_Y(DG2_Y,DG1_Y,PLAN_HC2R_Y)
   
    CALL TRANSPOSE_Y_TO_X(DG1_Y,DFI)
    

    
  END SUBROUTINE EVAL_DISCR_FOURIER_DFT_Y


  SUBROUTINE eval_discr_fourier_DFT_z(this,fi,dfi)
    implicit none
    CLASS(t_discr_fourier_DFT) :: this
    real(kind=8),dimension(:,:,:),allocatable :: fi,dfi
    integer :: i,j,k,nj,jm
    TYPE(DECOMP_INFO) :: PH

    print*,'Pas de dealiazing en z, à coder.'

    
  END SUBROUTINE EVAL_DISCR_FOURIER_DFT_Z

  SUBROUTINE eval_discr_fourier_DFT_x_clean(this,fi,dfi)
    implicit none
    CLASS(T_DISCR_FOURIER_DFT) :: THIS
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: FI,DFI
    TYPE(DECOMP_INFO) :: PH
    INTEGER ::I,J,K,n,im,K_max,cutoff

    
    CALL GET_DECOMP_INFO(PH)
    !CALL FFTW_EXECUTE_R2R(PLAN_R2HC_X,FI,DG1_X)
    CALL R2R_1M_X(FI,DG1_X,PLAN_R2HC_X)
    
    n = PH%XEN(1)-PH%XST(1)+1
    
    K_max = floor(1._DP/3._DP * N)
    

    if (THIS%order_of_derivative==1) then

       FORALL(I=1:1,J=PH%XST(2):PH%XEN(2),K=PH%XST(3):PH%XEN(3))
          DG2_X(I,J,K) = 0d0
       END FORALL
       FORALL(I=2:N/2+1,J=PH%XST(2):PH%XEN(2),K=PH%XST(3):PH%XEN(3))
          DG2_X(    I,J,K) = -DG1_X(N+2-I,J,K)*THIS%W(I)
          DG2_X(N+2-I,J,K) =  DG1_X(    I,J,K)*THIS%W(I)
       END FORALL
       
    elseif (THIS%order_of_derivative==2) then

       FORALL(I=1:N,J=PH%XST(2):PH%XEN(2),K=PH%XST(3):PH%XEN(3))
          DG2_X( I,J,K) = DG1_X(I,J,K)*THIS%W(I)
       END FORALL

    elseif (THIS%order_of_derivative==0) then

       FORALL(I=1:N,J=PH%XST(2):PH%XEN(2),K=PH%XST(3):PH%XEN(3))
          DG2_X( I,J,K) = DG1_X(I,J,K)*THIS%W(I)
       END FORALL

    end if
    
    !Nettoyage du spectre avec la règle des 2/3
    DO i = 1,N-2*k_max+1
       DG2_X(k_max+i,:,:) = 0._DP
    END DO

!    CALL FFTW_EXECUTE_R2R(PLAN_HC2R_X,DG2_X,DFI)
    CALL R2R_1M_X(DG2_x,DFI,PLAN_HC2R_X)    
  END SUBROUTINE EVAL_DISCR_FOURIER_DFT_X_CLEAN

  SUBROUTINE eval_discr_fourier_DFT_y_clean(this,fi,dfi)
    implicit none
    CLASS(t_discr_fourier_DFT) :: this
    real(kind=8),dimension(:,:,:),allocatable :: fi,dfi
    integer :: i,j,k,nj,jm,k_max
    TYPE(DECOMP_INFO) :: PH

    CALL GET_DECOMP_INFO(PH)
    CALL TRANSPOSE_X_TO_Y(FI,DG2_Y)
    
    CALL R2R_1M_Y(DG2_Y,DG1_Y,PLAN_R2HC_Y)

    NJ = PH%YEN(2)-PH%YST(2)+1
    
    k_max = floor(1._DP/3._DP*NJ)
    
    if (THIS%order_of_derivative==1) then

       FORALL(J=1:1,I=PH%YST(1):PH%YEN(1),K=PH%YST(3):PH%YEN(3))
          DG2_Y(I,J,K) = 0d0
       END FORALL
       FORALL(J=2:NJ/2+1,I=PH%YST(1):PH%YEN(1),K=PH%YST(3):PH%YEN(3))
          DG2_Y(I,     J,K) = -DG1_Y(I,NJ+2-J,K)*THIS%W(J)
          DG2_Y(I,NJ+2-J,K) =  DG1_Y(I,     J,K)*THIS%W(J)
       END FORALL

    else if (THIS%order_of_derivative==2) then

       FORALL(J=1:NJ,I=PH%YST(1):PH%YEN(1),K=PH%YST(3):PH%YEN(3))
          DG2_Y(I,     J,K) = DG1_Y(I,J,K)*THIS%W(J)
       END FORALL
    end if

    !Nettoyage du spectre avec la règle des 2/3
    DO j = 1,NJ-2*k_max+1
       DG2_Y(:,k_max+j,:) = 0._DP
    END DO

       
    CALL R2R_1M_Y(DG2_Y,DG1_Y,PLAN_HC2R_Y)
   
    CALL TRANSPOSE_Y_TO_X(DG1_Y,DFI)
    

    
  END SUBROUTINE EVAL_DISCR_FOURIER_DFT_Y_CLEAN

  SUBROUTINE eval_discr_fourier_DFT_z_clean(this,fi,dfi)
    implicit none
    CLASS(t_discr_fourier_DFT) :: this
    real(kind=8),dimension(:,:,:),allocatable :: fi,dfi
    integer :: i,j,k,nj,jm
    TYPE(DECOMP_INFO) :: PH

    print*,'Pas de dealiazing en z, à coder.'

    
  END SUBROUTINE EVAL_DISCR_FOURIER_DFT_Z_CLEAN

      
end module m_discr_fourier_DFT
