module m_discr_tcheby_dct
  use m_mesh_base
  use m_fourier_transform
  
  implicit none
  

  type t_discr_tcheby_dct
     integer :: N(3)
     integer :: order_of_derivative
     real(kind=8) , allocatable :: w(:)
   contains
     procedure :: init_discr_tcheby_dct 
     procedure :: eval_discr_tcheby_dct_x,eval_discr_tcheby_dct_y
  end type t_discr_tcheby_dct
  
contains

  
  SUBROUTINE INIT_DISCR_TCHEBY_DCT(this,xmin,xmax,n,OPT)
    use m_fourier_transform
    implicit none
    !.. args
    class(t_discr_tcheby_dct ),target :: this
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
          THIS%W(I) = 0D0
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
    
  end subroutine init_discr_tcheby_dct
  
  SUBROUTINE eval_discr_tcheby_dct_x(this,fi,dfi)
    use m_tensor_product
    use m_fourier_transform
    use decomp_2d
    use decomp_2d_mpi
    implicit none
    CLASS(T_DISCR_TCHEBY_DCT) :: THIS
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



    
    
  END SUBROUTINE EVAL_DISCR_TCHEBY_DCT_X


  SUBROUTINE eval_discr_tcheby_dct_y(this,fi,dfi)
    use m_tensor_product
    use m_fourier_transform
    use decomp_2d
    use decomp_2d_mpi
    implicit none
    CLASS(t_discr_tcheby_dct) :: this
    real(kind=8),dimension(:,:,:),allocatable :: fi,dfi
    integer :: i,j,k,nj,jm
    TYPE(DECOMP_INFO) :: PH

    CALL GET_DECOMP_INFO(PH)
    CALL TRANSPOSE_X_TO_Y(FI,DG2_Y)
    
    CALL R2R_1M_Y(DG2_Y,DG1_Y,PLAN_R2HC_Y)

    NJ = PH%YEN(2)-PH%YST(2)+1
    
    if (THIS%order_of_derivative==1) then

!!$       do j=2,nj/2-1
!!$          JM = HALF_COMPLEX_WP(j,NJ)
!!$          print*,j,jm,DG1_Y(3,J,3),DG1_Y(4,NJ+2-j,3)
!!$       end do
       
!!$       FORALL(I=PH%YST(1):PH%YEN(1),K=PH%YST(3):PH%YEN(3))
!!$          DG2_Y(i,1,k)=0
!!$          DG2_Y(i,NJ,k)=0
!!$       end FORALL
       
       FORALL(J=1:1,I=PH%YST(1):PH%YEN(1),K=PH%YST(3):PH%YEN(3))
          DG2_Y(I,J,K) = 0d0
       END FORALL
       FORALL(J=2:NJ/2+1,I=PH%YST(1):PH%YEN(1),K=PH%YST(3):PH%YEN(3))
          DG2_Y(I,     J,K) = -DG1_Y(I,NJ+2-J,K)*THIS%W(J)
          DG2_Y(I,NJ+2-J,K) =  DG1_Y(I,     J,K)*THIS%W(J)
       END FORALL

       
!!$       print*,"---"
!!$       print*,1,0,DG2_Y(2,1,2),DG2_Y(2,1,2)
!!$       do j=2,nj/2+1
!!$          JM = HALF_COMPLEX_WP(j,NJ)
!!$          print*,j,jm,DG2_Y(2,J,2),DG2_Y(2,NJ+2-j,2)
!!$       end do
       
    else if (THIS%order_of_derivative==2) then

!!$       do j=2,nj/2-1
!!$          JM = HALF_COMPLEX_WP(j,NJ)
!!$          print*,j,jm,DG1_Y(3,J,3),DG1_Y(4,NJ+2-j,3)
!!$       end do
       
       FORALL(J=1:NJ,I=PH%YST(1):PH%YEN(1),K=PH%YST(3):PH%YEN(3))
          DG2_Y(I,     J,K) = DG1_Y(I,J,K)*THIS%W(J)
       END FORALL
    end if
       
    CALL R2R_1M_Y(DG2_Y,DG1_Y,PLAN_HC2R_Y)
   
    CALL TRANSPOSE_Y_TO_X(DG1_Y,DFI)
    

    
  END SUBROUTINE EVAL_DISCR_TCHEBY_DCT_Y

  SUBROUTINE eval_discr_tcheby_dct_z(this,fi,dfi)
    use decomp_2d
    implicit none
    CLASS(t_discr_tcheby_dct) :: this
    real(kind=8),dimension(:,:,:),allocatable :: fi,dfi
    
  END SUBROUTINE EVAL_DISCR_TCHEBY_DCT_Z
      
end module m_discr_tcheby_dct
