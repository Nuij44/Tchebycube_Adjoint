module m_tensor_product
  implicit none


  real(kind=8), allocatable,public :: dg1_x(:,:,:),dg2_x(:,:,:)
  real(kind=8), allocatable,public :: dg1_y(:,:,:),dg2_y(:,:,:)  
  real(kind=8), allocatable,public :: dg1_z(:,:,:),dg2_z(:,:,:)  
  
  
contains
  subroutine start_tensor_product(&
       xst,xen,nx,&
       yst,yen,ny,&
       zst,zen,nz )
    implicit none
    integer :: xst(3),xen(3),nx
    integer :: yst(3),yen(3),ny
    integer :: zst(3),zen(3),nz

    allocate(dg1_x(xst(1):xen(1),xst(2):xen(2),xst(3):xen(3)))
    allocate(dg2_x(xst(1):xen(1),xst(2):xen(2),xst(3):xen(3)))
    allocate(dg1_y(yst(1):yen(1),yst(2):yen(2),yst(3):yen(3)))
    allocate(dg2_y(yst(1):yen(1),yst(2):yen(2),yst(3):yen(3)))
    allocate(dg1_z(zst(1):zen(1),zst(2):zen(2),zst(3):zen(3)))
    allocate(dg2_z(zst(1):zen(1),zst(2):zen(2),zst(3):zen(3)))
    
  end subroutine start_tensor_product
  
  SUBROUTINE m_tensor_product_ut(xst,xen,nx,yst,yen,ny,zst,zen,nz )
    use decomp_2d
    use decomp_2d_mpi
    implicit none
    integer :: xst(3),xen(3),nx ! xen - xst + 1 = nx+1
    integer :: yst(3),yen(3),ny
    integer :: zst(3),zen(3),nz
    
    real(KIND=8), allocatable,dimension(:,:,:) :: DG1_X, DG2_X
    real(KIND=8), allocatable,dimension(:,:,:) :: DG1_Y, DG2_Y
    real(KIND=8), allocatable,dimension(:,:,:) :: DG1_Z, DG2_Z
    real(KIND=8), allocatable,dimension(:,:,:) :: WRK

    real(KIND=8), allocatable,dimension(:,:) :: PX, PXM1
    real(KIND=8), allocatable,dimension(:,:) :: PY, PYM1
    real(KIND=8), allocatable,dimension(:,:) :: PZ, PZM1

    real(kind=8) :: diff
    integer :: ierr
    
    ALLOCATE( DG1_X( XST(1):XEN(1) , XST(2):XEN(2) , XST(3):XEN(3) ) )
    ALLOCATE( DG2_X( XST(1):XEN(1) , XST(2):XEN(2) , XST(3):XEN(3) ) )
    
    ALLOCATE( DG1_Y( YST(1):YEN(1) , YST(2):YEN(2) , YST(3):YEN(3) ) )
    ALLOCATE( DG2_Y( YST(1):YEN(1) , YST(2):YEN(2) , YST(3):YEN(3) ) )

    ALLOCATE( DG1_Z( ZST(1):ZEN(1) , ZST(2):ZEN(2) , ZST(3):ZEN(3) ) )
    ALLOCATE( DG2_Z( ZST(1):ZEN(1) , ZST(2):ZEN(2) , ZST(3):ZEN(3) ) )

    

    
    CALL BUILD_OPERATOR(PX,PXM1,XST(1),XEN(1))
    CALL BUILD_OPERATOR(PY,PYM1,YST(2),YEN(2))
    CALL BUILD_OPERATOR(PZ,PZM1,ZST(3),ZEN(3))

    ALLOCATE( WRK( XST(1):XEN(1) , XST(2):XEN(2) , XST(3):XEN(3) ) )

    CALL RANDOM_NUMBER(DG1_X)
    WRK = DG1_X
    
    CALL TENSOR_PRODUCT_I( PXM1 , DG1_X, DG2_X, XST, XEN )
    
    CALL TRANSPOSE_X_TO_Y( DG2_X, DG1_Y )
    CALL TENSOR_PRODUCT_J( PYM1  , DG1_Y, DG2_Y, YST, YEN )
    
    CALL TRANSPOSE_Y_TO_Z( DG2_Y, DG1_Z )
    CALL TENSOR_PRODUCT_K( PZM1  , DG1_Z, DG2_Z, ZST, ZEN )

    
    CALL TENSOR_PRODUCT_K( PZ   , DG2_Z, DG1_Z, ZST, ZEN )
    CALL TRANSPOSE_Z_TO_Y( DG1_Z, DG1_Y )
    
    CALL TENSOR_PRODUCT_J( PY   , DG1_Y, DG2_Y, YST, YEN)
    CALL TRANSPOSE_Y_TO_X( DG2_Y, DG1_X )
    
    CALL TENSOR_PRODUCT_I( PX   , DG1_X, DG2_X, XST, XEN )

    diff = maxval(abs(DG2_X-WRK))

    CALL MPI_ALLREDUCE(MPI_IN_PLACE,diff,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,IERR)
    if (nrank==0) then
       print'("error on matrix base change",e15.8)',diff
    end if
    
    deallocate(DG1_X,DG2_X)
    deallocate(DG1_Y,DG2_Y)
    deallocate(DG1_Z,DG2_Z)
    deallocate(PX,PXM1)
    deallocate(PY,PYM1)
    deallocate(PZ,PZM1)
    deallocate(WRK)
    
  contains
    
    subroutine build_operator(P,PM1,st,en)
      implicit none
      real(kind=8),dimension(:,:),allocatable :: P,PM1
      real(kind=8),dimension(:),allocatable :: VP
      integer :: st,en
      integer :: i,ii
      integer :: j,jj
      real(kind=8) :: pi = acos(dble(-1.))
      integer :: n
      real(kind=8),dimension(:,:),allocatable :: D
      
      allocate (P  (st:en,st:en))
      allocate (PM1(st:en,st:en))

      allocate (D  (st:en,st:en))
      n = en-st+1
      do i=st,en
         do j=st,en
            ii = (i-st)
            jj = (j-st)
            if (ii==jj) then
               D(i,j) = 0
            else
               D(i,j) = (0.5)*(-1)**(ii+jj)*sin(dble(ii-jj)*pi/(n))**(-1)
            end if
         end do
      end do
      
      D=matmul(D,D)
      allocate (VP(st:en))
      n = en-st+1
      CALL DIAGONALISE_REAL(D,P,PM1,VP,n)
      deallocate(vp)
      
    end subroutine build_operator
    
    
    SUBROUTINE DIAGONALISE_REAL(A,P,PM1,VP,N)
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: N
      REAL(KIND=8),INTENT(IN) , DIMENSION(N,N)::A
      REAL(KIND=8),INTENT(OUT), DIMENSION(N,N)::P,PM1
      REAL(KIND=8),INTENT(OUT), DIMENSION(N)  ::VP
      REAL(KIND=8),DIMENSION(N,N)   :: AWORK
      REAL(KIND=8),DIMENSION(N*(N+6))   :: WORK
      REAL(KIND=8),DIMENSION(N)   :: VPI
      INTEGER                     :: INFO,I
      INTEGER     ,DIMENSION(N)   :: IPIV
      external  ::  DGEEV, DGETRF, DGETRI
      
      AWORK = A
      CALL DGEEV('N','V',N,AWORK,N,VP,VPI,PM1,N,P,N,WORK,4*N,INFO)
      
    
      IF (INFO.NE.0) THEN
         PRINT*,'ECHEC DE  REAL_DIAG'
         PRINT*,INFO
      END IF
      
      DO I=1,N
         IF (ABS(VPI(I)).GT.1D-8) THEN
            PRINT*,'PRESENCE DE VALEURS PROPRES COMPLEXES'
            PRINT*,I,VP(I),VPI(I)
            PRINT*,'ECHEC DE  REAL_DIAG'
            STOP
       ENDIF
    ENDDO
    
    PM1=P
    CALL DGETRF(N,N,PM1,N,IPIV,INFO)
    CALL DGETRI(N,PM1,N,IPIV,WORK,N,INFO)
    IF (INFO.NE.0) THEN
       PRINT*,'ECHEC DE  REAL_DIAG DANS L INVERSION '
       PRINT*,INFO
       STOP
    END IF
  END SUBROUTINE DIAGONALISE_REAL




    
  END SUBROUTINE M_TENSOR_PRODUCT_UT
  
  SUBROUTINE TENSOR_PRODUCT_CL_I(CC,FI,DFI,XST,XEN)
    IMPLICIT NONE
    INTEGER      :: XST(3),XEN(3)
    real(KIND=8) :: CC(1:2,XST(1):XEN(1))
    real(KIND=8) :: FI(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
    real(KIND=8) :: DFI(1:2,XST(2):XEN(2),XST(3):XEN(3))
    
    real(KIND=8) ::  ALPHA,BETA
    INTEGER      :: N(3),NI,NJ,NK
    
    N = XEN-XST + 1
    NI = N(1)
    NJ = N(2)
    NK = N(3)
    ALPHA = dble(1.0)
    BETA = dble(0.0)
    
    CALL DGEMM( 'N', 'N', 2, NJ*NK, NI, ALPHA, CC, 2, FI, NI, BETA, DFI, 2 )
    
  END SUBROUTINE TENSOR_PRODUCT_CL_I
  
  SUBROUTINE TENSOR_PRODUCT_CL_J(CC,FI,DFI,YST,YEN)
    IMPLICIT NONE
    INTEGER      :: YST(3),YEN(3)
    real(kind=8) :: CC(1:2,YST(2):YEN(2))
    real(kind=8) ::  FI(YST(1):YEN(1),YST(2):YEN(2),YST(3):YEN(3))
    real(kind=8) :: DFI(YST(1):YEN(1),1:2,YST(3):YEN(3))
    real(kind=8) ::  ALPHA,BETA
    INTEGER      :: K,N(3),NI,NJ
    
    ALPHA = dble(1.)
    BETA = dble(0.)
    N = YEN - YST + 1
    NI = N(1)
    NJ = N(2)
    DO K=YST(3),YEN(3)
       CALL DGEMM('N','T',NI,2,NJ,ALPHA,FI(YST(1),YST(2),K),NI,CC,2,BETA,DFI(YST(1),1,K),NI)
    END DO
  END SUBROUTINE TENSOR_PRODUCT_CL_J
  
  
  SUBROUTINE TENSOR_PRODUCT_CL_K(CC,FI,DFI,ZST,ZEN)
    IMPLICIT NONE
    INTEGER      :: ZST(3),ZEN(3)
    REAL(KIND=8) :: CC(1:2,ZST(3):ZEN(3))
    REAL(KIND=8) ::  FI(ZST(1):ZEN(1),ZST(2):ZEN(2),ZST(3):ZEN(3))
    REAL(KIND=8) :: DFI(ZST(1):ZEN(1),ZST(2):ZEN(2),1:2)
    
    INTEGER      :: N(3),NI,NJ,NK
    real(KIND=8) :: ALPHA,BETA
    
    ALPHA = dble(1.)
    BETA = dble(0.)
    N = ZEN - ZST + 1
    NI = N(1)
    NJ = N(2)
    NK = N(3) 
    CALL DGEMM('N','T',NI*NJ,2,NK, ALPHA, FI, NI*NJ , CC, 2 , BETA , DFI, NI*NJ )
        
  END SUBROUTINE TENSOR_PRODUCT_CL_K
  
      
  subroutine tensor_product_i(DX,FI,DFI,XST,XEN)
    implicit none
    INTEGER :: XST(3),XEN(3)
    real(kind=8) ::  DX(XST(1):XEN(1),XST(1):XEN(1))
    real(kind=8) ::  FI(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
    real(kind=8) :: DFI(XST(1):XEN(1),XST(2):XEN(2),XST(3):XEN(3))
        
    real(kind=8) ::  ALPHA,BETA
    INTEGER :: N(3),NI,NJ,NK
    
    N = XEN-XST + 1
    NI = N(1)
    NJ = N(2)
    NK = N(3)
    ALPHA = dble(1.)
    BETA = dble(0.0)
    
    CALL DGEMM( 'N', 'N', NI, NJ*NK, NI, ALPHA, DX, NI, FI, NI, BETA, DFI, NI )
    
  end subroutine tensor_product_i
  
  SUBROUTINE TENSOR_PRODUCT_J(DY,FI,DFI,YST,YEN)
    IMPLICIT NONE
    INTEGER :: YST(3),YEN(3)
    real(kind=8) :: DY(YST(2):YEN(2),YST(2):YEN(2))
    real(kind=8) ::  FI(YST(1):YEN(1),YST(2):YEN(2),YST(3):YEN(3))
    real(kind=8) :: DFI(YST(1):YEN(1),YST(2):YEN(2),YST(3):YEN(3))
    real(kind=8) ::  ALPHA,BETA
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
       CALL DGEMM('N','T',NI,NJ,NJ, ALPHA ,FI(IS,JS,K), NI , DY, NJ , BETA , DFI(IS,JS,K), NI )
    END DO


  END SUBROUTINE TENSOR_PRODUCT_J
  
  SUBROUTINE TENSOR_PRODUCT_K(DZ,FI,DFI,ZST,ZEN)
    IMPLICIT NONE
    INTEGER      :: ZST(3),ZEN(3)
    real(kind=8) ::  DZ(ZST(3):ZEN(3),ZST(3):ZEN(3))
    real(kind=8) ::  FI(ZST(1):ZEN(1),ZST(2):ZEN(2),ZST(3):ZEN(3))
    real(kind=8) :: DFI(ZST(1):ZEN(1),ZST(2):ZEN(2),ZST(3):ZEN(3))
    
    INTEGER      :: N(3),NI,NJ,NK
    real(kind=8) ::  ALPHA,BETA
        
    ALPHA = dble(1.)
    BETA = dble(0.)
    N = ZEN - ZST + 1
    NI = N(1)
    NJ = N(2)
    NK = N(3)
    CALL DGEMM('N','T',NI*NJ,NK,NK, ALPHA, FI, NI*NJ , DZ, NK , BETA , DFI, NI*NJ )
    
  END SUBROUTINE TENSOR_PRODUCT_K
  
end module m_tensor_product
  
