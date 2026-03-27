module m_numerics
    use iso_fortran_env
    implicit none
    INTEGER, parameter  :: dp= REAL64
    REAL(dp), parameter ::  pi = acos(-1._dp)
    REAL(dp), parameter ::  one=1.0_dp, zero = 0.0_dp

    
    interface 
       real(kind=8) function udf_space(x,y,z)
         implicit none
         real(kind=8),intent(in):: x,y,z
       end function udf_space
       
       pure function udf_timespace(t,x,y,z) result(res)
         import dp
         implicit none
         real(kind=dp),intent(in):: t,x,y,z
         real(kind=dp) :: res
       end function udf_timespace
    

  function udf_timespace_npure(t,x,y,z) result(res)
         import dp
         implicit none
         real(kind=dp),intent(in):: t,x,y,z
         real(kind=dp) :: res
       end function udf_timespace_npure


        end interface

    
contains


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

  SUBROUTINE DIAGONALISE_cplx(A,P,PM1,VP,N)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: N
    REAL(KIND=8),INTENT(IN) , DIMENSION(N,N)::A
    REAL(KIND=8),INTENT(OUT), DIMENSION(N,N)::P,PM1
    complex(KIND=8),INTENT(OUT), DIMENSION(N)  ::VP
    REAL(KIND=8),DIMENSION(N,N)   :: AWORK
    REAL(KIND=8),DIMENSION(N*(N+6))   :: WORK
    REAL(KIND=8),DIMENSION(N)   :: VPI,VPR
    INTEGER                     :: INFO,I
    INTEGER     ,DIMENSION(N)   :: IPIV
    external  ::  DGEEV, DGETRF, DGETRI
    
    AWORK = A
    CALL DGEEV('N','V',N,AWORK,N,VPR,VPI,PM1,N,P,N,WORK,4*N,INFO)
    
    
    IF (INFO.NE.0) THEN
       PRINT*,'ECHEC DE  REAL_DIAG'
       PRINT*,INFO
    END IF
    
    DO I=1,N
       VP(I)=DCMPLX(VPR(I),VPI(I))
       !print*, VP(I)
       IF (ABS(VPI(I)).GT.1D-8) THEN
          PRINT*,'PRESENCE DE VALEURS PROPRES COMPLEXES'
          PRINT*,I,VP(I),VPI(I)
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
  END SUBROUTINE DIAGONALISE_CPLX

  
  SUBROUTINE command_line_read_input(ifile)
    use mpi
    implicit none
    character(len=*) :: ifile
    integer i
    character(len=4096) :: arg
    integer :: stat,rank
    integer :: ierr
    integer :: ll
    
    
    call mpi_comm_rank(mpi_comm_world,rank,ierr)

    if (rank==0) then
       i = 0
       DO
          CALL get_command_argument(i, arg)
          i = i+1
          IF (LEN_TRIM(arg) == 0) EXIT
          select case (arg)
          case ('--input')
             CALL get_command_argument(i, arg) ; i = i+1
             print *,trim(arg)
             ifile=trim(arg)
          case default
          end select
       END DO
       
    end if
    
    ll = len(ifile)
    CALL MPI_BCAST(ifile(1:ll) ,ll,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
    
  end SUBROUTINE command_line_read_input

  SUBROUTINE command_line_prod_scal(input_file,A,B)
    use mpi
    implicit none
    character(len=*) :: input_file,A,B
    integer i
    character(len=4096) :: arg
    integer :: stat,rank
    integer :: ierr
    integer :: ll
    
    
    call mpi_comm_rank(mpi_comm_world,rank,ierr)

    if (rank==0) then
       i = 0
       DO
          CALL get_command_argument(i, arg)
          i = i+1
          IF (LEN_TRIM(arg) == 0) EXIT
          select case (arg)
          case ('--input')
             CALL get_command_argument(i, arg) ; i = i+1
             print *,trim(arg)
             input_file=trim(arg)
          case ('--a-file')
             CALL get_command_argument(i, arg) ; i = i+1
             print *,trim(arg)
             A=trim(arg)
          case ('--b-file')
             CALL get_command_argument(i, arg) ; i = i+1
             print *,trim(arg)
             B=trim(arg)
          case default
          end select
       END DO
       
    end if
    
    ll = len(input_file)
    CALL MPI_BCAST(input_file(1:ll) ,ll,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
    ll = len(A)
    CALL MPI_BCAST(A(1:ll) ,ll,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
    ll = len(B)
    CALL MPI_BCAST(B(1:ll) ,ll,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
    
  end SUBROUTINE command_line_prod_scal


   SUBROUTINE STR2INT(STR,INT,STAT)
    IMPLICIT NONE
    ! ARGUMENTS
    CHARACTER(LEN=*),INTENT(IN) :: STR
    INTEGER,INTENT(OUT)         :: INT
    INTEGER,INTENT(OUT)         :: STAT

    READ(STR,*,IOSTAT=STAT)  INT
  END SUBROUTINE STR2INT

  SUBROUTINE STR2R8(STR,R8,STAT)
    IMPLICIT NONE
    ! ARGUMENTS
    CHARACTER(LEN=*),INTENT(IN) :: STR
    REAL(KIND=8),INTENT(OUT)         :: R8
    INTEGER,INTENT(OUT)         :: STAT

    READ(STR,*,IOSTAT=STAT)  R8
  END SUBROUTINE STR2R8

END module m_numerics
