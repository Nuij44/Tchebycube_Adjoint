module m_lap1d_fourier
  use decomp_2d
  use m_boundary_conditions
  implicit none
  
  type t_lap1d_fourier
     integer :: N 
     real(kind=8),allocatable :: P(:,:),PM1(:,:),lbd(:),D2(:,:)
   contains
     procedure :: initialize => init_lap1d_fourier
     procedure :: free => free_lap1d_fourier
  end type t_lap1d_fourier
  
contains
  
  subroutine init_lap1d_fourier(this,op,nu,n)
    use m_operator_fourier
    implicit none
    !.. Args
    class(t_lap1d_fourier ),target,intent(out)   ::  this
    type(t_operator_fourier), intent(in) ::  op
    real(dp),intent(in)  ::  nu
    integer, intent(in) :: n
    !.. Variables
    integer                   ::  ierr
    
    
    this%n = n

    
    if (op%n.ne.this%n) then
       print*,'size incompatibility '
       call mpi_finalize(ierr)
       stop
    end if


    
    ALLOCATE( THIS%D2(1:N+1,1:N+1) , SOURCE = nu*OP%GET_D2() )
    
    !... reduction of the operator
    allocate( this%p  (1:n+1,1:n+1) , source=0.0_dp )
    allocate( this%pm1(1:n+1,1:n+1) , source=0.0_dp )
    allocate( this%lbd(1:n+1)     , source=0.0_dp )
    
    CALL REDUCTION_OPERATOR_PERIODIC( &
         THIS%D2, THIS%P, THIS%PM1, THIS%LBD, THIS%N)
    
  end subroutine init_lap1d_fourier
  
  subroutine free_lap1d_fourier(this)
    implicit none
    class(t_lap1d_fourier) :: this
    deallocate(this%P, this%PM1, this%LBD)
    deallocate(this%D2)
  end subroutine free_lap1d_fourier
  
  subroutine reduction_operator_PERIODIC(D2,P,PM1,LBD,N)
    implicit none
    !... args in
    integer :: n
    real(kind=dp),intent(in) :: d2(1:n+1,1:n+1)
    !... args out
    real(kind=dp),intent(out) :: P(1:N+1,1:N+1),PM1(1:N+1,1:N+1)
    real(kind=dp),intent(out) :: LBD(1:N+1)
    !... variables
    real(kind=dp) :: D2_(1:N+1,1:N+1)
    real(kind=dp) :: P_(1:N+1,1:N+1),PM1_(1:N+1,1:N+1)
    real(kind=dp) :: LBD_(1:N+1)
    
    integer :: i,j
    
    forall(i=1:n+1,j=1:n+1)
       d2_(i,j) = d2(i,j) 
    end forall
    call diagonalise_real( d2_, p_, pm1_, lbd_, n+1 )
    forall(i=1:n+1,j=1:n+1)
       p(i,j) = p_(i,j)
       pm1(i,j) = pm1_(i,j)
    end forall
    forall(i=1:n+1) lbd(i) = lbd_(i)
    
  end subroutine reduction_operator_PERIODIC

end module m_lap1d_fourier
