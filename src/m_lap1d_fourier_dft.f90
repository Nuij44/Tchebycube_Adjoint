module m_lap1d_fourier_DFT
  use decomp_2d
  use m_boundary_conditions
  use m_operator_fourier_DFT
  implicit none
  
  type t_lap1d_fourier_DFT
     integer :: N 
     real(dp),allocatable :: wave(:)
   contains
     procedure :: initialize => init_lap1d_fourier_DFT
     procedure :: free => free_lap1d_fourier_DFT
  end type t_lap1d_fourier_DFT
  
contains
  
  subroutine init_lap1d_fourier_DFT(this,op,nu,n)

    implicit none
    !.. Args
    class(t_lap1d_fourier_DFT ),target,intent(out)   ::  this
    type(t_operator_fourier_DFT), intent(in) ::  op
    real(dp),intent(in)  ::  nu
    integer, intent(in)  ::  n
    !.. Variables
    integer                   ::  ierr, i
    REAL(dp) :: L
    
    
    this%n = n
    
    L = OP%DISCR_D1%L
    
    if (op%n.ne.this%n) then
       print*,'size incompatibility '
       call mpi_finalize(ierr)
       stop
    end if

    allocate( this%wave(1:n+1)     , source=0.0_dp )

    do i=1,(n+1)/2
       this%wave(i+1) = 2._dp*pi/L*i
       this%wave(n-i+2) = -2._dp*pi/L*i
    end do
    this%wave(1)=0._dp
       
    
  end subroutine init_lap1d_fourier_DFT
  
  subroutine free_lap1d_fourier_DFT(this)
    implicit none
    class(t_lap1d_fourier_DFT) :: this

  end subroutine free_lap1d_fourier_DFT
  
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

end module m_lap1d_fourier_DFT
