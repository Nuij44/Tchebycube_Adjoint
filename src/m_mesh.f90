module m_mesh_base
  use m_numerics
  implicit none
  type t_mesh_base
     integer :: N ! nb o
     real(kind=8) :: xmin,xmax
     real(kind=8),allocatable :: x(:)
     real(kind=8),allocatable :: g(:) ! $x^{-1}_X(X)$
     real(kind=8),allocatable :: h(:) ! $ to do $
     
   contains
     procedure :: free, print, initialize
  end type t_mesh_base
  
  private :: free,print, initialize
  
contains
  subroutine free(this)
    implicit none
    class(t_mesh_base ) :: this
    
    deallocate(this%x)
    deallocate(this%g)
    deallocate(this%h)
    
    this%N = 0
    this%xmin = 0
    this%xmax = 0
    
  end subroutine free

  subroutine print(this)
    implicit none
    !.. args
    class(t_mesh_base ) :: this
    !.. args
    integer :: i
    
    do i=0,this%N
       print*,i,this%x(i)
    end do
    
  end subroutine print
  
  subroutine initialize(this,xmin,xmax,N,periodic)
    implicit none
    !.. args
    class(t_mesh_base ) :: this
    real(kind=8) :: xmin,xmax
    integer :: N
    !.. variables
    integer :: i
    real(kind=8) :: xx
    logical, optional:: periodic
    
    this%n = n
    this%xmin = xmin
    this%xmax = xmax
    
    allocate( this%x(1:n+1) , this%g(1:n+1) , this%h(1:n+1) )
    
    do i=1,N+1
       xx = cos( (i-1)*pi/n )
       this%x(i) =  0.5*(1-xx)*(xmax-xmin) + xmin
       this%g(i) = -2._dp/(xmax-xmin)
    end do

    if (present(periodic)) then
       do i=1,N+1
          xx = 2*PI*(i-1)/(n+1)
          this%x(i) = (xx/(2*pi))*(xmax-xmin) + xmin
       end do
    end if
    
  end subroutine initialize
end module m_mesh_base
