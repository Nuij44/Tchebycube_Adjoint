module m_lap1d_tcheby
  use decomp_2d
  use m_boundary_conditions
  implicit none
  
  type t_lap1d_tcheby
     integer :: N 
     real(kind=8)             :: cl(2,2)
     real(kind=8),allocatable :: P(:,:),PM1(:,:),lbd(:)
     real(kind=8),allocatable :: FL(:),FR(:)
     real(kind=8),allocatable :: BL(:,:),BR(:,:),D2(:,:)
     real(kind=8),allocatable :: CC(:,:)
   contains
     procedure :: initialize => init_lap1d_tcheby
     procedure :: free => free_lap1d_tcheby
  end type t_lap1d_tcheby
  
contains
  
  subroutine init_lap1d_tcheby(this,op,nu,alpha_L,beta_L,alpha_R,beta_R,n)
    use m_operator_tcheby
    implicit none
    !.. Args
    class(t_lap1d_tcheby ),target,intent(out)   ::  this
    type(t_operator_tcheby), intent(in) ::  op
    real(dp),intent(in)  ::  nu
    real(dp),intent(in)  ::  alpha_L, beta_L, alpha_R, beta_R
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
    ALLOCATE( THIS%BL(1:N+1,1:N+1) , SOURCE = ( ALPHA_L*OP%GET_ID() + BETA_L*OP%GET_D1() ) )
    ALLOCATE( THIS%BR(1:N+1,1:N+1) , SOURCE = ( ALPHA_R*OP%GET_ID() + BETA_R*OP%GET_D1() ) )
    
    !... reduction of the operator
    allocate( this%p  (1:n+1,1:n+1) , source=0.0_dp )
    allocate( this%pm1(1:n+1,1:n+1) , source=0.0_dp )
    allocate( this%lbd(1:n+1)     , source=0.0_dp )
    allocate( this%fl (1:n+1)     , source=0.0_dp )
    allocate( this%fr (1:n+1)     , source=0.0_dp )
    allocate( this%cc (1:2,1:n+1) , source=0.0_dp )
    
    call reduction_operator(       &
         this%d2,this%BL,this%BR,   & ! problem definition
         this%P,this%PM1,this%LBD, & ! diagonalization
         this%CL,this%cc,this%FL,this%FR,  & ! aux variables
         this%N)
    
  end subroutine init_lap1d_tcheby
  
  subroutine free_lap1d_tcheby(this)
    implicit none
    class(t_lap1d_tcheby) :: this
    deallocate(this%P, this%PM1, this%LBD, this%FL, this%FR)
    deallocate(this%BL, this%BR,this%D2)
  end subroutine free_lap1d_tcheby
  
  subroutine reduction_operator(D2,BL,BR,P,PM1,LBD,CL,CC,FL,FR,N)
    implicit none
    !... args in
    integer :: n
    real(kind=dp),intent(in) :: d2(1:n+1,1:n+1)
    real(kind=dp),intent(in) :: bl(1:n+1,1:n+1)
    real(kind=dp),intent(in) :: br(1:n+1,1:n+1)
    !... args out
    real(kind=dp),intent(out) :: P(1:N+1,1:N+1),PM1(1:N+1,1:N+1)
    real(kind=dp),intent(out) :: LBD(1:N+1)
    real(kind=dp),intent(out) :: FL(1:N+1),FR(1:N+1)
    real(kind=dp),intent(out) :: CL(1:2,1:2),CC(1:2,1:N+1)
    !... variables
    real(kind=dp) :: D2_(2:N,2:N)
    real(kind=dp) :: P_(2:N,2:N),PM1_(2:N,2:N)
    real(kind=dp) :: LBD_(2:N)
    
    integer :: i,j
        
    ! inverse of cl attention au det
    ! CL(1,1:2) = BL(0,0) , BL(0,N)
    ! CL(2,1:2) = BR(N,0) , BR(N,N)
    
    cl(1,1:2) = [   br(n+1,n+1) , - bl(1,n+1) ]
    cl(2,1:2) = [ - br(n+1,  1) ,   bl(1,  1) ] 
    cl = cl / ( bl(1,1)*br(n+1,n+1) - br(n+1,1)*bl(1,n+1)  )
    
    forall(i=2:n,j=2:n)
       d2_(i,j) = d2(i,j) &
            - d2(i,1  )*( cl(1,1)*bl(1,j) + cl(1,2)*br(n+1,j)  ) &
            - d2(i,n+1)*( cl(2,1)*bl(1,j) + cl(2,2)*br(n+1,j)  ) 
    end forall
    
    call diagonalise_real( d2_, p_, pm1_, lbd_, n-1 )
    
    !> 
    forall(i=1:n+1) p(i,i) = 1._dp
    forall(i=1:n+1) pm1(i,i) = 1._dp
    
    forall(i=2:n,j=2:n)
       p(i,j) = p_(i,j)
       pm1(i,j) = pm1_(i,j)
    end forall
    
    CC(1,  1) = 1
    CC(1,2:n) = -bl(1,2:n)
    CC(1,n+1) = 0
    
    CC(2,  1) = 0
    CC(2,2:n) = -br(n+1,2:n)
    CC(2,n+1) = 1
    
    
    forall(i=1:n+1) lbd(i) = huge(dp)
    forall(i=2:n  ) lbd(i) = lbd_(i)
    
    forall(i=2:n)
       fl(i) = d2(i,1)*cl(1,1) + d2(i,n+1)*cl(2,1)
       fr(i) = d2(i,1)*cl(1,2) + d2(i,n+1)*cl(2,2)
    end forall
    
  end subroutine reduction_operator
    
end module m_lap1d_tcheby
