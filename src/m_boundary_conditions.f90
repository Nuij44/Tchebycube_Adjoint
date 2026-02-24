module m_boundary_conditions
  use m_numerics
  use m_mesh_base
  use decomp_2d
  implicit none
  
  type t_boundary_conditions
     integer :: N(3)
     real(kind=8) :: alpha_l,beta_l
     real(kind=8) :: alpha_r,beta_r
     real(kind=8),allocatable :: bv_l(:,:),bv_r(:,:)
   contains
     procedure :: initialize , free,set_bcs!, check
  end type t_boundary_conditions


!!$  type t_boundary_conditions_cplx
!!$     integer :: N(3)
!!$     complex(kind=8) :: alpha_l,beta_l
!!$     complex(kind=8) :: alpha_r,beta_r
!!$     complex(kind=8),allocatable :: bv_l(:,:),bv_r(:,:)
!!$   contains
!!$     procedure :: initialize , free,set_bcs!, check
!!$  end type t_boundary_conditions_cplx



  
contains

  subroutine set_bcs(this,axe,side,time,msh_xyz,udf,st,en,n)
    implicit none
    class(t_boundary_conditions ),intent(inout)   ::  this
    integer :: axe, side
    real(kind=8) :: time
    type(t_mesh_base) ::  msh_xyz(3)
    procedure(udf_timespace) :: udf
    integer :: st(3),en(3),n(3),i,j,k
    
    n = this%n
    
    if (side==-1) then
       if (axe==1) then

          do k=st(3),en(3)
             do j=st(2),en(2)
                this%bv_l(j,k) = udf(time,msh_xyz(1)%x(1),msh_xyz(2)%x(j),msh_xyz(3)%x(k))
             end do
          end do
          
          
       else if (axe==2) then
          do k=st(3),en(3)
             do i=st(1),en(1)
                this%bv_l(i,k) = udf(time,msh_xyz(1)%x(i),msh_xyz(2)%x(1),msh_xyz(3)%x(k))
             end do
          end do
       else if (axe==3) then
          do j=st(2),en(2)
             do i=st(1),en(1)
                this%bv_l(i,j) = udf(time,msh_xyz(1)%x(i),msh_xyz(2)%x(j),msh_xyz(3)%x(1))
             end do
          end do
       end if
       
    else if (side==+1) then
       
       if (axe==1) then
          do k=st(3),en(3)
             do j=st(2),en(2)
                this%bv_r(j,k) = udf(time,msh_xyz(1)%x(n(1)+1),msh_xyz(2)%x(j),msh_xyz(3)%x(k))
             end do
          end do
       else if (axe==2) then
          do k=st(3),en(3)
             do i=st(1),en(1)
                this%bv_r(i,k) = udf(time,msh_xyz(1)%x(i),msh_xyz(2)%x(n(2)+1),msh_xyz(3)%x(k))
             end do
          end do
       else if (axe==3) then
          do j=st(2),en(2)
             do i=st(1),en(1)
                this%bv_r(i,j) = udf(time,msh_xyz(1)%x(i),msh_xyz(2)%x(j),msh_xyz(3)%x(n(3)+1))
             end do
          end do
       end if
       
       
    end if
    
  end subroutine set_bcs
  
  
  subroutine initialize(this,axe,alpha_L,beta_L,alpha_R,beta_R,st,en,n)
    implicit none
    class(t_boundary_conditions ),intent(out)   ::  this
    integer ,intent(in) ::  axe
    real(dp),intent(in) ::  alpha_L, beta_L, alpha_R, beta_R
    integer ,intent(in) ::  st(3),en(3),n(3)
    
    this%n = n

    if (axe==1) then
       allocate(  this%bv_l(st(2):en(2),st(3):en(3)), source = 0._dp )
       allocate(  this%bv_r(st(2):en(2),st(3):en(3)), source = 0._dp ) 
    end if


    if (axe==2) then
       allocate(  this%bv_l(st(1):en(1),st(3):en(3)), source = 0._dp )
       allocate(  this%bv_r(st(1):en(1),st(3):en(3)), source = 0._dp ) 
    end if
    
    if (axe==3) then
       allocate(  this%bv_l(st(1):en(1),st(2):en(2)), source = 0._dp )
       allocate(  this%bv_r(st(1):en(1),st(2):en(2)), source = 0._dp ) 
    end if
    
    this%alpha_L = alpha_L
    this%beta_L = beta_L
    
    this%alpha_R = alpha_R
    this%beta_R = beta_R
    
  end subroutine initialize
  
  subroutine free(this)
    implicit none
    class(t_boundary_conditions ),target,intent(out)   ::  this
    deallocate(  this%bv_l )
    deallocate(  this%bv_r )
    
  end subroutine free

  
end module m_boundary_conditions


  
