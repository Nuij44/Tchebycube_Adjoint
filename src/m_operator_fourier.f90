module m_operator_fourier
    use m_operator_base
    use m_discr_fourier
    
    implicit none
    
    type, extends(t_operator_base)    ::  t_operator_fourier
       type(t_discr_fourier) :: discr_d1 
       type(t_discr_fourier) :: discr_d2 
       type(t_discr_fourier) :: discr_id 
     contains 
       procedure :: init_operator_fourier => init_operator_fourier
       procedure :: initialize => init_operator_fourier
       procedure :: d1 => t_operator_fourier_d1
       procedure :: d2 => t_operator_fourier_d2
       procedure :: id => t_operator_fourier_id
       procedure :: d1_clean => t_operator_fourier_d1_clean
       procedure :: d2_clean => t_operator_fourier_d2_clean
       procedure :: id_clean => t_operator_fourier_id_clean
       procedure :: Get_D1,Get_Id,Get_D2
    end type t_operator_fourier

    
    
  contains
    
    Function Get_D1(this)
      implicit none
      class(t_operator_fourier) :: this
      real(kind=8) Get_D1(1:this%n+1,1:this%n+1)
      Get_D1 = this%discr_d1%A
    End Function Get_D1

    Function Get_D2(this)
      implicit none
      class(t_operator_fourier) :: this
      real(kind=8) Get_D2(1:this%n+1,1:this%n+1)
      Get_D2 = this%discr_d2%A
    End Function Get_D2

    Function Get_Id(this)
      implicit none
      class(t_operator_fourier) :: this
      real(kind=8) Get_Id(1:this%n+1,1:this%n+1)
      Get_Id = this%discr_Id%A
    End Function Get_Id

    

    
    subroutine check_operator_fourier(this,mesh)
      use decomp_2d
      implicit none
      class(t_operator_fourier) :: this
      type(t_mesh_base)        :: mesh(3)
      
      real(kind=8),dimension(:,:,:),allocatable :: fi,dfi
      real(kind=8),dimension(:,:,:),allocatable :: dg1y,dg2y
      real(kind=8),dimension(:,:,:),allocatable :: dg1z,dg2z

      real(kind=8) :: x,y,z

      integer ::i,j,k
      
      TYPE(DECOMP_INFO)  :: dd 

      CALL GET_DECOMP_INFO(DD)

      call alloc_x(  fi , opt_global=.true. )
      call alloc_x( dfi , opt_global=.true. )

      call alloc_y( dg1y , opt_global=.true. )
      call alloc_y( dg2y , opt_global=.true. )
      
      call alloc_z( dg1z , opt_global=.true. )
      call alloc_z( dg2z , opt_global=.true. )
      
      

      do k=dd%xst(3),dd%xen(3)
         do j=dd%xst(2),dd%xen(2)
            do i=dd%xst(1),dd%xen(1)
               x = mesh(1)%x(i)
               y = mesh(2)%x(j)
               z = mesh(3)%x(k)
               fi(i,j,k) = x+2*y+3*z
            end do
         end do
      end do
      
      call this%d1(fi,dfi,dg1y,dg2y,dg1z,dg2z)
      !print*,dfi
      

      deallocate(fi,dfi,dg2y,dg1z,dg2z)
    end subroutine check_operator_fourier

    
    subroutine init_operator_fourier(this, mesh, axis)
      implicit none
      class(t_operator_fourier) :: this
      type(t_mesh_base)        :: mesh
      integer                  :: axis
      
      integer :: N
      
      
      N = mesh%N
      this%n = mesh%N
      this%axis = axis
      
      !!call this%discr%initialize_discret(mesh)
      call this%discr_id%init_discr_fourier(msh=mesh,opt=100)
      call this%discr_d1%init_discr_fourier(msh=mesh,opt=101)
      call this%discr_d2%init_discr_fourier(msh=mesh,opt=102)
      
    end subroutine init_operator_fourier
    
    subroutine t_operator_fourier_d1(this,fi,dfi,dg1y,dg2y,dg1z,dg2z)
      implicit none
      class(t_operator_fourier) :: this
      real(kind=8),dimension(:,:,:),allocatable :: fi,dfi
      real(kind=8),dimension(:,:,:),allocatable,optional :: dg1y,dg2y
      real(kind=8),dimension(:,:,:),allocatable,optional :: dg1z,dg2z
       

      if (this%axis==1) then
         call this%discr_d1%eval_discr_fourier_x(fi,dfi)
      else if (this%axis==2) then
         call this%discr_d1%eval_discr_fourier_y(fi,dfi,dg1y=dg1y,dg2y=dg2y)
      else if (this%axis==3) then
         call this%discr_d1%eval_discr_fourier_z(fi,dfi,dg1y=dg1y,dg2y=dg2y,dg1z=dg1z,dg2z=dg2z)
      end if
      
      
    end subroutine t_operator_fourier_d1
    
    subroutine t_operator_fourier_d2(this,fi,dfi,dg1y,dg2y,dg1z,dg2z)
      implicit none
      class(t_operator_fourier) :: this
      real(kind=8),dimension(:,:,:),allocatable :: fi,dfi
      real(kind=8),dimension(:,:,:),allocatable,optional :: dg1y,dg2y
      real(kind=8),dimension(:,:,:),allocatable,optional :: dg1z,dg2z

      if (this%axis==1) then
         call this%discr_d2%eval_discr_fourier_x(fi,dfi)
      else if (this%axis==2) then
         call this%discr_d2%eval_discr_fourier_y(fi,dfi,dg1y=dg1y,dg2y=dg2y)
      else if (this%axis==3) then
         call this%discr_d2%eval_discr_fourier_z(fi,dfi,dg1y=dg1y,dg2y=dg2y,dg1z=dg1z,dg2z=dg2z)
      end if

      
    end subroutine t_operator_fourier_d2

    
    subroutine t_operator_fourier_id(this,fi,dfi,dg1y,dg2y,dg1z,dg2z)
      implicit none
      class(t_operator_fourier) ::  this
      real(kind=8),dimension(:,:,:),allocatable :: fi,dfi
      real(kind=8),dimension(:,:,:),allocatable,optional :: dg1y,dg2y
      real(kind=8),dimension(:,:,:),allocatable,optional :: dg1z,dg2z
      
    end subroutine t_operator_fourier_id

    subroutine t_operator_fourier_d1_clean(this,fi,dfi,dg1y,dg2y,dg1z,dg2z)
      implicit none
      class(t_operator_fourier) :: this
      real(kind=8),dimension(:,:,:),allocatable :: fi,dfi
      real(kind=8),dimension(:,:,:),allocatable,optional :: dg1y,dg2y
      real(kind=8),dimension(:,:,:),allocatable,optional :: dg1z,dg2z
       

      if (this%axis==1) then
         call this%discr_d1%eval_discr_fourier_x(fi,dfi)
      else if (this%axis==2) then
         call this%discr_d1%eval_discr_fourier_y(fi,dfi,dg1y=dg1y,dg2y=dg2y)
      else if (this%axis==3) then
         call this%discr_d1%eval_discr_fourier_z(fi,dfi,dg1y=dg1y,dg2y=dg2y,dg1z=dg1z,dg2z=dg2z)
      end if
      
      
    end subroutine t_operator_fourier_d1_clean
    
    subroutine t_operator_fourier_d2_clean(this,fi,dfi,dg1y,dg2y,dg1z,dg2z)
      implicit none
      class(t_operator_fourier) :: this
      real(kind=8),dimension(:,:,:),allocatable :: fi,dfi
      real(kind=8),dimension(:,:,:),allocatable,optional :: dg1y,dg2y
      real(kind=8),dimension(:,:,:),allocatable,optional :: dg1z,dg2z

      if (this%axis==1) then
         call this%discr_d2%eval_discr_fourier_x(fi,dfi)
      else if (this%axis==2) then
         call this%discr_d2%eval_discr_fourier_y(fi,dfi,dg1y=dg1y,dg2y=dg2y)
      else if (this%axis==3) then
         call this%discr_d2%eval_discr_fourier_z(fi,dfi,dg1y=dg1y,dg2y=dg2y,dg1z=dg1z,dg2z=dg2z)
      end if

      
    end subroutine t_operator_fourier_d2_clean

    
    subroutine t_operator_fourier_id_clean(this,fi,dfi,dg1y,dg2y,dg1z,dg2z)
      implicit none
      class(t_operator_fourier) ::  this
      real(kind=8),dimension(:,:,:),allocatable :: fi,dfi
      real(kind=8),dimension(:,:,:),allocatable,optional :: dg1y,dg2y
      real(kind=8),dimension(:,:,:),allocatable,optional :: dg1z,dg2z
      
    end subroutine t_operator_fourier_id_clean

    
  end module m_operator_fourier
  
