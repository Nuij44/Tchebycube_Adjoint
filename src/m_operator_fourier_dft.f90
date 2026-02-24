module m_operator_fourier_DFT
    use m_operator_base
    use m_discr_fourier_DFT
    use m_mesh_base
    
    implicit none
    
    type, extends(t_operator_base)    ::  t_operator_fourier_DFT
       type(t_discr_fourier_DFT) :: discr_d1 
       type(t_discr_fourier_DFT) :: discr_d2 
       type(t_discr_fourier_DFT) :: discr_id 
     contains 
       procedure :: init_operator_fourier_DFT => init_operator_fourier_DFT
       procedure :: initialize => init_operator_fourier_DFT
       procedure :: d1 => t_operator_fourier_DFT_d1
       procedure :: d2 => t_operator_fourier_DFT_d2
       procedure :: id => t_operator_fourier_DFT_id
       procedure :: d1_clean => t_operator_fourier_DFT_d1_clean
       procedure :: d2_clean => t_operator_fourier_DFT_d2_clean
       procedure :: id_clean => t_operator_fourier_DFT_id_clean
    end type t_operator_fourier_DFT

    
    
  contains
    
    subroutine check_operator_fourier_DFT(this,mesh)
      use decomp_2d
      implicit none
      class(t_operator_fourier_DFT) :: this
      type(t_mesh_base)        :: mesh(3)
      
      real(dp),dimension(:,:,:),allocatable :: fi,dfi
      real(dp),dimension(:,:,:),allocatable :: dg1y,dg2y
      real(dp),dimension(:,:,:),allocatable :: dg1z,dg2z

      real(dp) :: x,y,z

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
      

      deallocate(fi,dfi,dg2y,dg1z,dg2z)
    end subroutine check_operator_fourier_DFT

    
    subroutine init_operator_fourier_DFT(this, mesh, axis)
      implicit none
      class(t_operator_fourier_DFT) :: this
      type(t_mesh_base)        :: mesh
      integer                  :: axis
      integer, dimension(3)    :: NDIM
      
      integer :: N
      real(kind=8) :: xmin,xmax
      
      N = mesh%N
      !this%n = mesh%N
      this%axis = axis
      xmin = mesh%x(1)
      xmax = (mesh%x(2)- mesh%x(1))*(n+1) 
      
      call this%discr_id%init_discr_fourier_DFT(xmin=xmin,xmax=xmax,n=n,opt=100)
      call this%discr_d1%init_discr_fourier_DFT(xmin=xmin,xmax=xmax,n=n,opt=101)
      call this%discr_d2%init_discr_fourier_DFT(xmin=xmin,xmax=xmax,n=n,opt=102)
      
      
    end subroutine init_operator_fourier_DFT
    
    subroutine t_operator_fourier_DFT_d1(this,fi,dfi,dg1y,dg2y,dg1z,dg2z)
      implicit none
      class(t_operator_fourier_DFT) :: this
      real(dp),dimension(:,:,:),allocatable :: fi,dfi
      real(dp),dimension(:,:,:),allocatable,optional :: dg1y,dg2y
      real(dp),dimension(:,:,:),allocatable,optional :: dg1z,dg2z


      if (this%axis==1) then
         call this%discr_d1%eval_discr_fourier_DFT_x(fi,dfi)
      else if (this%axis==2) then
         call this%discr_d1%eval_discr_fourier_DFT_y(fi,dfi)
      else if (this%axis==3) then
         call this%discr_d1%eval_discr_fourier_DFT_z(fi,dfi)
      end if
      
      
    end subroutine t_operator_fourier_DFT_d1

    subroutine t_operator_fourier_DFT_d1_clean(this,fi,dfi,dg1y,dg2y,dg1z,dg2z)
      implicit none
      class(t_operator_fourier_DFT) :: this
      real(dp),dimension(:,:,:),allocatable :: fi,dfi
      real(dp),dimension(:,:,:),allocatable,optional :: dg1y,dg2y
      real(dp),dimension(:,:,:),allocatable,optional :: dg1z,dg2z


      if (this%axis==1) then
         call this%discr_d1%eval_discr_fourier_DFT_x_clean(fi,dfi)
      else if (this%axis==2) then
         call this%discr_d1%eval_discr_fourier_DFT_y_clean(fi,dfi)
      else if (this%axis==3) then
         call this%discr_d1%eval_discr_fourier_DFT_z_clean(fi,dfi)
      end if
      
      
    end subroutine t_operator_fourier_DFT_d1_clean

    
    subroutine t_operator_fourier_DFT_d2(this,fi,dfi,dg1y,dg2y,dg1z,dg2z)
      implicit none
      class(t_operator_fourier_DFT) :: this
      real(dp),dimension(:,:,:),allocatable :: fi,dfi
      real(dp),dimension(:,:,:),allocatable,optional :: dg1y,dg2y
      real(dp),dimension(:,:,:),allocatable,optional :: dg1z,dg2z

      if (this%axis==1) then
         call this%discr_d2%eval_discr_fourier_DFT_x(fi,dfi)
      else if (this%axis==2) then
         call this%discr_d2%eval_discr_fourier_DFT_y(fi,dfi)
      else if (this%axis==3) then
         call this%discr_d2%eval_discr_fourier_DFT_z(fi,dfi)
      end if

      
    end subroutine t_operator_fourier_DFT_d2

    subroutine t_operator_fourier_DFT_d2_clean(this,fi,dfi,dg1y,dg2y,dg1z,dg2z)
      implicit none
      class(t_operator_fourier_DFT) :: this
      real(dp),dimension(:,:,:),allocatable :: fi,dfi
      real(dp),dimension(:,:,:),allocatable,optional :: dg1y,dg2y
      real(dp),dimension(:,:,:),allocatable,optional :: dg1z,dg2z

      if (this%axis==1) then
         call this%discr_d2%eval_discr_fourier_DFT_x_clean(fi,dfi)
      else if (this%axis==2) then
         call this%discr_d2%eval_discr_fourier_DFT_y_clean(fi,dfi)
      else if (this%axis==3) then
         call this%discr_d2%eval_discr_fourier_DFT_z_clean(fi,dfi)
      end if

      
    end subroutine t_operator_fourier_DFT_d2_clean

    
    subroutine t_operator_fourier_DFT_id(this,fi,dfi,dg1y,dg2y,dg1z,dg2z)
      implicit none
      class(t_operator_fourier_DFT) ::  this
      real(dp),dimension(:,:,:),allocatable :: fi,dfi
      real(dp),dimension(:,:,:),allocatable,optional :: dg1y,dg2y
      real(dp),dimension(:,:,:),allocatable,optional :: dg1z,dg2z
      dfi = fi
    end subroutine t_operator_fourier_DFT_id

    subroutine t_operator_fourier_DFT_id_clean(this,fi,dfi,dg1y,dg2y,dg1z,dg2z)
      implicit none
      class(t_operator_fourier_DFT) ::  this
      real(dp),dimension(:,:,:),allocatable :: fi,dfi
      real(dp),dimension(:,:,:),allocatable,optional :: dg1y,dg2y
      real(dp),dimension(:,:,:),allocatable,optional :: dg1z,dg2z

      if (this%axis==1) then
         call this%discr_id%eval_discr_fourier_DFT_x_clean(fi,dfi)
      else if (this%axis==2) then
         call this%discr_id%eval_discr_fourier_DFT_y_clean(fi,dfi)
      else if (this%axis==3) then
         call this%discr_id%eval_discr_fourier_DFT_z_clean(fi,dfi)
      end if


    end subroutine t_operator_fourier_DFT_id_clean

    
  end module m_operator_fourier_DFT
  
