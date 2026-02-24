module m_operator_base
  use m_mesh_base
  use decomp_2d
  use m_numerics
  implicit none
  
  type , abstract :: t_operator_base
     integer :: axis
     integer :: N
   contains
     
     !procedure (interface_init), deferred, pass(this) :: initialize_op
     ! fonctionnalité pour l'utilisateur
     procedure(eval_operator), deferred, pass(this) ::  d1
     procedure(eval_operator), deferred, pass(this) ::  d1_clean
     procedure(eval_operator), deferred, pass(this) ::  d2
     procedure(eval_operator), deferred, pass(this) ::  d2_clean
     procedure(eval_operator), deferred, pass(this) ::  Id
     procedure(eval_operator), deferred, pass(this) ::  Id_clean
     
  end type t_operator_base
  
  abstract interface
     
     subroutine interface_init (this, mesh, axis)
       use m_numerics
       import :: t_operator_base
       import :: t_mesh_base
       class (t_operator_base)   ::  this
       type(t_mesh_base) :: mesh
       integer :: axis
     end subroutine interface_init
     
     SUBROUTINE eval_operator(this,fi,dfi,dg1y,dg2y,dg1z,dg2z)
       use m_numerics
       IMPORT :: t_operator_base
       CLASS(t_operator_base) :: this
       real(dp),dimension(:,:,:),allocatable :: fi,dfi
       real(dp),dimension(:,:,:),allocatable,optional :: dg1y,dg2y
       real(dp),dimension(:,:,:),allocatable,optional :: dg1z,dg2z
     END SUBROUTINE EVAL_OPERATOR
     
  end interface
  
  TYPE(DECOMP_INFO),save :: dd

  integer,save :: xst_(3),xen_(3)
  integer,save :: yst_(3),yen_(3)
  integer,save :: zst_(3),zen_(3)
  
contains

  SUBROUTINE START_OP_BASE ()
    IMPLICIT NONE       
    CALL GET_DECOMP_INFO(DD)

    XST_ = DD%XST
    XEN_ = DD%XEN

    YST_ = DD%YST
    YEN_ = DD%YEN

    ZST_ = DD%ZST
    ZEN_ = DD%ZEN
    
!    PRINT*,DD%XST,DD%XEN
    
  END SUBROUTINE START_OP_BASE
  
  integer function GetDim(this)
    implicit none
    class(t_operator_base) :: this
    GetDim = this%N
  end function GetDim
  
  integer function GetAxis(this)
    implicit none
    class(t_operator_base) :: this
    GetAxis = this%Axis
  end function GetAxis
  
  
end module m_operator_base
