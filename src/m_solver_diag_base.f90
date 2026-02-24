

module m_solver_diag_base
  use decomp_2d
  use decomp_2d_fft
  use decomp_2d_mpi
  use m_boundary_conditions
  implicit none
  
  type , abstract :: t_solver_diag_base
     integer :: N(3)
     type(t_boundary_conditions) ::  bcs_x,bcs_y,bcs_z
     real(kind=8) :: sigma,nu(3)
   contains
     
  end type t_solver_diag_base
  

  
contains

  
end module m_solver_diag_base


