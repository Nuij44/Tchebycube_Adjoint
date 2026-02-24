module m_adjoint_tool
  use m_numerics
  use m_mesh_base
  use m_operator_base
  use mpi
  use decomp_2d
  use decomp_2d_mpi
  use m_hdf5_ifce

  INTERFACE SAVE_HDF5
    MODULE PROCEDURE EXPORT_HDF5_VELOCITY_STRAT
    MODULE PROCEDURE EXPORT_HDF5_FIELD
  END INTERFACE SAVE_HDF5

  INTERFACE IMPORT_HDF5
      MODULE PROCEDURE IMPORT_HDF5_BASIC
      MODULE PROCEDURE IMPORT_HDF5_BASIC_STRAT
  END INTERFACE IMPORT_HDF5
  
  INTERFACE IMPORT_HDF5_INIT
     MODULE PROCEDURE IMPORT_HDF5_INIT_BASIC
     MODULE PROCEDURE IMPORT_HDF5_INIT_STRAT
  END INTERFACE IMPORT_HDF5_INIT

  contains

  SUBROUTINE command_line_read_input_dns(ifile,ofile,initfile,save_it)
    use mpi
    implicit none
    character(len=*) :: ifile,ofile,initfile
    logical,optional :: save_it
    integer i
    character(len=4096) :: arg
    integer :: stat,rank
    integer :: ierr
    integer :: ll
    
    
    call mpi_comm_rank(mpi_comm_world,rank,ierr)

    save_it = .FALSE.
    
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
          case ('--output-dir')
             CALL get_command_argument(i, arg) ; i = i+1
             print *,trim(arg)
             ofile=trim(arg)
          case ('--init-data')
             CALL get_command_argument(i, arg) ; i = i+1
             print *,trim(arg)
             initfile=trim(arg)
          case ('--save_it')
             save_it = .TRUE.
             print*,'Save iteration = ',save_it
          case default
          end select
       END DO
       
    end if
    
    ll = len(ifile)
    CALL MPI_BCAST(ifile(1:ll) ,ll,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
    ll = len(ofile)
    CALL MPI_BCAST(ofile(1:ll) ,ll,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
    ll = len(initfile)
    CALL MPI_BCAST(initfile(1:ll) ,ll,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
    CALL MPI_BCAST(save_it ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
    
  end SUBROUTINE command_line_read_input_dns

  SUBROUTINE command_line_read_input_adj(ifile,ofile,initfile,adjoint)
    use mpi
    implicit none
    character(len=*) :: ifile,ofile,initfile
    logical :: adjoint
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
          case ('--init-data')
             CALL get_command_argument(i, arg) ; i = i+1
             print *,trim(arg)
             initfile=trim(arg)
          case ('--input')
             CALL get_command_argument(i, arg) ; i = i+1
             print *,trim(arg)
             ifile=trim(arg)
          case ('--output-dir')
             CALL get_command_argument(i, arg) ; i = i+1
             print *,trim(arg)
             ofile=trim(arg)
          case ('--do-adjoint')
             CALL get_command_argument(i, arg) ; i = i+1
             adjoint=.TRUE.
             print *,adjoint
          case default
          end select
       END DO
       
    end if
    
    ll = len(ifile)
    CALL MPI_BCAST(ifile(1:ll) ,ll,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
    ll = len(ofile)
    CALL MPI_BCAST(ofile(1:ll) ,ll,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
    ll = len(initfile)
    CALL MPI_BCAST(initfile(1:ll) ,ll,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)
    CALL MPI_BCAST(adjoint ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
    
  end SUBROUTINE command_line_read_input_adj


  SUBROUTINE IMPORT_HDF5_INIT_BASIC(FILENAME,u1,u2,u3)
    use decomp_2d
    use m_mesh_base
    use mpi
    implicit none
    CHARACTER(len=*)                          :: FILENAME
    REAL(KIND=8),dimension(:,:,:),allocatable :: u1,u2,u3
    
    !> VARS...
    INTEGER :: INFO,IERR
    INTEGER(KIND=HID_T) :: P_ID, F_ID, X_ID
    INTEGER, DIMENSION(3) :: IS_GLB,IE_GLB
    INTEGER, DIMENSION(3) :: IS_LOC,IE_LOC
    TYPE(DECOMP_INFO) :: ph
    
    CALL GET_DECOMP_INFO(PH)
    
    is_glb(1:3) = [ ph%xst(1) ,ph%yst(2) , ph%zst(3) ] 
    ie_glb(1:3) = [ ph%xen(1) ,ph%yen(2) , ph%zen(3) ] 
    
    is_loc(1:3) = [ ph%xst(1) ,ph%xst(2) , ph%xst(3) ]
    ie_loc(1:3) = [ ph%xen(1) ,ph%xen(2) , ph%xen(3) ] 
    
    
    CALL MPI_INFO_CREATE( INFO, IERR )
    CALL H5PCREATE_F( H5P_FILE_ACCESS_F, P_ID, IERR )
    CALL H5PSET_FAPL_MPIO_F( P_ID, MPI_COMM_WORLD, INFO, IERR )
    
    CALL H5FOPEN_F(FILENAME, H5F_ACC_RDWR_F, F_ID, IERR,ACCESS_PRP = P_ID)
    CALL H5PCREATE_F(H5P_DATASET_XFER_F, X_ID, IERR)
    CALL H5PSET_DXPL_MPIO_F(X_ID, H5FD_MPIO_COLLECTIVE_F, IERR)
    
    
    CALL READ_ARRAY_RANK_3( F_ID, X_ID,'/dump/u1',u1  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL READ_ARRAY_RANK_3( F_ID, X_ID,'/dump/u2',u2  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL READ_ARRAY_RANK_3( F_ID, X_ID,'/dump/u3',u3  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    
    CALL H5PCLOSE_F( X_ID, IERR )
    CALL H5FCLOSE_F( F_ID, IERR )
    CALL H5PCLOSE_F( P_ID, IERR )
    
    
  END SUBROUTINE IMPORT_HDF5_INIT_BASIC

  SUBROUTINE IMPORT_HDF5_INIT_STRAT(FILENAME,u1,u2,u3,T)
    use decomp_2d
    use m_mesh_base
    use mpi
    implicit none
    CHARACTER(len=*)                          :: FILENAME
    REAL(KIND=8),dimension(:,:,:),allocatable :: u1,u2,u3,T
    
    !> VARS...
    INTEGER :: INFO,IERR
    INTEGER(KIND=HID_T) :: P_ID, F_ID, X_ID
    INTEGER, DIMENSION(3) :: IS_GLB,IE_GLB
    INTEGER, DIMENSION(3) :: IS_LOC,IE_LOC
    TYPE(DECOMP_INFO) :: ph
    
    CALL GET_DECOMP_INFO(PH)
    
    is_glb(1:3) = [ ph%xst(1) ,ph%yst(2) , ph%zst(3) ] 
    ie_glb(1:3) = [ ph%xen(1) ,ph%yen(2) , ph%zen(3) ] 
    
    is_loc(1:3) = [ ph%xst(1) ,ph%xst(2) , ph%xst(3) ]
    ie_loc(1:3) = [ ph%xen(1) ,ph%xen(2) , ph%xen(3) ] 
    
    
    CALL MPI_INFO_CREATE( INFO, IERR )
    CALL H5PCREATE_F( H5P_FILE_ACCESS_F, P_ID, IERR )
    CALL H5PSET_FAPL_MPIO_F( P_ID, MPI_COMM_WORLD, INFO, IERR )
    
    CALL H5FOPEN_F(FILENAME, H5F_ACC_RDWR_F, F_ID, IERR,ACCESS_PRP = P_ID)
    CALL H5PCREATE_F(H5P_DATASET_XFER_F, X_ID, IERR)
    CALL H5PSET_DXPL_MPIO_F(X_ID, H5FD_MPIO_COLLECTIVE_F, IERR)
    
    
    CALL READ_ARRAY_RANK_3( F_ID, X_ID,'/dump/u1',u1  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL READ_ARRAY_RANK_3( F_ID, X_ID,'/dump/u2',u2  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL READ_ARRAY_RANK_3( F_ID, X_ID,'/dump/u3',u3  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL READ_ARRAY_RANK_3( F_ID, X_ID,'/dump/T' ,T   , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    
    CALL H5PCLOSE_F( X_ID, IERR )
    CALL H5FCLOSE_F( F_ID, IERR )
    CALL H5PCLOSE_F( P_ID, IERR )
    
    
  END SUBROUTINE IMPORT_HDF5_INIT_STRAT
  
  SUBROUTINE EXPORT_HDF5_VELOCITY_STRAT(FILENAME,MESH,u1,u2,u3,T)
    implicit none
    CHARACTER(len=*)                          :: FILENAME
    REAL(KIND=8),dimension(:,:,:),allocatable :: x1,x2,x3
    REAL(KIND=8),dimension(:,:,:),allocatable :: u1,u2,u3,T
    TYPE(T_MESH_BASE) :: MESH(3)
    
    !> VARS...
    INTEGER :: INFO,IERR
    INTEGER(KIND=HID_T) :: P_ID, F_ID, X_ID
    INTEGER, DIMENSION(3) :: IS_GLB,IE_GLB
    INTEGER, DIMENSION(3) :: IS_LOC,IE_LOC
    TYPE(DECOMP_INFO) :: ph
    
    CALL GET_DECOMP_INFO(PH)
    
    is_glb(1:3) = [ ph%xst(1) ,ph%yst(2) , ph%zst(3) ] 
    ie_glb(1:3) = [ ph%xen(1) ,ph%yen(2) , ph%zen(3) ] 
    
    is_loc(1:3) = [ ph%xst(1) ,ph%xst(2) , ph%xst(3) ]
    ie_loc(1:3) = [ ph%xen(1) ,ph%xen(2) , ph%xen(3) ] 
    
    
    
    !> BUILD RESTART FILE =============================================================
    !CALL H5OPEN_F( IERR )
    CALL MPI_INFO_CREATE( INFO, IERR )
    CALL H5PCREATE_F( H5P_FILE_ACCESS_F, P_ID, IERR )
    CALL H5PSET_FAPL_MPIO_F( P_ID,MPI_COMM_WORLD, INFO, IERR )
    CALL H5FCREATE_F(trim(FILENAME), H5F_ACC_TRUNC_F, F_ID, IERR, ACCESS_PRP = P_ID)
    CALL H5PCREATE_F(H5P_DATASET_XFER_F, X_ID, IERR)
    CALL H5PSET_DXPL_MPIO_F(X_ID, H5FD_MPIO_COLLECTIVE_F, IERR)
    
    !> ARBORECENCES DE MES DATA
    call CREATE_A_GROUP(F_ID,'/grid')
    call CREATE_A_GROUP(F_ID,'/save')
    
    !> X1-GRID ========================================================================================================
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/grid/x1', (/IS_GLB(1)/), (/IE_GLB(1)/), 1 )
    CALL WRITE_ARRAY_RANK_1( F_ID, X_ID,'/grid/x1', mesh(1)%x, IS_GLB(1), IE_GLB(1), IS_LOC(1), IE_LOC(1) )
    
    !> X2-GRID ========================================================================================================
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/grid/x2', (/IS_GLB(2)/), (/IE_GLB(2)/), 1 )
    CALL WRITE_ARRAY_RANK_1( F_ID, X_ID,'/grid/x2', mesh(2)%x, IS_GLB(2), IE_GLB(2), IS_LOC(2), IE_LOC(2) )
    
    !> ZF GRID ========================================================================================================
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/grid/x3', (/IS_GLB(3)/), (/IE_GLB(3)/), 1 )
    CALL WRITE_ARRAY_RANK_1( F_ID, X_ID,'/grid/x3', mesh(3)%x, IS_GLB(3), IE_GLB(3), IS_LOC(3), IE_LOC(3) )
    
    !> the data to store
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/save/u1', IS_GLB, IE_GLB, 3 )
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/save/u2', IS_GLB, IE_GLB, 3 )
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/save/u3', IS_GLB, IE_GLB, 3 )
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/save/T', IS_GLB, IE_GLB, 3 )
    
    
    !store the velocity fields
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/save/u1' ,u1  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/save/u2' ,u2  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/save/u3' ,u3  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/save/T'  ,T   , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    
    
    CALL H5PCLOSE_F( X_ID, IERR )
    CALL H5FCLOSE_F( F_ID, IERR )
    CALL H5PCLOSE_F( P_ID, IERR )
    
    
  END SUBROUTINE EXPORT_HDF5_VELOCITY_STRAT

  SUBROUTINE EXPORT_HDF5_VELOCITY(FILENAME,MESH,u1,u2,u3)
    implicit none
    CHARACTER(len=*)                          :: FILENAME
    REAL(KIND=8),dimension(:,:,:),allocatable :: x1,x2,x3
    REAL(KIND=8),dimension(:,:,:),allocatable :: u1,u2,u3
    TYPE(T_MESH_BASE) :: MESH(3)
    
    !> VARS...
    INTEGER :: INFO,IERR
    INTEGER(KIND=HID_T) :: P_ID, F_ID, X_ID
    INTEGER, DIMENSION(3) :: IS_GLB,IE_GLB
    INTEGER, DIMENSION(3) :: IS_LOC,IE_LOC
    TYPE(DECOMP_INFO) :: ph
    
    CALL GET_DECOMP_INFO(PH)
    
    is_glb(1:3) = [ ph%xst(1) ,ph%yst(2) , ph%zst(3) ] 
    ie_glb(1:3) = [ ph%xen(1) ,ph%yen(2) , ph%zen(3) ] 
    
    is_loc(1:3) = [ ph%xst(1) ,ph%xst(2) , ph%xst(3) ]
    ie_loc(1:3) = [ ph%xen(1) ,ph%xen(2) , ph%xen(3) ] 
    
    
    
    !> BUILD RESTART FILE =============================================================
    !CALL H5OPEN_F( IERR )
    CALL MPI_INFO_CREATE( INFO, IERR )
    CALL H5PCREATE_F( H5P_FILE_ACCESS_F, P_ID, IERR )
    CALL H5PSET_FAPL_MPIO_F( P_ID,MPI_COMM_WORLD, INFO, IERR )
    CALL H5FCREATE_F(trim(FILENAME), H5F_ACC_TRUNC_F, F_ID, IERR, ACCESS_PRP = P_ID)
    CALL H5PCREATE_F(H5P_DATASET_XFER_F, X_ID, IERR)
    CALL H5PSET_DXPL_MPIO_F(X_ID, H5FD_MPIO_COLLECTIVE_F, IERR)
    
    !> ARBORECENCES DE MES DATA
    call CREATE_A_GROUP(F_ID,'/grid')
    call CREATE_A_GROUP(F_ID,'/save')
    
    !> X1-GRID ========================================================================================================
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/grid/x1', (/IS_GLB(1)/), (/IE_GLB(1)/), 1 )
    CALL WRITE_ARRAY_RANK_1( F_ID, X_ID,'/grid/x1', mesh(1)%x, IS_GLB(1), IE_GLB(1), IS_LOC(1), IE_LOC(1) )
    
    !> X2-GRID ========================================================================================================
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/grid/x2', (/IS_GLB(2)/), (/IE_GLB(2)/), 1 )
    CALL WRITE_ARRAY_RANK_1( F_ID, X_ID,'/grid/x2', mesh(2)%x, IS_GLB(2), IE_GLB(2), IS_LOC(2), IE_LOC(2) )
    
    !> ZF GRID ========================================================================================================
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/grid/x3', (/IS_GLB(3)/), (/IE_GLB(3)/), 1 )
    CALL WRITE_ARRAY_RANK_1( F_ID, X_ID,'/grid/x3', mesh(3)%x, IS_GLB(3), IE_GLB(3), IS_LOC(3), IE_LOC(3) )
    
    !> the data to store
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/save/u1', IS_GLB, IE_GLB, 3 )
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/save/u2', IS_GLB, IE_GLB, 3 )
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/save/u3', IS_GLB, IE_GLB, 3 )
        
    !store the velocity fields
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/save/u1' ,u1  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/save/u2' ,u2  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/save/u3' ,u3  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )    
    
    CALL H5PCLOSE_F( X_ID, IERR )
    CALL H5FCLOSE_F( F_ID, IERR )
    CALL H5PCLOSE_F( P_ID, IERR )
    
    
  END SUBROUTINE EXPORT_HDF5_VELOCITY


  
  SUBROUTINE IMPORT_HDF5_BASIC(FILENAME,msh,u1,u2,u3)
    use decomp_2d
    use m_mesh_base
    use mpi
    implicit none
    CHARACTER(len=*)                          :: FILENAME
    TYPE(T_MESH_base)                         :: msh(3)
    REAL(KIND=8),dimension(:,:,:),allocatable :: u1,u2,u3
    
    !> VARS...
    INTEGER :: INFO,IERR
    INTEGER(KIND=HID_T) :: P_ID, F_ID, X_ID
    INTEGER, DIMENSION(3) :: IS_GLB,IE_GLB
    INTEGER, DIMENSION(3) :: IS_LOC,IE_LOC
    TYPE(DECOMP_INFO) :: ph
    
    CALL GET_DECOMP_INFO(PH)
    
    is_glb(1:3) = [ ph%xst(1) ,ph%yst(2) , ph%zst(3) ] 
    ie_glb(1:3) = [ ph%xen(1) ,ph%yen(2) , ph%zen(3) ] 
    
    is_loc(1:3) = [ ph%xst(1) ,ph%xst(2) , ph%xst(3) ]
    ie_loc(1:3) = [ ph%xen(1) ,ph%xen(2) , ph%xen(3) ] 
    
    
    CALL MPI_INFO_CREATE( INFO, IERR )
    CALL H5PCREATE_F( H5P_FILE_ACCESS_F, P_ID, IERR )
    CALL H5PSET_FAPL_MPIO_F( P_ID, MPI_COMM_WORLD, INFO, IERR )
    
    CALL H5FOPEN_F(FILENAME, H5F_ACC_RDWR_F, F_ID, IERR,ACCESS_PRP = P_ID)
    CALL H5PCREATE_F(H5P_DATASET_XFER_F, X_ID, IERR)
    CALL H5PSET_DXPL_MPIO_F(X_ID, H5FD_MPIO_COLLECTIVE_F, IERR)
    
    
    CALL READ_ARRAY_RANK_3( F_ID, X_ID,'/save/u1',u1  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL READ_ARRAY_RANK_3( F_ID, X_ID,'/save/u2',u2  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL READ_ARRAY_RANK_3( F_ID, X_ID,'/save/u3',u3  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    
    CALL H5PCLOSE_F( X_ID, IERR )
    CALL H5FCLOSE_F( F_ID, IERR )
    CALL H5PCLOSE_F( P_ID, IERR )
    
    
  END SUBROUTINE IMPORT_HDF5_BASIC
  
  SUBROUTINE IMPORT_HDF5_BASIC_STRAT(FILENAME,msh,u1,u2,u3,T)
    use decomp_2d
    use m_mesh_base
    use mpi
    implicit none
    CHARACTER(len=*)                          :: FILENAME
    TYPE(T_MESH_base)                         :: msh(3)
    REAL(KIND=8),dimension(:,:,:),allocatable :: u1,u2,u3,T
    
    !> VARS...
    INTEGER :: INFO,IERR
    INTEGER(KIND=HID_T) :: P_ID, F_ID, X_ID
    INTEGER, DIMENSION(3) :: IS_GLB,IE_GLB
    INTEGER, DIMENSION(3) :: IS_LOC,IE_LOC
    TYPE(DECOMP_INFO) :: ph
    
    CALL GET_DECOMP_INFO(PH)
    
    is_glb(1:3) = [ ph%xst(1) ,ph%yst(2) , ph%zst(3) ] 
    ie_glb(1:3) = [ ph%xen(1) ,ph%yen(2) , ph%zen(3) ] 
    
    is_loc(1:3) = [ ph%xst(1) ,ph%xst(2) , ph%xst(3) ]
    ie_loc(1:3) = [ ph%xen(1) ,ph%xen(2) , ph%xen(3) ] 
    
    
    CALL MPI_INFO_CREATE( INFO, IERR )
    CALL H5PCREATE_F( H5P_FILE_ACCESS_F, P_ID, IERR )
    CALL H5PSET_FAPL_MPIO_F( P_ID, MPI_COMM_WORLD, INFO, IERR )
    
    CALL H5FOPEN_F(FILENAME, H5F_ACC_RDWR_F, F_ID, IERR,ACCESS_PRP = P_ID)
    CALL H5PCREATE_F(H5P_DATASET_XFER_F, X_ID, IERR)
    CALL H5PSET_DXPL_MPIO_F(X_ID, H5FD_MPIO_COLLECTIVE_F, IERR)
    
    
    CALL READ_ARRAY_RANK_3( F_ID, X_ID,'/save/u1',u1  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL READ_ARRAY_RANK_3( F_ID, X_ID,'/save/u2',u2  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL READ_ARRAY_RANK_3( F_ID, X_ID,'/save/u3',u3  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL READ_ARRAY_RANK_3( F_ID, X_ID,'/save/T' ,T   , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    
    CALL H5PCLOSE_F( X_ID, IERR )
    CALL H5FCLOSE_F( F_ID, IERR )
    CALL H5PCLOSE_F( P_ID, IERR )
    
    
  END SUBROUTINE IMPORT_HDF5_BASIC_STRAT
  
  
  subroutine COMPUTE_ADJOINT_NON_LINEAR_TERMS(&
    op_x1,op_x2,op_x3, u, v, w, u_adj, v_adj, w_adj, hu, hv, hw,&
    dg1, dg2, dg3, dg4, dg5, dg6, dg7, dg8, dg9)

    use m_navier_stokes_cart
    implicit none
    
    class(t_operator_base)     :: op_x1, op_x2, op_x3
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: U,HU,U_ADJ
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: V,HV,V_ADJ
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: W,HW,W_ADJ
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: DG1,DG2,DG3
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: DG4,DG5,DG6
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: DG7,DG8,DG9
    
 
!D'après le papier de Eaves et Caulfield
    CALL CROSS(OP_X1, OP_X2, OP_X3, &
         U_ADJ, V_ADJ, W_ADJ, &
         U, V, W, &
         DG1,DG2,DG3)

    CALL CURL(OP_X1, OP_X2, OP_X3, &
         DG1, DG2, DG3, &
         DG7, DG8, DG9, &
         DG4, DG5, DG6)

    CALL CURL(OP_X1, OP_X2, OP_X3, &
         U, V, W, &
         DG1, DG2, DG3, &
         DG4, DG5, DG6)
    
    CALL CROSS(OP_X1, OP_X2, OP_X3, &
         U_ADJ, V_ADJ, W_ADJ, &
         DG1, DG2, DG3, &
         DG4,DG5,DG6)

    HU = DG7 - DG4
    HV = DG8 - DG5
    HW = DG9 - DG6


  END subroutine COMPUTE_ADJOINT_NON_LINEAR_TERMS
  
  subroutine COMPUTE_STRAT_ADJOINT_NON_LINEAR_TERMS(&
    op_x1,op_x2,op_x3, u, v, w, t_adj, nlt,dg01, dg02, dg03)
    
    implicit none
    
    class(t_operator_base)     :: op_x1, op_x2, op_x3
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: U,V,W,T_ADJ,NLT
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: DG01,DG02,DG03

    CALL OP_X1%D1(T_ADJ,DG01)
    CALL OP_X2%D1(T_ADJ,DG02)
    CALL OP_X3%D1(T_ADJ,DG03)

    NLT = U*DG01 + V*DG02 + W*DG03

  END SUBROUTINE COMPUTE_STRAT_ADJOINT_NON_LINEAR_TERMS


  subroutine COMPUTE_ADJOINT_NON_LINEAR_TERMS_CYL(&
              A,Z,R,OP_A,OP_Z,OP_R, UA, UZ, UR, UA_ADJ, UZ_ADJ, UR_ADJ,HA, HZ, HR,&
       DG1, DG2, DG3, DG4, DG5, DG6)

    use m_navier_stokes_cyl
    IMPLICIT NONE
    
    CLASS(T_OPERATOR_BASE)     :: OP_A, OP_Z, OP_R
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: A,UA,UA_ADJ,HA
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: Z,UZ,UZ_ADJ,HZ
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: R,UR,UR_ADJ,HR
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: DG1,DG2,DG3
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: DG4,DG5,DG6
 



    CALL CROSS(OP_R, OP_A, OP_Z, R, &
         UR_ADJ, UA_ADJ, UZ_ADJ, &
         UR, UA, UZ, &
         DG1,DG2,DG3)

    CALL CURL(OP_R, OP_A, OP_Z, R, &
         DG1, DG2, DG3, &
         HR, HA, HZ, &
         DG4, DG5, DG6)

    CALL CURL(OP_R, OP_A, OP_Z, R, &
         UR, UA, UZ, &
         DG1, DG2, DG3, &
         DG4, DG5, DG6)
    
    CALL CROSS(OP_R, OP_A, OP_Z, R, &
         UR_ADJ, UA_ADJ, UZ_ADJ, &
         DG1, DG2, DG3, &
         DG4,DG5,DG6)

    HR = HR - DG4
    HA = HA - DG5
    HZ = HZ - DG6
    

  END subroutine COMPUTE_ADJOINT_NON_LINEAR_TERMS_CYL

    

end module m_adjoint_tool
