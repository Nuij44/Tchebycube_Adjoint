module m_snapshots
  use hdf5
  use m_hdf5_ifce
  implicit none
  type t_snapshots
     character(len=1024) :: db_name
     integer :: id
     integer :: dn_export
     
  end type t_snapshots

contains

  SUBROUTINE update_id_and_time(FILENAME,snap_dt,id)
    use mpi
    implicit none
    real(kind=8) :: snap_dt
    CHARACTER(len=*)   :: FILENAME
    integer :: id
    !> VARS...
    INTEGER :: INFO,IERR
    INTEGER(KIND=HID_T) :: P_ID, F_ID, X_ID
    INTEGER, DIMENSION(3) :: IS_GLB,IE_GLB
    INTEGER, DIMENSION(3) :: IS_LOC,IE_LOC
    real(kind=8) :: tc

    
    !> BUILD RESTART FILE =============================================================
    !CALL H5OPEN_F( IERR )
    CALL MPI_INFO_CREATE( INFO, IERR )
    CALL H5PCREATE_F( H5P_FILE_ACCESS_F, P_ID, IERR )
    CALL H5PSET_FAPL_MPIO_F( P_ID, MPI_COMM_WORLD, INFO, IERR )
    !>
    CALL H5FOPEN_F(trim(FILENAME), H5F_ACC_RDWR_F, F_ID, IERR,ACCESS_PRP = P_ID)
    CALL H5PCREATE_F(H5P_DATASET_XFER_F, X_ID, IERR)
    CALL H5PSET_DXPL_MPIO_F(X_ID, H5FD_MPIO_COLLECTIVE_F, IERR)
    
    CALL READ_REAL8_ATTRIBUTE(F_ID,'/','tc',TC)
    tc = tc + snap_dt
    CALL WRITE_REAL8_ATTRIBUTE(F_ID,'/','tc',TC)
    snap_dt = 0
    
    CALL READ_INTEGER_ATTRIBUTE(F_ID,'/','id',id)
    id = id +1
    CALL WRITE_INTEGER_ATTRIBUTE(F_ID,'/','id',id)
    
    CALL H5PCLOSE_F( X_ID, IERR )
    CALL H5FCLOSE_F( F_ID, IERR )
    CALL H5PCLOSE_F( P_ID, IERR )
    
  end SUBROUTINE update_id_and_time


  
  SUBROUTINE snapshot_driver(this,it,tc,x1,x2,x3,u1,u2,u3,p,t)
    implicit none
    class(t_snapshots) :: this
    integer :: it
    real(kind=8) :: tc
    
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE,OPTIONAL :: X1,X2,X3
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE,OPTIONAL :: U1,U2,U3,P,T
    
    if (it==0) then
 !      call EXPORT_GRID(FILENAME,X1,X2,X3)
    end if
    
    if (mod(it,this%dn_export)==0) then
       
!       call EXPORT_GRID(FILENAME,u1,u2,u3,p,t)
    end if
    
  end SUBROUTINE snapshot_driver
  
  

  SUBROUTINE EXPORT_snapshot(FILENAME,u1,u2,u3,P,T)
    use decomp_2d
    use m_mesh_base
    use mpi
    implicit none
    CHARACTER(LEN=*)                          :: FILENAME
    TYPE(T_MESH_BASE)                         :: MSH(3)
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: X1,X2,X3
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: U1,U2,U3,P
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE  :: T

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
    
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/u1', IS_GLB, IE_GLB, 3 )
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/u2', IS_GLB, IE_GLB, 3 )
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/u3', IS_GLB, IE_GLB, 3 )
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/p' , IS_GLB, IE_GLB, 3 )
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/T' , IS_GLB, IE_GLB, 3 )
    
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/u1' ,u1  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/u2' ,u2  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/u3' ,u3  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/p'  ,P  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/T'  ,T  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    
    CALL H5PCLOSE_F( X_ID, IERR )
    CALL H5FCLOSE_F( F_ID, IERR )
    CALL H5PCLOSE_F( P_ID, IERR )
    
  END SUBROUTINE EXPORT_SNAPSHOT


  SUBROUTINE EXPORT_GRID(FILENAME,X1,X2,X3)
    use decomp_2d
    use m_mesh_base
    use mpi
    implicit none
    CHARACTER(len=*)   :: FILENAME
    TYPE(T_MESH_base)     :: msh(3)
    REAL(KIND=8),dimension(:,:,:),allocatable :: x1,x2,x3
    
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

    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/x1', IS_GLB, IE_GLB, 3 )
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/x2', IS_GLB, IE_GLB, 3 )
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/x3', IS_GLB, IE_GLB, 3 )
    
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/x1' ,x1  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/x2' ,x2  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/x3' ,x3  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )

    
    call CREATE_SCALAR_ATTRIBUTE(f_id,'/','tc',H5T_NATIVE_DOUBLE)
    call CREATE_SCALAR_ATTRIBUTE(f_id,'/','id',H5T_NATIVE_INTEGER)
    
!    CALL WRITE_REAL8_ATTRIBUTE(F_ID,'/','tc',TC)
!    CALL WRITE_INTEGER_ATTRIBUTE(F_ID,'/','id',id)

    
    CALL H5PCLOSE_F( X_ID, IERR )
    CALL H5FCLOSE_F( F_ID, IERR )
    CALL H5PCLOSE_F( P_ID, IERR )
    
  END SUBROUTINE EXPORT_GRID


  
  
end module m_snapshots
