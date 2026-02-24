module m_hdf5_ifce
  use hdf5

  implicit none
  
  REAL(KIND=8),dimension(:,:,:),allocatable,public :: H5_BUFF
  REAL(KIND=8),dimension(:,:,:),allocatable,public :: H5_BUFF_MEM
  
  REAL(KIND=8),dimension(:,:,:),allocatable,public :: H5_BUFF_FCC
  REAL(KIND=8),dimension(:,:,:),allocatable,public :: H5_BUFF_CFC
  REAL(KIND=8),dimension(:,:,:),allocatable,public :: H5_BUFF_CCF
  REAL(KIND=8),dimension(:,:,:),allocatable,public :: H5_BUFF_CCC
  
  REAL(KIND=8),dimension(:),allocatable,public :: H5_X1_C,H5_X1_F
  REAL(KIND=8),dimension(:),allocatable,public :: H5_X2_C,H5_X2_F
  REAL(KIND=8),dimension(:),allocatable,public :: H5_X3_C,H5_X3_F
  
contains


  
  SUBROUTINE EXPORT_HDF5_FIELD_I(FILENAME,X1,X2,X3,u1,u2)
    use decomp_2d
    use m_mesh_base
    use mpi
    implicit none
    CHARACTER(len=*)                          :: FILENAME
    REAL(KIND=8),dimension(:,:,:),allocatable :: x1,x2,x3
    REAL(KIND=8),dimension(:,:,:),allocatable :: u1,u2
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
    
    !> the data to store
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/x1', IS_GLB, IE_GLB, 3 )
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/x2', IS_GLB, IE_GLB, 3 )
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/x3', IS_GLB, IE_GLB, 3 )

    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/x1' ,x1  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/x2' ,x2  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/x3' ,x3  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/u1', IS_GLB, IE_GLB, 3 )
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/u2', IS_GLB, IE_GLB, 3 )
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/u1' ,u1  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/u2' ,u2  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    
    
    CALL H5PCLOSE_F( X_ID, IERR )
    CALL H5FCLOSE_F( F_ID, IERR )
    CALL H5PCLOSE_F( P_ID, IERR )
    
  END SUBROUTINE EXPORT_HDF5_FIELD_I



  
  SUBROUTINE EXPORT_HDF5_FIELD(FILENAME,X1,X2,X3,u1,u2,u3,P,T,C)
    use decomp_2d
    use m_mesh_base
    use mpi
    implicit none
    CHARACTER(len=*)                          :: FILENAME
    REAL(KIND=8),dimension(:,:,:),allocatable :: x1,x2,x3
    REAL(KIND=8),dimension(:,:,:),allocatable :: u1,u2,u3,p
    REAL(KIND=8),dimension(:,:,:),allocatable, optional :: T,C

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
    
    !> the data to store
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/x1', IS_GLB, IE_GLB, 3 )
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/x2', IS_GLB, IE_GLB, 3 )
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/x3', IS_GLB, IE_GLB, 3 )

    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/x1' ,x1  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/x2' ,x2  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/x3' ,x3  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/u1', IS_GLB, IE_GLB, 3 )
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/u2', IS_GLB, IE_GLB, 3 )
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/u3', IS_GLB, IE_GLB, 3 )
    CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/p' , IS_GLB, IE_GLB, 3 )
    
    
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/u1' ,u1  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/u2' ,u2  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/u3' ,u3  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/p'  ,P  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    
    IF (PRESENT(T)) THEN
       CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/T' , IS_GLB, IE_GLB, 3 )
       CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/T'  ,T  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
    END IF
    
    
    CALL H5PCLOSE_F( X_ID, IERR )
    CALL H5FCLOSE_F( F_ID, IERR )
    CALL H5PCLOSE_F( P_ID, IERR )
    
  END SUBROUTINE EXPORT_HDF5_FIELD


  
  SUBROUTINE DUMP_HDF5_BASIC(FILENAME,MODE,TC,DT,msh,u1,u2,u3,P,T,C)
    use decomp_2d
    use m_mesh_base
    use mpi
    implicit none
    CHARACTER(len=*)                          :: FILENAME
    CHARACTER(len=*)                          :: MODE
    TYPE(T_MESH_base)                         :: msh(3)
    REAL(KIND=8)                              :: TC,DT
    REAL(KIND=8),dimension(:,:,:),allocatable :: u1,u2,u3,p
    REAL(KIND=8),dimension(:,:,:),allocatable, optional :: T,C

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
    
    
    IF (TRIM(MODE)=='NEW') THEN
       
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
       call CREATE_A_GROUP(F_ID,'/dump')
       call CREATE_SCALAR_ATTRIBUTE(f_id,'/dump','tc',H5T_NATIVE_DOUBLE)
       call CREATE_SCALAR_ATTRIBUTE(f_id,'/dump','dt',H5T_NATIVE_DOUBLE)

       !> X1-GRID ========================================================================================================
       CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/grid/x1', (/IS_GLB(1)/), (/IE_GLB(1)/), 1 )
       CALL WRITE_ARRAY_RANK_1( F_ID, X_ID,'/grid/x1', msh(1)%x, IS_GLB(1), IE_GLB(1), IS_LOC(1), IE_LOC(1) )

       !> X2-GRID ========================================================================================================
       CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/grid/x2', (/IS_GLB(2)/), (/IE_GLB(2)/), 1 )
       CALL WRITE_ARRAY_RANK_1( F_ID, X_ID,'/grid/x2', msh(2)%x, IS_GLB(2), IE_GLB(2), IS_LOC(2), IE_LOC(2) )

       !> ZF GRID ========================================================================================================
       CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/grid/x3', (/IS_GLB(3)/), (/IE_GLB(3)/), 1 )
       CALL WRITE_ARRAY_RANK_1( F_ID, X_ID,'/grid/x3', msh(3)%x, IS_GLB(3), IE_GLB(3), IS_LOC(3), IE_LOC(3) )

       !> the data to store
       CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/dump/u1', IS_GLB, IE_GLB, 3 )
       CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/dump/u2', IS_GLB, IE_GLB, 3 )
       CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/dump/u3', IS_GLB, IE_GLB, 3 )
       CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/dump/p' , IS_GLB, IE_GLB, 3 )

       IF (PRESENT(T)) THEN
          CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/dump/T' , IS_GLB, IE_GLB, 3 )
       END IF
       IF (PRESENT(C)) THEN
          CALL CREATE_ARRAY_RANK_NDIMS( F_ID, '/dump/C' , IS_GLB, IE_GLB, 3 )
       END IF

       
       CALL H5PCLOSE_F( X_ID, IERR )
       CALL H5FCLOSE_F( F_ID, IERR )
       CALL H5PCLOSE_F( P_ID, IERR )

    END IF

    IF (TRIM(MODE)=='WRITE') THEN
       !> WRITE A RESTART FILE ============================================================

       CALL MPI_INFO_CREATE( INFO, IERR )
       CALL H5PCREATE_F( H5P_FILE_ACCESS_F, P_ID, IERR )
       CALL H5PSET_FAPL_MPIO_F( P_ID, MPI_COMM_WORLD, INFO, IERR )
       !>
       CALL H5FOPEN_F(FILENAME, H5F_ACC_RDWR_F, F_ID, IERR,ACCESS_PRP = P_ID)
       CALL H5PCREATE_F(H5P_DATASET_XFER_F, X_ID, IERR)
       CALL H5PSET_DXPL_MPIO_F(X_ID, H5FD_MPIO_COLLECTIVE_F, IERR)

       CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/dump/u1' ,u1  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
       CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/dump/u2' ,u2  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
       CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/dump/u3' ,u3  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
       CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/dump/p'  ,P  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
       CALL WRITE_REAL8_ATTRIBUTE(F_ID,'/dump','tc',TC)
       CALL WRITE_REAL8_ATTRIBUTE(F_ID,'/dump','dt',dt)

       IF (PRESENT(T)) THEN
          CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/dump/T'  ,T  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
       END IF
       IF (PRESENT(C)) THEN
          CALL WRITE_ARRAY_RANK_3( F_ID, X_ID,'/dump/C'  ,C  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
       END IF

       CALL H5PCLOSE_F( X_ID, IERR )
       CALL H5FCLOSE_F( F_ID, IERR )
       CALL H5PCLOSE_F( P_ID, IERR )
    END IF

    !> READ A RESTART FILE ============================================================

    IF (TRIM(MODE)=='READ') THEN
       !>
       CALL MPI_INFO_CREATE( INFO, IERR )
       CALL H5PCREATE_F( H5P_FILE_ACCESS_F, P_ID, IERR )
       CALL H5PSET_FAPL_MPIO_F( P_ID, MPI_COMM_WORLD, INFO, IERR )
       !>
       CALL H5FOPEN_F(FILENAME, H5F_ACC_RDWR_F, F_ID, IERR,ACCESS_PRP = P_ID)
       CALL H5PCREATE_F(H5P_DATASET_XFER_F, X_ID, IERR)
       CALL H5PSET_DXPL_MPIO_F(X_ID, H5FD_MPIO_COLLECTIVE_F, IERR)

       
       CALL READ_ARRAY_RANK_3( F_ID, X_ID,'/dump/u1',u1  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
       CALL READ_ARRAY_RANK_3( F_ID, X_ID,'/dump/u2',u2  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
       CALL READ_ARRAY_RANK_3( F_ID, X_ID,'/dump/u3',u3  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
       CALL READ_ARRAY_RANK_3( F_ID, X_ID,'/dump/p' ,p  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )

       IF (PRESENT(T)) THEN
          CALL READ_ARRAY_RANK_3( F_ID, X_ID,'/dump/T' ,T  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
       END IF
       IF (PRESENT(C)) THEN
          CALL READ_ARRAY_RANK_3( F_ID, X_ID,'/dump/C' ,C  , IS_GLB, IE_GLB, IS_LOC, IE_LOC )
       END IF
       
       CALL READ_REAL8_ATTRIBUTE(F_ID,'/dump','tc',TC)
       CALL READ_REAL8_ATTRIBUTE(F_ID,'/dump','dt',dt)
       CALL H5PCLOSE_F( X_ID, IERR )
       CALL H5FCLOSE_F( F_ID, IERR )
       CALL H5PCLOSE_F( P_ID, IERR )
    END IF
    
  END SUBROUTINE DUMP_HDF5_BASIC



  
  SUBROUTINE CREATE_ARRAY_RANK_NDIMS(GROUP_ID,DATASET_NAME,IS_GLB,IE_GLB,NDIMS)
    IMPLICIT NONE
    INTEGER(KIND=HID_T)       :: GROUP_ID
    character(len=*)          :: DATASET_NAME
    INTEGER                   :: NDIMS
    INTEGER, dimension(NDIMS) :: IS_GLB
    INTEGER, dimension(NDIMS) :: IE_GLB
    INTEGER(KIND=HSIZE_T) :: COUNT(NDIMS)                                          ! GLOBAL HYPER SLAB INFO
    INTEGER(KIND=HID_T)   :: GLOBAL_DATASPACE_ID, DATASET_ID
    INTEGER               :: IERR
    
    !> CREATE GLOBAL DATA SET
    COUNT = IE_GLB-IS_GLB+1
    CALL H5SCREATE_SIMPLE_F ( NDIMS, COUNT, GLOBAL_DATASPACE_ID, IERR )
    ! CREATE THE DATASPACE 
    CALL H5DCREATE_F( GROUP_ID, trim(DATASET_NAME), H5T_NATIVE_DOUBLE, GLOBAL_DATASPACE_ID, DATASET_ID, IERR )
    
    CALL H5DCLOSE_F( DATASET_ID, IERR )
    CALL H5SCLOSE_F( GLOBAL_DATASPACE_ID, IERR )
    
  end subroutine CREATE_ARRAY_RANK_NDIMS
  

  SUBROUTINE WRITE_ARRAY_RANK_3( F_ID, X_ID, DSET_NAME, FI, IS_GLB, IE_GLB, IS_LOC, IE_LOC  )
    IMPLICIT NONE
    INTEGER(KIND=HID_T)                         :: F_ID
    INTEGER(KIND=HID_T)                         :: X_ID
    CHARACTER(LEN=*)                            :: DSET_NAME
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: FI
    INTEGER, DIMENSION(3)                       :: IS_GLB
    INTEGER, DIMENSION(3)                       :: IE_GLB
    INTEGER, DIMENSION(3)                       :: IS_LOC
    INTEGER, DIMENSION(3)                       :: IE_LOC
    
    
    INTEGER, PARAMETER                      :: NDIMS=3
    INTEGER(KIND=HSIZE_T), DIMENSION(NDIMS) :: DATA_DIMS
    INTEGER(KIND=HSIZE_T), DIMENSION(NDIMS) :: OFFSET 
    INTEGER(KIND=HID_T)                     :: DSET_ID
    INTEGER(KIND=HID_T)                     :: LOCAL_DATASPACE_ID
    INTEGER(KIND=HID_T)                     :: GLOBAL_DATASPACE_ID
    INTEGER                                 :: IERR
    
    !> LOCAL 
    DATA_DIMS = IE_LOC - IS_LOC + 1 !> dimensions locales
    OFFSET    = 0                   !> pas d'offset sur les données locales
    CALL H5SCREATE_SIMPLE_F( NDIMS, DATA_DIMS, LOCAL_DATASPACE_ID, IERR )
    
    CALL H5SSELECT_HYPERSLAB_F(         &
         SPACE_ID = LOCAL_DATASPACE_ID, &
         OPERATOR = H5S_SELECT_SET_F,   &
         START    = OFFSET,             &
         COUNT    = DATA_DIMS,          &
         HDFERR   = IERR                )
    
    !> file or globale space
    DATA_DIMS = IE_GLB - IS_GLB + 1
    CALL H5SCREATE_SIMPLE_F ( NDIMS, DATA_DIMS, GLOBAL_DATASPACE_ID, IERR )
    
    DATA_DIMS = IE_LOC - IS_LOC + 1
    OFFSET = IS_LOC - IS_GLB
    CALL H5SSELECT_HYPERSLAB_F(           &
         SPACE_ID = GLOBAL_DATASPACE_ID,  &
         OPERATOR = H5S_SELECT_SET_F,     &
         START    = OFFSET,               &
         COUNT    = DATA_DIMS,            &
         HDFERR   = IERR                  )
    
    CALL H5DOPEN_F(F_ID, DSET_NAME, DSET_ID, IERR)
    DATA_DIMS = IE_LOC-IS_LOC+1
    CALL H5DWRITE_F(DSET_ID, H5T_NATIVE_DOUBLE,&
         FI(IS_LOC(1):IE_LOC(1),IS_LOC(2):IE_LOC(2),IS_LOC(3):IE_LOC(3)), DATA_DIMS, IERR,&
         FILE_SPACE_ID=GLOBAL_DATASPACE_ID, MEM_SPACE_ID=LOCAL_DATASPACE_ID, XFER_PRP=X_ID )
    
    CALL H5DCLOSE_F( DSET_ID, IERR )
    
  end subroutine WRITE_ARRAY_RANK_3
  
  SUBROUTINE READ_ARRAY_RANK_3( F_ID, X_ID, DSET_NAME, FI, IS_GLB, IE_GLB, IS_LOC, IE_LOC  )
    IMPLICIT NONE
    INTEGER(KIND=HID_T)                         :: F_ID
    INTEGER(KIND=HID_T)                         :: X_ID
    CHARACTER(LEN=*)                            :: DSET_NAME
    REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: FI
    INTEGER, DIMENSION(3)                       :: IS_GLB
    INTEGER, DIMENSION(3)                       :: IE_GLB
    INTEGER, DIMENSION(3)                       :: IS_LOC
    INTEGER, DIMENSION(3)                       :: IE_LOC
    

    INTEGER, PARAMETER                      :: NDIMS=3
    INTEGER(KIND=HSIZE_T), DIMENSION(NDIMS) :: DATA_DIMS
    INTEGER(KIND=HSIZE_T), DIMENSION(NDIMS) :: OFFSET 
    INTEGER(KIND=HID_T)                     :: DSET_ID
    INTEGER(KIND=HID_T)                     :: LOCAL_DATASPACE_ID
    INTEGER(KIND=HID_T)                     :: GLOBAL_DATASPACE_ID
    INTEGER                                 :: IERR
    integer(kind=hsize_t)                   :: maxdims(ndims),ndims_(ndims)
    
    !> LOCAL 
    DATA_DIMS = IE_LOC - IS_LOC + 1 !> dimensions locales
    OFFSET    = 0                   !> pas d'offset sur les données locales
    CALL H5SCREATE_SIMPLE_F( NDIMS, DATA_DIMS, LOCAL_DATASPACE_ID, IERR )
    
    CALL H5SSELECT_HYPERSLAB_F(         &
         SPACE_ID = LOCAL_DATASPACE_ID, &
         OPERATOR = H5S_SELECT_SET_F,   &
         START    = OFFSET,             &
         COUNT    = DATA_DIMS,          &
         HDFERR   = IERR                )
    
    !> file or globale space
    DATA_DIMS = IE_GLB - IS_GLB + 1
    CALL H5SCREATE_SIMPLE_F ( NDIMS, DATA_DIMS, GLOBAL_DATASPACE_ID, IERR )
    
    DATA_DIMS = IE_LOC - IS_LOC + 1
    OFFSET = IS_LOC - IS_GLB
      CALL H5SSELECT_HYPERSLAB_F(           &
           SPACE_ID = GLOBAL_DATASPACE_ID,  &
           OPERATOR = H5S_SELECT_SET_F,     &
           START    = OFFSET,               &
           COUNT    = DATA_DIMS,            &
           HDFERR   = IERR                  )
      
      CALL H5DOPEN_F(F_ID, DSET_NAME, DSET_ID, IERR)
      CALL H5SGET_SIMPLE_EXTENT_DIMS_F( GLOBAL_DATASPACE_ID, NDIMS_, MAXDIMS, IERR) 
      DATA_DIMS = IE_LOC-IS_LOC+1
      CALL H5DREAD_F(DSET_ID, H5T_NATIVE_DOUBLE, FI(IS_LOC(1):IE_LOC(1),IS_LOC(2):IE_LOC(2),IS_LOC(3):IE_LOC(3)),&
           DATA_DIMS, IERR,FILE_SPACE_ID=GLOBAL_DATASPACE_ID, MEM_SPACE_ID=LOCAL_DATASPACE_ID, XFER_PRP=X_ID )
      
      CALL H5DCLOSE_F( DSET_ID, IERR )
      
    end subroutine READ_ARRAY_RANK_3
    
    
    SUBROUTINE WRITE_ARRAY_RANK_2( F_ID, X_ID, DSET_NAME, FI, IS_GLB, IE_GLB, IS_LOC, IE_LOC  )
      IMPLICIT NONE
      INTEGER(KIND=HID_T)                       :: F_ID
      INTEGER(KIND=HID_T)                       :: X_ID
      CHARACTER(LEN=*)                          :: DSET_NAME
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: FI
      INTEGER, DIMENSION(2)                     :: IS_GLB
      INTEGER, DIMENSION(2)                     :: IE_GLB
      INTEGER, DIMENSION(2)                     :: IS_LOC
      INTEGER, DIMENSION(2)                     :: IE_LOC
      

      INTEGER, PARAMETER                      :: NDIMS=2
      INTEGER(KIND=HSIZE_T), DIMENSION(NDIMS) :: DATA_DIMS
      INTEGER(KIND=HSIZE_T), DIMENSION(NDIMS) :: OFFSET 
      INTEGER(KIND=HID_T)                     :: DSET_ID
      INTEGER(KIND=HID_T)                     :: LOCAL_DATASPACE_ID
      INTEGER(KIND=HID_T)                     :: GLOBAL_DATASPACE_ID
      INTEGER                                 :: IERR

      !> LOCAL 
      DATA_DIMS = IE_LOC - IS_LOC + 1 !> dimensions locales
      OFFSET    = 0                   !> pas d'offset sur les données locales
      CALL H5SCREATE_SIMPLE_F( NDIMS, DATA_DIMS, LOCAL_DATASPACE_ID, IERR )
      
      CALL H5SSELECT_HYPERSLAB_F(         &
           SPACE_ID = LOCAL_DATASPACE_ID, &
           OPERATOR = H5S_SELECT_SET_F,   &
           START    = OFFSET,             &
           COUNT    = DATA_DIMS,          &
           HDFERR   = IERR                )
      
      !> file or globale space
      DATA_DIMS = IE_GLB - IS_GLB + 1
      CALL H5SCREATE_SIMPLE_F ( NDIMS, DATA_DIMS, GLOBAL_DATASPACE_ID, IERR )
      
      DATA_DIMS = IE_LOC - IS_LOC + 1
      OFFSET = IS_LOC - IS_GLB
      CALL H5SSELECT_HYPERSLAB_F(           &
           SPACE_ID = GLOBAL_DATASPACE_ID,  &
           OPERATOR = H5S_SELECT_SET_F,     &
           START    = OFFSET,               &
           COUNT    = DATA_DIMS,            &
           HDFERR   = IERR                  )
      
      CALL H5DOPEN_F(F_ID, DSET_NAME, DSET_ID, IERR)
      DATA_DIMS = IE_LOC-IS_LOC+1
      CALL H5DWRITE_F(DSET_ID, H5T_NATIVE_DOUBLE, FI(IS_LOC(1):IE_LOC(1),IS_LOC(2):IE_LOC(2)), DATA_DIMS, IERR,&
           FILE_SPACE_ID=GLOBAL_DATASPACE_ID, MEM_SPACE_ID=LOCAL_DATASPACE_ID, XFER_PRP=X_ID )
      
      CALL H5DCLOSE_F( DSET_ID, IERR )
      
    end subroutine WRITE_ARRAY_RANK_2

    SUBROUTINE READ_ARRAY_RANK_2( F_ID, X_ID, DSET_NAME, FI, IS_GLB, IE_GLB, IS_LOC, IE_LOC  )
      IMPLICIT NONE
      INTEGER(KIND=HID_T)                       :: F_ID
      INTEGER(KIND=HID_T)                       :: X_ID
      CHARACTER(LEN=*)                          :: DSET_NAME
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: FI
      INTEGER, DIMENSION(2)                     :: IS_GLB
      INTEGER, DIMENSION(2)                     :: IE_GLB
      INTEGER, DIMENSION(2)                     :: IS_LOC
      INTEGER, DIMENSION(2)                     :: IE_LOC
      

      INTEGER, PARAMETER                      :: NDIMS=2
      INTEGER(KIND=HSIZE_T), DIMENSION(NDIMS) :: DATA_DIMS
      INTEGER(KIND=HSIZE_T), DIMENSION(NDIMS) :: OFFSET 
      INTEGER(KIND=HID_T)                     :: DSET_ID
      INTEGER(KIND=HID_T)                     :: LOCAL_DATASPACE_ID
      INTEGER(KIND=HID_T)                     :: GLOBAL_DATASPACE_ID
      INTEGER                                 :: IERR
      integer(kind=hsize_t)                   :: maxdims(ndims),ndims_(ndims)
      
      !> LOCAL 
      DATA_DIMS = IE_LOC - IS_LOC + 1 !> dimensions locales
      OFFSET    = 0                   !> pas d'offset sur les données locales
      CALL H5SCREATE_SIMPLE_F( NDIMS, DATA_DIMS, LOCAL_DATASPACE_ID, IERR )
      
      CALL H5SSELECT_HYPERSLAB_F(         &
           SPACE_ID = LOCAL_DATASPACE_ID, &
           OPERATOR = H5S_SELECT_SET_F,   &
           START    = OFFSET,             &
           COUNT    = DATA_DIMS,          &
           HDFERR   = IERR                )
      
      !> file or globale space
      DATA_DIMS = IE_GLB - IS_GLB + 1
      CALL H5SCREATE_SIMPLE_F ( NDIMS, DATA_DIMS, GLOBAL_DATASPACE_ID, IERR )
      
      DATA_DIMS = IE_LOC - IS_LOC + 1
      OFFSET = IS_LOC - IS_GLB
      CALL H5SSELECT_HYPERSLAB_F(           &
           SPACE_ID = GLOBAL_DATASPACE_ID,  &
           OPERATOR = H5S_SELECT_SET_F,     &
           START    = OFFSET,               &
           COUNT    = DATA_DIMS,            &
           HDFERR   = IERR                  )
      
      CALL H5DOPEN_F(F_ID, DSET_NAME, DSET_ID, IERR)
      CALL H5SGET_SIMPLE_EXTENT_DIMS_F( GLOBAL_DATASPACE_ID, NDIMS_, MAXDIMS, IERR) 
      DATA_DIMS = IE_LOC-IS_LOC+1
      CALL H5DREAD_F(DSET_ID, H5T_NATIVE_DOUBLE, FI(IS_LOC(1):IE_LOC(1),IS_LOC(2):IE_LOC(2)), DATA_DIMS, IERR,&
           FILE_SPACE_ID=GLOBAL_DATASPACE_ID, MEM_SPACE_ID=LOCAL_DATASPACE_ID, XFER_PRP=X_ID )
      
     CALL H5DCLOSE_F( DSET_ID, IERR )
      
    end subroutine READ_ARRAY_RANK_2


    SUBROUTINE WRITE_ARRAY_RANK_1( F_ID, X_ID, DSET_NAME, FI, IS_GLB, IE_GLB, IS_LOC, IE_LOC  )
      IMPLICIT NONE
      INTEGER(KIND=HID_T)                       :: F_ID
      INTEGER(KIND=HID_T)                       :: X_ID
      CHARACTER(LEN=*)                          :: DSET_NAME
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE   :: FI
      INTEGER                                   :: IS_GLB
      INTEGER                                   :: IE_GLB
      INTEGER                                   :: IS_LOC
      INTEGER                                   :: IE_LOC
      

      INTEGER, PARAMETER                      :: NDIMS=1
      INTEGER(KIND=HSIZE_T), DIMENSION(NDIMS) :: DATA_DIMS
      INTEGER(KIND=HSIZE_T), DIMENSION(NDIMS) :: OFFSET 
      INTEGER(KIND=HID_T)                     :: DSET_ID
      INTEGER(KIND=HID_T)                     :: LOCAL_DATASPACE_ID
      INTEGER(KIND=HID_T)                     :: GLOBAL_DATASPACE_ID
      INTEGER                                 :: IERR

      !> LOCAL 
      DATA_DIMS = IE_LOC - IS_LOC + 1 !> dimensions locales
      OFFSET    = 0                   !> pas d'offset sur les données locales
      CALL H5SCREATE_SIMPLE_F( NDIMS, DATA_DIMS, LOCAL_DATASPACE_ID, IERR )
      
      CALL H5SSELECT_HYPERSLAB_F(         &
           SPACE_ID = LOCAL_DATASPACE_ID, &
           OPERATOR = H5S_SELECT_SET_F,   &
           START    = OFFSET,             &
           COUNT    = DATA_DIMS,          &
           HDFERR   = IERR                )
      
      !> file or globale space
      DATA_DIMS = IE_GLB - IS_GLB + 1
      CALL H5SCREATE_SIMPLE_F ( NDIMS, DATA_DIMS, GLOBAL_DATASPACE_ID, IERR )
      
      DATA_DIMS = IE_LOC - IS_LOC + 1
      OFFSET = IS_LOC - IS_GLB 
      CALL H5SSELECT_HYPERSLAB_F(           &
           SPACE_ID = GLOBAL_DATASPACE_ID,  &
           OPERATOR = H5S_SELECT_SET_F,     &
           START    = OFFSET,               &
           COUNT    = DATA_DIMS,            &
           HDFERR   = IERR                  )
      CALL H5DOPEN_F(F_ID, TRIM(DSET_NAME), DSET_ID, IERR)
      DATA_DIMS = IE_LOC-IS_LOC+1
      CALL H5DWRITE_F(DSET_ID, H5T_NATIVE_DOUBLE, FI(IS_LOC:IE_LOC), DATA_DIMS, IERR,&
           FILE_SPACE_ID=GLOBAL_DATASPACE_ID, MEM_SPACE_ID=LOCAL_DATASPACE_ID, XFER_PRP=X_ID )
      
      CALL H5DCLOSE_F( DSET_ID, IERR )
      
    end subroutine WRITE_ARRAY_RANK_1

    SUBROUTINE READ_ARRAY_RANK_1( F_ID, X_ID, DSET_NAME, FI, IS_GLB, IE_GLB, IS_LOC, IE_LOC  )
      IMPLICIT NONE
      INTEGER(KIND=HID_T)                       :: F_ID
      INTEGER(KIND=HID_T)                       :: X_ID
      CHARACTER(LEN=*)                          :: DSET_NAME
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE   :: FI
      INTEGER                                   :: IS_GLB
      INTEGER                                   :: IE_GLB
      INTEGER                                   :: IS_LOC
      INTEGER                                   :: IE_LOC

      INTEGER, PARAMETER                      :: NDIMS=1
      INTEGER(KIND=HSIZE_T), DIMENSION(NDIMS) :: DATA_DIMS
      INTEGER(KIND=HSIZE_T), DIMENSION(NDIMS) :: OFFSET 
      INTEGER(KIND=HID_T)                     :: DSET_ID
      INTEGER(KIND=HID_T)                     :: LOCAL_DATASPACE_ID
      INTEGER(KIND=HID_T)                     :: GLOBAL_DATASPACE_ID
      INTEGER                                 :: IERR
      integer(kind=hsize_t)                   :: maxdims(ndims),ndims_(ndims)
      
      !> LOCAL 
      DATA_DIMS = IE_LOC - IS_LOC + 1 !> dimensions locales
      OFFSET    = 0                   !> pas d'offset sur les données locales
      CALL H5SCREATE_SIMPLE_F( NDIMS, DATA_DIMS, LOCAL_DATASPACE_ID, IERR )
      
      CALL H5SSELECT_HYPERSLAB_F(         &
           SPACE_ID = LOCAL_DATASPACE_ID, &
           OPERATOR = H5S_SELECT_SET_F,   &
           START    = OFFSET,             &
           COUNT    = DATA_DIMS,          &
           HDFERR   = IERR                )
      
      !> file or globale space
      DATA_DIMS = IE_GLB - IS_GLB + 1
      CALL H5SCREATE_SIMPLE_F ( NDIMS, DATA_DIMS, GLOBAL_DATASPACE_ID, IERR )
      
      DATA_DIMS = IE_LOC - IS_LOC + 1
      OFFSET = IS_LOC - IS_GLB
      CALL H5SSELECT_HYPERSLAB_F(           &
           SPACE_ID = GLOBAL_DATASPACE_ID,  &
           OPERATOR = H5S_SELECT_SET_F,     &
           START    = OFFSET,               &
           COUNT    = DATA_DIMS,            &
           HDFERR   = IERR                  )
      
      CALL H5DOPEN_F(F_ID, DSET_NAME, DSET_ID, IERR)
      CALL H5SGET_SIMPLE_EXTENT_DIMS_F( GLOBAL_DATASPACE_ID, NDIMS_, MAXDIMS, IERR) 
      DATA_DIMS = IE_LOC-IS_LOC+1
      CALL H5DREAD_F(DSET_ID, H5T_NATIVE_DOUBLE, FI(IS_LOC:IE_LOC), DATA_DIMS, IERR,&
           FILE_SPACE_ID=GLOBAL_DATASPACE_ID, MEM_SPACE_ID=LOCAL_DATASPACE_ID, XFER_PRP=X_ID )
      
     CALL H5DCLOSE_F( DSET_ID, IERR )
      
   end subroutine READ_ARRAY_RANK_1


   SUBROUTINE CREATE_VECTOR_ATTRIBUTE(F_ID,GROUP_NAME,ATTRIBUTE_NAME,ATTRIBUTE_TYPE,DIM)
      IMPLICIT NONE
      INTEGER(KIND=HID_T)  :: F_ID
      CHARACTER(LEN=*)     :: GROUP_NAME
      CHARACTER(LEN=*)     :: ATTRIBUTE_NAME
      INTEGER(KIND=HID_T)  :: ATTRIBUTE_TYPE
      INTEGER              :: DIM
      INTEGER(KIND=HID_T)  :: GROUP_ID, SPACE_ID, ATTR_ID
      INTEGER,                           PARAMETER :: RANK=1
      INTEGER(HSIZE_T), DIMENSION(RANK)            :: DIMS
      INTEGER :: IERR
      
      DIMS = DIM
      CALL H5GOPEN_F( F_ID, GROUP_NAME, GROUP_ID, IERR )
      CALL H5SCREATE_SIMPLE_F(RANK, DIMS, SPACE_ID, IERR)
      CALL H5ACREATE_F( GROUP_ID,TRIM(ATTRIBUTE_NAME), ATTRIBUTE_TYPE, SPACE_ID, ATTR_ID, IERR)
      CALL H5ACLOSE_F(  ATTR_ID, IERR )
      CALL H5SCLOSE_F( SPACE_ID, IERR )
      CALL H5GCLOSE_F( GROUP_ID, IERR )
    END SUBROUTINE CREATE_VECTOR_ATTRIBUTE
    
    

    SUBROUTINE WRITE_VECTOR_INTEGER_ATTRIBUTE(F_ID,GROUP_NAME,ATTRIBUTE_NAME,VAL,DIM)
      IMPLICIT NONE
      INTEGER(KIND=HID_T)    :: F_ID
      CHARACTER(LEN=*)       :: GROUP_NAME
      CHARACTER(LEN=*)       :: ATTRIBUTE_NAME
      INTEGER                :: DIM
      INTEGER,dimension(DIM) :: VAL
      INTEGER(KIND=HID_T)  :: GROUP_ID, SPACE_ID, ATTR_ID
      INTEGER,                        PARAMETER :: RANK=1
      INTEGER(HSIZE_T), DIMENSION(1), PARAMETER :: DIMS=1
      INTEGER, DIMENSION(DIM) :: BUFF
      INTEGER :: IERR

      CALL H5GOPEN_F( F_ID, GROUP_NAME, GROUP_ID, IERR )
      CALL H5AOPEN_F(GROUP_ID,ATTRIBUTE_NAME, ATTR_ID, IERR)
      
      CALL H5SCREATE_SIMPLE_F(RANK, DIMS, SPACE_ID, IERR)
      BUFF = VAL
      CALL H5AWRITE_F( ATTR_ID, H5T_NATIVE_INTEGER, BUFF, DIMS , IERR)

      CALL H5SCLOSE_F( SPACE_ID, IERR )
      CALL H5ACLOSE_F( ATTR_ID, IERR )
      CALL H5GCLOSE_F( GROUP_ID, IERR )
      
    END SUBROUTINE WRITE_VECTOR_INTEGER_ATTRIBUTE
    
    SUBROUTINE READ_VECTOR_INTEGER_ATTRIBUTE(F_ID,GROUP_NAME,ATTRIBUTE_NAME,VAL,DIM)
      IMPLICIT NONE
      INTEGER(KIND=HID_T)    :: F_ID
      CHARACTER(LEN=*)       :: GROUP_NAME
      CHARACTER(LEN=*)       :: ATTRIBUTE_NAME
      INTEGER                :: DIM
      INTEGER,dimension(DIM) :: VAL
      INTEGER(KIND=HID_T)  :: GROUP_ID, SPACE_ID, ATTR_ID
      INTEGER,                        PARAMETER :: RANK=1
      INTEGER(HSIZE_T), DIMENSION(1), PARAMETER :: DIMS=1
      INTEGER, DIMENSION(1) :: BUFF
      INTEGER :: IERR

      CALL H5GOPEN_F( F_ID, GROUP_NAME, GROUP_ID, IERR )
      CALL H5AOPEN_F(GROUP_ID,ATTRIBUTE_NAME, ATTR_ID, IERR)
      
      CALL H5SCREATE_SIMPLE_F(RANK, DIMS, SPACE_ID, IERR)
      CALL H5AREAD_F( ATTR_ID, H5T_NATIVE_INTEGER, BUFF, DIMS , IERR)
      VAL = BUFF(DIM)
      CALL H5SCLOSE_F( SPACE_ID, IERR )
      CALL H5ACLOSE_F( ATTR_ID, IERR )
      CALL H5GCLOSE_F( GROUP_ID, IERR )
      
    END SUBROUTINE READ_VECTOR_INTEGER_ATTRIBUTE




    
    SUBROUTINE CREATE_SCALAR_ATTRIBUTE(F_ID,GROUP_NAME,ATTRIBUTE_NAME,ATTRIBUTE_TYPE)
      IMPLICIT NONE
      INTEGER(KIND=HID_T)  :: F_ID
      CHARACTER(LEN=*)     :: GROUP_NAME
      CHARACTER(LEN=*)     :: ATTRIBUTE_NAME
      INTEGER(KIND=HID_T)  :: ATTRIBUTE_TYPE
      INTEGER(KIND=HID_T)  :: GROUP_ID, SPACE_ID, ATTR_ID
      INTEGER,                        PARAMETER :: RANK=1
      INTEGER(HSIZE_T), DIMENSION(1), PARAMETER :: DIMS=1
      INTEGER :: IERR
      
      CALL H5GOPEN_F( F_ID, GROUP_NAME, GROUP_ID, IERR )
      CALL H5SCREATE_SIMPLE_F(RANK, DIMS, SPACE_ID, IERR)
      CALL H5ACREATE_F( GROUP_ID,TRIM(ATTRIBUTE_NAME), ATTRIBUTE_TYPE, SPACE_ID, ATTR_ID, IERR)
      CALL H5ACLOSE_F(  ATTR_ID, IERR )
      CALL H5SCLOSE_F( SPACE_ID, IERR )
      CALL H5GCLOSE_F( GROUP_ID, IERR )
    END SUBROUTINE CREATE_SCALAR_ATTRIBUTE
    
    
    SUBROUTINE WRITE_REAL8_ATTRIBUTE(F_ID,GROUP_NAME,ATTRIBUTE_NAME,VAL)
      IMPLICIT NONE
      INTEGER(KIND=HID_T)  :: F_ID
      CHARACTER(LEN=*)     :: GROUP_NAME
      CHARACTER(LEN=*)     :: ATTRIBUTE_NAME
      REAL(KIND=8)         :: VAL
      INTEGER(KIND=HID_T)  :: GROUP_ID, SPACE_ID, ATTR_ID
      INTEGER,                        PARAMETER :: RANK=1
      INTEGER(HSIZE_T), DIMENSION(1), PARAMETER :: DIMS=1
      REAL(KIND=8), DIMENSION(1) :: BUFF
      INTEGER :: IERR

      CALL H5GOPEN_F( F_ID, GROUP_NAME, GROUP_ID, IERR )
      CALL H5AOPEN_F(GROUP_ID,ATTRIBUTE_NAME, ATTR_ID, IERR)
      
      CALL H5SCREATE_SIMPLE_F(RANK, DIMS, SPACE_ID, IERR)
      BUFF = VAL
      CALL H5AWRITE_F( ATTR_ID, H5T_NATIVE_DOUBLE, BUFF, DIMS , IERR)

      CALL H5SCLOSE_F( SPACE_ID, IERR )
      CALL H5ACLOSE_F( ATTR_ID, IERR )
      CALL H5GCLOSE_F( GROUP_ID, IERR )

    END SUBROUTINE WRITE_REAL8_ATTRIBUTE

    SUBROUTINE READ_REAL8_ATTRIBUTE(F_ID,GROUP_NAME,ATTRIBUTE_NAME,VAL)
      IMPLICIT NONE
      INTEGER(KIND=HID_T)  :: F_ID
      CHARACTER(LEN=*)     :: GROUP_NAME
      CHARACTER(LEN=*)     :: ATTRIBUTE_NAME
      REAL(KIND=8)         :: VAL
      INTEGER(KIND=HID_T)  :: GROUP_ID, SPACE_ID, ATTR_ID
      INTEGER,                        PARAMETER :: RANK=1
      INTEGER(HSIZE_T), DIMENSION(1), PARAMETER :: DIMS=1
      REAL(KIND=8), DIMENSION(1) :: BUFF
      INTEGER :: IERR

      CALL H5GOPEN_F( F_ID, GROUP_NAME, GROUP_ID, IERR )
      CALL H5AOPEN_F(GROUP_ID,ATTRIBUTE_NAME, ATTR_ID, IERR)
      
      CALL H5SCREATE_SIMPLE_F(RANK, DIMS, SPACE_ID, IERR)
      
      CALL H5AREAD_F( ATTR_ID, H5T_NATIVE_DOUBLE, BUFF, DIMS , IERR )
      VAL = BUFF(1)
      CALL H5SCLOSE_F( SPACE_ID, IERR )
      CALL H5ACLOSE_F( ATTR_ID, IERR )
      CALL H5GCLOSE_F( GROUP_ID, IERR )

    END SUBROUTINE READ_REAL8_ATTRIBUTE
    
    SUBROUTINE WRITE_INTEGER_ATTRIBUTE(F_ID,GROUP_NAME,ATTRIBUTE_NAME,VAL)
      IMPLICIT NONE
      INTEGER(KIND=HID_T)  :: F_ID
      CHARACTER(LEN=*)     :: GROUP_NAME
      CHARACTER(LEN=*)     :: ATTRIBUTE_NAME
      INTEGER              :: VAL
      INTEGER(KIND=HID_T)  :: GROUP_ID, SPACE_ID, ATTR_ID
      INTEGER,                        PARAMETER :: RANK=1
      INTEGER(HSIZE_T), DIMENSION(1), PARAMETER :: DIMS=1
      INTEGER, DIMENSION(1) :: BUFF
      INTEGER :: IERR


      CALL H5GOPEN_F( F_ID, GROUP_NAME, GROUP_ID, IERR )
      CALL H5AOPEN_F(GROUP_ID,ATTRIBUTE_NAME, ATTR_ID, IERR)
      
      CALL H5SCREATE_SIMPLE_F(RANK, DIMS, SPACE_ID, IERR)
      BUFF = VAL
      CALL H5AWRITE_F( ATTR_ID, H5T_NATIVE_INTEGER, BUFF, DIMS , IERR)

      CALL H5SCLOSE_F( SPACE_ID, IERR )
      CALL H5ACLOSE_F( ATTR_ID, IERR )
      CALL H5GCLOSE_F( GROUP_ID, IERR )
      
    END SUBROUTINE WRITE_INTEGER_ATTRIBUTE
    
    SUBROUTINE READ_INTEGER_ATTRIBUTE(F_ID,GROUP_NAME,ATTRIBUTE_NAME,VAL)
      IMPLICIT NONE
      INTEGER(KIND=HID_T)  :: F_ID
      CHARACTER(LEN=*)     :: GROUP_NAME
      CHARACTER(LEN=*)     :: ATTRIBUTE_NAME
      INTEGER              :: VAL
      INTEGER(KIND=HID_T)  :: GROUP_ID, SPACE_ID, ATTR_ID
      INTEGER,                        PARAMETER :: RANK=1
      INTEGER(HSIZE_T), DIMENSION(1), PARAMETER :: DIMS=1
      INTEGER, DIMENSION(1) :: BUFF
      INTEGER :: IERR

      CALL H5GOPEN_F( F_ID, GROUP_NAME, GROUP_ID, IERR )
      CALL H5AOPEN_F(GROUP_ID,ATTRIBUTE_NAME, ATTR_ID, IERR)
      
      CALL H5SCREATE_SIMPLE_F(RANK, DIMS, SPACE_ID, IERR)
      CALL H5AREAD_F( ATTR_ID, H5T_NATIVE_INTEGER, BUFF, DIMS , IERR)
      VAL = BUFF(1)
      CALL H5SCLOSE_F( SPACE_ID, IERR )
      CALL H5ACLOSE_F( ATTR_ID, IERR )
      CALL H5GCLOSE_F( GROUP_ID, IERR )
      
    END SUBROUTINE READ_INTEGER_ATTRIBUTE
    
    SUBROUTINE CREATE_A_GROUP(F_ID,PATH)
      IMPLICIT NONE
      CHARACTER(LEN=*)     :: PATH
      INTEGER(KIND=HID_T)  :: F_ID
      INTEGER(KIND=HID_T)  :: G_ID
      INTEGER :: IERR
      
      CALL H5GCREATE_F (F_ID,PATH, G_ID, IERR)
      CALL H5GCLOSE_F( G_ID, IERR )
      
    END SUBROUTINE CREATE_A_GROUP



    
  

  end module m_hdf5_ifce
  
