program tcheby_1d
  use m_mesh_base
  use m_numerics
  use m_operator_base
  use m_operator_tcheby
  use m_operator_fourier_dft
  use m_adjoint_tool_cyl
  use m_fourier_transform
  use m_tensor_product

  use m_hdf5_ifce
  use mpi
  use decomp_2d
!  use m_udf_mms
  !>
  use m_solver_diag_cyl_hhi
  
  !>
  use m_navier_stokes_cyl
  use m_snapshots
  use m_quadrature
  use ifport
  implicit none
  
  type(t_mesh_base)       ::  msh(3)
  !>
  type(t_solver_diag_cyl_hhi) :: eqn_ua, eqn_uz, eqn_ur, eqn_ua_ad, eqn_uz_ad, eqn_ur_ad
  type(t_solver_diag_cyl_hhi) :: eqn_fi,eqn_fi_ad
  
  !>
  type(t_operator_fourier_dft) ::  opa,opz
  type(t_operator_tcheby) :: opr
  type(t_quadrature) :: quad
  !> 
  real(kind=8),allocatable     ::  dg01(:,:,:),dg02(:,:,:),dg03(:,:,:)
  real(kind=8),allocatable     ::  dg04(:,:,:),dg05(:,:,:),dg06(:,:,:)
  real(kind=8),allocatable     ::  dg07(:,:,:),dg08(:,:,:),dg09(:,:,:),DG10(:,:,:)
  real(kind=8),allocatable     ::  dg11(:,:,:),dg12(:,:,:),dg13(:,:,:),DG14(:,:,:)

  integer                  ::  i,j,k

  real(kind=8) :: dirichl(2),neumann(2)
  integer :: n(3),cpu_grid(2),ierr
  REAL(kind=8) ::  xmin(3), xmax(3),err
  REAL(kind=8) ::  tc

  TYPE(DECOMP_INFO) :: ph
  
  real(kind=8),allocatable ,dimension(:,:,:) ::  UA_T,UZ_T,UR_T
  real(kind=8),allocatable ,dimension(:,:,:) ::  UAM1,UZM1,URM1
  real(kind=8),allocatable ,dimension(:,:,:) ::  NLAM1,NLZM1,NLRM1
  real(kind=8),allocatable ,dimension(:,:,:) ::  UA,UZ,UR,PRES,FI
  real(kind=8),allocatable ,dimension(:,:,:) ::  SA,SZ,SR,SFI
  real(kind=8),allocatable ,dimension(:,:,:) ::  NLA,NLZ,NLR
  real(kind=8),allocatable ,dimension(:,:,:) ::  A,Z,R,A_tot,Z_tot,R_tot
  real(kind=8) :: w2,t1,h
  integer :: it_time,is(3),ie(3)
  integer , parameter :: OX=1,OY=2,OZ=3
  REAL(DP) :: prm_A,prm_B,eta,OMEGA_i,OMEGA_o,integ
  
  ! attention au type à lire
  REAL(kind=8) :: a_min,a_max,z_min,z_max,r_min,r_max
  integer :: na,nz,nr,nb_cpu_y,nb_cpu_z
  namelist /parameters_cube/ na,nz,nr,z_min,z_max,r_min,r_max,nb_cpu_y,nb_cpu_z
  
  REAL(KIND=8) :: NU,prm_K,RE
  INTEGER :: nb_iter
  namelist /parameters_physical/RE

  character(len=1024) ::  root_dir,io_nrj,output_dir
  REAL(kind=8) ::  dt,cfl,tmax
  integer dn_dump
  logical :: resume , explicit,correction_pres
  namelist /parameters_timescheme/ tmax,dt, cfl, dn_dump,root_dir,io_nrj
  
  logical :: resume_snap=.false.
  integer :: dn_snap
  character(len=1024) ::  snap_dir
  real(kind=8) :: snap_dt,starttime,endtime
  integer :: snap_id
  namelist /parameters_diagnostics/ dn_snap,snap_dir

  character(len=1024) ::  input_file,init_file
  integer :: rank,iostat,seed_size
  integer, allocatable :: myseed(:)
  character(len=1024) ::  file_dump, base_snap, filename, base_save
  character(len=6) :: num
  REAL(kind=8) ::  DIV_MAX,tmp,DJ_ADJ,J_DU0,J_U
  
  REAL(KIND=8),DIMENSION(3) :: NU_MOMENTUM
  REAL(KIND=8),DIMENSION(3) :: NU_ENERGY
  REAL(KIND=8),DIMENSION(3) :: NU_POISSON
  REAL(kind=8) :: sigma,omega_tilde,noise,res,alpha,beta
  
  REAL(DP),DIMENSION(:,:,:,:),ALLOCATABLE :: SAVE_UA,SAVE_UZ,SAVE_UR

  COMPLEX(DP),ALLOCATABLE,DIMENSION(:,:,:) :: DGA,DGZ,DGR
  COMPLEX(DP),ALLOCATABLE,DIMENSION(:,:,:) :: DGA_Y,DGZ_Y,DGR_Y
  
  LOGICAL :: memoire = .TRUE.
  LOGICAL :: do_adj = .FALSE.
  call mpi_init(ierr)
  CALL H5OPEN_F(IERR)

  call mpi_comm_rank(mpi_comm_world,rank,ierr)
  
  call command_line_read_input_adj(input_file,output_dir,init_file,do_adj)

  if (rank==0) then
     OPEN (UNIT=24, FILE=TRIM(input_file),status='old', action='read')
     read(24, nml=parameters_cube, IOSTAT=iostat)
     read(24, nml=parameters_physical, IOSTAT=iostat)
     read(24, nml=parameters_timescheme, IOSTAT=iostat)
     read(24, nml=parameters_diagnostics, IOSTAT=iostat)
     close(24)
     ierr = SYSTEM('mkdir -p '//trim(root_dir)//'/dump' )
     ierr = SYSTEM('mkdir -p '//trim(root_dir)//'/save' )
     ierr = SYSTEM('mkdir -p '//trim(root_dir)//trim(snap_dir) )
     write(*,parameters_cube)
     write(*,parameters_physical)
     write(6,*)TRIM(TRIM(root_dir)//'timevar')
     OPEN(UNIT=42, FILE=TRIM(TRIM(root_dir)//'timevar'))
  end if

  
  CALL MPI_BCAST( na      , 1,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
  CALL MPI_BCAST( nz      , 1,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
  CALL MPI_BCAST( nr      , 1,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
  CALL MPI_BCAST( nb_iter , 1,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)


  CALL MPI_BCAST( a_min   , 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  CALL MPI_BCAST( a_max   , 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  CALL MPI_BCAST( z_min   , 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  CALL MPI_BCAST( z_max   , 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  CALL MPI_BCAST( r_min   , 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  CALL MPI_BCAST( r_max   , 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  
  CALL MPI_BCAST( nb_cpu_y, 1,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
  CALL MPI_BCAST( nb_cpu_z, 1,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
  
  
  CALL MPI_BCAST( RE, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR ) 
  
  CALL MPI_BCAST( resume  , 1   , MPI_LOGICAL         , 0, MPI_COMM_WORLD, IERR ) 
  CALL MPI_BCAST( root_dir, 1024, MPI_character       , 0, MPI_COMM_WORLD, IERR )
  CALL MPI_BCAST( dn_dump , 1   , MPI_INTEGER         , 0, MPI_COMM_WORLD, IERR )
  CALL MPI_BCAST( dt      , 1   , MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR ) 
  CALL MPI_BCAST( cfl     , 1   , MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR )
  CALL MPI_BCAST( tmax     , 1   , MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR )
  
  CALL MPI_BCAST( dn_snap , 1   , MPI_INTEGER         , 0, MPI_COMM_WORLD, IERR )
  CALL MPI_BCAST( snap_dir, 1024, MPI_character       , 0, MPI_COMM_WORLD, IERR )
  
  N = [NA,NZ,NR]  
  CPU_GRID = [NB_CPU_Y,NB_CPU_Z]

  ALPHA = 1._DP
  BETA = 5._DP
  
  xmin(1:3) = [0._DP,z_min,r_min]
  xmax(1:3) = [2._DP*pi,z_max,r_max]

  eta = r_min/r_max

  OMEGA_I = 1._DP
  OMEGA_O = (r_max)**(-1.5_DP)
  
  prm_A = (r_max/(1._DP - eta**2)) * (OMEGA_o - OMEGA_i*eta**2)
  prm_B =(r_min**2/(1._DP - eta**2)) * (OMEGA_i - OMEGA_o) 
  
  nb_iter = floor(tmax/dt)

  if (rank == 0) print*,nb_iter
  
  file_dump = trim(root_dir)//'dump/dump.h5'
  base_snap = trim(root_dir)//trim(snap_dir)//'snap_'
  base_save = trim(root_dir)//'save/save_'

   ! GRILLE 2D
  CALL DECOMP_2D_INIT( &
       NX = N(1)+1, NY = N(2)+1, NZ = N(3)+1 , P_ROW= CPU_GRID(1) , P_COL= CPU_GRID(2) )

  CALL GET_DECOMP_INFO(PH)

  CALL START_TENSOR_PRODUCT(&
       PH%XST,PH%XEN,N(1)+1,PH%YST,PH%YEN,N(2)+1,PH%ZST,PH%ZEN,N(3)+1)
  CALL START_FOURIER_TRANSFORMS( &
       PH%XST, PH%XEN, N(1)+1, PH%YST, PH%YEN,N(2)+1, PH%ZST, PH%ZEN, N(3)+1 )

  call test_dct(PH%XST,PH%XEN,NA+1,PH%YST,PH%YEN,NZ+1,PH%ZST,PH%ZEN,NR+1) ! attention 
  

  CALL MSH(1)%INITIALIZE(XMIN(1),XMAX(1),N(1),.TRUE. ) 
  CALL MSH(2)%INITIALIZE(XMIN(2),XMAX(2),N(2),.TRUE. ) 
  CALL MSH(3)%INITIALIZE(XMIN(3),XMAX(3),N(3) ) 


  CALL START_OP_BASE()
  CALL OPA%INIT_OPERATOR_FOURIER_DFT( MESH=MSH(OX), AXIS=OX )
  CALL OPZ%INIT_OPERATOR_FOURIER_DFT( MESH=MSH(OY), AXIS=OY )
  CALL OPR%INIT_OPERATOR_TCHEBY( MESH=MSH(OZ), AXIS=OZ )
  

  CALL INIT_quadrature_hhi(quad,xmin,xmax,&
       ph%xst,ph%xen,na,ph%yst,ph%yen,nz,ph%zst,ph%zen,nr)
  
  
  dirichl =  [1._DP,0._DP]
  neumann =  [0._DP,1._DP]

  NU = RE**(-1)
  
  NU_MOMENTUM = -  nu*[1,1,1]!*0.5
  NU_POISSON  =  1.

  CALL EQN_UA%SET_PARAMS( NU=NU_MOMENTUM , SIGMA=0D0 , AXIS=1 )
  CALL EQN_UA%SET_BCS( AXIS=3 , BCS_MINUS=DIRICHL , BCS_PLUS=DIRICHL )
  CALL EQN_UA%INITIALISE_SVV( MSH, OPA, OPZ, OPR, ALPHA, BETA,PH )
  CALL EQN_UA%SET_BVS( MESH = MSH , AXIS=3 , UDF_MINUS=UDF_NULL , UDF_PLUS=UDF_NULL )
  
  CALL EQN_UZ%SET_PARAMS( NU=NU_MOMENTUM , SIGMA=0D0 , AXIS=2 )
  CALL EQN_UZ%SET_BCS( AXIS=3 , BCS_MINUS=DIRICHL , BCS_PLUS=DIRICHL )
  CALL EQN_UZ%INITIALISE_SVV(MSH,OPA,OPZ,OPR,ALPHA,BETA,PH)
  CALL EQN_UZ%SET_BVS( MESH = MSH , AXIS=3 , UDF_MINUS=UDF_NULL , UDF_PLUS=UDF_NULL )
  
  CALL EQN_UR%SET_PARAMS( NU=NU_MOMENTUM , SIGMA=0D0 , AXIS=3 )
  CALL EQN_UR%SET_BCS( AXIS=3 , BCS_MINUS=DIRICHL , BCS_PLUS=DIRICHL )
  CALL EQN_UR%INITIALISE_SVV(MSH,OPA,OPZ,OPR,ALPHA,BETA,PH)
  CALL EQN_UR%SET_BVS( MESH = MSH , AXIS=3 , UDF_MINUS=UDF_NULL , UDF_PLUS=UDF_NULL )

  CALL EQN_FI%SET_PARAMS( NU=NU_POISSON , SIGMA=0D0 , AXIS=0 )
  CALL EQN_FI%SET_BCS( AXIS=3 , BCS_MINUS=NEUMANN , BCS_PLUS=NEUMANN )
  CALL EQN_FI%INITIALISE_SVV(MSH,OPA,OPZ,OPR,ALPHA,BETA,PH)
  CALL EQN_FI%SET_BVS( MESH = MSH , AXIS=3 , UDF_MINUS=UDF_NULL , UDF_PLUS=UDF_NULL )
  
  
  call preproc()

  !Lecture de la condition initiale
  CALL IMPORT_HDF5_INIT(TRIM(init_file),UA,UZ,UR)

!  CALL NORMALIZATION(UA,UZ,UR,1E-4)

  DG04 = UA*UA*R
  DG05 = UZ*UZ*R
  DG06 = UR*UR*R

  CALL integrate_spec(quad,DG04,J_U,PH,NA,NZ,NR,xmax,xmin)
  CALL integrate_spec(quad,DG05,integ,PH,NA,NZ,NR,xmax,xmin)
  J_U = integ + J_U
  CALL integrate_spec(quad,DG06,integ,PH,NA,NZ,NR,xmax,xmin)
  J_U = integ + J_U
  
  if (nrank == 0)write(6,*)"E0 : ",J_U
  
  
  !Allocattion des tableaux pour sauver les iterations   
  if (do_adj) then
     if (memoire) then
        ALLOCATE(SAVE_UA(0:nb_iter/2,PH%XST(1):PH%XEN(1),PH%XST(2):PH%XEN(2),PH%XST(3):PH%XEN(3))); SAVE_UA = 0._DP
        ALLOCATE(SAVE_UZ(0:nb_iter/2,PH%XST(1):PH%XEN(1),PH%XST(2):PH%XEN(2),PH%XST(3):PH%XEN(3))); SAVE_UZ = 0._DP
        ALLOCATE(SAVE_UR(0:nb_iter/2,PH%XST(1):PH%XEN(1),PH%XST(2):PH%XEN(2),PH%XST(3):PH%XEN(3))); SAVE_UR = 0._DP
  
        SAVE_UA(0,:,:,:) = UA(:,:,:)
        SAVE_UZ(0,:,:,:) = UZ(:,:,:)
        SAVE_UR(0,:,:,:) = UR(:,:,:)
     else
        write(num,'(I6.6)')0
        filename = trim(base_save)//num//'.h5'
        call save_hdf5(trim(FILENAME),MSH,UA,UZ,UR)
     end if
  end if

  it_time = 1
  
  SA = UA/DT
  SZ = UZ/DT
  SR = UR/DT

  
  CALL DUMP_HDF5_BASIC(FILE_DUMP,'NEW',TC,DT,MSH,UA,UZ,UR,PRES,DG01,DG09)

  UAM1=UA
  UZM1=UZ
  URM1=UR
  
  CALL dealiazing(ua,uz,ur)
  
  DG10 = UA + PRM_A*R + PRM_B/R
  
  CALL COMPUTE_NON_LINEAR_TERMS(&
       A, Z, R, OPA, OPZ, OPR, DG10, UZ, UR , NLAM1, NLZM1, NLRM1,&
       DG01, DG02, DG03, DG04, DG05, DG06, DG07, DG08, DG09)

  DG10 = PRM_A*R + PRM_B/R
  DG14 = 0._DP
  
  CALL COMPUTE_NON_LINEAR_TERMS(&
       A, Z, R, OPA, OPZ, OPR, DG10, DG14, DG14 , DG11, DG12, DG13,&
       DG01, DG02, DG03, DG04, DG05, DG06, DG07, DG08, DG09)

  NLA = NLAM1 - DG11
  NLZ = NLZM1 - DG12
  NLR = NLRM1 - DG13

  TC = dt

  CALL OPA%D1(UA,SFI)
  DG01 = NU*(-2._DP/R**2)*SFI
  CALL OPA%D1(UR,SFI)
  DG02 = NU*(+2._DP/R**2)*SFI

  NLR = NLR - DG01 
  NLA = NLA - DG02 
  NLZ = NLZ 
     
  SA = SA - NLA 
  SZ = SZ - NLZ 
  SR = SR - NLR
     
  NLAM1 = NLA
  NLZM1 = NLZ
  NLRM1 = NLR

  SIGMA = 1._DP/DT

  CALL EQN_UA%SOLVE(UA, SA, SIGMA ,PH)
  CALL EQN_UZ%SOLVE(UZ, SZ, SIGMA ,PH)
  CALL EQN_UR%SOLVE(UR, SR, SIGMA ,PH)
        
        
  call DIV( A, Z, R, OPA, OPZ, OPR, UA, UZ, UR, SFI, dg01, dg02 , dg03 )
  SFI = SFI/DT

  call EQN_FI%SOLVE(FI,SFI,0._DP,PH)
        
  CALL GRAD(A, Z, R, OPA, OPZ, OPR, FI, DG01, DG02, DG03)
        
  PRES = FI
        
        
  is = get_is_b([0,0,0])
  ie = get_ie_b([0,0,0])
  FORALL(I=IS(1):IE(1),J=IS(2):IE(2),K=IS(3):IE(3))
     UA(I,J,K) = UA(I,J,K) - DG01(I,J,K) * DT
  END FORALL
  is = get_is_b([0,0,0])
  ie = get_ie_b([0,0,0])
  FORALL(I=IS(1):IE(1),J=IS(2):IE(2),K=IS(3):IE(3))
     UZ(I,J,K) = UZ(I,J,K) - DG02(I,J,K) * DT 
  END FORALL
  is = get_is_b([0,0,0])
  ie = get_ie_b([0,0,0])
  FORALL(I=IS(1):IE(1),J=IS(2):IE(2),K=IS(3):IE(3))
     UR(I,J,K) = UR(I,J,K) - DG03(I,J,K) * DT 
  END FORALL

  if (do_adj .AND. MOD(IT_TIME,2)==0) then
     if (memoire) then
        SAVE_UA(it_time/2,:,:,:) = UA(:,:,:)
        SAVE_UZ(it_time/2,:,:,:) = UZ(:,:,:)
        SAVE_UR(it_time/2,:,:,:) = UR(:,:,:)
     else
        write(num,'(I6.6)')it_time
        filename = trim(base_save)//num//'.h5'
        call save_hdf5(trim(FILENAME),MSH,UA,UZ,UR)
     end if
  end if

  
  !Calcul de J(u) = int_domaine <u,u> + <t,t>
  DG04 = UA*UA*R
  DG05 = UZ*UZ*R
  DG06 = UR*UR*R

  CALL integrate_spec(quad,DG04,J_U,PH,NA,NZ,NR,xmax,xmin)
  CALL integrate_spec(quad,DG05,integ,PH,NA,NZ,NR,xmax,xmin)
  J_U = integ + J_U
  CALL integrate_spec(quad,DG06,integ,PH,NA,NZ,NR,xmax,xmin)
  J_U = integ + J_U
  
  if (nrank==0) write(42,*)dt,-J_U
  
  ! check divergence 
  call DIV( A, Z ,R, OPA, OPZ, OPR, UA, UZ, UR, DG09, dg01, dg02 , dg03 )
  is = get_is()
  ie = get_ie()
  DIV_MAX = MAXVAL(ABS(DG09(IS(1):IE(1),IS(2):IE(2),IS(3):IE(3))))
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,DIV_MAX,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,IERR)
  
  endtime   = MPI_Wtime();
  endtime =  endtime-starttime
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,endtime,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,IERR)
     
  
  call GetCFL(msh(1),msh(2),msh(3), UA, UZ, UR, dt, cfl)
  it_time = 1
  if (rank==0) print'(i9,11(1x,e15.8))',it_time,tc,dt,cfl,DIV_MAX,endtime,J_U
  
    
  DO IT_time=2,nb_iter
     
     tc = tc + dt
     snap_dt = snap_dt + dt
     
     starttime = MPI_Wtime();
        
     CALL GRAD( A,Z,R, &
          OPA, OPZ, OPR, PRES, DG01,DG02, DG03)
     
     SA = (2._DP*UA-0.5_DP*UAM1)/DT - DG01 
     SZ = (2._DP*UZ-0.5_DP*UZM1)/DT - DG02 
     SR = (2._DP*UR-0.5_DP*URM1)/DT - DG03
     
     UAM1=UA
     UZM1=UZ
     URM1=UR

     CALL dealiazing(ua,uz,ur)
     
     DG10 = UA + PRM_A*R + PRM_B/R

     CALL COMPUTE_NON_LINEAR_TERMS(&
          A, Z, R, OPA, OPZ, OPR, DG10, UZ, UR, NLA, NLZ, NLR,&
          DG01, DG02, DG03, DG04, DG05, DG06, DG07, DG08, DG09)

     DG10 = PRM_A*R + PRM_B/R
     DG14 = 0._DP
  
     CALL COMPUTE_NON_LINEAR_TERMS(&
          A, Z, R, OPA, OPZ, OPR, DG10, DG14, DG14 , DG11, DG12, DG13,&
          DG01, DG02, DG03, DG04, DG05, DG06, DG07, DG08, DG09)

     CALL OPA%D1(UA,SFI)
     DG01 = NU*(-2._DP/R**2)*SFI
     CALL OPA%D1(UR,SFI)
     DG02 = NU*(+2._DP/R**2)*SFI

     NLR = NLR - DG01 - DG13
     NLA = NLA - DG02 - DG11
     NLZ = NLZ - DG12
     
     SA = SA - 2._DP*NLA + NLAM1 
     SZ = SZ - 2._DP*NLZ + NLZM1 
     SR = SR - 2._DP*NLR + NLRM1
     
     NLAM1 = NLA
     NLZM1 = NLZ
     NLRM1 = NLR

     SIGMA = 1.5_DP/DT

     CALL EQN_UA%SOLVE(UA, SA, SIGMA ,PH)
     CALL EQN_UZ%SOLVE(UZ, SZ, SIGMA ,PH)
     CALL EQN_UR%SOLVE(UR, SR, SIGMA ,PH)
        
        
     call DIV( A, Z, R, OPA, OPZ, OPR, UA, UZ, UR, SFI, dg01, dg02 , dg03 )
     SFI = SFI*1.5_DP/DT

     call EQN_FI%SOLVE(FI,SFI,0._DP,PH)
        
     CALL GRAD(A, Z, R, OPA, OPZ, OPR, FI, DG01, DG02, DG03)
        
     PRES = PRES + FI
        
        
     is = get_is_b([0,0,0])
     ie = get_ie_b([0,0,0])
     FORALL(I=IS(1):IE(1),J=IS(2):IE(2),K=IS(3):IE(3))
        UA(I,J,K) = UA(I,J,K) - DG01(I,J,K) * 2._DP*DT/3._DP
     END FORALL
     is = get_is_b([0,0,0])
     ie = get_ie_b([0,0,0])
     FORALL(I=IS(1):IE(1),J=IS(2):IE(2),K=IS(3):IE(3))
        UZ(I,J,K) = UZ(I,J,K) - DG02(I,J,K) * 2._DP*DT/3._DP 
     END FORALL
     is = get_is_b([0,0,0])
     ie = get_ie_b([0,0,0])
     FORALL(I=IS(1):IE(1),J=IS(2):IE(2),K=IS(3):IE(3))
        UR(I,J,K) = UR(I,J,K) - DG03(I,J,K) * 2._DP*DT/3._DP 
     END FORALL
  
     
     ! check divergence 
     call DIV( A, Z ,R, OPA, OPZ, OPR, UA, UZ, UR, DG09, dg01, dg02 , dg03 )
     is = get_is()
     ie = get_ie()
     DIV_MAX = MAXVAL(ABS(DG09(IS(1):IE(1),IS(2):IE(2),IS(3):IE(3))))
     CALL MPI_ALLREDUCE(MPI_IN_PLACE,DIV_MAX,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,IERR)
     
     endtime   = MPI_Wtime();
     endtime =  endtime-starttime
     CALL MPI_ALLREDUCE(MPI_IN_PLACE,endtime,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,IERR)
     
     !Calcul de J(u) = int_domaine <u,u> + <t,t>
     DG04 = UA*UA*R
     DG05 = UZ*UZ*R
     DG06 = UR*UR*R


     CALL integrate_spec(quad,DG04,J_U,PH,NA,NZ,NR,xmax,xmin)
     CALL integrate_spec(quad,DG05,integ,PH,NA,NZ,NR,xmax,xmin)
     J_U = integ + J_U
     CALL integrate_spec(quad,DG06,integ,PH,NA,NZ,NR,xmax,xmin)
     J_U = integ + J_U

     
     call GetCFL(msh(1),msh(2),msh(3), UA, UZ, UR, dt, cfl)
     if (rank==0) print'(i9,11(1x,e15.8))',it_time,tc,dt,cfl,DIV_MAX,endtime,J_U

     if (do_adj .AND. MOD(IT_TIME,2)==0) then
        if (memoire) then
           SAVE_UA(it_time/2,:,:,:) = UA(:,:,:)
           SAVE_UZ(it_time/2,:,:,:) = UZ(:,:,:)
           SAVE_UR(it_time/2,:,:,:) = UR(:,:,:)
        else
           write(num,'(I6.6)')it_time
           filename = trim(base_save)//num//'.h5'
           call save_hdf5(trim(FILENAME),MSH,UA,UZ,UR)
        end if
     end if

     if (nrank==0) write(42,*)dt,-J_U


     if (cfl .GT. 10.) then
        if (rank == 0) print'("CFL TOO BIG.")'
        call MPI_FINALIZE(ierr)
        stop
     end if

     if (tc>=tmax) exit
     
  end DO

  CALL DUMP_HDF5_BASIC(file_dump,'WRITE',TC,DT,msh,UA,UZ,UR,Pres,DG01,dg09)
  
  UA_T = UA
  UZ_T = UZ
  UR_T = UR

  if (do_adj) then

     if (rank == 0) then
        write(6,'("--------------------------------------")')
        write(6,'("          Adjoint looping")')
        write(6,'("--------------------------------------")')
     end if

     NU_MOMENTUM = nu*[1,1,1]!*0.5   ! changment de signe pour le problème adjoint
     NU_POISSON  = 1.
     dt = - dt                      ! DT < 0 pour remonter vers U0_tilde


     CALL EQN_UA_AD%SET_PARAMS( NU=NU_MOMENTUM , SIGMA=0D0 , AXIS=1 )
     CALL EQN_UA_AD%SET_BCS( AXIS=3 , BCS_MINUS=DIRICHL , BCS_PLUS=DIRICHL )
     CALL EQN_UA_AD%INITIALISE_SVV( MSH, OPA, OPZ, OPR, ALPHA, BETA, PH )
     CALL EQN_UA_AD%SET_BVS( MESH = MSH , AXIS=3 , UDF_MINUS=UDF_NULL , UDF_PLUS=UDF_NULL )
     
     CALL EQN_UZ_AD%SET_PARAMS( NU=NU_MOMENTUM , SIGMA=0D0 , AXIS=2 )
     CALL EQN_UZ_AD%SET_BCS( AXIS=3 , BCS_MINUS=DIRICHL , BCS_PLUS=DIRICHL )
     CALL EQN_UZ_AD%INITIALISE_SVV(MSH,OPA,OPZ,OPR,ALPHA,BETA,PH)
     CALL EQN_UZ_AD%SET_BVS( MESH = MSH , AXIS=3 , UDF_MINUS=UDF_NULL , UDF_PLUS=UDF_NULL )
  
     CALL EQN_UR_AD%SET_PARAMS( NU=NU_MOMENTUM , SIGMA=0D0 , AXIS=3)
     CALL EQN_UR_AD%SET_BCS( AXIS=3 , BCS_MINUS=DIRICHL , BCS_PLUS=DIRICHL )
     CALL EQN_UR_AD%INITIALISE_SVV(MSH,OPA,OPZ,OPR,ALPHA,BETA,PH)
     CALL EQN_UR_AD%SET_BVS( MESH = MSH , AXIS=3 , UDF_MINUS=UDF_NULL , UDF_PLUS=UDF_NULL )

     CALL EQN_FI_AD%SET_PARAMS( NU=NU_POISSON , SIGMA=0D0 , AXIS=0)
     CALL EQN_FI_AD%SET_BCS( AXIS=3 , BCS_MINUS=NEUMANN , BCS_PLUS=NEUMANN )
     CALL EQN_FI_AD%INITIALISE_SVV(MSH,OPA,OPZ,OPR,ALPHA,BETA,PH)
     CALL EQN_FI_AD%SET_BVS( MESH = MSH , AXIS=3 , UDF_MINUS=UDF_NULL , UDF_PLUS=UDF_NULL )

  
     !Restart
     
     !  TC = 0._DP

     UA = 0._DP
     UZ = 0._DP
     UR = 0._DP

     UAM1 = UA
     UZM1 = UZ
     URM1 = UR

     NLAM1 = 0._DP
     NLZM1 = 0._DP
     NLRM1 = 0._DP
     
     PRES = 0._DP
     if (memoire) then
        DG04(:,:,:) = SAVE_UA(nb_iter/2,:,:,:)
        DG05(:,:,:) = SAVE_UZ(nb_iter/2,:,:,:)
        DG06(:,:,:) = SAVE_UR(nb_iter/2,:,:,:)
     else
        write(num,'(I6.6)')nb_iter
        filename = trim(base_save)//num//'.h5'
        CALL IMPORT_HDF5(FILENAME,msh,DG04,DG05,DG06)
     end if

     UAM1=UA
     UZM1=UZ
     URM1=UR
     
     CALL dealiazing(ua,uz,ur)
     CALL dealiazing(DG04,DG05,DG06)
     
     DG10 = DG04 + PRM_A*R + PRM_B/R
     
     CALL COMPUTE_ADJOINT_NON_LINEAR_TERMS_CYL(&
          A, Z, R, OPA, OPZ, OPR, DG10, DG05, DG06, UA, UZ, UR, NLAM1, NLZM1, NLRM1,&
          DG01, DG02, DG03, DG07, DG08, DG09)

     DG10 = PRM_A*R + PRM_B/R
     DG14 = 0._DP
     
     CALL COMPUTE_ADJOINT_NON_LINEAR_TERMS_CYL(&
          A, Z, R, OPA, OPZ, OPR, DG10, DG14, DG14, UA, UZ, UR, DG11, DG12, DG13,&
          DG01, DG02, DG03, DG07, DG08, DG09)

     NLAM1 = NLAM1 - DG11
     NLZM1 = NLZM1 - DG12
     NLRM1 = NLRM1 - DG13
     

     DO IT_time=1,nb_iter
        
        tc = tc + dt
        snap_dt = snap_dt + dt
        
        starttime = MPI_Wtime();
        
        ! Récuperation de U(t),V(t) et W(t)
        if (memoire) then
           CALL INTERPOLATION(2,SAVE_UA(FLOOR((nb_iter-it_time)/2._DP)+1,:,:,:), SAVE_UA(CEILING((nb_iter-it_time)/2._DP)+1,:,:,:), MOD(it_time+1,2), DG04)
           CALL INTERPOLATION(2,SAVE_UZ(FLOOR((nb_iter-it_time)/2._DP)+1,:,:,:), SAVE_UZ(CEILING((nb_iter-it_time)/2._DP)+1,:,:,:), MOD(it_time+1,2), DG05)
           CALL INTERPOLATION(2,SAVE_UR(FLOOR((nb_iter-it_time)/2._DP)+1,:,:,:), SAVE_UR(CEILING((nb_iter-it_time)/2._DP)+1,:,:,:), MOD(it_time+1,2), DG06)
        end if
        
        
        CALL GRAD(A,Z,R, &
             OPA, OPZ, OPR, PRES, DG01,DG02, DG03)
        
        SA = (2._DP*UA-0.5_DP*UAM1)/DT - DG01 + 2._DP*DG04
        SZ = (2._DP*UZ-0.5_DP*UZM1)/DT - DG02 + 2._DP*DG05
        SR = (2._DP*UR-0.5_DP*URM1)/DT - DG03 + 2._DP*DG06

        UAM1=UA
        UZM1=UZ
        URM1=UR
        
        CALL dealiazing(ua,uz,ur)
        CALL dealiazing(DG04,DG05,DG06)
        
        DG10 = DG04 + PRM_A*R + PRM_B/R
        
        CALL COMPUTE_ADJOINT_NON_LINEAR_TERMS_CYL(&
             A, Z, R, OPA, OPZ, OPR, DG10, DG05, DG06, UA, UZ, UR , NLA, NLZ, NLR,&
             DG01, DG02, DG03, DG07, DG08, DG09)

        DG10 = PRM_A*R + PRM_B/R
        DG14 = 0._DP
        
        CALL COMPUTE_ADJOINT_NON_LINEAR_TERMS_CYL(&
             A, Z, R, OPA, OPZ, OPR, DG10, DG14, DG14, UA, UZ, UR , DG11, DG12, DG13,&
             DG01, DG02, DG03, DG07, DG08, DG09)


        
 !       DG10 = PRM_A*R + PRM_B/R
     
 !       NLA = NLA + 2._DP*DG10*UA/R
    
        CALL OPA%D1(UA,SFI)
        DG01 = -NU*(-2._DP/R**2)*SFI
        CALL OPA%D1(UR,SFI)
        DG02 = -NU*(+2._DP/R**2)*SFI
     
        NLR = NLR - DG01 - DG13
        NLA = NLA - DG02 - DG11
        NLZ = NLZ - DG12
        
        SA = SA - 2._DP*NLA + NLAM1 
        SZ = SZ - 2._DP*NLZ + NLZM1 
        SR = SR - 2._DP*NLR + NLRM1
     
        NLAM1 = NLA
        NLZM1 = NLZ
        NLRM1 = NLR
                       
        SIGMA = 1.5_DP/DT
        
        CALL EQN_UA_AD%SOLVE(UA, SA, SIGMA ,PH)
        CALL EQN_UZ_AD%SOLVE(UZ, SZ, SIGMA ,PH)
        CALL EQN_UR_AD%SOLVE(UR, SR, SIGMA ,PH)
        
        
        call DIV( A, Z, R, OPA, OPZ, OPR, UA, UZ, UR, SFI, dg01, dg02 , dg03 )
        SFI = SFI*1.5_DP/DT
        
        call EQN_FI_AD%SOLVE(FI,SFI,0._DP,PH)
        
        CALL GRAD(A,Z,R,OPA, OPZ, OPR, FI, DG01, DG02, DG03)
        
        PRES = PRES + FI
        
        
        is = get_is_b([0,0,0])
        ie = get_ie_b([0,0,0])
        FORALL(I=IS(1):IE(1),J=IS(2):IE(2),K=IS(3):IE(3))
           UA(I,J,K) = UA(I,J,K) - DG01(I,J,K) * 2._DP*DT/3._DP 
        END FORALL
        is = get_is_b([0,0,0])
        ie = get_ie_b([0,0,0])
        FORALL(I=IS(1):IE(1),J=IS(2):IE(2),K=IS(3):IE(3))
           UZ(I,J,K) = UZ(I,J,K) - DG02(I,J,K) * 2._DP*DT/3._DP 
        END FORALL
        is = get_is_b([0,0,0])
        ie = get_ie_b([0,0,0])
        FORALL(I=IS(1):IE(1),J=IS(2):IE(2),K=IS(3):IE(3))
           UR(I,J,K) = UR(I,J,K) - DG03(I,J,K) * 2._DP*DT/3._DP 
        END FORALL
     
     
        ! check divergence 
        call DIV( A,Z,R,OPA, OPZ, OPR, UA, UZ, UR, DG09, dg01, dg02 , dg03 )
        is = get_is()
        ie = get_ie()
        DIV_MAX = MAXVAL(ABS(DG09(IS(1):IE(1),IS(2):IE(2),IS(3):IE(3))))     
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,DIV_MAX,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,IERR)
     
        endtime   = MPI_Wtime();
        endtime =  endtime-starttime
        CALL MPI_ALLREDUCE(MPI_IN_PLACE,endtime,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,IERR)
        

        call GetCFL(msh(1),msh(2),msh(3), UA, UZ, UR, dt, cfl)
        if (rank==0) print'(i9,10(1x,e15.8))',it_time,tc,dt,cfl,DIV_MAX,endtime

        if (cfl .GT. 10.) exit

     
     end DO

     !projection du gradient pour divegence nulle.
     call DIV( A, Z, R, OPA, OPZ, OPR, UA, UZ, UR, SFI, dg01, dg02 , dg03 )
     
     call EQN_FI_AD%SOLVE(FI,SFI,0._DP,PH)
     
     CALL GRAD(A,Z,R,OPA, OPZ, OPR, FI, DG01, DG02, DG03)

     UA = UA - DG01
     UZ = UZ - DG02
     UR = UR - DG03

     call DIV( A,Z,R,OPA, OPZ, OPR, UA, UZ, UR, DG09, dg01, dg02 , dg03 )
     is = get_is()
     ie = get_ie()
     DIV_MAX = MAXVAL(ABS(DG09(IS(1):IE(1),IS(2):IE(2),IS(3):IE(3))))     
     CALL MPI_ALLREDUCE(MPI_IN_PLACE,DIV_MAX,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,IERR)
     if (nrank==0) print*,"Div Grad J = ",div_max 

     !Ecriture du gradient
     FILENAME = TRIM(root_dir)//'grad.h5'
     call save_hdf5(trim(FILENAME),MSH,UA,UZ,UR)
     
     CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

  END if

  CALL MPI_FINALIZE(ierr)

  STOP

  
contains

  subroutine interpolation(order,U_min,U_sup,reste,RES)
    implicit none
    REAL(DP),DIMENSION(:,:,:),intent(IN)  :: U_min, U_sup
    INTEGER :: reste,order
    REAL(DP),DIMENSION(:,:,:),intent(OUT) :: RES

    
    RES = ((order-reste)*U_min+reste*U_sup)/REAL(order)

  END subroutine interpolation

  
  function get_is_b(halo) result(res)
        implicit none
        integer halo(3)
        integer res(3)
        
        res = ph%xst

        if (ph%xst(1) == 1) res(1) = res(1)+halo(1)
        if (ph%xst(2) == 1) res(2) = res(2)+halo(2)
        if (ph%xst(3) == 1) res(3) = res(3)+halo(3)
                
      end function get_is_b

      function get_ie_b(halo) result(res)
        implicit none
        integer halo(3)
        integer res(3)
        res = ph%xen
        if (ph%xen(1) == n(1)+1) res(1) = res(1)-halo(1)
        if (ph%xen(2) == n(2)+1) res(2) = res(2)-halo(2)
        if (ph%xen(3) == n(3)+1) res(3) = res(3)-halo(3)
        
      end function get_ie_b
      

  
        subroutine GetCFL_VAR(msh_x,msh_y,msh_z, Ux, Uy, Uz, dt, cfl, cfl_target)
        implicit none
        TYPE(T_MESH_base) :: msh_x,msh_y,msh_z
        REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: Ux, Uy, Uz
        REAL(KIND=8) :: dt,CFL,aa,dx,dy,dz,cfl_target,dt_update
        INTEGER :: I,J,K,IS(3),IE(3),ierr
        !is = get_is()
        !ie = get_ie()
        is = get_is_b([1,1,1])
        ie = get_ie_b([1,1,1])
        CFL = 0
        DO K=is(3),ie(3)
           DO J=is(2),ie(2)
              DO I=is(1),ie(1)
                 DX = (MSH_X%X(I+1)-MSH_X%X(I-1))*0.5
                 DY = (MSH_Y%X(J+1)-MSH_Y%X(J-1))*0.5
                 DZ = (MSH_Z%X(K+1)-MSH_Z%X(K-1))*0.5
                 AA = 0
                 AA = AA + ABS(UX(I,J,K))/DX
                 AA = AA + ABS(UY(I,J,K))/DY
                 AA = AA + ABS(UZ(I,J,K))/DZ
                 CFL = MAX(CFL,AA)
              END DO
           END DO
        END DO
        CFL = CFL * DT

        CALL MPI_ALLREDUCE(MPI_IN_PLACE,CFL,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,IERR)
        ! > update de la cfl

        dt_update = dt*CFL_target/CFL
        
        if (dt_update>2*dt) then
           dt = 1.1*dt
        else 
           if ( abs(dt-dt_update)/dt > 1e-1 )  then
              dt = dt_update
           else
              
           end if
        end if
        


        
      END subroutine GetCFL_VAR

      

      subroutine update_bcs(tc)
        implicit none
        real(kind=8) ::tc
        
!!$        CALL EQN_U%SET_BVS_TIME(TC, MSH=MSH , AXIS=3 , UDF_MINUS=UDF_NULL ,  UDF_PLUS=UDF_NULL )
!!$        CALL EQN_V%SET_BVS_TIME(TC, MSH=MSH , AXIS=3 , UDF_MINUS=UDF_NULL ,  UDF_PLUS=UDF_NULL )
!!$        CALL EQN_W%SET_BVS_TIME(TC, MSH=MSH , AXIS=3 , UDF_MINUS=UDF_NULL ,  UDF_PLUS=UDF_NULL )
        
      end subroutine update_bcs
      
  
  function get_is() result(res)
        implicit none
        integer res(3)
        res = ph%xst
        if (ph%xst(1) == 1) res(1) = res(1)+1
        if (ph%xst(2) == 1) res(2) = res(2)+1
        if (ph%xst(3) == 1) res(3) = res(3)+1
 !       return res
      end function get_is

      function get_ie() result(res)
        implicit none
        integer res(3)
        res = ph%xen
        if (ph%xen(1) == n(1)+1) res(1) = res(1)-1
        if (ph%xen(2) == n(2)+1) res(2) = res(2)-1
        if (ph%xen(3) == n(3)+1) res(3) = res(3)-1
!        return res
      end function get_ie


   subroutine GetCFL(msh_x,msh_y,msh_z, Ux, Uy, Uz, dt, cfl)
    implicit none
    TYPE(T_MESH_base) :: msh_x,msh_y,msh_z
    REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: Ux, Uy, Uz
    REAL(KIND=8) :: dt,CFL,aa,dx,dy,dz
    INTEGER :: I,J,K,IS(3),IE(3),ierr
    is = get_is()
    ie = get_ie()
    CFL = 0
    DO K=is(3),ie(3)
       DO J=is(2),ie(2)
          DO I=is(1),ie(1)
!             DX = (MSH_X%X(is(1)+1)-MSH_X%X(is(1)))
!             DY = (MSH_Y%X(is(2)+1)-MSH_Y%X(is(2)))
             DX = (MSH_X%X(i+1)-MSH_X%X(i))
             DY = (MSH_Y%X(j+1)-MSH_Y%X(j))
             DZ = (MSH_Z%X(K+1)-MSH_Z%X(K))
             AA = 0
             AA = AA + ABS(UX(I,J,K))/DX
             AA = AA + ABS(UY(I,J,K))/DY
             AA = AA + ABS(UZ(I,J,K))/DZ
             AA = abs(AA*dt)
             CFL = MAX(CFL,AA)
          END DO
       END DO
    END DO
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,CFL,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,IERR)

  END subroutine GetCFL


  
  subroutine preproc()
    implicit none
    real(kind=8) ::coeff_p(3)

    call alloc_x(DGA , OPT_GLOBAL=.TRUE.) ; DGA = 0._DP
    call alloc_x(DGZ , OPT_GLOBAL=.TRUE.) ; DGZ = 0._DP
    call alloc_x(DGR , OPT_GLOBAL=.TRUE.) ; DGR = 0._DP

    call alloc_y(DGA_Y , OPT_GLOBAL=.TRUE.) ; DGA_Y = 0._DP
    call alloc_y(DGZ_Y , OPT_GLOBAL=.TRUE.) ; DGZ_Y = 0._DP
    call alloc_y(DGR_Y , OPT_GLOBAL=.TRUE.) ; DGR_Y = 0._DP
   
    call alloc_x(DG01 , OPT_GLOBAL=.TRUE.) ; DG01 = 0
    call alloc_x(DG02 , OPT_GLOBAL=.TRUE.) ; DG02 = 0
    call alloc_x(DG03 , OPT_GLOBAL=.TRUE.) ; DG03 = 0
    call alloc_x(DG04 , OPT_GLOBAL=.TRUE.) ; DG04 = 0
    call alloc_x(DG05 , OPT_GLOBAL=.TRUE.) ; DG05 = 0
    call alloc_x(DG06 , OPT_GLOBAL=.TRUE.) ; DG06 = 0
    call alloc_x(DG07 , OPT_GLOBAL=.TRUE.) ; DG07 = 0
    call alloc_x(DG08 , OPT_GLOBAL=.TRUE.) ; DG08 = 0
    call alloc_x(DG09 , OPT_GLOBAL=.TRUE.) ; DG09 = 0
    call alloc_x(DG10 , OPT_GLOBAL=.TRUE.) ; DG10 = 0
    call alloc_x(DG11 , OPT_GLOBAL=.TRUE.) ; DG11 = 0
    call alloc_x(DG12 , OPT_GLOBAL=.TRUE.) ; DG12 = 0
    call alloc_x(DG13 , OPT_GLOBAL=.TRUE.) ; DG13 = 0
    call alloc_x(DG14 , OPT_GLOBAL=.TRUE.) ; DG14 = 0 
    
    call alloc_x(UA_T   , OPT_GLOBAL=.TRUE.) ; UA_T = 0 
    call alloc_x(UZ_T   , OPT_GLOBAL=.TRUE.) ; UZ_T = 0
    call alloc_x(UR_T   , OPT_GLOBAL=.TRUE.) ; UR_T = 0
    
    call alloc_x(PRES , OPT_GLOBAL=.TRUE.) ; PRES = 0
    call alloc_x(FI   , OPT_GLOBAL=.TRUE.) ; FI = 0
    call alloc_x(SFI  , OPT_GLOBAL=.TRUE.) ; SFI = 0
    !
    call alloc_x(UA   , OPT_GLOBAL=.TRUE.) ; UA = 0
    call alloc_x(SA  , OPT_GLOBAL=.TRUE.) ; SA = 0
    call alloc_x(NLA , OPT_GLOBAL=.TRUE.) ; NLA = 0
    !
    call alloc_x(UZ   , OPT_GLOBAL=.TRUE.) ; UZ = 0
    call alloc_x(SZ  , OPT_GLOBAL=.TRUE.) ; SZ = 0
    call alloc_x(NLZ , OPT_GLOBAL=.TRUE.) ; NLZ = 0
    !
    call alloc_x(UR   , OPT_GLOBAL=.TRUE.) ; UR = 0
    call alloc_x(SR  , OPT_GLOBAL=.TRUE.) ; SR = 0
    call alloc_x(NLR , OPT_GLOBAL=.TRUE.) ; NLR = 0
  
    call alloc_x(A  , OPT_GLOBAL=.TRUE.) ; A = 0
    call alloc_x(Z  , OPT_GLOBAL=.TRUE.) ; Z = 0
    call alloc_x(R  , OPT_GLOBAL=.TRUE.) ; R = 0

    call alloc_x(UAM1  , OPT_GLOBAL=.TRUE.) ; UAM1 = 0._DP
    call alloc_x(UZM1  , OPT_GLOBAL=.TRUE.) ; UZM1 = 0._DP
    call alloc_x(URM1  , OPT_GLOBAL=.TRUE.) ; URM1 = 0._DP

    call alloc_x(NLAM1  , OPT_GLOBAL=.TRUE.) ; NLAM1 = 0._DP
    call alloc_x(NLZM1  , OPT_GLOBAL=.TRUE.) ; NLZM1 = 0._DP
    call alloc_x(NLRM1  , OPT_GLOBAL=.TRUE.) ; NLRM1 = 0._DP

    ALLOCATE(A_tot(1:NA+1,1:NZ+1,1:NR+1))
    ALLOCATE(Z_tot(1:NA+1,1:NZ+1,1:NR+1))
    ALLOCATE(R_tot(1:NA+1,1:NZ+1,1:NR+1))


    
    FORALL(I=ph%XST(1):ph%XEN(1),J=ph%XST(2):ph%XEN(2),K=ph%XST(3):ph%XEN(3))
       A(I,J,K) = MSH(1)%X(I)
       Z(I,J,K) = MSH(2)%X(J)
       R(I,J,K) = MSH(3)%X(K)
    END FORALL

    FORALL(I=1:NA+1,J=1:NZ+1,K=1:NR+1)
       A_tot(i,j,k) = MSH(1)%X(I)
       Z_tot(i,j,k) = MSH(2)%X(J)
       R_tot(i,j,k) = MSH(3)%X(K)
    END FORALL

    
  end subroutine preproc

  subroutine source_term(A,Z,R,PRES,nu,prm_K,F_A,F_Z,F_R)
    REAL(DP),DIMENSION(:,:,:),ALLOCATABLE,INTENT(IN) :: A,Z,R,PRES
    REAL(DP),INTENT(IN) :: nu,prm_k

    REAL(DP),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT) :: F_A, F_Z, F_R

    CALL GRAD(A,Z,R,OPA,OPZ,OPR,PRES,F_A,F_Z,F_R)

!    F_A = 0._DP
!    F_Z = 0._DP
!    F_R = 0._DP
    
    FORALL(I=PH%XST(1):PH%XEN(1),J=PH%XST(2):PH%XEN(2),K=PH%XST(3):PH%XEN(3))
       F_R(I,J,K) = F_R(I,J,K) - (prm_K/R(I,J,K))**2
       F_A(I,J,K) = F_A(I,J,K) + nu*0.75_DP*prm_K*R(I,J,K)**(-2.5_DP)
    END FORALL

  END subroutine source_term
   
  real(kind=8) pure function udf_one(t,x,y,z)
    real(kind=8),intent(in):: t,x,y,z
    udf_one = 1.0
  end function udf_one
  
  
  real(kind=8) pure function udf_null(t,x,y,z)
    real(kind=8),intent(in):: t,x,y,z
    udf_null = 0.
  end function udf_null

  subroutine export_tecplot(FILENAME,X,Y,Z,FI)
    implicit none
    character(len=*) :: filename
    real(kind=8),dimension(:),allocatable :: x,y,z
    real(kind=8),dimension(:,:,:),allocatable :: fi
    integer :: i,j,k,id,rank
    integer :: cpt_line
    integer :: line_max=100
    integer :: is(3),ie(3),n(3)
    CHARACTER(LEN=3) :: RX
    
    call mpi_comm_rank(mpi_comm_world,rank,ierr)
    
    is = lbound(FI)
    ie = ubound(FI)
    n = ie-is+1
    ID = 24
    WRITE (RX,'(I3.3)')RANK
    id=1000+rank
    OPEN (UNIT=ID,FILE=TRIM(FILENAME)//TRIM(RX)//'.dat')
    WRITE(ID,*)'VARIABLES= X,Y,Z,FI'
    WRITE(ID,*)'ZONE  I=',N(1),', J=',N(2),',K=',N(3),&
         ', DATAPACKING=BLOCK'
    cpt_line=0
    do k=is(3),ie(3)
       do j=is(2),ie(2)
          do i=is(1),ie(1)
             write(ID,'(1x,e15.8,1x,$)')x(I)
             cpt_line = cpt_line + 1
             if (mod(cpt_line,line_max)==0)  write(ID,*)
          end do
       end do
    end do
    cpt_line=0
    do k=is(3),ie(3)
       do j=is(2),ie(2)
          do i=is(1),ie(1)
             write(ID,'(1x,e15.8,1x,$)')Y(J)
             cpt_line = cpt_line + 1
             if (mod(cpt_line,line_max)==0)  write(ID,*)
          end do
       end do
    end do
    cpt_line=0
    do k=is(3),ie(3)
       do j=is(2),ie(2)
          do i=is(1),ie(1)
             write(ID,'(1x,e15.8,1x,$)')Z(K)
             cpt_line = cpt_line + 1
             if (mod(cpt_line,line_max)==0)  write(ID,*)
          end do
       end do
    end do

    cpt_line=0
    do k=is(3),ie(3)
       do j=is(2),ie(2)
          do i=is(1),ie(1)
             write(ID,'(1x,e15.8,1x,$)')FI(I,J,K)
             cpt_line = cpt_line + 1
             if (mod(cpt_line,line_max)==0)  write(ID,*)
          end do
       end do
    end do
    
    close(id)
  end subroutine export_tecplot

  subroutine dealiazing(ua,uz,ur)
    REAL(DP),ALLOCATABLE,DIMENSION(:,:,:) :: UA,UZ,UR
    integer :: i,k_max,k_cut_p,k_cut_m

    k_max = NA/3
    k_cut_p = k_max +2
    k_cut_m = NA - k_max + 1


    DGA = UA
    DGZ = UZ
    DGR = UR

    
    call c2c_1m_x(DGA,plan_fwd_x)
    call c2c_1m_x(DGZ,plan_fwd_x)
    call c2c_1m_x(DGR,plan_fwd_x)

    DGA(k_cut_p : k_cut_m,:,:) = CMPLX(0._DP,0._DP)
    DGZ(k_cut_p : k_cut_m,:,:) = CMPLX(0._DP,0._DP)
    DGR(k_cut_p : k_cut_m,:,:) = CMPLX(0._DP,0._DP)

    call transpose_x_to_y(DGA, DGA_Y)
    call transpose_x_to_y(DGZ, DGZ_Y)
    call transpose_x_to_y(DGR, DGR_Y)

    k_max = NZ/3
    k_cut_p = k_max +2
    k_cut_m = NZ - k_max + 1

    
    call c2c_1m_y(DGA_Y,plan_fwd_y)
    call c2c_1m_y(DGZ_Y,plan_fwd_y)
    call c2c_1m_y(DGR_Y,plan_fwd_y)

    DGA_Y(:,k_cut_p : k_cut_m,:) = CMPLX(0._DP,0._DP)
    DGZ_Y(:,k_cut_p : k_cut_m,:) = CMPLX(0._DP,0._DP)
    DGR_Y(:,k_cut_p : k_cut_m,:) = CMPLX(0._DP,0._DP)
    
    call c2c_1m_y(DGA_Y,plan_bck_y)
    call c2c_1m_y(DGZ_Y,plan_bck_y)
    call c2c_1m_y(DGR_Y,plan_bck_y)

    call transpose_y_to_x(DGA_Y, DGA)
    call transpose_y_to_x(DGZ_Y, DGZ)
    call transpose_y_to_x(DGR_Y, DGR)

    call c2c_1m_x(DGA,plan_bck_x)
    call c2c_1m_x(DGZ,plan_bck_x)
    call c2c_1m_x(DGR,plan_bck_x)
    
    UA = DGA
    UZ = DGZ
    UR = DGR
    
    
  end subroutine dealiazing

  subroutine NORMALIZATION(UA,UZ,UR,E0)
    REAL(DP),ALLOCATABLE,DIMENSION(:,:,:) :: UA,UZ,UR
    REAL(DP) :: E0, NRJ, vol

    DG01 = UA * UA * R
    CALL integrate_spec(quad,DG01,NRJ,PH,NA,NZ,NR,xmax,xmin)
    DG01 = UZ * UZ * R
    CALL integrate_spec(quad,DG01,VOL,PH,NA,NZ,NR,xmax,xmin)
    NRJ = NRJ + VOL
    DG01 = UR * UR * R
    CALL integrate_spec(quad,DG01,VOL,PH,NA,NZ,NR,xmax,xmin)
    NRJ = NRJ + VOL  
    UA = E0 * UA / NRJ
    UZ = E0 * UZ / NRJ
    UR = E0 * UR / NRJ
  end subroutine NORMALIZATION
  
end program tcheby_1d
