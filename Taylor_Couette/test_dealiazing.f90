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
  type(t_solver_diag_cyl_hhi) :: eqn_ua, eqn_uz, eqn_ur
  type(t_solver_diag_cyl_hhi) :: eqn_fi
  
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

  real(kind=8),allocatable ,dimension(:,:,:) ::  UAM1,UZM1,URM1
  real(kind=8),allocatable ,dimension(:,:,:) ::  NLAM1,NLZM1,NLRM1
  real(kind=8),allocatable ,dimension(:,:,:) ::  UA,UZ,UR,PRES,FI
  real(kind=8),allocatable ,dimension(:,:,:) ::  SA,SZ,SR,SFI
  real(kind=8),allocatable ,dimension(:,:,:) ::  NLA,NLZ,NLR
  real(kind=8),allocatable ,dimension(:,:,:) ::  A,Z,R
  real(kind=8) :: w2,t1,h
  integer :: it_time,is(3),ie(3)
  integer , parameter :: OX=1,OY=2,OZ=3
  REAL(DP) :: prm_A,prm_B,eta,OMEGA_i,OMEGA_o
  
  ! attention au type à lire
  REAL(kind=8) :: a_min,a_max,z_min,z_max,r_min,r_max
  integer :: na,nz,nr,nb_cpu_y,nb_cpu_z
  namelist /parameters_cube/ na,nz,nr,z_min,z_max,r_min,r_max,nb_cpu_y,nb_cpu_z
  
  REAL(KIND=8) :: NU,RE
  INTEGER :: nb_iter
  namelist /parameters_physical/RE

  character(len=1024) ::  root_dir,io_nrj,output_dir
  REAL(kind=8) ::  dt,cfl,tmax
  integer dn_dump
  logical :: resume , explicit,correction_pres
  namelist /parameters_timescheme/ tmax,dt, cfl, dn_dump,root_dir,io_nrj
  
  logical :: resume_snap=.false.
  integer :: dn_snap
  character(len=1024) ::  snap_dir,init_file
  real(kind=8) :: snap_dt,starttime,endtime
  integer :: snap_id
  namelist /parameters_diagnostics/ dn_snap,snap_dir,init_file

  integer :: na0,nz0,nr0
  REAL(KIND=8), dimension(:,:,:),ALLOCATABLE :: UA0, UZ0, UR0
  namelist /parameters_initial_condition/ na0,nz0,nr0
  
  character(len=1024) ::  input_file
  integer :: rank,iostat,seed_size
  integer, allocatable :: myseed(:)
  character(len=1024) ::  file_dump, base_snap, filename, base_save, vort_snap
  character(len=6) :: num
  character(len=6) :: form
  REAL(kind=8) ::  DIV_MAX,tmp,DJ_ADJ,J_DU0,J_U
  
  REAL(KIND=8),DIMENSION(3) :: NU_MOMENTUM
  REAL(KIND=8),DIMENSION(3) :: NU_ENERGY
  REAL(KIND=8),DIMENSION(3) :: NU_POISSON
  REAL(kind=8) :: sigma,omega_tilde,noise,res,theta,signal_attendu

  LOGICAL :: memoire = .TRUE.
  LOGICAL :: do_adj = .FALSE.

  COMPLEX(DP),ALLOCATABLE,DIMENSION(:,:,:) :: DGA,DGZ,DGR
  COMPLEX(DP),ALLOCATABLE,DIMENSION(:,:,:) :: DGA_Y,DGZ_Y,DGR_Y

  REAL(DP),ALLOCATABLE,DIMENSION(:) :: E_azi, E_ver, sp_a, sp_z

  INTEGER,PARAMETER :: max_mod=33

  call mpi_init(ierr)
  CALL H5OPEN_F(IERR)

  call mpi_comm_rank(mpi_comm_world,rank,ierr)
  
  call command_line_read_input(input_file)
  
  if (rank==0) then
     OPEN (UNIT=24, FILE=TRIM(input_file),status='old', action='read')
     read(24, nml=parameters_cube, IOSTAT=iostat)
     read(24, nml=parameters_physical, IOSTAT=iostat)
     read(24, nml=parameters_timescheme, IOSTAT=iostat)
     read(24, nml=parameters_diagnostics, IOSTAT=iostat)
 !    read(24, nml=parameters_initial_condition, IOSTAT=iostat)
     close(24)
     ierr = SYSTEM('mkdir -p '//trim(root_dir)//'/dump' )
     ierr = SYSTEM('mkdir -p '//trim(root_dir)//'/save' )
     ierr = SYSTEM('mkdir -p '//trim(root_dir)//trim(snap_dir) )
     ierr = SYSTEM('mkdir -p '//trim(root_dir)//'/vorticity' )
     ierr = SYSTEM('mkdir -p '//trim(root_dir)//'/spectre' )
     ierr = SYSTEM('mkdir -p '//trim(root_dir)//'/modal_nrj_azimutal' )
     ierr = SYSTEM('mkdir -p '//trim(root_dir)//'/modal_nrj_vertical' )
     write(*,parameters_cube)
     write(*,parameters_physical)
     write(*,parameters_diagnostics)
     OPEN(UNIT=42, FILE=trim(root_dir)//'U_vect.dat')
     OPEN(UNIT=50, FILE=trim(root_dir)//'a_vect.dat')
     OPEN(UNIT=60, FILE=trim(root_dir)//'z_vect.dat')
     OPEN(UNIT=70, FILE=trim(root_dir)//'r_vect.dat')
     OPEN(UNIT=80, FILE=trim(root_dir)//'/spectre/azimut.dat')
     OPEN(UNIT=90, FILE=trim(root_dir)//'/spectre/vertic.dat')
     OPEN(UNIT=101, FILE=trim(root_dir)//'/modal_nrj_azimutal/ua.dat')
     OPEN(UNIT=102, FILE=trim(root_dir)//'/modal_nrj_azimutal/uz.dat')
     OPEN(UNIT=103, FILE=trim(root_dir)//'/modal_nrj_azimutal/ur.dat')
     OPEN(UNIT=111, FILE=trim(root_dir)//'/modal_nrj_vertical/ua.dat')
     OPEN(UNIT=112, FILE=trim(root_dir)//'/modal_nrj_vertical/uz.dat')
     OPEN(UNIT=113, FILE=trim(root_dir)//'/modal_nrj_vertical/ur.dat')

  end if

  
  CALL MPI_BCAST( na      , 1,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
  CALL MPI_BCAST( nz      , 1,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
  CALL MPI_BCAST( nr      , 1,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
!  CALL MPI_BCAST( na0     , 1,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
!  CALL MPI_BCAST( nz0     , 1,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
!  CALL MPI_BCAST( nr0     , 1,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
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
  
  CALL MPI_BCAST( dn_snap  , 1   , MPI_INTEGER         , 0, MPI_COMM_WORLD, IERR )
  CALL MPI_BCAST( snap_dir , 1024, MPI_character       , 0, MPI_COMM_WORLD, IERR )
  CALL MPI_BCAST( init_file, 1024, MPI_character       , 0, MPI_COMM_WORLD, IERR )

  
  N = [NA,NZ,NR]  
  CPU_GRID = [NB_CPU_Y,NB_CPU_Z]

  
  xmin(1:3) = [0._DP,z_min,r_min]
  xmax(1:3) = [2._DP*pi,z_max,r_max]

  eta = r_min/r_max

  OMEGA_I = 1._DP
  OMEGA_O = (r_max)**(-1.5_DP)
  
  prm_A = (r_max/(1._DP - eta**2)) * (OMEGA_o - OMEGA_i*eta**2)
  prm_B =(r_min**2/(1._DP - eta**2)) * (OMEGA_i - OMEGA_o) 

  
  nb_iter = floor(tmax/dt)
  
  file_dump = trim(root_dir)//'dump/dump.h5'
  base_snap = trim(root_dir)//trim(snap_dir)//'snap_'
  vort_snap = trim(root_dir)//'vorticity/snap_'
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
  CALL EQN_UA%INITIALISE( MSH, OPA, OPZ, OPR, PH )
  CALL EQN_UA%SET_BVS( MESH = MSH , AXIS=3 , UDF_MINUS=UDF_NULL , UDF_PLUS=UDF_NULL )
  
  CALL EQN_UZ%SET_PARAMS( NU=NU_MOMENTUM , SIGMA=0D0 , AXIS=2 )
  CALL EQN_UZ%SET_BCS( AXIS=3 , BCS_MINUS=DIRICHL , BCS_PLUS=DIRICHL )
  CALL EQN_UZ%INITIALISE(MSH,OPA,OPZ,OPR,PH)
  CALL EQN_UZ%SET_BVS( MESH = MSH , AXIS=3 , UDF_MINUS=UDF_NULL , UDF_PLUS=UDF_NULL )
  
  CALL EQN_UR%SET_PARAMS( NU=NU_MOMENTUM , SIGMA=0D0 , AXIS=3 )
  CALL EQN_UR%SET_BCS( AXIS=3 , BCS_MINUS=DIRICHL , BCS_PLUS=DIRICHL )
  CALL EQN_UR%INITIALISE(MSH,OPA,OPZ,OPR,PH)
  CALL EQN_UR%SET_BVS( MESH = MSH , AXIS=3 , UDF_MINUS=UDF_NULL , UDF_PLUS=UDF_NULL )

  CALL EQN_FI%SET_PARAMS( NU=NU_POISSON , SIGMA=0D0 , AXIS=0 )
  CALL EQN_FI%SET_BCS( AXIS=3 , BCS_MINUS=NEUMANN , BCS_PLUS=NEUMANN )
  CALL EQN_FI%INITIALISE(MSH,OPA,OPZ,OPR,PH)
  CALL EQN_FI%SET_BVS( MESH = MSH , AXIS=3 , UDF_MINUS=UDF_NULL , UDF_PLUS=UDF_NULL )
  
  
  call preproc()

  DGA = UA
  DGZ = UZ
  DGR = UR
  is = get_is_b([0,0,0])
  ie = get_ie_b([0,0,0])

  CALL RANDOM_SEED()

  do i = 2, (NA+1)/2 +1
     CALL RANDOM_NUMBER(err)
     DG01(i,:,:) = err
     DG01((NA+1)/2 +i ,:,:) = err

     CALL RANDOM_NUMBER(err)
     DG02(i,:,:) = err
     DG02((NA+1)/2 +i ,:,:) = err
  end do
!     CALL Random_Number(DG01(IS(1):IE(1),IS(2):IE(2),IS(3):IE(3)))
!     CALL Random_Number(DG02(IS(1):IE(1),IS(2):IE(2),IS(3):IE(3)))

!  DG01 = 1._DP
  DG02 = 1._DP
  
  
  DGA = CMPLX(DG01,0._DP)
  DGZ = CMPLX(DG02,0._DP)

  !aliasing
  call c2c_1m_x(DGA, plan_bck_x)
  call c2c_1m_x(DGZ, plan_bck_x)
  DGR = DGA*DGZ
  call c2c_1m_x(DGR, plan_fwd_x)
  call c2c_1m_x(DGA, plan_fwd_x)
  if (nrank == 0) then
     open(UNIT=300, FILE="Prod_base.dat")
     open(UNIT=301, FILE="U_base.dat")
     do i=is(1),ie(1)
        write(300,"(e15.8)")REAL(DGR(i,is(2),is(3)))
        write(301,"(e15.8)")REAL(DGA(i,is(2),is(3)))
     end do
  end if

  !Coupure 1/2
  DGA = CMPLX(DG01,0._DP)
  DGZ = CMPLX(DG02,0._DP)

  call dealiazing_spec(DGA,DGA,DGA,NA/4)
  call dealiazing_spec(DGZ,DGZ,DGZ,NA/4)


  
  call c2c_1m_x(DGA, plan_bck_x)
  call c2c_1m_x(DGZ, plan_bck_x)
  DGR = DGA*DGZ
  call c2c_1m_x(DGR, plan_fwd_x)
  call c2c_1m_x(DGA, plan_fwd_x)
  if (nrank == 0) then
     open(UNIT=400, FILE="Prod_05.dat")
     open(UNIT=401, FILE="U_05.dat")
     do i=is(1),ie(1)
        write(400,"(e15.8)")REAL(DGR(i,is(2),is(3)))
        write(401,"(e15.8)")REAL(DGA(i,is(2),is(3)))
     end do
  end if

  !Coupure 1/3
  DGA = CMPLX(DG01,0._DP)
  DGZ = CMPLX(DG02,0._DP)

  call dealiazing_spec(DGA,DGA,DGA,NA/3)
  call dealiazing_spec(DGZ,DGZ,DGZ,NA/3)

  call c2c_1m_x(DGA, plan_bck_x)
  call c2c_1m_x(DGZ, plan_bck_x)
  DGR = DGA*DGZ
  call c2c_1m_x(DGR, plan_fwd_x)
  call c2c_1m_x(DGA, plan_fwd_x)
  if (nrank == 0) then
     open(UNIT=500, FILE="Prod_03.dat")
     open(UNIT=501, FILE="U_03.dat")
     do i=is(1),ie(1)
        write(500,"(e15.8)")REAL(DGR(i,is(2),is(3)))
        write(501,"(e15.8)")REAL(DGA(i,is(2),is(3)))
     end do
  end if

  !Coupure 1/4
  DGA = CMPLX(DG01,0._DP)
  DGZ = CMPLX(DG02,0._DP)

  call dealiazing_spec(DGA,DGA,DGA,NA/8)
  call dealiazing_spec(DGZ,DGZ,DGZ,NA/8)

  call c2c_1m_x(DGA, plan_bck_x)
  call c2c_1m_x(DGZ, plan_bck_x)
  DGR = DGA*DGZ
  call c2c_1m_x(DGR, plan_fwd_x)
  if (nrank == 0) then
     open(UNIT=600, FILE="Prod_025.dat")
     do i=is(1),ie(1)
        write(600,"(e15.8)")REAL(DGR(i,is(2),is(3)))
     end do
  end if

  !Coupure 1/5
  DGA = CMPLX(DG01,0._DP)
  DGZ = CMPLX(DG02,0._DP)

  call dealiazing_spec(DGA,DGA,DGA,NA/10)
  call dealiazing_spec(DGZ,DGZ,DGZ,NA/10)

  call c2c_1m_x(DGA, plan_bck_x)
  call c2c_1m_x(DGZ, plan_bck_x)
  DGR = DGA*DGZ
  call c2c_1m_x(DGR, plan_fwd_x)
  if (nrank == 0) then
     open(UNIT=700, FILE="Prod_02.dat")
     do i=is(1),ie(1)
        write(700,"(e15.8)")REAL(DGR(i,is(2),is(3)))
     end do
  end if

  !Coupure 3/4
  DGA = CMPLX(DG01,0._DP)
  DGZ = CMPLX(DG02,0._DP)

  call dealiazing_spec(DGA,DGA,DGA,3*NA/8)
  call dealiazing_spec(DGZ,DGZ,DGZ,3*NA/8)

  call c2c_1m_x(DGA, plan_bck_x)
  call c2c_1m_x(DGZ, plan_bck_x)
  DGR = DGA*DGZ
  call c2c_1m_x(DGR, plan_fwd_x)
  call c2c_1m_x(DGA, plan_fwd_x)
  if (nrank == 0) then
     open(UNIT=800, FILE="Prod_075.dat")
     open(UNIT=801, FILE="U_075.dat")
     do i=is(1),ie(1)
        write(800,"(e15.8)")REAL(DGR(i,is(2),is(3)))
        write(801,"(e15.8)")REAL(DGA(i,is(2),is(3)))
     end do
  end if

  !Coupure 9/10
  DGA = CMPLX(DG01,0._DP)
  DGZ = CMPLX(DG02,0._DP)

  call dealiazing_spec(DGA,DGA,DGA,9*NA/20)
  call dealiazing_spec(DGZ,DGZ,DGZ,9*NA/20)

  call c2c_1m_x(DGA, plan_bck_x)
  call c2c_1m_x(DGZ, plan_bck_x)
  DGR = DGA*DGZ
  call c2c_1m_x(DGR, plan_fwd_x)
  call c2c_1m_x(DGA, plan_fwd_x)
  if (nrank == 0) then
     open(UNIT=900, FILE="Prod_09.dat")
     open(UNIT=901, FILE="U_09.dat")
     do i=is(1),ie(1)
        write(900,"(e15.8)")REAL(DGR(i,is(2),is(3)))
        write(901,"(e15.8)")REAL(DGA(i,is(2),is(3)))
     end do
  end if


  
  CALL MPI_FINALIZE(ierr)
  
  stop

  
  do i = is(1),ie(1)!1, (Na+1)/2 +1
     do j = is(2),ie(2)
        do k = is(3),ie(3)
           do it_time = 1,30
              !    DGA(i,:,:) = ((33._DP - i)/31._DP)!*(1 - 5E-1) + 5E-1
              !     DGA((NA+1)/2 + i,:,:) = ((NA+1)/2 + i - 34._DP)/31._DP !*(1 - 5E-1) + 5E-1
              DGA(I,J,K) = DGA(I,J,K) + 1._DP/it_time * cos(it_time*(A(I,J,K))) !+ cos(it_time*(Z(I,J,K))) 
           end do
        end do
     end do
  end do
  write(6,*)"U base = ",REAL(DGA(:,is(2),is(3)))
  !DGA (1,:,:) = 1._DP
!  DGA = CMPLX(0._DP,0._DP)
!  DGA(23:44,:,:) = DGA(23:44,:,:)*4._DP

  call c2c_1m_x(DGA, plan_bck_x)

  UA = DGA
  
  if (nrank == 0) then
!     write(6,*)"U base = ",REAL(DGA(:,is(2),is(3)))
     open(UNIT=200, FILE="Phy_nndea.dat")
     do i=is(1),ie(1)
        write(200,"(e15.8)")REAL(DGA(i,is(2),is(3)))
!        if (nrank == 0) write(6,*)i," ",DGA(i,is(2),is(3))
     end do
  end if

  DG01 = DGA
  
  CALL DEALIAZING(DG01,DG02,DG03)
  UZ = DG01
  if (nrank == 0) then
!     write(6,*)"U dea = ",DG01(:,is(2),is(3))
     open(UNIT=201, FILE="Phy_dea.dat")
     do i=is(1),ie(1)
        write(201,"(e15.8)")REAL(DG01(i,is(2),is(3)))
!        if (nrank == 0) write(6,*)i," ",DG01(i,is(2),is(3))
     end do
  end if


  CALL COMPUTE_NON_LINEAR_TERMS(&
       A, Z, R, OPA, OPZ, OPR, UA, UA, UA , NLAM1, NLZM1, NLRM1,&
       DG01, DG02, DG03, DG04, DG05, DG06, DG07, DG08, DG09)

  CALL COMPUTE_NON_LINEAR_TERMS(&
       A, Z, R, OPA, OPZ, OPR, UZ, UZ, UZ , NLA, NLZ, NLR,&
       DG01, DG02, DG03, DG04, DG05, DG06, DG07, DG08, DG09)

  if (nrank == 0) then
     open(UNIT=202, FILE="NL_base.dat")
     do i=is(1),ie(1)
        write(202,"(e15.8)")NLAM1(i,is(2),is(3))
     end do
  end if
  if (nrank == 0) then
     open(UNIT=203, FILE="NL_dea.dat")
     do i=is(1),ie(1)
        write(203,"(e15.8)")NLA(i,is(2),is(3))
     end do
  end if


  

  
  CALL MPI_FINALIZE(ierr)
  stop
  
  ! === TEST 2 : avec troncature ===
  do i = 1, NA+1
     div_max = 2._DP * PI * real(i-1, DP) / real(NA+1, DP)
     DGA(i,:,:) = cmplx( cos(2*div_max) + cos((NA/3+2)*div_max), 0._DP, DP)
  end do
  DG01 = DGA
  call c2c_1m_x(DGA, plan_fwd_x)
  DGA(NA/3+2 : NA/2+1,   :, :) = 0._DP
  DGA(NA/2+2 : NA-NA/3+1,:, :) = 0._DP
  call c2c_1m_x(DGA, plan_bck_x)
  
  if (nrank==0) then
     write(*,'(A)') "=== Apres troncature ==="
     write(*,'(A6,A14,A14,A14)') "i","obtenu","attendu","erreur"
     do i = 1, NA+1
        div_max = 2._DP * PI * real(i-1, DP) / real(NA+1, DP)
        signal_attendu = cos(2*div_max)   ! seulement le mode k=2
        write(*,'(I6,3F14.8)') i, real(DGA(i,1,1),DP),signal_attendu,real(DGA(i,1,1),DP) - signal_attendu
     end do
  end if

  DG03 = 0._DP
  FORALL(I=IS(1):IE(1),J=IS(2):IE(2),K=IS(3):IE(3))
     DG03(I,J,K) = ABS(DG01(I,J,K) - DGA(I,J,K))
  END FORALL

    print*,nrank," error : ",MAXVAL(DG03),MAXLOC(DG03)
    stop
  

  
  !Lecture de la condition initiale

  if (nrank==0) print*,"Lecture de la condition initiale"
  if (nrank==0) print*,TRIM(init_file)
  
  CALL IMPORT_HDF5_INIT(TRIM(init_file),DG01,DG02,DG03)

  if (nrank==0) print*,"Fin de lecture de la condition initiale"

  UA = DG01
  UZ = DG02
  UR = DG03

  CALL DEALIAZING(UA,UZ,UR)
  
  filename = trim(base_snap)//'grid.h5'
  call update_id_and_time(trim(filename),snap_dt,snap_id)
  WRITE(num,'(I6.6)')snap_id
  filename = trim(base_snap)//num//'.h5'
  call EXPORT_snapshot(trim(FILENAME),UA,UZ,UR,pres,PRES)

  CALL COMPUTE_VORTICITY(A,Z,R,OPA,OPZ,OPR,UA,UZ,UR,DG01,DG02,DG03,DG07,DG08,DG09)
  filename = trim(vort_snap)//num//'.h5'
  call EXPORT_snapshot(trim(FILENAME),DG01,DG02,DG03,pres,PRES)
  
  UAM1 = UA
  UZM1 = UZ
  URM1 = UR

  CALL DUMP_HDF5_BASIC(FILE_DUMP,'NEW',TC,DT,MSH,UA,UZ,UR,PRES,DG01,DG09)
  
  DG10 = UA + PRM_A*R + PRM_B/R
  
  CALL COMPUTE_NON_LINEAR_TERMS(&
       A, Z, R, OPA, OPZ, OPR, DG10, UZ, UR , NLAM1, NLZM1, NLRM1,&
       DG01, DG02, DG03, DG04, DG05, DG06, DG07, DG08, DG09)
  
  DG10 = PRM_A*R + PRM_B/R
  DG14 = 0._DP
  
  CALL COMPUTE_NON_LINEAR_TERMS(&
       A, Z, R, OPA, OPZ, OPR, DG10, DG14, DG14 , DG11, DG12, DG13,&
       DG01, DG02, DG03, DG04, DG05, DG06, DG07, DG08, DG09)

  NLAM1 = NLAM1 - DG11
  NLZM1 = NLZM1 - DG12
  NLRM1 = NLRM1 - DG13
  
    
  DO IT_time=1,nb_iter
     
     tc = tc + dt
     snap_dt = snap_dt + dt
     
     starttime = MPI_Wtime();
        
     CALL GRAD( A,Z,R, &
          OPA, OPZ, OPR, PRES, DG01,DG02, DG03)
     
     SA = (2._DP*UA-0.5_DP*UAM1)/DT - DG01 
     SZ = (2._DP*UZ-0.5_DP*UZM1)/DT - DG02 
     SR = (2._DP*UR-0.5_DP*URM1)/DT - DG03 

     DG10 = PRM_A*R + PRM_B/R  + UA
     
     CALL COMPUTE_NON_LINEAR_TERMS(&
          A, Z, R, OPA, OPZ, OPR, DG10, UZ, UR , NLA, NLZ, NLR,&
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

     UAM1=UA
     UZM1=UZ
     URM1=UR

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
     
       
     call GetCFL(msh(1),msh(2),msh(3), UA, UZ, UR, dt, cfl)


     !Calcul de J(u) = int_domaine <u,u> + <t,t>
     DG04 = UA*UA*R
     DG05 = UZ*UZ*R
     DG06 = UR*UR*R

     CALL get_quadrature_hhi(quad,DG04,DG01,ph%xst,ph%xen,n(1),ph%yst,ph%yen,n(2),ph%zst,ph%zen,n(3))
     CALL get_quadrature_hhi(quad,DG05,DG02,ph%xst,ph%xen,n(1),ph%yst,ph%yen,n(2),ph%zst,ph%zen,n(3))
     CALL get_quadrature_hhi(quad,DG06,DG03,ph%xst,ph%xen,n(1),ph%yst,ph%yen,n(2),ph%zst,ph%zen,n(3))
     
     J_U = (DG01(PH%XST(1),PH%XST(2),PH%XST(3)) + DG02(PH%XST(1),PH%XST(2),PH%XST(3)) + DG03(PH%XST(1),PH%XST(2),PH%XST(3)))


     if (rank==0) print'(i9,11(1x,e15.8))',it_time,tc,dt,cfl,DIV_MAX,endtime,J_U
     
     if (cfl .GT. 10.) then
        if (rank == 0) print'("CFL TOO BIG.")'
        call MPI_FINALIZE(ierr)
        stop
     end if

     if (mod(it_time,1000)==1) then

        CALL AZIMUT_MOD(UA,UZ,UR,E_AZI,max_mod)
        CALL VERTIC_MOD(UA,UZ,UR,E_VER,max_mod)
        
        if (nrank==0) then
           write(42,*)TC,J_U
           write(50,'(15(e15.8))')TC,DG01(PH%XST(1),PH%XST(2),PH%XST(3)),E_AZI
           write(60,'(15(e15.8))')TC,DG02(PH%XST(1),PH%XST(2),PH%XST(3)),E_VER
           write(70,'(15(e15.8))')TC,DG03(PH%XST(1),PH%XST(2),PH%XST(3))
        end if
        
        CALL AZIMUT_MOD_SCAL(UA,E_AZI,max_mod)
        CALL VERTIC_MOD_SCAL(UA,E_VER,max_mod)
     
        if (nrank == 0) then
           write(101,'(15(e15.8))')TC,E_AZI
           write(111,'(15(e15.8))')TC,E_VER
        end if
        
        CALL AZIMUT_MOD_SCAL(UZ,E_AZI,max_mod)
        CALL VERTIC_MOD_SCAL(UZ,E_VER,max_mod)
        
        if (nrank == 0) then
           write(102,'(15(e15.8))')TC,E_AZI
           write(112,'(15(e15.8))')TC,E_VER
        end if

        CALL AZIMUT_MOD_SCAL(UR,E_AZI,max_mod)
        CALL VERTIC_MOD_SCAL(UR,E_VER,max_mod)
        
        if (nrank == 0) then
           write(103,'(15(e15.8))')TC,E_AZI
           write(113,'(15(e15.8))')TC,E_VER
        end if

        filename = trim(base_snap)//'grid.h5'
        call update_id_and_time(trim(filename),snap_dt,snap_id)
        WRITE(num,'(I6.6)')snap_id
        filename = trim(base_snap)//num//'.h5'

        call EXPORT_snapshot(trim(FILENAME),UA,UZ,UR,pres,PRES)

        CALL COMPUTE_VORTICITY(A,Z,R,OPA,OPZ,OPR,UA,UZ,UR,DG01,DG02,DG03,DG07,DG08,DG09)
        filename = trim(vort_snap)//num//'.h5'
        call EXPORT_snapshot(trim(FILENAME),DG01,DG02,DG03,pres,PRES)

     end if

     if (mod(it_time,1000)==1) then
        CALL AZIMUT_MOD(UA,UZ,UR,SP_A,NA/2)
        CALL VERTIC_MOD(UA,UZ,UR,SP_Z,NZ/2)

        if (nrank==0) then
           write(form,'(i3)')na/2
           write(80,'(e15.8,'//TRIM(form)//'(e15.8))') tc,sp_a
               
           write(form,'(i3)')nz/2
           write(90,'(e15.8,'//TRIM(form)//'(e15.8))') tc,sp_z
        end if
     end if
     
     if (tc>=tmax) exit
     
  end DO

  CALL DUMP_HDF5_BASIC(file_dump,'WRITE',TC,DT,msh,UA,UZ,UR,Pres,DG01,dg09)
  
  CALL MPI_FINALIZE(ierr)

  STOP

  
contains

  SUBROUTINE azimut_mod(UA,UZ,UR,E_azi,max_mod)
    REAL(DP),ALLOCATABLE,DIMENSION(:,:,:),intent(in) :: UA,UZ,UR
    REAL(DP),DIMENSION(0:max_mod-1),intent(out) :: E_azi
    INTEGER,intent(in) :: max_mod

    type(C_PTR) :: plan
    INTEGER :: m

    DGA = UA
    DGZ = UZ
    DGR = UR
    
    call c2c_1m_x(DGA,plan_fwd_x)
    call c2c_1m_x(DGZ,plan_fwd_x)
    call c2c_1m_x(DGR,plan_fwd_x)

    E_azi = 0._DP
    do m = 0,max_mod-1
       do k = PH%XST(3), PH%XEN(3)
          do j = PH%XST(2), PH%XEN(2)
             E_azi(m) = E_azi(m) + ABS(DGA(m+1,j,k))**2 + ABS(DGZ(m+1,j,k))**2 + ABS(DGR(m+1,j,k))**2
          end do
       end do
    end do

    CALL MPI_ALLREDUCE(MPI_IN_PLACE,E_azi,max_mod,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)

  end SUBROUTINE azimut_mod

    SUBROUTINE azimut_mod_scal(U,E_azi,max_mod)
    REAL(DP),ALLOCATABLE,DIMENSION(:,:,:),intent(in) :: U
    REAL(DP),DIMENSION(0:max_mod-1),intent(out) :: E_azi
    INTEGER,intent(in) :: max_mod

    type(C_PTR) :: plan
    INTEGER :: m

    DGA = U
    call c2c_1m_x(DGA,plan_fwd_x)
    
    E_azi = 0._DP
    do m = 0,max_mod-1
       do k = PH%XST(3), PH%XEN(3)
          do j = PH%XST(2), PH%XEN(2)
             E_azi(m) = E_azi(m) + ABS(DGA(m+1,j,k))**2
          end do
       end do
    end do

    CALL MPI_ALLREDUCE(MPI_IN_PLACE,E_azi,max_mod,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)

  end SUBROUTINE azimut_mod_scal


  SUBROUTINE vertic_mod(UA,UZ,UR,E_ver,max_mod)
    REAL(DP),ALLOCATABLE,DIMENSION(:,:,:),intent(in) :: UA,UZ,UR
    REAL(DP),DIMENSION(0:max_mod-1),intent(out) :: E_ver
    INTEGER,intent(in) :: max_mod

    type(C_PTR) :: plan
    INTEGER :: m
    
    DGA = UA
    DGZ = UZ
    DGR = UR

    call transpose_x_to_y(DGA, DGA_Y)
    call transpose_x_to_y(DGZ, DGZ_Y)
    call transpose_x_to_y(DGR, DGR_Y)
    
    call c2c_1m_y(DGA_Y,plan_fwd_y)
    call c2c_1m_y(DGZ_Y,plan_fwd_y)
    call c2c_1m_y(DGR_Y,plan_fwd_y)

    E_ver = 0._DP
    do m = 0,max_mod-1
       do k = PH%YST(3), PH%YEN(3)
          do i = PH%YST(1), PH%YEN(1)
             E_ver(m) = E_ver(m) + ABS(DGA_Y(i,m+1,k))**2 + ABS(DGZ_Y(i,m+1,k))**2 + ABS(DGR_Y(i,m+1,k))**2
          end do
       end do
    end do

    CALL MPI_ALLREDUCE(MPI_IN_PLACE,E_ver,max_mod,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)

  end SUBROUTINE vertic_mod

  SUBROUTINE vertic_mod_scal(U,E_ver,max_mod)
    REAL(DP),ALLOCATABLE,DIMENSION(:,:,:),intent(in) :: U
    REAL(DP),DIMENSION(0:max_mod-1),intent(out) :: E_ver
    INTEGER,intent(in) :: max_mod

    type(C_PTR) :: plan
    INTEGER :: m
    
    DGA = U

    call transpose_x_to_y(DGA, DGA_Y)
    
    call c2c_1m_y(DGA_Y,plan_fwd_y)

    E_ver = 0._DP
    do m = 0,max_mod-1
       do k = PH%YST(3), PH%YEN(3)
          do i = PH%YST(1), PH%YEN(1)
             E_ver(m) = E_ver(m) + ABS(DGA_Y(i,m+1,k))**2
          end do
       end do
    end do

    CALL MPI_ALLREDUCE(MPI_IN_PLACE,E_ver,max_mod,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
    
  end SUBROUTINE vertic_mod_scal

  
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

    call alloc_x(DGA , OPT_GLOBAL=.TRUE.) ; DGA = 0._DP
    call alloc_x(DGZ , OPT_GLOBAL=.TRUE.) ; DGZ = 0._DP
    call alloc_x(DGR , OPT_GLOBAL=.TRUE.) ; DGR = 0._DP

    call alloc_y(DGA_Y , OPT_GLOBAL=.TRUE.) ; DGA_Y = 0._DP
    call alloc_y(DGZ_Y , OPT_GLOBAL=.TRUE.) ; DGZ_Y = 0._DP
    call alloc_y(DGR_Y , OPT_GLOBAL=.TRUE.) ; DGR_Y = 0._DP

    ALLOCATE(E_AZI(0:MAX_MOD-1))
    ALLOCATE(E_VER(0:MAX_MOD-1))

    ALLOCATE(SP_A(0:NA/2 -1))
    ALLOCATE(SP_Z(0:NZ/2 -1))
    
    FORALL(I=ph%XST(1):ph%XEN(1),J=ph%XST(2):ph%XEN(2),K=ph%XST(3):ph%XEN(3))
       A(I,J,K) = MSH(1)%X(I)
       Z(I,J,K) = MSH(2)%X(J)
       R(I,J,K) = MSH(3)%X(K)
    END FORALL

    CALL EXPORT_GRID(trim(base_snap)//'grid.h5',A,Z,R)
    
  end subroutine preproc

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
    REAL(DP),ALLOCATABLE,DIMENSION(:,:,:) :: UA,uz,ur
    integer :: i,k_max,k_cut_p,k_cut_m

    k_max = NA/3
    k_cut_p = k_max +2
    k_cut_m = NA - k_max + 1


    DGA = UA

    call c2c_1m_x(DGA, plan_fwd_x)

    DGA(k_cut_p : k_cut_m,:,:) = CMPLX(0._DP,0._DP)

    if (nrank == 0)print*,"range : ",k_cut_p," : ", k_cut_m
    
    call c2c_1m_x(DGA, plan_bck_x)

    DG01 = DGA
    
    
  end subroutine dealiazing


   subroutine dealiazing_spec(ua,uz,ur,k_max)
    COMPLEX(DP),ALLOCATABLE,DIMENSION(:,:,:) :: UA,uz,ur
    integer :: i,k_max,k_cut_p,k_cut_m

!    k_max = NA/3
    k_cut_p = k_max +2
    k_cut_m = NA - k_max + 1


    UA(k_cut_p : k_cut_m,:,:) = CMPLX(0._DP,0._DP)

    if (nrank == 0)print*,"K_max : ",k_max," range : ",k_cut_p," : ", k_cut_m
    
    
  end subroutine dealiazing_spec

    
end program tcheby_1d
