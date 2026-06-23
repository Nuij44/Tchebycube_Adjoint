program tcheby_1d
   use m_mesh_base
   use m_numerics
   use m_operator_base
   use m_operator_tcheby
   use m_operator_fourier_dft
   use m_adjoint_tool
   use m_fourier_transform
   use m_tensor_product
   
   use m_hdf5_ifce
   use mpi
   use decomp_2d
   use m_udf_mms
   !>
   use m_solver_diag_cart_hhi
   use m_solver_diag_cart_hhi_cplx
   
   !>
   use m_navier_stokes_cart
   use m_snapshots
   use m_quadrature
   use ifport
   implicit none
   
   type(t_mesh_base)       ::  msh(3)
   !>
   type(t_solver_diag_cart_hhi) :: eqn_u, eqn_v, eqn_w, eqn_u_ad, eqn_v_ad, eqn_w_ad
   type(t_solver_diag_cart_hhi) :: eqn_p,eqn_p_ad
   type(t_solver_diag_cart_hhi) :: eqn_T,eqn_T_ad
   
   !>
   type(t_operator_fourier_dft) ::  opx,opy
   type(t_operator_tcheby) :: opz
   type(t_quadrature) :: quad
   
   !> 
   real(kind=8),allocatable     ::  dg01(:,:,:),dg02(:,:,:),dg03(:,:,:)
   real(kind=8),allocatable     ::  dg04(:,:,:),dg05(:,:,:),dg06(:,:,:)
   real(kind=8),allocatable     ::  dg07(:,:,:),dg08(:,:,:),dg09(:,:,:)
   real(kind=8),allocatable     ::  dg10(:,:,:),dg11(:,:,:),dg12(:,:,:)
   
   integer                  ::  i,j,k
   
   real(kind=8) :: dirichl(2),neumann(2)  
   integer :: n(3),cpu_grid(2),ierr
   REAL(kind=8) ::  xmin(3), xmax(3),err,integ
   REAL(kind=8) ::  tc
   
   TYPE(DECOMP_INFO) :: ph
   
   real(kind=8),allocatable ,dimension(:,:,:) ::  U_Tot,T_tot
   real(kind=8),allocatable ,dimension(:,:,:) ::  U_T,V_T,W_T
   real(kind=8),allocatable ,dimension(:,:,:) ::  U_0,V_0,W_0,T_0
   real(kind=8),allocatable ,dimension(:,:,:) ::  UM1,VM1,WM1,TM1
   real(kind=8),allocatable ,dimension(:,:,:) ::  NLUM1,NLVM1,NLWM1,NLTM1
   real(kind=8),allocatable ,dimension(:,:,:) ::  DU_0,DV_0,DW_0
   real(kind=8),allocatable ,dimension(:,:,:) ::  DJ_U_ADJ,DJ_V_ADJ,DJ_W_ADJ
   real(kind=8),allocatable ,dimension(:,:,:) ::  U,V,W,PRES,T,FI
   real(kind=8),allocatable ,dimension(:,:,:) ::  SU,SV,SW,SFI,ST
   real(kind=8),allocatable ,dimension(:,:,:) ::  NLU,NLV,NLW,NLT
   real(kind=8),allocatable ,dimension(:,:,:) ::  LU,LV,LW,LT,X,Y,Z
   real(dp),allocatable     ::  NOISE_U(:,:,:),NOISE_V(:,:,:),NOISE_W(:,:,:),NOISE_T(:,:,:)
   real(kind=8) :: w2,t1,h
   integer :: it_time,is(3),ie(3)
   integer , parameter :: OX=1,OY=2,OZ=3
   
   
   ! attention au type à lire
   REAL(kind=8) :: x_min,x_max,y_min,y_max,z_min,z_max
   integer :: nx,ny,nz,nb_cpu_y,nb_cpu_z
   namelist /parameters_cube/ nx,ny,nz,x_max,y_max,z_max,nb_cpu_y,nb_cpu_z
   
   REAL(KIND=8) :: Re,Pr,Ri
   INTEGER :: nb_iter
   namelist /parameters_physical/ Re,Pr,Ri
   
   character(len=1024) ::  root_dir,io_nrj
   REAL(kind=8) ::  dt,cfl,tmax
   integer dn_dump
   logical :: resume , explicit,correction_pres
   namelist /parameters_timescheme/ tmax,dt, cfl, dn_dump,root_dir,resume,explicit,io_nrj,correction_pres
   
   logical :: resume_snap=.false.
   integer :: dn_snap
   character(len=1024) ::  snap_dir
   real(kind=8) :: snap_dt,starttime,endtime
   integer :: snap_id
   namelist /parameters_diagnostics/ dn_snap,snap_dir
   
   character(len=1024) ::  input_file,output_dir,init_file
   integer :: rank,iostat,seed_size
   integer, allocatable :: myseed(:)
   character(len=1024) ::  file_dump, base_snap, filename, base_save
   character(len=6) :: num
   REAL(kind=8) ::  DIV_MAX,tmp,DJ_ADJ,J_DU0,J_U
   
   REAL(KIND=8),DIMENSION(3) :: NU_MOMENTUM
   REAL(KIND=8),DIMENSION(3) :: NU_ENERGY
   REAL(KIND=8),DIMENSION(3) :: NU_POISSON
   REAL(kind=8) :: sigma,alpha_buoy,noise,res
   
   LOGICAL :: do_adj = .FALSE.
   
   REAL(DP),DIMENSION(:,:,:,:),ALLOCATABLE :: SAVE_U,SAVE_V,SAVE_W,SAVE_T

   REAL(DP),DIMENSION(:),allocatable :: plot_X,plot_Y
   integer :: UNIT

   COMPLEX(DP),ALLOCATABLE,DIMENSION(:,:,:) :: DGA,DGZ,DGR
   COMPLEX(DP),ALLOCATABLE,DIMENSION(:,:,:) :: DGA_Y,DGZ_Y,DGR_Y

   
   call mpi_init(ierr)
   CALL H5OPEN_F(IERR)
   
   call mpi_comm_rank(mpi_comm_world,rank,ierr)
   
   call command_line_read_input(input_file)
  


!   root_dir = 'data_fun/'
   init_file = '../Result/init_opt.h5'
!   snap_dir = 'snapshots/'
   if (rank==0) then
   print*,'input_file',input_file
      write(6,*)TRIM('mkdir -p '//trim(root_dir)//trim(snap_dir))
      OPEN (UNIT=24, FILE=TRIM(input_file),status='old', action='read', iostat=ierr)
      print*,'ierr',ierr
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
      write(*,parameters_timescheme)
      write(*,parameters_diagnostics)
      OPEN(UNIT=42, FILE='nrj_SSP.dat', action='write')
    end if
   
   
   
   x_min = 0._DP
   y_min = 0._DP
   z_min = -1._DP
   
   
   CALL MPI_BCAST( nx      , 1,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST( ny      , 1,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST( nz      , 1,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST( nb_iter , 1,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
   
   CALL MPI_BCAST( x_min   , 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST( x_max   , 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   
   CALL MPI_BCAST( y_min   , 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST( y_max   , 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST( z_min   , 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST( z_max   , 1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
   
   CALL MPI_BCAST( nb_cpu_y, 1,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST( nb_cpu_z, 1,MPI_INTEGER         ,0,MPI_COMM_WORLD,IERR)
   
   CALL MPI_BCAST( Re, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR ) 
   CALL MPI_BCAST( Pr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR ) 
   CALL MPI_BCAST( Ri, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR ) 
   
   CALL MPI_BCAST( resume  , 1   , MPI_LOGICAL         , 0, MPI_COMM_WORLD, IERR ) 
   CALL MPI_BCAST( root_dir, 1024, MPI_character       , 0, MPI_COMM_WORLD, IERR )
   CALL MPI_BCAST( dn_dump , 1   , MPI_INTEGER         , 0, MPI_COMM_WORLD, IERR )
   CALL MPI_BCAST( dt      , 1   , MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR ) 
   CALL MPI_BCAST( cfl     , 1   , MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR )
   CALL MPI_BCAST( explicit, 1   , MPI_LOGICAL         , 0, MPI_COMM_WORLD, IERR )
   CALL MPI_BCAST( correction_pres, 1   , MPI_LOGICAL         , 0, MPI_COMM_WORLD, IERR )
   CALL MPI_BCAST( tmax     , 1   , MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR )
   
   CALL MPI_BCAST( dn_snap , 1   , MPI_INTEGER         , 0, MPI_COMM_WORLD, IERR )
   CALL MPI_BCAST( snap_dir, 1024, MPI_character       , 0, MPI_COMM_WORLD, IERR )
   

   
   xmin(1:3) = [x_min,y_min,z_min]
   xmax(1:3) = [4.35*pi,1.05*pi,1.]
 


   N = [NX,NY,NZ]

   CPU_GRID = [nb_cpu_y,nb_cpu_z]
  
   nb_iter = floor(tmax/dt)

   file_dump = trim(root_dir)//'dump/dump.h5'
   base_snap = trim(root_dir)//trim(snap_dir)//'snap_'
   base_save = trim(root_dir)//'save/save_'

   if (rank==0) then
      print*,'file_dump:',trim(file_dump)
      print*,'base_snap:',trim(base_snap)
      print*,'base_save:',trim(base_save)
   end if
      


   ! GRILLE 2D
   CALL DECOMP_2D_INIT( &
   NX = N(1)+1, NY = N(2)+1, NZ = N(3)+1 , P_ROW= CPU_GRID(1) , P_COL= CPU_GRID(2) )
   
   CALL GET_DECOMP_INFO(PH)
   
   CALL START_TENSOR_PRODUCT(&
   PH%XST,PH%XEN,N(1)+1,PH%YST,PH%YEN,N(2)+1,PH%ZST,PH%ZEN,N(3)+1)
   CALL START_FOURIER_TRANSFORMS( &
   PH%XST, PH%XEN, N(1)+1, PH%YST, PH%YEN,N(2)+1, PH%ZST, PH%ZEN, N(3)+1 )
   
   call test_dct(PH%XST,PH%XEN,NX+1,PH%YST,PH%YEN,NY+1,PH%ZST,PH%ZEN,NZ+1) ! attention 
   
   CALL MSH(1)%INITIALIZE(XMIN(1),XMAX(1),N(1),.TRUE. ) 
   CALL MSH(2)%INITIALIZE(XMIN(2),XMAX(2),N(2),.TRUE. ) 
   CALL MSH(3)%INITIALIZE(XMIN(3),XMAX(3),N(3) ) 
   
   
   CALL START_OP_BASE()
   CALL OPX%INIT_OPERATOR_FOURIER_DFT( MESH=MSH(OX), AXIS=OX )
   CALL OPY%INIT_OPERATOR_FOURIER_DFT( MESH=MSH(OY), AXIS=OY )
   CALL OPZ%INIT_OPERATOR_TCHEBY( MESH=MSH(OZ), AXIS=OZ )
   
   
   CALL INIT_quadrature_hhi(quad,xmin,xmax,&
   ph%xst,ph%xen,nx,ph%yst,ph%yen,ny,ph%zst,ph%zen,nz)
   
   
   dirichl =  [1.,0.]
   neumann =  [0.,1.]
   
   NU_MOMENTUM = -  1._DP/Re*[1,1,1]
   NU_POISSON  =  1.
   NU_ENERGY = - 1._DP/(Re*Pr)
   ALPHA_BUOY = Ri
   
   CALL EQN_U%SET_PARAMS( NU=NU_MOMENTUM , SIGMA=0D0 )
   CALL EQN_U%SET_BCS( AXIS=3 , BCS_MINUS=DIRICHL , BCS_PLUS=DIRICHL )
   CALL EQN_U%INITIALISE( MSH, OPX, OPY, OPZ, PH )
   CALL EQN_U%SET_BVS( MESH = MSH , AXIS=3 , UDF_MINUS=UDF_NULL , UDF_PLUS=UDF_NULL )
   
   CALL EQN_V%SET_PARAMS( NU=NU_MOMENTUM , SIGMA=0D0 )
   CALL EQN_V%SET_BCS( AXIS=3 , BCS_MINUS=DIRICHL , BCS_PLUS=DIRICHL )
   CALL EQN_V%INITIALISE(MSH,OPX,OPY,OPZ,PH)
   CALL EQN_V%SET_BVS( MESH = MSH , AXIS=3 , UDF_MINUS=UDF_NULL , UDF_PLUS=UDF_NULL )
   
   CALL EQN_W%SET_PARAMS( NU=NU_MOMENTUM , SIGMA=0D0 )
   CALL EQN_W%SET_BCS( AXIS=3 , BCS_MINUS=DIRICHL , BCS_PLUS=DIRICHL )
   CALL EQN_W%INITIALISE(MSH,OPX,OPY,OPZ,PH)
   CALL EQN_W%SET_BVS( MESH = MSH , AXIS=3 , UDF_MINUS=UDF_NULL , UDF_PLUS=UDF_NULL )
   
   CALL EQN_T%SET_PARAMS( NU=NU_ENERGY , SIGMA=0D0 )
   CALL EQN_T%SET_BCS( AXIS=3 , BCS_MINUS=DIRICHL , BCS_PLUS=DIRICHL )
   CALL EQN_T%INITIALISE( MSH , OPX , OPY , OPZ , PH )
   CALL EQN_T%SET_BVS( MESH = MSH , AXIS=3 , UDF_MINUS=UDF_NULL , UDF_PLUS=UDF_NULL )
   
   CALL EQN_P%SET_PARAMS( NU=NU_POISSON , SIGMA=0D0 )
   CALL EQN_P%SET_BCS( AXIS=3 , BCS_MINUS=NEUMANN , BCS_PLUS=NEUMANN )
   CALL EQN_P%INITIALISE(MSH,OPX,OPY,OPZ,PH)
   CALL EQN_P%SET_BVS( MESH = MSH , AXIS=3 , UDF_MINUS=UDF_NULL , UDF_PLUS=UDF_NULL )
   
   
   call preproc()

   !Allocattion des tableaux pour sauver les iterations   
   if (do_adj) then
      ALLOCATE(SAVE_U(0:nb_iter,PH%XST(1):PH%XEN(1),PH%XST(2):PH%XEN(2),PH%XST(3):PH%XEN(3))); SAVE_U = 0._DP
      ALLOCATE(SAVE_V(0:nb_iter,PH%XST(1):PH%XEN(1),PH%XST(2):PH%XEN(2),PH%XST(3):PH%XEN(3))); SAVE_V = 0._DP
      ALLOCATE(SAVE_W(0:nb_iter,PH%XST(1):PH%XEN(1),PH%XST(2):PH%XEN(2),PH%XST(3):PH%XEN(3))); SAVE_W = 0._DP
      ALLOCATE(SAVE_T(0:nb_iter,PH%XST(1):PH%XEN(1),PH%XST(2):PH%XEN(2),PH%XST(3):PH%XEN(3))); SAVE_T = 0._DP
   end if

!   init_file = '../Result/init_opt.h5'
!   init_file = 'output_data_interpolated.h5'

!   if (rank == 0) print*,'TRIM(init_file):',TRIM(init_file)
!   CALL IMPORT_HDF5_INIT_STRAT(TRIM(init_file),U,V,W,T)

!   U = UM1*0.1
!   V = VM1*0.1
!   W = WM1*0.1
!   T = TM1*0.1

   !Condition initiale laminaire
!   U = 0._DP
!   V = 0._DP
!   W = 0._DP
   
   !condition initiale random
   if (rank == 0) print*,'Initialisation Random'
   CALL RANDOM_SEED()
   CALL Random_Number(U)
   CALL Random_Number(V)
   CALL Random_Number(W)
   U = 2._DP*U - 1._DP
   V = 2._DP*V - 1._DP
   W = 2._DP*W - 1._DP

   CALL normalization(U,V,W,1e-2)

!   NOISE = 1E-2

!   U = NOISE*U_0
!   V = NOISE*V_0
!   W = NOISE*W_0
!   T = NOISE*T_0

   if (rank ==0) then
      print*,'U:',maxval(U),minval(U)
      print*,'V:',maxval(V),minval(V)
      print*,'W:',maxval(W),minval(W)
   end if


  CALL DUMP_HDF5_BASIC(FILE_DUMP,'NEW',TC,DT,MSH,U,V,W,PRES,T,DG09)

   
   UM1 = U
   VM1 = V
   WM1 = W
   
   tc = 0
   snap_dt = 0
   J_U = 0._DP


   CALL EXPORT_GRID(trim(base_snap)//'grid.h5',X,Y,Z)


   filename = trim(base_snap)//'grid.h5'
   call update_id_and_time(trim(filename),snap_dt,snap_id)
   WRITE(num,'(I6.6)')snap_id

   NLUM1=NLU
   NLVM1=NLV
   NLWM1=NLW

   CALL dealiazing(u,v,w)

   FORALL(I=PH%XST(1):PH%XEN(1),J=PH%XST(2):PH%XEN(2),K=PH%XST(3):PH%XEN(3))
      U_tot(I,J,K) = U(I,J,K) + Z(I,J,K)
   END FORALL


   
   CALL COMPUTE_NON_LINEAR_TERMS( OPX,OPY,OPZ, U_tot, V, W, NLU, NLV, NLW, dg01, dg02, dg03, dg04, dg05, dg06, dg07, dg08, dg09)
   


   DO IT_time=1,nb_iter
      
      tc = tc + dt
      snap_dt = snap_dt + dt
      
      starttime = MPI_Wtime();
      
      FORALL(I=PH%XST(1):PH%XEN(1),J=PH%XST(2):PH%XEN(2),K=PH%XST(3):PH%XEN(3))
         U_tot(I,J,K) = U(I,J,K) + Z(I,J,K)
      END FORALL  
      
      CALL GRAD(&
      OPX, OPY, OPZ, PRES, DG01,DG02, DG03)
      
      SU = (2._DP*U-0.5_DP*UM1)/DT - DG01
      SV = (2._DP*V-0.5_DP*VM1)/DT - DG02
      SW = (2._DP*W-0.5_DP*WM1)/DT - DG03

      UM1=U
      VM1=V
      WM1=W
            
      CALL dealiazing(u,v,w)
      
      CALL COMPUTE_NON_LINEAR_TERMS( OPX,OPY,OPZ, U_tot, V, W, NLU, NLV, NLW, dg01, dg02, dg03, dg04, dg05, dg06, dg07, dg08, dg09)
      
      SU = SU - 2._DP*NLU + NLUM1
      SV = SV - 2._DP*NLV + NLVM1
      SW = SW - 2._DP*NLW + NLWM1
      
      
      NLUM1 = NLU
      NLVM1 = NLV
      NLWM1 = NLW
      
      SIGMA = 1.5_DP/DT
      
      CALL EQN_U%SOLVE(U, SU, SIGMA ,PH)
      CALL EQN_V%SOLVE(V, SV, SIGMA ,PH)
      CALL EQN_W%SOLVE(W, SW, SIGMA ,PH)
            
      call DIV( OPX, OPY, OPZ, U, V, W, SFI, dg01, dg02 , dg03 )
      SFI = SFI*1.5_DP/DT
      
      call EQN_P%SOLVE_POISSON(FI,SFI,NU_POISSON,PH)
      
      CALL GRAD(OPX, OPY, OPZ, FI, DG01, DG02, DG03)
      
      PRES = PRES + FI
      
      
      is = get_is_b([0,0,0])
      ie = get_ie_b([0,0,0])
      FORALL(I=IS(1):IE(1),J=IS(2):IE(2),K=IS(3):IE(3))
      U(I,J,K) = U(I,J,K) - DG01(I,J,K) * 2._DP*DT/3._DP
      END FORALL
      is = get_is_b([0,0,0])
      ie = get_ie_b([0,0,0])
      FORALL(I=IS(1):IE(1),J=IS(2):IE(2),K=IS(3):IE(3))
      V(I,J,K) = V(I,J,K) - DG02(I,J,K) * 2._DP*DT/3._DP
      END FORALL
      is = get_is_b([0,0,0])
      ie = get_ie_b([0,0,0])
      FORALL(I=IS(1):IE(1),J=IS(2):IE(2),K=IS(3):IE(3))
      W(I,J,K) = W(I,J,K) - DG03(I,J,K) * 2._DP*DT/3._DP 
      END FORALL
      
      
      ! check divergence 
      call DIV( OPX, OPY, OPZ, U, V, W, DG09, dg01, dg02 , dg03 )
      is = get_is()
      ie = get_ie()
      DIV_MAX = MAXVAL(ABS(DG09(IS(1):IE(1),IS(2):IE(2),IS(3):IE(3))))     
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,DIV_MAX,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,IERR)
      
      endtime   = MPI_Wtime();
      endtime =  endtime-starttime
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,endtime,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,IERR)
         
      call GetCFL(msh(1),msh(2),msh(3), U, V, W, dt, cfl)

      
      !Calcul de J(u) = int_domaine <u,u> + <t,t>
      DG04 = U**2
      DG05 = V**2
      DG06 = W**2

      CALL integrate_spec(quad,DG04,J_U,PH,N(1),N(2),N(3),xmax,xmin)
      CALL integrate_spec(quad,DG05,integ,PH,N(1),N(2),N(3),xmax,xmin)
      J_U = integ + J_U
      CALL integrate_spec(quad,DG06,integ,PH,N(1),N(2),N(3),xmax,xmin)
      J_U = integ + J_U


     if (rank==0) then
        write(42,'(e15.8,1x,e15.8)') tc, J_U*0.5_DP
     end if
      
     if (rank==0) print'(i9,10(1x,e15.8))',it_time,tc,dt,cfl,J_U,DIV_MAX,endtime,J_U*0.5_DP
      
      if (cfl .GT. 10.) then
         if (rank == 0) print'("CFL TOO BIG.")'
         call MPI_FINALIZE(ierr)
         stop
      end if
 
     if (mod(it_time,10000)==0) then
        CALL DUMP_HDF5_BASIC('dump/dump.h5','WRITE',TC,DT,msh,u,v,w,Pres,T,dg09)
     end if

     if (mod(it_time,100000)==0) then
        filename = trim(base_snap)//'grid.h5'
        call update_id_and_time(trim(filename),snap_dt,snap_id)
        WRITE(num,'(I6.6)')snap_id
        filename = trim(base_snap)//num//'.h5'

        FORALL(I=PH%XST(1):PH%XEN(1),J=PH%XST(2):PH%XEN(2),K=PH%XST(3):PH%XEN(3))
           U_tot(I,J,K) = U(I,J,K) + Z(I,J,K)
        END FORALL
        call EXPORT_snapshot(trim(FILENAME),u_tot,v,w,pres,pres)
     end if
     


     if (tc>=tmax) exit


   end DO
   
!   call stop_gnuplot()
   CALL MPI_FINALIZE(ierr)
   STOP
   


   contains

     subroutine start_gnuplot(nb_iter)
       integer :: nb_iter, rc
       character(7) :: pipe_name = 'GNUPLOT'

       ALLOCATE(plot_X(nb_iter),plot_Y(nb_iter))
       plot_X = 0._DP
       plot_Y = 0._DP
       
       ! Create the named pipe
       rc = SYSTEM("mkfifo " // trim(pipe_name))
       if (rc /= 0) then
          write(*, *) "Failed to create the named pipe."
          stop
       end if
       write(*, *) "Named pipe created: ", pipe_name

       !  Open the named pipe for writing
       open(newunit=unit, file=trim(pipe_name), status="old", access="stream", form="formatted", iostat=rc)
       if (rc /= 0) then
          write(*, *) "Failed to open the named pipe for writing."
          stop
       end if
       write(*, *) "Pipe opened writing writing."
       
       ! Start gnuplot in the background and read from the named pipe
       rc = system('gnuplot -p < ' // trim(pipe_name) // ' &')
       
       ! Define plot attributes
       write(unit, "(a)") "set term x11; &
            set title 'Fortran Realtime Plot'; &
            set grid; &
            set xlabel 'X Axis'; &
            set ylabel 'Y Axis'; &
            set logscale y;"

     end subroutine start_gnuplot

   
     subroutine plot_gnuplot(iter,T,E)
       integer :: iter,i
       REAL(DP) :: T,E

       ! Keep previous values
       plot_x(iter) = T
       plot_y(iter) = E
       
       ! Set the x and y ranges dynamically
!       write(unit, "(a,F10.3,a)") "set xrange [0:", maxval(plot_x(1:iter)), "]"
!       write(unit, "(a,F10.3,a)") "set yrange [0:", maxval(plot_y(1:iter)), "]"

       write(unit, '(a)') "plot '-' with linespoints"
       do i = 1, iter
          write(unit, '(F10.3,a,F10.3)') plot_x(i), " ", plot_y(i)
          call flush(unit)
       enddo
       write(unit, '(a)') 'e'
     end subroutine plot_gnuplot

     subroutine stop_gnuplot()
       character(7) :: pipe_name='GNUPLOT'
       integer :: rc
       ! Close the named pipe
       close(unit)
       
       ! Remove the named pipe  
       call execute_command_line("rm "// pipe_name, exitstat=rc)
       if (rc /= 0) then
          print *, "Failed to remove named pipe"
          stop
       end if
     end subroutine stop_gnuplot

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
      
      call alloc_x(U_T   , OPT_GLOBAL=.TRUE.) ; U_T = 0 
      call alloc_x(V_T   , OPT_GLOBAL=.TRUE.) ; V_T = 0
      call alloc_x(W_T   , OPT_GLOBAL=.TRUE.) ; W_T = 0

      call alloc_x(U_tot , OPT_GLOBAL=.TRUE.) ; U_tot = 0 
      call alloc_x(T_tot , OPT_GLOBAL=.TRUE.) ; T_tot = 0
      
      call alloc_x(PRES , OPT_GLOBAL=.TRUE.) ; PRES = 0
      call alloc_x(FI   , OPT_GLOBAL=.TRUE.) ; FI = 0
      call alloc_x(SFI  , OPT_GLOBAL=.TRUE.) ; SFI = 0
      !
      call alloc_x(U   , OPT_GLOBAL=.TRUE.) ; U = 0
      call alloc_x(SU  , OPT_GLOBAL=.TRUE.) ; SU = 0
      call alloc_x(NLU , OPT_GLOBAL=.TRUE.) ; NLU = 0
      call alloc_x(LU  , OPT_GLOBAL=.TRUE.) ; LU = 0
      !
      call alloc_x(V   , OPT_GLOBAL=.TRUE.) ; V = 0
      call alloc_x(SV  , OPT_GLOBAL=.TRUE.) ; SV = 0
      call alloc_x(NLV , OPT_GLOBAL=.TRUE.) ; NLV = 0
      call alloc_x(LV  , OPT_GLOBAL=.TRUE.) ; LV = 0
      !
      call alloc_x(W   , OPT_GLOBAL=.TRUE.) ; W = 0
      call alloc_x(SW  , OPT_GLOBAL=.TRUE.) ; SW = 0
      call alloc_x(NLW , OPT_GLOBAL=.TRUE.) ; NLW = 0
      call alloc_x(LW  , OPT_GLOBAL=.TRUE.) ; LW = 0
      
      call alloc_x(T   , OPT_GLOBAL=.TRUE.) ; T = 0
      call alloc_x(NLT   , OPT_GLOBAL=.TRUE.) ; NLT = 0
      call alloc_x(ST   , OPT_GLOBAL=.TRUE.) ; ST = 0
      
      call alloc_x(X  , OPT_GLOBAL=.TRUE.) ; X = 0
      call alloc_x(Y  , OPT_GLOBAL=.TRUE.) ; Y = 0
      call alloc_x(Z  , OPT_GLOBAL=.TRUE.) ; Z = 0
      
      call alloc_x(DJ_U_ADJ  , OPT_GLOBAL=.TRUE.) ; DJ_U_ADJ = 0
      call alloc_x(DJ_V_ADJ  , OPT_GLOBAL=.TRUE.) ; DJ_V_ADJ = 0
      call alloc_x(DJ_W_ADJ  , OPT_GLOBAL=.TRUE.) ; DJ_W_ADJ = 0
      
      call alloc_x(NOISE_U  , OPT_GLOBAL=.TRUE.) 
      call alloc_x(NOISE_V  , OPT_GLOBAL=.TRUE.) 
      call alloc_x(NOISE_W  , OPT_GLOBAL=.TRUE.) 
      call alloc_x(NOISE_T  , OPT_GLOBAL=.TRUE.) 
      
      call alloc_x(UM1  , OPT_GLOBAL=.TRUE.) 
      call alloc_x(VM1  , OPT_GLOBAL=.TRUE.) 
      call alloc_x(WM1  , OPT_GLOBAL=.TRUE.) 
      call alloc_x(TM1  , OPT_GLOBAL=.TRUE.) 
      
      call alloc_x(NLUM1  , OPT_GLOBAL=.TRUE.) 
      call alloc_x(NLVM1  , OPT_GLOBAL=.TRUE.) 
      call alloc_x(NLWM1  , OPT_GLOBAL=.TRUE.)
      call alloc_x(NLTM1  , OPT_GLOBAL=.TRUE.) 
      
      call alloc_x(U_0  , OPT_GLOBAL=.TRUE.) 
      call alloc_x(V_0  , OPT_GLOBAL=.TRUE.) 
      call alloc_x(W_0  , OPT_GLOBAL=.TRUE.) 
      call alloc_x(T_0  , OPT_GLOBAL=.TRUE.) 
      
      call alloc_x(DU_0  , OPT_GLOBAL=.TRUE.) 
      call alloc_x(DV_0  , OPT_GLOBAL=.TRUE.) 
      call alloc_x(DW_0  , OPT_GLOBAL=.TRUE.) 
      
      call alloc_x(DGA , OPT_GLOBAL=.TRUE.) ; DGA = 0._DP
      call alloc_x(DGZ , OPT_GLOBAL=.TRUE.) ; DGZ = 0._DP
      call alloc_x(DGR , OPT_GLOBAL=.TRUE.) ; DGR = 0._DP

      call alloc_y(DGA_Y , OPT_GLOBAL=.TRUE.) ; DGA_Y = 0._DP
      call alloc_y(DGZ_Y , OPT_GLOBAL=.TRUE.) ; DGZ_Y = 0._DP
      call alloc_y(DGR_Y , OPT_GLOBAL=.TRUE.) ; DGR_Y = 0._DP


      
      FORALL(I=ph%XST(1):ph%XEN(1),J=ph%XST(2):ph%XEN(2),K=ph%XST(3):ph%XEN(3))
      X(I,J,K) = MSH(1)%X(I)
      Y(I,J,K) = MSH(2)%X(J)
      Z(I,J,K) = MSH(3)%X(K)
      END FORALL
      
      
      
      
   end subroutine preproc
   
   
   real(kind=8) pure function udf_one(t,x,y,z)
   real(kind=8),intent(in):: t,x,y,z
   udf_one = 1.0
end function udf_one



real(kind=8) pure function udf_sol(t,x,y,z)
real(kind=8),intent(in):: t,x,y,z
udf_sol = sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)    
end function udf_sol

real(kind=8) pure function udf_sol_u(t,x,y,z)
real(kind=8),intent(in):: t,x,y,z
udf_sol_u = sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z)    
end function udf_sol_u

real(kind=8) pure function udf_sol_v(t,x,y,z)
real(kind=8),intent(in):: t,x,y,z
udf_sol_v = sin(4*pi*x)*sin(4*pi*y)*sin(4*pi*z)    
end function udf_sol_v

real(kind=8) pure function udf_scm_u(t,x,y,z)
real(kind=8),intent(in):: t,x,y,z
udf_scm_u = x
end function udf_scm_u
real(kind=8) pure function udf_scm_v(t,x,y,z)
real(kind=8),intent(in):: t,x,y,z
udf_scm_v = y
end function udf_scm_v

real(kind=8) pure function udf_perturb(t,x,y,z)
real(kind=8),intent(in):: t,x,y,z
udf_perturb = 1 + sin(4*2*pi*x/x_max)*sin(4*2*pi*y/y_max)
end function udf_perturb


real(kind=8) pure function udf_scm(t,x,y,z)
real(kind=8),intent(in):: t,x,y,z
udf_scm = x
end function udf_scm

real(kind=8) pure function udf_poiseuille(t,x,y,z)
real(kind=8),intent(in):: t,x,y,z
udf_poiseuille = 1-z**2
end function udf_poiseuille

real(kind=8) pure function udf_linear(t,x,y,z)
real(kind=8),intent(in):: t,x,y,z
udf_linear = 1-x
end function udf_linear


real(kind=8) pure function udf_null(t,x,y,z)
real(kind=8),intent(in):: t,x,y,z
udf_null = 0.
end function udf_null

real(kind=8) pure function udf_unit(t,x,y,z)
real(kind=8),intent(in):: t,x,y,z
udf_unit = 1.
end function udf_unit


real(kind=8) pure function udf_minus(t,x,y,z)
real(kind=8),intent(in):: t,x,y,z
udf_minus = 0.
end function udf_minus

real(kind=8) pure function udf_plus(t,x,y,z)
real(kind=8),intent(in):: t,x,y,z
udf_plus = 1.
end function udf_plus


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

    k_max = N(1)/3
    k_cut_p = k_max +2
    k_cut_m = N(1) - k_max + 1


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

    k_max = N(2)/3
    k_cut_p = k_max +2
    k_cut_m = N(2) - k_max + 1

    
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

    DG01 = UA * UA 
    CALL integrate_spec(quad,DG01,NRJ,PH,N(1),N(2),N(3),xmax,xmin)
    DG01 = UZ * UZ 
    CALL integrate_spec(quad,DG01,VOL,PH,N(1),N(2),N(3),xmax,xmin)
    NRJ = NRJ + VOL
    DG01 = UR * UR 
    CALL integrate_spec(quad,DG01,VOL,PH,N(1),N(2),N(3),xmax,xmin)
    NRJ = NRJ + VOL  
    UA = E0 * UA / NRJ
    UZ = E0 * UZ / NRJ
    UR = E0 * UR / NRJ
  end subroutine NORMALIZATION

end program tcheby_1d
