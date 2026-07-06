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
   !>
   use m_solver_diag_cart_hhi
   use m_solver_diag_cart_hhi_cplx
   
   !>
   use m_navier_stokes_cart
   use m_snapshots
   use m_quadrature

   implicit none
   
   type(t_mesh_base)       ::  msh(3)
   !>
   type(t_solver_diag_cart_hhi) :: eqn_u, eqn_v, eqn_w, eqn_u_ad, eqn_v_ad, eqn_w_ad
   type(t_solver_diag_cart_hhi) :: eqn_p,eqn_p_ad
   type(t_solver_diag_cart_hhi) :: eqn_T,eqn_T_ad
   type(t_solver_diag_cart_hhi_cplx) :: eqn_uv,eqn_uv_ad
   
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
   REAL(kind=8) ::  xmin(3), xmax(3),err
   REAL(kind=8) ::  tc
   
   TYPE(DECOMP_INFO) :: ph
   
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
   integer :: l
   REAL(DP), allocatable :: carrot(:,:)
   
   real(kind=8), parameter, dimension(3) :: bet= (/ 4./15. ,  1./15. , 1./6. /) ! ck
   real(kind=8), parameter, dimension(3) :: gam= (/ 0.     ,-17./60. ,-5./12./) ! bk
   real(kind=8), parameter, dimension(3) :: sig= (/ 8./15. ,  5./12. , 3./4. /) ! ak
   
   
   real(kind=8), parameter, dimension(3) :: rho_l = (/ 0.     ,-17./60. ,-5./12./) ! bk
   real(kind=8), parameter, dimension(3) :: gam_l = (/ 8./15. ,  5./12. , 3./4. /) ! ak
   real(kind=8), parameter, dimension(3) :: alf_l = rho_l + gam_l ! (/ 4./15. ,  1./15. , 1./6. /) ! ck
   
   
   ! attention au type à lire
   REAL(kind=8) :: x_min,x_max,y_min,y_max,z_min,z_max,kx,ky
   integer :: nx,ny,nz,nb_cpu_y,nb_cpu_z
   namelist /parameters_cube/ nx,ny,nz,x_max,y_max,z_max,nb_cpu_y,nb_cpu_z
   
   REAL(KIND=8) :: NU,EPS,OMEGA,KAPPA
   INTEGER :: nb_iter
   namelist /parameters_physical/ NU,EPS,nb_iter,omega,kappa
   
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
   namelist /parameters_diagnostics/ dn_snap,snap_dir,resume_snap
   
   character(len=1024) ::  input_file,output_dir,init_file
   integer :: rank,iostat,seed_size
   integer, allocatable :: myseed(:)
   character(len=1024) ::  file_dump, base_snap, filename, base_save
   character(len=6) :: num
   REAL(kind=8) ::  DIV_MAX,tmp,DJ_ADJ,J_DU0,J_U
   
   REAL(KIND=8),DIMENSION(3) :: NU_MOMENTUM
   REAL(KIND=8),DIMENSION(3) :: NU_ENERGY
   REAL(KIND=8),DIMENSION(3) :: NU_POISSON
   REAL(kind=8) :: sigma,omega_tilde,noise,res
   
   LOGICAL :: do_adj = .FALSE.
   
   REAL(DP),DIMENSION(:,:,:,:),ALLOCATABLE :: SAVE_U,SAVE_V,SAVE_W
   

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
   end if
   
   
   
   x_min = 0
   y_min = 0
   z_min = 0
   
   
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
   
   
   CALL MPI_BCAST( NU, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR ) 
   CALL MPI_BCAST( EPS, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR ) 
   CALL MPI_BCAST( OMEGA, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR ) 
   CALL MPI_BCAST( KAPPA, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERR ) 
   
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
   
   N = [NX,NY,NZ]  
   CPU_GRID = [NB_CPU_Y,NB_CPU_Z]
   
   
   xmin(1:3) = [x_min,y_min,z_min]
   xmax(1:3) = [x_max,y_max,z_max]
   
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
   
   call test_dct(PH%XST,PH%XEN,NX+1,PH%YST,PH%YEN,NY+1,PH%ZST,PH%ZEN,NZ+1) ! attention 
   
   if (rank==0) then
      open(unit=20,file=TRIM(io_nrj))
   end if
   
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
   
   NU_MOMENTUM = -  nu*[1,1,1]!*0.5
   NU_POISSON  =  1.
   NU_ENERGY = -kappa
   
   CALL PREPROC()

   ! Lecture de la condition initiale
   CALL IMPORT_HDF5_INIT_BASIC("data_verif/grad.h5",U,V,W)
   CALL IMPORT_HDF5_INIT_BASIC("data_verif/Variation.h5",DU_0,DV_0,DW_0)

   !Saving the gradient computed by adjoint method
   
   FORALL(I=PH%XST(1):PH%XEN(1),J=PH%XST(2):PH%XEN(2),K=PH%XST(3):PH%XEN(3))
      DJ_U_ADJ(I,J,K) = U(I,J,K)*DU_0(I,J,K)
      DJ_V_ADJ(I,J,K) = V(I,J,K)*DV_0(I,J,K)
      DJ_W_ADJ(I,J,K) = W(I,J,K)*DW_0(I,J,K)
   END FORALL
   
   
   CALL get_quadrature_hhi(quad,DJ_U_ADJ,DG01,ph%xst,ph%xen,n(1),ph%yst,ph%yen,n(2),ph%zst,ph%zen,n(3))
   CALL get_quadrature_hhi(quad,DJ_V_ADJ,DG02,ph%xst,ph%xen,n(1),ph%yst,ph%yen,n(2),ph%zst,ph%zen,n(3))
   CALL get_quadrature_hhi(quad,DJ_W_ADJ,DG03,ph%xst,ph%xen,n(1),ph%yst,ph%yen,n(2),ph%zst,ph%zen,n(3))
   
   DJ_ADJ =  (DG01(PH%XST(1),PH%XST(2),PH%XST(3)) + DG02(PH%XST(1),PH%XST(2),PH%XST(3)) + DG03(PH%XST(1),PH%XST(2),PH%XST(3)))!/N(3)

   if (rank==0) then
      open(unit=42,file='data_verif/DJ_adj',action="write",status='new')
      write(42,*)DJ_ADJ
      close(42)
   end if
   
   call mpi_finalize(ierr)
   stop
   
   
   
   contains
   
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
      
      ALLOCATE(CARROT(N(3)+1,2))
      
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

end program tcheby_1d
