program tcheby_1d
  use m_mesh_base
  use m_numerics
  use m_operator_base
  use m_operator_tcheby
  use m_operator_fourier_dft
  use m_adjoint_tool_cyl
  use m_fourier_transform
  use m_tensor_product
  use m_udf_mms
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
  !>
  type(t_operator_fourier_dft) ::  opa,opz
  type(t_operator_tcheby) :: opr
  type(t_quadrature) :: quad
  !> 
  real(kind=8),allocatable     ::  dg01(:,:,:),dg02(:,:,:),dg03(:,:,:)
  real(kind=8),allocatable     ::  dg04(:,:,:),dg05(:,:,:),dg06(:,:,:)
    real(kind=8),allocatable     ::  dg07(:,:,:),dg08(:,:,:),dg09(:,:,:)

  integer                  ::  i,j,k

  real(kind=8) :: dirichl(2),neumann(2)
  integer :: n(3),cpu_grid(2),ierr
  REAL(kind=8) ::  xmin(3), xmax(3),err
  REAL(kind=8) ::  tc

  TYPE(DECOMP_INFO) :: ph
  real(kind=8),allocatable ,dimension(:,:,:) ::  A,Z,R,A_tot,Z_tot,R_tot  
  integer :: it_time,is(3),ie(3)
  integer , parameter :: OX=1,OY=2,OZ=3
   
  ! attention au type à lire
  REAL(kind=8) :: a_min,a_max,z_min,z_max,r_min,r_max
  integer :: na,nz,nr,nb_cpu_y,nb_cpu_z,nprocs
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
  character(len=1024) ::  file_dump, base_snap, filename, base_save
  character(len=6) :: num
  REAL(kind=8) ::  DIV_MAX,tmp,DJ_ADJ,J_DU0,J_U
  
  REAL(KIND=8),DIMENSION(3) :: NU_MOMENTUM
  REAL(KIND=8),DIMENSION(3) :: NU_ENERGY
  REAL(KIND=8),DIMENSION(3) :: NU_POISSON,err_mms
  REAL(kind=8) :: sigma,omega_tilde,noise,res,alpha,beta,vol,err_grad,err_lap,err_nl
  
  COMPLEX(DP),ALLOCATABLE,DIMENSION(:,:,:) :: DGA,DGZ,DGR
  COMPLEX(DP),ALLOCATABLE,DIMENSION(:,:,:) :: DGA_Y,DGZ_Y,DGR_Y
  
  LOGICAL :: memoire = .TRUE.
  LOGICAL :: do_adj = .FALSE.

  INTEGER, parameter :: inter_order = 1
  
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
  
  
  
  call preproc()

  CALL IMPORT_HDF5_INIT("prod_scal/x1.h5",DG01,DG02,DG03)
  CALL IMPORT_HDF5_INIT("prod_scal/x2.h5",DG04,DG05,DG06)

  !Calcul de J(u) = int_domaine <u,u> + <t,t>
  DG07 = DG01*DG04*R
  DG08 = DG02*DG05*R
  DG09 = DG03*DG06*R

  CALL integrate_spec(quad,DG07,J_U,PH,NA,NZ,NR,xmax,xmin)
  CALL integrate_spec(quad,DG08,RES,PH,NA,NZ,NR,xmax,xmin)
  J_U = RES + J_U
  CALL integrate_spec(quad,DG09,RES,PH,NA,NZ,NR,xmax,xmin)
  J_U = RES + J_U
  
  if (nrank == 0) then
     OPEN (UNIT=100, STATUS="new", FILE="prod_scal/Prod.dat", action='write')
     write(100,'(e15.8)')J_U
  end if
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
    
    call alloc_x(A  , OPT_GLOBAL=.TRUE.) ; A = 0
    call alloc_x(Z  , OPT_GLOBAL=.TRUE.) ; Z = 0
    call alloc_x(R  , OPT_GLOBAL=.TRUE.) ; R = 0

    
    FORALL(I=ph%XST(1):ph%XEN(1),J=ph%XST(2):ph%XEN(2),K=ph%XST(3):ph%XEN(3))
       A(I,J,K) = MSH(1)%X(I)
       Z(I,J,K) = MSH(2)%X(J)
       R(I,J,K) = MSH(3)%X(K)      
    END FORALL

    ALLOCATE(A_tot(1:NA+1,1:NZ+1,1:NR+1))
    ALLOCATE(Z_tot(1:NA+1,1:NZ+1,1:NR+1))
    ALLOCATE(R_tot(1:NA+1,1:NZ+1,1:NR+1))
    
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

  subroutine integrate_volume(U, X, Y, Z, integral)
    
    implicit none
    real(DP),ALLOCATABLE,DIMENSION(:,:,:),intent(in)  :: U,X,Y,Z
    real(DP), allocatable  :: dx(:), dy(:), dz(:)
    integer               :: i, j, k
    REAL(DP) :: integral
    
    allocate(dx(1:N(1)+1), dy(1:N(2)+1), dz(1:N(3)+1))
    
    ! --- Pas locaux : différences centrées + demi-cellule aux bords ---

    ! Direction X

    dx(1)  = 0.5_DP *  (X(2,1,1) - X(1,1,1))
    dx(N(1)+1) = 0.5_DP * (X(N(1)+1,1,1) - X(N(1),1,1))
    
    do i = 2, N(1)
       dx(i) = 0.5d0 * (X(i+1,1,1) - X(i-1,1,1))
    end do

    ! Direction Y

    dy(1)  = 0.5_DP * (Y(1,2,1) - Y(1,1,1))
    dy(N(2)+1) = 0.5_DP * (Y(1,N(2)+1,1)  - Y(1,N(2),1) )
    
    do j = 2, N(2)
       dy(j) = 0.5d0 * (Y(1,j+1,1) - Y(1,j-1,1))
    end do

    ! Direction Z

    dz(1)  = 0.5_DP * (Z(1,1,2) - Z(1,1,1))
    dz(N(3)+1) = 0.5_DP * (Z(1,1,N(3)+1) - Z(1,1,N(3)))
    
    do k = 2, N(3)
       dz(k) = 0.5d0 * (Z(1,1,k+1) - Z(1,1,k-1))
    end do

    ! --- Sommation pondérée (boucle i en interne : mémoire contiguë) ---
    integral = 0.0d0
    do k = PH%XST(3), PH%XEN(3)
       do j = PH%XST(2), PH%XEN(2)
          do i = PH%XST(1), PH%XEN(1)
             integral = integral + U(i,j,k) * dx(i) * dy(j) * dz(k)
          end do
       end do
    end do

    integral = integral * (N(2)+2)/(N(2)+1)
    
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,integral,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
    
    deallocate(dx, dy, dz)
    
  end subroutine integrate_volume
  
end program tcheby_1d
