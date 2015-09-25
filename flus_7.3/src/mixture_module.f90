module mixture_module

use mpi
use config_module
use grid_module
use numerics_module
use data_module
implicit none; private
public t_Mixt
public DelMixt
public SetMixtFromInitial
public SetMixtFromFile
public SetMixtFromBackup
public Postprocessing
public ComputeTimeStep
public ComputeError
public ComputeResidual
public ComputeDiagonalTerm
public ComputeDeltaFlux
public SumConservative

! Global variable
integer :: nDim ! dimension
integer :: nDimP1 ! dimension+1
integer :: nDimP2 ! dimension+2
integer :: nDimP3 ! dimension+3
integer :: iLamMu ! laminar viscosity index in pv
integer :: iEdyMu ! eddy viscosity index in pv
integer :: iBlend ! blending function index in pv

! Global parameter
real(8), parameter :: Eps = 1d-10 ! epsilon for limiter
integer, parameter :: C1 = 10 ! (1~10) ! for reference Omega
integer, parameter :: C2 = 5 ! (2~5) ! for reference EddyMu
real(8), parameter :: pi = 4d0*datan(1d0)

! SendRecv data type
type t_SendRecv
  real(8), allocatable :: pv(:) ! primitive variable for send / recv
  real(8), allocatable :: gpv(:) ! gradient primitive variable for send / recv
  real(8), allocatable :: limiter(:) ! limiter variable for send /recv
end type

! Mixture data type
type t_Mixt
  ! Mixture constant
  real(8) :: Gamma ! heat capacity ratio
  real(8) :: Gamm1 ! gamma - 1
  real(8) :: GasConstant ! specific gas constant
  real(8) :: Prandtl_Lam ! laminar Prandtl number
  real(8) :: Prandtl_Turb ! turbulent Prandtl number
  ! Time data
  real(8) :: Time ! mixture solution time
  real(8) :: TimeStepMin ! Minimum time step
  real(8) :: TImeStepDualTime ! time step for dual time step
  real(8), allocatable :: TimeStep(:) ! time step
  ! Conservative variable data
  integer :: nCVar ! number of conservative variable
  real(8), allocatable :: cv(:,:) ! conservative variable multiplied by cell volume
  real(8), allocatable :: cv0(:,:) ! old conservative variable muliplied by cell volume
  real(8), allocatable :: Residual(:,:) ! residual at cell
  ! Primitive variable data
  integer :: nPVar ! number of primitive variable
  real(8), allocatable :: pv(:,:) ! primitive variable
  ! Gradient variable data
  integer :: nGVar ! number of gradient variable
  real(8), allocatable :: gpv(:,:,:) ! gradient primitive variable
  real(8), allocatable :: limiter(:,:) ! limiter variable
  real(8), allocatable :: pvmin(:,:) ! min primitive variable
  real(8), allocatable :: pvmax(:,:) ! max primitive variable
  ! Connected domain send/recv mixture data
  type(t_SendRecv), allocatable :: ConnSend(:) ! connected domain send mixture data
  type(t_SendRecv), allocatable :: ConnRecv(:) ! connected domain recv mixture data
  integer, allocatable :: Send_Req(:), Send_Stat(:,:) ! send request / status
  integer, allocatable :: Recv_Req(:), Recv_Stat(:,:) ! recv request / status
  ! Implicit time scheme variable
  real(8) :: L2Norm ! L2Norm error for checking convergence
  real(8), allocatable :: diag(:,:) ! diagonal term for LUSGS
  real(8), allocatable :: ram(:,:) ! overrelaxation factor 'r_A' for LUSGS
  real(8), allocatable :: dwst(:,:)
  real(8), allocatable :: dw(:,:)
  ! FSI simulation variable
  logical :: MembraneBroken = .false. ! for checking membrane
  integer, allocatable :: bFlag(:) ! boundary face state flag
  real(8), allocatable :: mFlux(:) ! boundary face mass flux
  real(8), allocatable :: fTemp(:) ! boundary face adiabatic flame temperature
end type

contains

subroutine DelMixt(Mixt)
  type(t_Mixt), intent(inout) :: Mixt

  integer :: iConn

  ! Time data
  if( allocated( Mixt%TimeStep ) ) deallocate( Mixt%TimeStep )
  ! Conservative variable data
  if( allocated( Mixt%cv ) ) deallocate( Mixt%cv )
  if( allocated( Mixt%cv0 ) ) deallocate( Mixt%cv0 )
  if( allocated( Mixt%Residual ) ) deallocate( Mixt%Residual )
  ! Primitive variable data
  if( allocated( Mixt%pv ) ) deallocate( Mixt%pv )
  ! Gradient variable data
  if( allocated( Mixt%gpv ) ) deallocate( Mixt%gpv )
  if( allocated( Mixt%limiter ) ) deallocate( Mixt%limiter )
  if( allocated( Mixt%pvmin ) ) deallocate( Mixt%pvmin )
  if( allocated( Mixt%pvmax ) ) deallocate( Mixt%pvmax )
  ! SendRecv mixture data
  if( allocated( Mixt%ConnSend ) ) then
    do  iConn = 1, size(Mixt%ConnSend)
      if( allocated( Mixt%ConnSend(iConn)%pv ) ) deallocate( Mixt%ConnSend(iConn)%pv )
      if( allocated( Mixt%ConnSend(iConn)%gpv ) ) deallocate( Mixt%ConnSend(iConn)%gpv )
      if( allocated( Mixt%ConnSend(iConn)%limiter ) ) deallocate( Mixt%ConnSend(iConn)%limiter )
    end do
    deallocate( Mixt%ConnSend )
  end if
  if( allocated( Mixt%ConnRecv ) ) then
    do  iConn = 1, size(Mixt%ConnRecv)
      if( allocated( Mixt%ConnRecv(iConn)%pv ) ) deallocate( Mixt%ConnRecv(iConn)%pv )
      if( allocated( Mixt%ConnRecv(iConn)%gpv ) ) deallocate( Mixt%ConnRecv(iConn)%gpv )
      if( allocated( Mixt%ConnRecv(iConn)%limiter ) ) deallocate( Mixt%ConnRecv(iConn)%limiter )
    end do
    deallocate( Mixt%ConnRecv )
  end if
  if( allocated( Mixt%Send_Req ) ) deallocate( Mixt%Send_Req )
  if( allocated( Mixt%Send_Stat ) ) deallocate( Mixt%Send_Stat )
  if( allocated( Mixt%Recv_Req ) ) deallocate( Mixt%Recv_Req )
  if( allocated( Mixt%Recv_Stat ) ) deallocate( Mixt%Recv_Stat )
  ! Implicit time scheme data
  if( allocated( Mixt%diag ) ) deallocate( Mixt%diag )
  if( allocated( Mixt%ram )  ) deallocate( Mixt%ram )
  if( allocated( Mixt%dwst ) ) deallocate( Mixt%dwst )
  if( allocated( Mixt%dw )   ) deallocate( Mixt%dw )
  ! FSBI simulation data
  if( allocated( Mixt%bFlag ) ) deallocate( Mixt%bFlag )
  if( allocated( Mixt%mFlux ) ) deallocate( Mixt%mFlux )
  if( allocated( Mixt%fTemp ) ) deallocate( Mixt%fTemp )


end subroutine

subroutine AllocMixtData(Mixt,Grid,Conf)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf

  integer :: iConn
  integer :: nSendCell
  integer :: nRecvCell

  ! Allocate / initialize time data
  allocate( Mixt%TimeStep(Grid%nCellWhole) )
  Mixt%TimeStep(:) = 0d0

  ! Allocate / initialize conservative data
  allocate( Mixt%cv(0:Mixt%nCVar-1,Grid%nCellWhole) )
  allocate( Mixt%cv0(0:Mixt%nCVar-1,Grid%nCellWhole) )
  allocate( Mixt%Residual(0:Mixt%nCVar-1,Grid%nCellWhole) )
  Mixt%cv(:,:) = 0d0
  Mixt%cv0(:,:) = 0d0
  Mixt%Residual(:,:) = 0d0

  ! Allocate / initialize primitive data
  allocate( Mixt%pv(0:Mixt%nPVar-1,Grid%nCellWhole) )
  Mixt%pv(:,:) = 0d0

  ! Allocate / initialize gradient data
  allocate( Mixt%gpv(nDim,0:Mixt%nGVar-1,Grid%nCellWhole) )
  allocate( Mixt%limiter(0:Mixt%nGVar-1,Grid%nCellWhole) )
  select case(Conf%iLimit)
  case(4,5) ! for node base limiter (MLP)
    allocate( Mixt%pvmin(0:Mixt%nGVar-1,Grid%nNodeTotal) )
    allocate( Mixt%pvmax(0:Mixt%nGVar-1,Grid%nNodeTotal) )
  case default ! for cell base limiter (Barth)
    allocate( Mixt%pvmin(0:Mixt%nGVar-1,Grid%nCellWhole) )
    allocate( Mixt%pvmax(0:Mixt%nGVar-1,Grid%nCellWhole) )
  end select
  Mixt%gpv(:,:,:) = 0d0
  Mixt%limiter(:,:) = 0d0
  Mixt%pvmin(:,:) = 0d0
  Mixt%pvmax(:,:) = 0d0

  ! Allocate send / recv data
  allocate( Mixt%ConnSend(Grid%nConn) )
  allocate( Mixt%ConnRecv(Grid%nConn) )

  ! Allocate send / recv variable data
  do iConn = 1, Grid%nConn
    nSendCell = Grid%ConnSend(iConn)%nCell
    allocate( Mixt%ConnSend(iConn)%pv(0:nSendCell*Mixt%nPVar-1) )
    allocate( Mixt%ConnSend(iConn)%gpv(0:nSendCell*Mixt%nGVar*nDim-1) )
    allocate( Mixt%ConnSend(iConn)%limiter(0:nSendCell*Mixt%nGVar-1) )
    Mixt%ConnSend(iConn)%pv(:) = 0d0
    Mixt%ConnSend(iConn)%gpv(:) = 0d0
    Mixt%ConnSend(iConn)%limiter(:) = 0d0
    nRecvCell = Grid%ConnRecv(iConn)%nCell
    allocate( Mixt%ConnRecv(iConn)%pv(0:nRecvCell*Mixt%nPVar-1) )
    allocate( Mixt%ConnRecv(iConn)%gpv(0:nRecvCell*Mixt%nGVar*nDim-1) )
    allocate( Mixt%ConnRecv(iConn)%limiter(0:nRecvCell*Mixt%nGVar-1) )
    Mixt%ConnRecv(iConn)%pv(:) = 0d0
    Mixt%ConnRecv(iConn)%gpv(:) = 0d0
    Mixt%ConnRecv(iConn)%limiter(:) = 0d0
  end do

  ! Allocate send / recv request & status
  allocate( Mixt%Send_Req(Grid%nConn), Mixt%Send_Stat(MPI_STATUS_SIZE,Grid%nConn) )
  allocate( Mixt%Recv_Req(Grid%nConn), Mixt%Recv_Stat(MPI_STATUS_SIZE,Grid%nConn) )

  ! Allocate FSI simulation data
  allocate( Mixt%bFlag(Grid%nFaceBound) )
  allocate( Mixt%mFlux(Grid%nFaceBound) )
  allocate( Mixt%fTemp(Grid%nFaceBound) )
  Mixt%bFlag(:) = 0
  Mixt%mFlux(:) = 0d0
  Mixt%fTemp(:) = 0d0

  ! Allocate implicit time scheme data
  allocate( Mixt%diag(0:Mixt%nCVar-1,Grid%nCellWhole) )
  allocate( Mixt%ram(0:Mixt%nCVar-1,Grid%nFace) )
  allocate( Mixt%dwst(0:Mixt%nCVar-1,Grid%nCellWhole) )
  allocate( Mixt%dw(0:Mixt%nCVar-1,Grid%nCellWhole) )
  Mixt%diag(:,:) = 0d0
  Mixt%ram(:,:) = 0d0
  Mixt%dwst(:,:) = 0d0
  Mixt%dw(:,:) = 0d0

end subroutine

subroutine SetMixtFromInitial(Mixt,Grid,Conf)
  type(t_Mixt), intent(out) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf

  integer :: iCell
  real(8) :: Temperature, Pressure
  real(8) :: Mach, AOA, SoundSpeed
  real(8) :: Density, Velocity(3), Energy
  real(8) :: Kine, Omega, LamMu, EddyMu

  ! Set dimension
  nDim = Grid%nDim
  nDimP1 = Grid%nDim+1
  nDimP2 = Grid%nDim+2
  nDimP3 = Grid%nDim+3

  ! Set gamma / gamma-1 / gas constant
  ! Laminar / turbulent Prandtl number
  Mixt%Gamma = Conf%Gamma
  Mixt%Gamm1 = Conf%Gamma-1d0
  Mixt%GasConstant = Conf%GasConstant
  Mixt%Prandtl_Lam = Conf%Prandtl_Lam
  Mixt%Prandtl_Turb = Conf%Prandtl_Turb

  ! Set number of variable
  select case(Conf%FlowModel)
  case(1) ! Euler equation
    Mixt%nCVar = Grid%nDim+2 ! rho,rho*u,rho*v,rho*w,rho*E
    Mixt%nPVar = Grid%nDim+2 ! T,u,v,w,p
    Mixt%nGVar = Grid%nDim+2 ! T,u,v,w,p
  case(2) ! N-S equation
    Mixt%nCVar = Grid%nDim+2 ! rho,rho*u,rho*v,rho*w,rho*E
    Mixt%nPVar = Grid%nDim+4 ! T,u,v,w,p,LamMu,EddyMu
    Mixt%nGVar = Grid%nDim+2 ! T,u,v,w,p
  case(3) ! RANS equation
    Mixt%nCVar = Grid%nDim+4 ! rho,rho*u,rho*v,rho*w,rho*E,rho*K,rho*omega
    Mixt%nPVar = Grid%nDim+7 ! T,u,v,w,p,K,omega,LamMu,EddyMu,Blend
    Mixt%nGVar = Grid%nDim+4 ! T,u,v,w,p,K,omega
  end select

  ! Set index of the primitive variable
  select case(Conf%FlowModel)
  case(2) ! N-S equation
    iLamMu = Grid%nDim+2 ! laminar viscosity index in pv
    iEdyMu = Grid%nDim+3 ! eddy viscosity index in pv
  case(3) ! RANS equation
    iLamMu = Grid%nDim+4 ! laminar viscosity index in pv
    iEdyMu = Grid%nDim+5 ! eddy viscosity index in pv
    iBlend = Grid%nDim+6 ! blending function index in pv
  end select

  ! Allocate mixture / sendrecv variable
  call AllocMixtData(Mixt,Grid,Conf)

  ! Set initial time
  Mixt%Time = 0d0

  ! Initial primitive variables
  Temperature = Conf%iniTemp
  Pressure = Conf%iniPres
  Mach = Conf%iniMach
  AOA = Conf%iniAngl
  SoundSpeed = dsqrt( Mixt%Gamma * Mixt%GasConstant * Temperature )
  Velocity(1) = Mach * SoundSpeed * dcos( pi * AOA / 180d0 )
  Velocity(2) = Mach * SoundSpeed * dsin( pi * AOA / 180d0 )
  Velocity(3) = 0d0
  Energy = Mixt%GasConstant * Temperature / Mixt%Gamm1 + 0.5d0 * sum( Velocity(1:nDim)**2 )
  Density = Pressure / ( Mixt%GasConstant * Temperature )

  ! Set conservative variable from initial condition
  do iCell = 1, Grid%nCell
    Mixt%cv(0,iCell) = Density * Grid%vol(iCell)
    Mixt%cv(1:nDim,iCell) = Density * Velocity(1:nDim) * Grid%vol(iCell)
    Mixt%cv(nDimP1,iCell) = Density * Energy * Grid%vol(iCell)
  end do

  ! Turbulence transport equation
  if( Conf%FlowModel == 3 ) then
    ! Initial primitive variables
    Omega = C1 * dsqrt( sum( Velocity(1:nDim)**2 ) ) / Conf%Length_Ref
    call SetSutherlandLaw(Conf,Temperature,LamMu)
    EddyMu = LamMu * 10d0**(-C2)
    Kine = EddyMu * Omega / Density
    ! Set conservative variable from initial condition
    do iCell = 1, Grid%nCell
      Mixt%cv(nDimP2,iCell) = Density * Kine * Grid%vol(iCell)
      Mixt%cv(nDimP3,iCell) = Density * Omega * Grid%vol(iCell)
    end do
  end if

#ifdef test
  select case(Conf%FlowModel)
  case(1); call SetSODShockTube(Mixt,Grid)
  case(2); call SetVisShockTube(Mixt,Grid)
  end select
#endif

  ! Postprocess for mixture
  call Postprocessing(Mixt,Grid,Conf)

end subroutine

subroutine SetSODShockTube(Mixt,Grid)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid

  integer :: iCell
  real(8) :: Density, Velocity(nDim), Pressure, Energy

  ! Set conservative variable ( modified SOD's shock tube problem )
  do iCell = 1, Grid%nCell

    ! Initial left / right states
    if( Grid%cen(1,iCell) .lt. 0.3d0 ) then
      Density = 1d0
      Velocity(:) = 0d0
      Velocity(1) = 0.75d0
      Pressure = 1d0
    else
      Density = 0.125d0
      Velocity(:) = 0d0
      Pressure = 0.1d0
    end if
    Energy = Pressure / ( Density * Mixt%Gamm1 ) + 0.5d0 * sum( Velocity(1:nDim)**2 )

    ! Set conservative variables
    Mixt%cv(0,iCell) = Density * Grid%vol(iCell)
    Mixt%cv(1:nDim,iCell) = Density * Velocity(1:nDim) * Grid%vol(iCell)
    Mixt%cv(nDimP1,iCell) = Density * Energy * Grid%vol(iCell)

  end do

end subroutine

subroutine SetVisShockTube(Mixt,Grid)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid

  integer :: iCell
  real(8) :: Density, Velocity(nDim), Pressure, Energy

  ! Set conservative variable ( viscous shock tube problem )
  do iCell = 1, Grid%nCell

    ! Initial left / right states
    if( Grid%cen(1,iCell) .lt. 0.5d0 ) then
      Density = 120d0
      Velocity(:) = 0d0
      Pressure = Density / Mixt%Gamma
    else
      Density = 1.2d0
      Velocity(:) = 0d0
      Pressure = Density / Mixt%Gamma
    end if
    Energy = Pressure / ( Density * Mixt%Gamm1 ) + 0.5d0 * sum( Velocity(1:nDim)**2 )

    ! Set conservative variables
    Mixt%cv(0,iCell) = Density * Grid%vol(iCell)
    Mixt%cv(1:nDim,iCell) = Density * Velocity(1:nDim) * Grid%vol(iCell)
    Mixt%cv(nDimP1,iCell) = Density * Energy * Grid%vol(iCell)

  end do

end subroutine

subroutine SetMixtFromFile(Mixt,Grid,Conf)
  type(t_Mixt), intent(out) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf

  integer :: io, ier
  integer :: iVar, iCell, iCellGlobal, iCellLocal
  real(8) :: Temperature, Velocity(3), Pressure
  real(8) :: Density, Energy, Kine, Omega
  integer :: Header(5)
  character(32) :: VarName(0:11)
  integer :: intsize, real8size, CellReal8Type
  integer(kind=MPI_OFFSET_KIND) :: disp
  integer, allocatable :: CellFileMap(:)
  real(8), allocatable :: DataReal8(:)

  ! Set dimension
  nDim = Grid%nDim
  nDimP1 = Grid%nDim+1
  nDimP2 = Grid%nDim+2
  nDimP3 = Grid%nDim+3

  ! Set gamma / gamma-1 / gas constant
  ! Laminar / turbulent Prandtl number
  Mixt%Gamma = Conf%Gamma
  Mixt%Gamm1 = Conf%Gamma-1d0
  Mixt%GasConstant = Conf%GasConstant
  Mixt%Prandtl_Lam = Conf%Prandtl_Lam
  Mixt%Prandtl_Turb = Conf%Prandtl_Turb

  ! Set number of variable
  select case(Conf%FlowModel)
  case(1) ! Euler equation
    Mixt%nCVar = Grid%nDim+2 ! rho,rho*u,rho*v,rho*w,rho*E
    Mixt%nPVar = Grid%nDim+2 ! T,u,v,w,p
    Mixt%nGVar = Grid%nDim+2 ! T,u,v,w,p
  case(2) ! N-S equation
    Mixt%nCVar = Grid%nDim+2 ! rho,rho*u,rho*v,rho*w,rho*E
    Mixt%nPVar = Grid%nDim+4 ! T,u,v,w,p,LamMu,EddyMu
    Mixt%nGVar = Grid%nDim+2 ! T,u,v,w,p
  case(3) ! RANS equation
    Mixt%nCVar = Grid%nDim+4 ! rho,rho*u,rho*v,rho*w,rho*E,rho*K,rho*omega
    Mixt%nPVar = Grid%nDim+7 ! T,u,v,w,p,K,omega,LamMu,EddyMu,Blend
    Mixt%nGVar = Grid%nDim+4 ! T,u,v,w,p,K,omega
  end select

  ! Set index of the primitive variable
  select case(Conf%FlowModel)
  case(2) ! N-S equation
    iLamMu = Grid%nDim+2 ! laminar viscosity index in pv
    iEdyMu = Grid%nDim+3 ! eddy viscosity index in pv
  case(3) ! RANS equation
    iLamMu = Grid%nDim+4 ! laminar viscosity index in pv
    iEdyMu = Grid%nDim+5 ! eddy viscosity index in pv
    iBlend = Grid%nDim+6 ! blending function index in pv
  end select

  ! Allocate mixture / sendrecv variable
  call AllocMixtData(Mixt,Grid,Conf)

  ! Create cell file map
  allocate( CellFileMap(Grid%nCell) )
  iCell = 0
  do iCellGlobal = 1, size(Grid%Global_to_Local_Cell)
    iCellLocal = Grid%Global_to_Local_Cell(iCellGlobal)
    if( iCellLocal == 0 .or. iCellLocal > Grid%nCell ) cycle
    iCell = iCell + 1
    CellFileMap(iCell) = iCellGlobal-1
  end do

  ! Set type size
  call MPI_Type_size(MPI_INTEGER,intsize,ier)
  call MPI_Type_size(MPI_REAL8,real8size,ier)

  ! Create cell real8 type for reading real8 data in cell
  call MPI_Type_create_indexed_block(Grid%nCell,1,CellFileMap,MPI_REAL8,CellReal8Type,ier)
  call MPI_Type_commit(CellReal8Type,ier)

  ! Open FLUS file using MPI-IO
  call MPI_File_open(MPI_COMM_WORLD,trim(Conf%GridFileName),MPI_MODE_RDONLY,MPI_INFO_NULL,io,ier)

  ! Read header / variable name
  call MPI_File_read_all(io,Header,5,MPI_INTEGER,MPI_STATUS_IGNORE,ier)
  call MPI_File_read_all(io,VarName,32*Header(5),MPI_CHARACTER,MPI_STATUS_IGNORE,ier)

  ! Read time data
  call MPI_File_read_all(io,Mixt%Time,1,MPI_REAL8,MPI_STATUS_IGNORE,ier)

  ! Set displacement
  disp = 5 * intsize + 32 * Header(5) + real8size ! header / variable / time
  disp = disp + Grid%nDim * Grid%nNodeGlobal * real8size ! node coordinates
  disp = disp + Grid%nNodeCell * Grid%nCellGlobal * intsize ! node connectivity of the cell
  disp = disp + Grid%nNodeFace * Grid%nFaceGlobal * intsize ! node connectivity of the face

  ! Allocate / initialize real8 data array
  allocate( DataReal8(Grid%nCell) )
  DataReal8(:) = 0d0

  ! Loop over primitive variable
  do iVar = 0, Mixt%nPVar-1

    ! Set displacement & view for reading primitive variable of the cell
    if( iVar == 0 ) disp = disp + Grid%nFaceGlobal * intsize
    if( iVar /= 0 ) disp = disp + Grid%nCellGlobal * real8size
    call MPI_File_set_view(io,disp,MPI_REAL8,CellReal8Type,'native',MPI_INFO_NULL,ier)

    ! Read & set primitive variable
    call MPI_File_read_all(io,DataReal8,Grid%nCell,MPI_REAL8,MPI_STATUS_IGNORE,ier)
    do iCell = 1, Grid%nCell
      iCellGlobal = CellFileMap(iCell)+1
      iCellLocal = Grid%Global_to_Local_Cell(iCellGlobal)
      Mixt%pv(iVar,iCellLocal) = DataReal8(iCell)
    end do

  end do

  ! Deallocate real8 data array
  deallocate( DataReal8 )

  ! Deallocate cell file map
  deallocate( CellFileMap )

  ! Close FLUS file
  call MPI_File_close(io,ier)

  ! Set conservative variable from primitive variable
  do iCell = 1, Grid%nCell
    Temperature = Mixt%pv(0,iCell)
    Velocity(1:nDim) = Mixt%pv(1:nDim,iCell)
    Pressure = Mixt%pv(nDimP1,iCell)
    Energy = Mixt%GasConstant * Temperature / Mixt%Gamm1 + 0.5d0 * sum( Velocity(1:nDim)**2 )
    Density = Pressure / ( Mixt%GasConstant * Temperature )
    Mixt%cv(0,iCell) = Density * Grid%vol(iCell)
    Mixt%cv(1:nDim,iCell) = Density * Velocity(1:nDim) * Grid%vol(iCell)
    Mixt%cv(nDimP1,iCell) = Density * Energy * Grid%vol(iCell)
  end do

  ! Turbulence transport equation
  if( Conf%FlowModel == 3 ) then
    ! Initial primitive variables
    Temperature = Mixt%pv(0,iCell)
    Pressure = Mixt%pv(nDimP1,iCell)
    Density = Pressure / ( Mixt%GasConstant * Temperature )
    Kine = Mixt%pv(nDimP2,iCell)
    Omega = Mixt%pv(nDimP3,iCell)
    ! Set conservative variable from initial condition
    do iCell = 1, Grid%nCell
      Mixt%cv(nDimP2,iCell) = Density * Kine * Grid%vol(iCell)
      Mixt%cv(nDimP3,iCell) = Density * Omega * Grid%vol(iCell)
    end do
  end if

  ! Postprocess for mixture
  call Postprocessing(Mixt,Grid,Conf)

end subroutine

subroutine SetMixtFromBackup(Mixt,Grid,Conf,MixtBackup,GridBackup)
  type(t_Mixt), intent(out) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf
  type(t_Mixt), intent(in) :: MixtBackup
  type(t_Grid), intent(in) :: GridBackup

  integer :: iCell, error, ier
  real(8) :: Temperature, Pressure
  real(8) :: Density, Velocity(3), Energy

  ! Check flow model
  if( Conf%FlowModel /= 1 ) then ! this routine assume that the flow model is euler.
    write(*,*) 'FLUS >> Flow model must be euler in SetMixtFromBackup routine'
    call MPI_Abort(MPI_COMM_WORLD,error,ier)
  end if

  ! Set gamma / gamma-1 / gas constant
  Mixt%Gamma = MixtBackup%Gamma
  Mixt%Gamm1 = MixtBackup%Gamm1
  Mixt%GasConstant = MixtBackup%GasConstant

  ! Set number of variable
  Mixt%nCVar = Grid%nDim+2 ! rho,rho*u,rho*v,rho*w,rho*E
  Mixt%nPVar = Grid%nDim+2 ! T,u,v,w,p
  Mixt%nGVar = Grid%nDim+2 ! T,u,v,w,p

  ! Allocate mixture / sendrecv variable
  call AllocMixtData(Mixt,Grid,Conf)

  ! Set initial time
  Mixt%Time = MixtBackup%Time

  ! Transfer primitive variable from backup
  call DataTransfer(GridBackup,Grid,Mixt%nPVar,MixtBackup%pv,Mixt%pv)

  ! Set conservative variable from primitive variable
  do iCell = 1, Grid%nCell
    Temperature = Mixt%pv(0,iCell)
    Velocity(1:nDim) = Mixt%pv(1:nDim,iCell)
    Pressure = Mixt%pv(nDimP1,iCell)
    Energy = Mixt%GasConstant * Temperature / Mixt%Gamm1 + 0.5d0 * sum( Velocity(1:nDim)**2 )
    Density = Pressure / ( Mixt%GasConstant * Temperature )
    Mixt%cv(0,iCell) = Density * Grid%vol(iCell)
    Mixt%cv(1:nDim,iCell) = Density * Velocity(1:nDim) * Grid%vol(iCell)
    Mixt%cv(nDimP1,iCell) = Density * Energy * Grid%vol(iCell)
  end do

  ! Transfer membrane status
  Mixt%MembraneBroken = MixtBackup%MembraneBroken

  ! Postprocess for mixture
  call Postprocessing(Mixt,Grid,Conf)

end subroutine

subroutine Postprocessing(Mixt,Grid,Conf)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf

  !integer :: iCell

  ! Set first primitive variable (T,u,v,w,p,K,Omega)
  call SetPrimitive1(Mixt,Grid,Conf)

  ! Set gradient for 2nd-order
  if( Conf%iLimit > 1 ) call SetGradient_LS(Mixt,Grid)

  ! Set second primitive variable (LamMu,EddyMu,Blend)
  call SetPrimitive2(Mixt,Grid,Conf)

  ! Set primitive variable in ghost cell using MPI
  call SetMPIPrimitive(Mixt,Grid)

end subroutine

subroutine SetPrimitive1(Mixt,Grid,Conf)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf

  integer :: iCell
  real(8) :: Coeff1, Coeff2
  real(8) :: Density

  ! Set mean primitive variables
  do iCell = 1, Grid%nCell
    Coeff1 = 1d0 / Grid%vol(iCell)
    Coeff2 = 1d0 / Mixt%cv(0,iCell)
    Density = Mixt%cv(0,iCell) * Coeff1
    Mixt%pv(1:nDim,iCell) = Mixt%cv(1:nDim,iCell) * Coeff2
    Mixt%pv(nDimP1,iCell) = Mixt%Gamm1 * ( Mixt%cv(nDimP1,iCell) * Coeff1 - 0.5d0 * Density * sum( Mixt%pv(1:nDim,iCell)**2 ) )
    Mixt%pv(0,iCell) = Mixt%pv(nDimP1,iCell) / ( Density * Mixt%GasConstant )
  end do

  ! Set turb primitivie variables with lower bound
  ! And, reset conservative variables
  if( Conf%FlowModel == 3 ) then
    do iCell = 1, Grid%nCell
      Coeff2 = 1d0 / Mixt%cv(0,iCell)
      Mixt%pv(nDimP2,iCell) = dmax1( Mixt%cv(nDimP2,iCell) * Coeff2, 1d-10 )
      Mixt%pv(nDimP3,iCell) = dmax1( Mixt%cv(nDimP3,iCell) * Coeff2, 1d-10 )
      Mixt%cv(nDimP2,iCell) = Mixt%cv(0,iCell) * Mixt%pv(nDimP2,iCell)
      Mixt%cv(nDimP3,iCell) = Mixt%cv(0,iCell) * Mixt%pv(nDimP3,iCell)
    end do
  end if

end subroutine

subroutine SetPrimitive2(Mixt,Grid,Conf)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf

  integer :: iCell
  real(8) :: Temperature, Pressure, Density
  real(8) :: Kine, Omega, Dist
  real(8) :: LamMu, EddyMu, Blend

  select case(Conf%FlowModel)
  case(2) ! N-S equation
    do iCell = 1, Grid%nCell
      Temperature = Mixt%pv(0,iCell)
      call SetSutherlandLaw(Conf,Temperature,LamMu)
      Mixt%pv(iLamMu,iCell) = LamMu
      Mixt%pv(iEdyMu,iCell) = 0d0
    end do
  case(3) ! RANS equation
    do iCell = 1, Grid%nCell
      Temperature = Mixt%pv(0,iCell)
      Pressure = Mixt%pv(nDimP1,iCell)
      Density = Pressure / ( Mixt%GasConstant * Temperature )
      Kine = Mixt%pv(nDimP2,iCell)
      Omega = Mixt%pv(nDimP3,iCell)
      Dist = Grid%WallDist(iCell)
      call SetSutherlandLaw(Conf,Temperature,LamMu)
      call SetTurbSSTModel(Density,Kine,Omega,LamMu,Dist,Mixt%gpv(:,:,iCell),EddyMu,Blend)
      Mixt%pv(iLamMu,iCell) = LamMu
      Mixt%pv(iEdyMu,iCell) = EddyMu
      Mixt%pv(iBlend,iCell) = Blend
    end do
  end select

end subroutine

subroutine SetSutherlandLaw(Conf,Temperature,Laminar_Viscosity)
  type(t_Conf), intent(in) :: Conf
  real(8), intent(in) :: Temperature
  real(8), intent(out) :: Laminar_Viscosity

  real(8), parameter :: Coeff = 3d0 / 2d0

  ! Compute laminar viscosity from Sutherland law
  Laminar_Viscosity = Conf%Viscosity_Ref * ( Temperature / Conf%Temperature_Ref )**Coeff &
                    * ( Conf%Temperature_Ref + Conf%Temperature_Eff ) / ( Temperature + Conf%Temperature_Eff )

#ifdef test
  ! Constant laminar viscosity
  Laminar_Viscosity = Conf%Viscosity_Ref
#endif

end subroutine

subroutine SetTurbSSTModel(Density,Kine,Omega,LamMu,Dist,Grad,EddyMu,Blend)
  real(8), intent(in) :: Density
  real(8), intent(in) :: Kine
  real(8), intent(in) :: Omega
  real(8), intent(in) :: LamMu
  real(8), intent(in) :: Dist
  real(8), intent(in) :: Grad(:,0:)
  real(8), intent(out) :: EddyMu
  real(8), intent(out) :: Blend

  integer :: iDim
  real(8) :: arg2A, arg2B, arg2, arg1
  real(8) :: CDkw, StrainMag, zeta, F2
  real(8), parameter :: inv3 = 1d0 / 3d0
  real(8), parameter :: Five9 = 5d0 / 9d0

  ! Compute CDkw (positive part of the cross-diffusion)
  CDkw = sum( Grad(1:nDim,nDimP2)*Grad(1:nDim,nDimP3) )
  CDkw = CDkw*2d0*Density*0.856d0/Omega
  CDkw = dmax1(CDkw,1d-10)

  ! Compute auxiliary function F1(Blend)
  arg2A = dsqrt(Kine)/(0.09d0*Omega*Dist)
  arg2B = 500d0*LamMu/(Density*Omega*Dist**2)
  arg2 = dmax1(arg2A,arg2B)
  arg1 = dmin1(arg2,4d0*Density*0.856d0*Kine/(CDkw*Dist**2))
  Blend = dtanh(arg1**4)

  ! Compute auxiliary function F2
  arg2 = dmax1(2d0*arg2A,arg2B)
  F2 = dtanh(arg2**2)

  StrainMag = 0d0
  ! Compute strain magnitude (add diagonal part)
  do iDim = 1, nDim
    StrainMag = StrainMag + 2d0*Grad(iDim,iDim)**2
  end do

  ! Compute strain magnitude (add off diagonals)
  StrainMag = StrainMag + ( Grad(1,2) + Grad(2,1) )**2
  if( nDim == 3 ) then
    StrainMag = StrainMag + ( Grad(1,3) + Grad(3,1) )**2
    StrainMag = StrainMag + ( Grad(2,3) + Grad(3,2) )**2
  end if

  ! Compute strain magnitude
  StrainMag = dsqrt( StrainMag )

  ! Compute eddy viscosity
  zeta = dmax1(Five9*Omega,F2*StrainMag)
  EddyMu = Five9*Density*Kine/zeta

end subroutine

subroutine SetMPIPrimitive(Mixt,Grid)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid

  integer :: rank, srcs, dest, Tag, ier
  integer :: iVar, iCell, iCell_List
  integer :: iLocal, iConn, iDomain, nBuffer

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Non-blocking recv comm
  do iConn = 1, Grid%nConn
    iDomain = Grid%ConnList(iConn); srcs = iDomain-1; Tag = iDomain
    nBuffer = Grid%ConnRecv(iConn)%nCell * Mixt%nPVar
    call MPI_Irecv(Mixt%ConnRecv(iConn)%pv,nBuffer,MPI_REAL8,srcs,Tag,MPI_COMM_WORLD,Mixt%Recv_Req(iConn),ier)
  end do

  ! Set buffer send data
  do iConn = 1, Grid%nConn
    do iCell_List = 1, Grid%ConnSend(iConn)%nCell
      iCell = Grid%ConnSend(iConn)%Cell(iCell_List)
      do iVar = 0, Mixt%nPVar-1
        iLocal = Mixt%nPVar*(iCell_List-1)+iVar
        Mixt%ConnSend(iConn)%pv(iLocal) = Mixt%pv(iVar,iCell)
      end do
    end do
  end do

  ! Non-blocking send comm
  do iConn = 1, Grid%nConn
    iDomain = Grid%ConnList(iConn); dest = iDomain-1; Tag = rank+1
    nBuffer = Grid%ConnSend(iConn)%nCell * Mixt%nPVar
    call MPI_Isend(Mixt%ConnSend(iConn)%pv,nBuffer,MPI_REAL8,dest,Tag,MPI_COMM_WORLD,Mixt%Send_Req(iConn),ier)
  end do

  ! Wait for the set of non-blocking comm
  call MPI_Waitall(Grid%nConn,Mixt%Send_Req,Mixt%Send_Stat,ier)
  call MPI_Waitall(Grid%nConn,Mixt%Recv_Req,Mixt%Recv_Stat,ier)

  ! Set primitive variables from buffer recv
  do iConn = 1, Grid%nConn
    do iCell_List = 1, Grid%ConnRecv(iConn)%nCell
      iCell = Grid%ConnRecv(iConn)%Cell(iCell_List)
      do iVar = 0, Mixt%nPVar-1
        iLocal = Mixt%nPVar*(iCell_List-1)+iVar
        Mixt%pv(iVar,iCell) = Mixt%ConnRecv(iConn)%pv(iLocal)
      end do
    end do
  end do

end subroutine

subroutine ComputeTimeStep(Mixt,Grid,Conf,TargetTime)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf
  real(8), intent(in), optional :: TargetTime

  integer :: iDim, iLocal, iCell, iFace, ier
  real(8) :: Density, SoundSpeed, Visco
  real(8) :: DeltaArea, TimeStepMin
  real(8), parameter :: Coeff = 1d0
  real(8), parameter :: Four3 = 4d0 / 3d0

  ! Compute time step (Ref: J. Blazek, Computational Fluid Dynamics: Principles and Applications)
  do iCell = 1, Grid%nCell
    SoundSpeed = dsqrt( Mixt%Gamma * Mixt%GasConstant * Mixt%pv(0,iCell) )
    Mixt%TimeStep(iCell) = 0d0
    if( Conf%FlowModel > 1 ) then
      Density = Mixt%pv(nDimP1,iCell) / ( Mixt%GasConstant * Mixt%pv(0,iCell) )
      Visco = dmax1(Four3,Mixt%Gamma) / Density * ( Mixt%pv(iLamMu,iCell) / Mixt%Prandtl_Lam + Mixt%pv(iEdyMu,iCell) / Mixt%Prandtl_Turb )
    end if
    do iDim = 1, nDim
      DeltaArea = 0d0
      do iLocal = 1, Grid%nFaceCell
        iFace = Grid%c2f(iLocal,iCell)
        DeltaArea = DeltaArea + dabs( Grid%fn(iDim,iFace) * Grid%fa(iFace) )
      end do
      DeltaArea = 0.5d0 * DeltaArea
      Mixt%TimeStep(iCell) = Mixt%TimeStep(iCell) + ( dabs(Mixt%pv(iDim,iCell)) + SoundSpeed ) * DeltaArea
      if( Conf%FlowModel > 1 ) Mixt%TimeStep(iCell) = Mixt%TimeStep(iCell) + Coeff * Visco * DeltaArea**2 / Grid%Vol(iCell)
    end do
    Mixt%TimeStep(iCell) = Conf%CFL * Grid%vol(iCell) / Mixt%TimeStep(iCell)
  end do

  ! Set minimum time step and MPI comm
  TimeStepMin = minval( Mixt%TimeStep(1:Grid%nCell) )
  call MPI_Allreduce(TimeStepMin,Mixt%TimeStepMin,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ier)

  ! Set time step for unsteady simulation
  if( Conf%Unsteady .and. present(TargetTime) ) then
    if( Mixt%Time + Mixt%TimeStepMin > TargetTime ) Mixt%TimeStepMin = TargetTime - Mixt%Time
    Mixt%Time = Mixt%Time + Mixt%TimeStepMin
  end if

  ! Set global time step
  if( .not. Conf%LocalTimeStep ) then
    Mixt%TimeStep(:) = Mixt%TimeStepMin
  end if

end subroutine

subroutine ComputeError(Mixt,Grid)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid

  integer :: iCell, ier
  real(8) :: L2Norm

  L2Norm = 0d0
  do iCell = 1, Grid%nCell
    ! Check error with pressure
    L2Norm = L2Norm + dabs( Mixt%dw(nDimP1,iCell) )**2
  end do
  call MPI_Allreduce(L2Norm,Mixt%L2Norm,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
  Mixt%L2Norm = dsqrt( Mixt%L2Norm / Grid%nCellGlobal )

end subroutine

subroutine ComputeResidual(Mixt,Grid,Conf)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf

  integer :: iVar, iFace, iCell, iCell1, iCell2
  real(8) :: PRateMin, PMin, PMinRate
  real(8) :: drl(nDim), pvl(0:Mixt%nPVar-1)
  real(8) :: drr(nDim), pvr(0:Mixt%nPVar-1)
  real(8) :: Flux(0:Mixt%nCVar-1)
  real(8) :: Srcs(0:Mixt%nCVar-1)

  ! Set primitive boundary condition
  call SetBoundaryCondition(Mixt,Grid,Conf)

  if( Conf%iLimit > 1 ) then
    ! Set limiter barth or mlp
    select case(Conf%iLimit)
    case(2:3); call SetLimiter_Barth(Mixt,Grid,Conf)
    case(4:5); call SetLimiter_MLP(Mixt,Grid,Conf)
    end select
    ! Set wall gradient / limiter boundary condition
    call SetBoundaryGradient(Mixt,Grid,Conf)
  end if

  ! Initialize residual
  Mixt%Residual(:,:) = 0d0

  ! Compute convective flux
  do iFace = 1, Grid%nFace
    iCell1 = Grid%f2c(1,iFace)
    iCell2 = Grid%f2c(2,iFace)

    ! Set vector cell center to face center
    drl(:) = Grid%fc(:,iFace) - Grid%cen(:,iCell1)
    drr(:) = Grid%fc(:,iFace) - Grid%cen(:,iCell2)

    ! Set left / right primitive
    pvl(:) = Mixt%pv(:,iCell1)
    pvr(:) = Mixt%pv(:,iCell2)
    do iVar = 0, nDimP1 ! do not reconstruct K / Omega
      pvl(iVar) = pvl(iVar) + Mixt%limiter(iVar,iCell1) * sum( Mixt%gpv(:,iVar,iCell1) * drl(:) )
      pvr(iVar) = pvr(iVar) + Mixt%limiter(iVar,iCell2) * sum( Mixt%gpv(:,iVar,iCell2) * drr(:) )
    end do

    ! Mean flow convective flux
    select case(Conf%iFlux)
    case(1) ! Roe convective flux
      call ComputeConvRoe(pvl(:),pvr(:),Grid%fn(:,iFace),Grid%fv(iFace),Flux(:))
    case(2) ! RoeM2 convective flux
      PRateMin = 1d0
      call SetRateMinPressure(Mixt,Grid,iCell1,PRateMin)
      call SetRateMinPressure(Mixt,Grid,iCell2,PRateMin)
      call ComputeConvRoeM2(pvl(:),pvr(:),Grid%fn(:,iFace),Grid%fv(iFace),PRateMin,Flux(:))
    case(3) ! AUSM+ convective flux
      call ComputeConvAUSMP(pvl(:),pvr(:),Grid%fn(:,iFace),Grid%fv(iFace),Flux(:))
    case(4) ! AUSMPW+ convective flux
      PMin = 2d0**1023
      call SetMinPressure(Mixt,Grid,iCell1,iCell2,PMin)
      call SetMinPressure(Mixt,Grid,iCell1,iCell2,PMin)
      PMinRate = PMin / dmin1( Mixt%pv(nDimP1,iCell1), Mixt%pv(nDimP1,iCell2) )
      call ComputeConvAUSMPWP(pvl(:),pvr(:),Grid%fn(:,iFace),Grid%fv(iFace),PMinRate,Flux(:))
    end select
    Mixt%Residual(0:nDimP1,iCell1) = Mixt%Residual(0:nDimP1,iCell1) + Flux(0:nDimP1) * Grid%fa(iFace)
    Mixt%Residual(0:nDimP1,iCell2) = Mixt%Residual(0:nDimP1,iCell2) - Flux(0:nDimP1) * Grid%fa(iFace)

    ! Turbulence convective flux
    if( Conf%FlowModel == 3 ) then
      call ComputeTurbConvSca(pvl(:),pvr(:),Grid%fn(:,iFace),Grid%fv(iFace),Flux(:))
      Mixt%Residual(nDimP2:nDimP3,iCell1) = Mixt%Residual(nDimP2:nDimP3,iCell1) + Flux(nDimP2:nDimP3) * Grid%fa(iFace)
      Mixt%Residual(nDimP2:nDimP3,iCell2) = Mixt%Residual(nDimP2:nDimP3,iCell2) - Flux(nDimP2:nDimP3) * Grid%fa(iFace)
    end if

  end do

  ! Compute diffusive flux
  if( Conf%FlowModel > 1 ) then
    do iFace = 1, Grid%nFace
      iCell1 = Grid%f2c(1,iFace)
      iCell2 = Grid%f2c(2,iFace)

      ! Mean flow diffusive flux
      call ComputeDiffAvgGrad(Mixt%pv(:,iCell1),Mixt%pv(:,iCell2),Mixt%gpv(:,:,iCell1),Mixt%gpv(:,:,iCell2),Mixt%limiter(:,iCell1),Mixt%limiter(:,iCell2),&
                              Grid%cen(:,iCell1),Grid%cen(:,iCell2),Grid%fn(:,iFace),Flux(:))
      Mixt%Residual(0:nDimP1,iCell1) = Mixt%Residual(0:nDimP1,iCell1) + Flux(:) * Grid%fa(iFace)
      Mixt%Residual(0:nDimP1,iCell2) = Mixt%Residual(0:nDimP1,iCell2) - Flux(:) * Grid%fa(iFace)

      ! Turbulence diffusive flux
      if( Conf%FlowModel == 3 ) then
        call ComputeTurbDiffAvgGrad(Mixt%pv(:,iCell1),Mixt%pv(:,iCell2),Mixt%gpv(:,:,iCell1),Mixt%gpv(:,:,iCell2),Mixt%limiter(:,iCell1),Mixt%limiter(:,iCell2),&
                                    Grid%cen(:,iCell1),Grid%cen(:,iCell2),Grid%fn(:,iFace),Flux(:))
        Mixt%Residual(nDimP2:nDimP3,iCell1) = Mixt%Residual(nDimP2:nDimP3,iCell1) + Flux(nDimP2:nDimP3) * Grid%fa(iFace)
        Mixt%Residual(nDimP2:nDimP3,iCell2) = Mixt%Residual(nDimP2:nDimP3,iCell2) - Flux(nDimP2:nDimP3) * Grid%fa(iFace)
      end if
    end do
  end if

  ! Compute turbulence source residual
  if( Conf%FlowModel == 3 ) then
    do iCell = 1, Grid%nCell
      call ComputeTurbSrcs(Mixt%pv(:,iCell),Mixt%gpv(:,:,iCell),Srcs(:))
      Mixt%Residual(nDimP2:nDimP3,iCell) = Mixt%Residual(nDimP2:nDimP3,iCell) + Srcs(nDimP2:nDimP3) * Grid%vol(iCell)
    end do
  end if

  ! Compute axisymmetric source
  if( Conf%Axisymmetric ) then
    do iCell = 1, Grid%nCell
      call ComputeAxiSrcs(Mixt%pv(:,iCell),Grid%cen(:,iCell),Srcs(:))
      Mixt%Residual(0:nDimP1,iCell) = Mixt%Residual(0:nDimP1,iCell) + Srcs(0:nDimP1) * Grid%vol(iCell)
    end do
  end if

end subroutine

subroutine SetRateMinPressure(Mixt,Grid,iCell1,PRateMin)
  type(t_Mixt), intent(in) :: Mixt
  type(t_Grid), intent(in) :: Grid
  integer, intent(in) :: iCell1
  real(8), intent(inout) :: PRateMin

  integer :: iLocal, iCell2
  real(8) :: Prate

  if( iCell1 <= Grid%nCell ) then ! this have to be changed
    do iLocal = 1, Grid%nFaceCell
      iCell2 = Grid%c2c(iLocal,iCell1)
      Prate = Mixt%pv(nDimP1,iCell2) / Mixt%pv(nDimP1,iCell1)
      PRateMin = dmin1( PRateMin, Prate, 1d0 / Prate )
    end do
  end if

end subroutine

subroutine SetMinPressure(Mixt,Grid,iCell1,iCell2,PMin)
  type(t_Mixt), intent(in) :: Mixt
  type(t_Grid), intent(in) :: Grid
  integer, intent(in) :: iCell1
  integer, intent(in) :: iCell2
  real(8), intent(inout) :: PMin

  integer :: iLocal, iCell3

  if( iCell1 <= Grid%nCell ) then ! this have to be changed
    do iLocal = 1, Grid%nFaceCell
      iCell3 = Grid%c2c(iLocal,iCell1)
      if( iCell3 /= iCell2 ) then
        PMin = dmin1( PMin, Mixt%pv(nDimP1,iCell3) )
      end if
    end do
  end if

end subroutine

subroutine ComputeDiagonalTerm(Mixt,Grid,Conf)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf

  integer :: iCell, iFace
  integer :: iCell1, iCell2
  real(8) :: AvgTemperature, AvgSoundSpeed
  real(8) :: ProjVelocity(2), AvgProjVelocity
  real(8) :: Density, Visco, Length
  real(8) :: F1, beta_Blend
  real(8), parameter :: Four3 = 4d0 / 3d0

  ! Compute mean flow diagonal term
  do iCell = 1, Grid%nCell
    Mixt%diag(0:nDimP1,iCell) = Grid%Vol(iCell) / Mixt%TimeStep(iCell)
  end do

  ! Compute turbulent diagonal term
  if( Conf%FlowModel == 3 ) then
    do iCell = 1, Grid%nCell
      F1 = Mixt%pv(iBlend,iCell); beta_Blend = F1*3d0/40d0 + (1d0-F1)*0.0828d0
      Mixt%diag(nDimP2,iCell) = Grid%Vol(iCell) / Mixt%TimeStep(iCell) + 0.09d0 * Mixt%pv(nDimP3,iCell) * Grid%Vol(iCell)
      Mixt%diag(nDimP3,iCell) = Grid%Vol(iCell) / Mixt%TimeStep(iCell) + 2d0 * beta_Blend * Mixt%pv(nDimP3,iCell) * Grid%Vol(iCell)
    end do
  end if

  do iFace = 1, Grid%nFace
    iCell1 = Grid%f2c(1,iFace)
    iCell2 = Grid%f2c(2,iFace)

    ! Arithmetic Averaging for Face Values
    AvgTemperature = 0.5d0 * ( Mixt%pv(0,iCell1) + Mixt%pv(0,iCell2) )
    AvgSoundSpeed = dsqrt( Mixt%gamma * Mixt%GasConstant * AvgTemperature )
    ProjVelocity(1) = sum( Grid%fn(:,iFace) * Mixt%pv(1:nDim,iCell1) )
    ProjVelocity(2) = sum( Grid%fn(:,iFace) * Mixt%pv(1:nDim,iCell2) )
    AvgProjVelocity = 0.5d0 * sum( ProjVelocity(:) ) - Grid%fv(iFace)

    ! Compute mean flow convective spectral radii (Over Relaxation Factor Range 1 < w <= 2)
    Mixt%ram(0:nDimP1,iFace) = ( dabs(AvgProjVelocity) + AvgSoundSpeed ) * Grid%fa(iFace) * 1.01d0

    ! Compute mean flow viscous spectral radii
    if( Conf%FlowModel > 1 ) then
      Density = Mixt%pv(nDimP1,iCell) / ( Mixt%GasConstant * Mixt%pv(0,iCell) )
      Visco = dmax1(Four3,Mixt%Gamma) / Density * ( Mixt%pv(iLamMu,iCell) / Mixt%Prandtl_Lam + Mixt%pv(iEdyMu,iCell) / Mixt%Prandtl_Turb )
      Length = dsqrt( sum( ( Grid%cen(:,iCell2) - Grid%cen(:,iCell1) )**2 ) )
      Mixt%ram(0:nDimP1,iFace) = Mixt%ram(0:nDimP1,iFace) + Visco * Grid%fa(iFace) / Length
    end if

    if( Conf%FlowModel == 3 ) then
      ! Compute turbulent convective spectral radii
      Mixt%ram(nDimP2:nDimP3,iFace) = dabs(AvgProjVelocity) * Grid%fa(iFace) * 1.01d0
      ! Compute turbulent viscous spectral radii
      Mixt%ram(nDimP2:nDimP3,iFace) = Mixt%ram(nDimP2:nDimP3,iFace) + Visco * Grid%fa(iFace) / Length
    end if

    ! Update diagonal term
    Mixt%diag(:,iCell1) = Mixt%diag(:,iCell1) + 0.5d0 * Mixt%ram(:,iFace)
    Mixt%diag(:,iCell2) = Mixt%diag(:,iCell2) + 0.5d0 * Mixt%ram(:,iFace)
  end do

end subroutine

subroutine ComputeDeltaFlux(Mixt,Grid,iCell,iFace,dwst,df )
  type(t_Mixt), intent(in) :: Mixt
  type(t_Grid), intent(in) :: Grid
  integer, intent(in) :: iCell
  integer, intent(in) :: iFace
  real(8), intent(in) :: dwst(0:Mixt%nCVar-1)
  real(8), intent(out) :: df(0:Mixt%nCVar-1)
  real(8) :: vol, UnitNormal(nDim), Area, FaceVelocity
  real(8) :: Temperature_l, Velocity_l(nDim), Pressure_l, Density_l, Energy_l, Enthalpy_l, ProjVelocity_l, RelVelocity_l
  real(8) :: Temperature_r, Velocity_r(nDim), Pressure_r, Density_r, Energy_r, Enthalpy_r, ProjVelocity_r, RelVelocity_r
  real(8) :: Flux_l(0:Mixt%nCVar-1), Flux_r(0:Mixt%nCVar-1)
  real(8) :: cv_l(0:Mixt%nCVar-1), cv_r(0:Mixt%nCVar-1)

  ! Set initial grid value
  UnitNormal(:) = Grid%fn(:,iFace)
  Area = Grid%fa(iFace)
  FaceVelocity = Grid%fv(iFace)

  ! Set former step primitive variables -> this comment needs to be changed
  Temperature_l = Mixt%pv(0,iCell)
  Velocity_l(1:nDim) = Mixt%pv(1:nDim,iCell)
  Pressure_l = Mixt%pv(nDimP1,iCell)
  Density_l = Pressure_l / ( Mixt%GasConstant * Temperature_l )
  Energy_l = Mixt%GasConstant * Temperature_l / Mixt%Gamm1 + 0.5d0 * sum( Velocity_l(:)**2 )
  Enthalpy_l = Energy_l + Mixt%GasConstant * Temperature_l
  ProjVelocity_l = sum( Velocity_l(:) * UnitNormal(:) )
  RelVelocity_l = ProjVelocity_l - FaceVelocity

  ! Set volume of the dummy cell
  if ( iCell > Grid%nCellTotal ) then
    vol = Grid%vol(Grid%f2c(1,iFace))
  else
    vol = Grid%vol(iCell)
  end if

  ! Compute left conservative
  cv_l(0) = Density_l
  cv_l(1:nDim) = Density_l * Velocity_l(1:nDim)
  cv_l(nDimP1) = Density_l * Energy_l

  ! Compute right conservative
  cv_r(0) = cv_l(0) + dwst(0)
  cv_r(1:nDim) = cv_l(1:nDim) + dwst(1:nDim)
  cv_r(nDimP1) = cv_l(nDimP1) + dwst(nDimP1)

  ! Compute right primitive
  Density_r = cv_r(0)
  Velocity_r(1:nDim) = cv_r(1:nDim) / cv_r(0)
  Energy_r = cv_r(nDimP1) / cv_r(0)
  Temperature_r = Mixt%Gamm1 / Mixt%GasConstant * ( Energy_r - 0.5d0 * sum( Velocity_r(:)**2 ) )
  Pressure_r = Density_r * Mixt%GasConstant * Temperature_r
  Enthalpy_r = Energy_r + Mixt%GasConstant * Temperature_r
  ProjVelocity_r = sum( Velocity_r(:) * UnitNormal(:) )
  RelVelocity_r = ProjVelocity_r - FaceVelocity

  ! Compute left flux
  FLux_l(0) = Density_l * RelVelocity_l
  Flux_l(1:nDim) = Density_l * RelVelocity_l * Velocity_l(:) + Pressure_l * UnitNormal(:)
  Flux_l(nDimP1) = Density_l * RelVelocity_l * Enthalpy_l + Pressure_l * FaceVelocity

  ! Compute right flux
  FLux_r(0) = Density_r * RelVelocity_r
  Flux_r(1:nDim) = Density_r * RelVelocity_r * Velocity_r(:) + Pressure_r * UnitNormal(:)
  Flux_r(nDimP1) = Density_r * RelVelocity_r * Enthalpy_r + Pressure_r * FaceVelocity

  ! Compute delta flux
  df(:) = Flux_r(:) - Flux_l(:)
  df(:) = df(:) * Area

end subroutine

subroutine SetGradient_LS(Mixt,Grid)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid

  integer :: iCell, jCell, kCell
  integer :: iLocal, jLocal, iVar, counter
  real(8) :: r11, r12, r22, r13, r23, r33
  real(8) :: alpha1, alpha2, alpha3, beta, Weight
  real(8) :: DeltaPrim, DeltaCG(nDim), Coeff(nDim)
  logical :: isBoundary

  ! Set gradient as zero
  Mixt%gpv(:,:,:) = 0d0

  ! Compute gradient using Least-Square
  do iCell = 1, Grid%nCell

    ! Initialize coefficient
    r11 = 0d0; r12 = 0d0; r22 = 0d0
    r13 = 0d0; r23 = 0d0; r33 = 0d0

    counter = 0
    do iLocal = 1, Grid%nFaceCell
      jCell = Grid%c2c(iLocal,iCell)
      if( jCell > Grid%nCellTotal ) cycle
      counter = counter + 1

      ! Vector icell center to jcell center
      DeltaCG(:) = Grid%cen(:,jCell) - Grid%cen(:,iCell)

      ! summation for coefficient
      r11 = r11 + DeltaCG(1) * DeltaCG(1)
      r12 = r12 + DeltaCG(1) * DeltaCG(2)
      r22 = r22 + DeltaCG(2) * DeltaCG(2)
      if( nDim == 3 ) then
        r13 = r13 + DeltaCG(1) * DeltaCG(3)
        r23 = r23 + DeltaCG(2) * DeltaCG(3)
        r33 = r33 + DeltaCG(3) * DeltaCG(3)
      end if
    end do

    isBoundary = .false.
    if( counter <= nDim ) then
      isBoundary = .true.
      do iLocal = 1, Grid%nFaceCell
        jCell = Grid%c2c(iLocal,iCell)
        if( jCell > Grid%nCellTotal ) cycle
        do jLocal = 1, Grid%nFaceCell
          kCell = Grid%c2c(jLocal,jCell)
          if( kCell > Grid%nCellTotal .or. kCell < 1 .or. kCell == iCell ) cycle

          ! Vector icell center to jcell center
          DeltaCG(:) = Grid%cen(:,kCell) - Grid%cen(:,iCell)

          ! summation for coefficient
          r11 = r11 + DeltaCG(1) * DeltaCG(1)
          r12 = r12 + DeltaCG(1) * DeltaCG(2)
          r22 = r22 + DeltaCG(2) * DeltaCG(2)
          if( nDim == 3 ) then
            r13 = r13 + DeltaCG(1) * DeltaCG(3)
            r23 = r23 + DeltaCG(2) * DeltaCG(3)
            r33 = r33 + DeltaCG(3) * DeltaCG(3)
          end if
        end do
      end do
    end if

    ! Compute coefficient
    r11 = dsqrt( r11 )
    r12 = r12 / r11
    r22 = dsqrt( r22 - r12 * r12 )
    if( nDim == 3 ) then
      r13 = r13 / r11
      r23 = ( r23 - r12 * r13 ) / r22
      r33 = dsqrt( r33 - r13 * r13 - r23 * r23 )
      beta = ( r12 * r23 - r13 * r22 ) / ( r11 * r22 )
    end if

    do iLocal = 1, Grid%nFaceCell
      jCell = Grid%c2c(iLocal,iCell)
      if( jCell > Grid%nCellTotal ) cycle

      ! Vector icell center to jcell center
      DeltaCG(:) = Grid%cen(:,jCell) - Grid%cen(:,iCell)

      ! Weight of jCell
      alpha1 = DeltaCG(1) / r11**2
      alpha2 = ( DeltaCG(2) - r12 / r11 * DeltaCG(1) ) / r22**2
      if( nDim == 3 ) then
        alpha3 = ( DeltaCG(3) - r23 / r22 * DeltaCG(2) + beta * DeltaCG(1) ) / r33**2
      end if

      ! Distance weight
      !Weight = 1d0 / dsqrt( sum( DeltaCG(:)**2 ) ) ! ocsillation occurs when MLP used
      Weight = 1d0

      ! Compute total coefficient
      if( nDim == 2 ) then
        Coeff(1) = Weight * ( alpha1 - r12 / r11 * alpha2 )
        Coeff(2) = Weight * ( alpha2 )
      else if( nDim == 3 ) then
        Coeff(1) = Weight * ( alpha1 - r12 / r11 * alpha2 + beta * alpha3 )
        Coeff(2) = Weight * ( alpha2 - r23 / r22 * alpha3 )
        Coeff(3) = Weight * ( alpha3 )
      end if

      ! Compute gradient
      do iVar = 0, Mixt%nGVar-1
        DeltaPrim = Mixt%pv(iVar,jCell) - Mixt%pv(iVar,iCell)
        Mixt%gpv(:,iVar,iCell) = Mixt%gpv(:,iVar,iCell) + Coeff(:) * DeltaPrim
      end do

      if( isBoundary ) then

        do jLocal = 1, Grid%nFaceCell
          kCell = Grid%c2c(jLocal,jCell)
          if( kCell > Grid%nCellTotal .or. kCell < 1 .or. kCell == iCell ) cycle

          ! Vector icell center to jcell center
          DeltaCG(:) = Grid%cen(:,kCell) - Grid%cen(:,iCell)

          ! Weight of kCell
          alpha1 = DeltaCG(1) / r11**2
          alpha2 = ( DeltaCG(2) - r12 / r11 * DeltaCG(1) ) / r22**2
          if( nDim == 3 ) then
            alpha3 = ( DeltaCG(3) - r23 / r22 * DeltaCG(2) + beta * DeltaCG(1) ) / r33**2
          end if

          ! Distance weight
          !Weight = 1d0 / dsqrt( sum( DeltaCG(:)**2 ) ) ! ocsillation occurs when MLP used
          Weight = 1d0

          ! Compute total coefficient
          if( nDim == 2 ) then
            Coeff(1) = Weight * ( alpha1 - r12 / r11 * alpha2 )
            Coeff(2) = Weight * ( alpha2 )
          else if( nDim == 3 ) then
            Coeff(1) = Weight * ( alpha1 - r12 / r11 * alpha2 + beta * alpha3 )
            Coeff(2) = Weight * ( alpha2 - r23 / r22 * alpha3 )
            Coeff(3) = Weight * ( alpha3 )
          end if

          ! Compute gradient
          do iVar = 0, Mixt%nGVar-1
            DeltaPrim = Mixt%pv(iVar,kCell) - Mixt%pv(iVar,iCell)
            Mixt%gpv(:,iVar,iCell) = Mixt%gpv(:,iVar,iCell) + Coeff(:) * DeltaPrim
          end do
        end do

      end if

    end do
  end do

  ! Set gradient in ghost cell using MPI
  call SetMPIGradient(Mixt,Grid)

end subroutine

subroutine SetMPIGradient(Mixt,Grid)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid

  integer :: rank, srcs, dest, Tag, ier
  integer :: iVar, iCell, iCell_List, iDim
  integer :: iLocal, iConn, iDomain, nBuffer

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Non-blocking recv comm
  do iConn = 1, Grid%nConn
    iDomain = Grid%ConnList(iConn); srcs = iDomain-1; Tag = iDomain
    nBuffer = Grid%ConnRecv(iConn)%nCell * Mixt%nGVar * nDim
    call MPI_Irecv(Mixt%ConnRecv(iConn)%gpv,nBuffer,MPI_REAL8,srcs,Tag,MPI_COMM_WORLD,Mixt%Recv_Req(iConn),ier)
  end do

  ! Set buffer send data
  do iConn = 1, Grid%nConn
    do iCell_List = 1, Grid%ConnSend(iConn)%nCell
      iCell = Grid%ConnSend(iConn)%Cell(iCell_List)
      do iVar = 0, Mixt%nGVar-1
        do iDim = 1, nDim
          iLocal = Mixt%nGVar*nDim*(iCell_List-1)+nDim*iVar+iDim-1
          Mixt%ConnSend(iConn)%gpv(iLocal) = Mixt%gpv(iDim,iVar,iCell)
        end do
      end do
    end do
  end do

  ! Non-blocking sent comm
  do iConn = 1, Grid%nConn
    iDomain = Grid%ConnList(iConn); dest = iDomain-1; Tag = rank+1
    nBuffer = Grid%ConnSend(iConn)%nCell * Mixt%nGVar * nDim
    call MPI_Isend(Mixt%ConnSend(iConn)%gpv,nBuffer,MPI_REAL8,dest,Tag,MPI_COMM_WORLD,Mixt%Send_Req(iConn),ier)
  end do

  ! Wait for the set of non-blocking comm
  call MPI_Waitall(Grid%nConn,Mixt%Send_Req,Mixt%Send_Stat,ier)
  call MPI_Waitall(Grid%nConn,Mixt%Recv_Req,Mixt%Recv_Stat,ier)

  ! Set gradient variables from buffer recv
  do iConn = 1, Grid%nConn
    do iCell_List = 1, Grid%ConnRecv(iConn)%nCell
      iCell = Grid%ConnRecv(iConn)%Cell(iCell_List)
      do iVar = 0, Mixt%nGVar-1
        do iDim = 1, nDim
          iLocal = Mixt%nGVar*nDim*(iCell_List-1)+nDim*iVar+iDim-1
          Mixt%gpv(iDim,iVar,iCell) = Mixt%ConnRecv(iConn)%gpv(iLocal)
        end do
      end do
    end do
  end do

end subroutine

subroutine SetLimiter_Barth(Mixt,Grid,Conf)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf

  integer :: iLocal, iVar
  integer :: iFace, iCell, iNode
  integer :: iCell1, iCell2
  real(8) :: DeltaPrim1, DeltaPrim2
  real(8) :: Limiter, Eps2

  ! Set limiter as one
  Mixt%limiter(:,:) = 1d0

  ! Set primitive mix / max
  Mixt%pvmin(:,:) = Mixt%pv(0:Mixt%nGVar-1,:)
  Mixt%pvmax(:,:) = Mixt%pv(0:Mixt%nGVar-1,:)
  do iFace = 1, Grid%nFace
    iCell1 = Grid%f2c(1,iFace)
    iCell2 = Grid%f2c(2,iFace)
    Mixt%pvmin(:,iCell1) = dmin1(Mixt%pvmin(:,iCell1),Mixt%pv(0:Mixt%nGVar-1,iCell2))
    Mixt%pvmax(:,iCell1) = dmax1(Mixt%pvmax(:,iCell1),Mixt%pv(0:Mixt%nGVar-1,iCell2))
    Mixt%pvmin(:,iCell2) = dmin1(Mixt%pvmin(:,iCell2),Mixt%pv(0:Mixt%nGVar-1,iCell1))
    Mixt%pvmax(:,iCell2) = dmax1(Mixt%pvmax(:,iCell2),Mixt%pv(0:Mixt%nGVar-1,iCell1))
  end do

  select case(Conf%iLimit)
  case(2) ! Barth limiter for unsteady
    do iCell = 1, Grid%nCell
      do iVar = 0, Mixt%nGVar-1
        do iLocal = 1, Grid%nNodeCell
          iNode = Grid%c2n(iLocal,iCell)
          DeltaPrim2 = sum( Mixt%gpv(:,iVar,iCell) * ( Grid%xyz(:,iNode) - Grid%cen(:,iCell) ) )
          if( DeltaPrim2 > Eps ) then
            DeltaPrim1 = Mixt%pvmax(iVar,iCell) - Mixt%pv(iVar,iCell)
            if( dabs( DeltaPrim1 ) < Eps ) DeltaPrim1 = 0d0
            Limiter = dmin1( 1d0, DeltaPrim1 / DeltaPrim2 )
          else if( DeltaPrim2 < -Eps ) then
            DeltaPrim1 = Mixt%pvmin(iVar,iCell) - Mixt%pv(iVar,iCell)
            if( dabs( DeltaPrim1 ) < Eps ) DeltaPrim1 = 0d0
            Limiter = dmin1( 1d0, DeltaPrim1 / DeltaPrim2 )
          else
            Limiter = 0d0
          end if
          Mixt%limiter(iVar,iCell) = dmin1( Mixt%limiter(iVar,iCell), Limiter )
        end do
      end do
    end do

  case(3) ! Venkatakrishnan limiter for steady
    do iCell = 1, Grid%nCell
      Eps2 = ( Conf%vkt * (Grid%vol(iCell))**(1d0/nDim) )**3
      do iVar = 0, Mixt%nGVar-1
        do iLocal = 1, Grid%nNodeCell
          iNode = Grid%c2n(iLocal,iCell)
          DeltaPrim2 = sum( Mixt%gpv(:,iVar,iCell) * ( Grid%xyz(:,iNode) - Grid%cen(:,iCell) ) )
          if( DeltaPrim2 > Eps ) then
            DeltaPrim1 = Mixt%pvmax(iVar,iCell) - Mixt%pv(iVar,iCell)
            if( dabs( DeltaPrim1 ) < Eps ) DeltaPrim1 = 0d0
            Limiter = dmin1( 1d0, ( DeltaPrim1**2 + Eps2 + 2d0 * DeltaPrim2 * DeltaPrim1 ) / ( DeltaPrim1**2 + 2d0 * DeltaPrim2**2 + DeltaPrim1 * DeltaPrim2 + Eps2 ) )
          else if( DeltaPrim2 < -Eps ) then
            DeltaPrim1 = Mixt%pvmin(iVar,iCell) - Mixt%pv(iVar,iCell)
            if( dabs( DeltaPrim1 ) < Eps ) DeltaPrim1 = 0d0
            Limiter = dmin1( 1d0, ( DeltaPrim1**2 + Eps2 + 2d0 * DeltaPrim2 * DeltaPrim1 ) / ( DeltaPrim1**2 + 2d0 * DeltaPrim2**2 + DeltaPrim1 * DeltaPrim2 + Eps2 ) )
          else
            Limiter = 0d0
          end if
          Mixt%limiter(iVar,iCell) = dmin1( Mixt%limiter(iVar,iCell), Limiter )
        end do
      end do
    end do

  end select

  ! Set limiter in ghost cell using MPI
  call SetMPILimiter(Mixt,Grid)

end subroutine

subroutine SetLimiter_MLP(Mixt,Grid,Conf)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf

  integer :: iLocal, iVar
  integer :: iCell, iNode
  real(8) :: DeltaPrim1, DeltaPrim2
  real(8) :: Limiter, Eps2

  ! Set limiter as one
  Mixt%limiter(:,:) = 1d0

  ! Set primitive min / max
  Mixt%pvmin(:,:) =  2d0**1023
  Mixt%pvmax(:,:) = -2d0**1023
  do iCell = 1, Grid%nCellWhole
    do iLocal = 1, Grid%nNodeCell
      iNode = Grid%c2n(iLocal,iCell)
      if( iNode < 1 ) cycle
      Mixt%pvmin(:,iNode) = dmin1(Mixt%pvmin(:,iNode),Mixt%pv(0:Mixt%nGVar-1,iCell))
      Mixt%pvmax(:,iNode) = dmax1(Mixt%pvmax(:,iNode),Mixt%pv(0:Mixt%nGVar-1,iCell))
    end do
  end do

  select case(Conf%iLimit)
  case(4) ! MLP-u1 limiter for unsteady
    do iCell = 1, Grid%nCell
      do iVar = 0, Mixt%nGVar-1
        do iLocal = 1, Grid%nNodeCell
          iNode = Grid%c2n(iLocal,iCell)
          DeltaPrim2 = sum( Mixt%gpv(:,iVar,iCell) * ( Grid%xyz(:,iNode) - Grid%cen(:,iCell) ) )
          if( DeltaPrim2 > Eps ) then
            DeltaPrim1 = Mixt%pvmax(iVar,iNode) - Mixt%pv(iVar,iCell)
            if( dabs( DeltaPrim1 ) < Eps ) DeltaPrim1 = 0d0
            Limiter = dmin1( 1d0, DeltaPrim1 / DeltaPrim2 )
          else if( DeltaPrim2 < -Eps ) then
            DeltaPrim1 = Mixt%pvmin(iVar,iNode) - Mixt%pv(iVar,iCell)
            if( dabs( DeltaPrim1 ) < Eps ) DeltaPrim1 = 0d0
            Limiter = dmin1( 1d0, DeltaPrim1 / DeltaPrim2 )
          else
            Limiter = 0d0
          end if
          Mixt%limiter(iVar,iCell) = dmin1( Mixt%limiter(iVar,iCell), Limiter )
        end do
      end do
    end do

  case(5) ! MLP-u2 limiter for steady
    do iCell = 1, Grid%nCell
      Eps2 = ( Conf%vkt * (Grid%vol(iCell))**(1d0/nDim) )**3
      do iVar = 0, Mixt%nGVar-1
        do iLocal = 1, Grid%nNodeCell
          iNode = Grid%c2n(iLocal,iCell)
          DeltaPrim2 = sum( Mixt%gpv(:,iVar,iCell) * ( Grid%xyz(:,iNode) - Grid%cen(:,iCell) ) )
          if( DeltaPrim2 > Eps ) then
            DeltaPrim1 = Mixt%pvmax(iVar,iNode) - Mixt%pv(iVar,iCell)
            if( dabs( DeltaPrim1 ) < Eps ) DeltaPrim1 = 0d0
            Limiter = dmin1( 1d0, ( DeltaPrim1**2 + Eps2 + 2d0 * DeltaPrim2 * DeltaPrim1 ) / ( DeltaPrim1**2 + 2d0 * DeltaPrim2**2 + DeltaPrim1 * DeltaPrim2 + Eps2 ) )
          else if( DeltaPrim2 < -Eps ) then
            DeltaPrim1 = Mixt%pvmin(iVar,iNode) - Mixt%pv(iVar,iCell)
            if( dabs( DeltaPrim1 ) < Eps ) DeltaPrim1 = 0d0
            Limiter = dmin1( 1d0, ( DeltaPrim1**2 + Eps2 + 2d0 * DeltaPrim2 * DeltaPrim1 ) / ( DeltaPrim1**2 + 2d0 * DeltaPrim2**2 + DeltaPrim1 * DeltaPrim2 + Eps2 ) )
          else
            Limiter = 0d0
          end if
          Mixt%limiter(iVar,iCell) = dmin1( Mixt%limiter(iVar,iCell), Limiter )
        end do
      end do
    end do

  end select

  ! Set limiter in ghost cell using MPI
  call SetMPILimiter(Mixt,Grid)

end subroutine

subroutine SetMPILimiter(Mixt,Grid)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid

  integer :: rank, srcs, dest, Tag, ier
  integer :: iVar, iCell, iCell_List
  integer :: iLocal, iConn, iDomain, nBuffer

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Non-blocking recv comm
  do iConn = 1, Grid%nConn
    iDomain = Grid%ConnList(iConn); srcs = iDomain-1; Tag = iDomain
    nBuffer = Grid%ConnRecv(iConn)%nCell * Mixt%nGVar
    call MPI_Irecv(Mixt%ConnRecv(iConn)%limiter,nBuffer,MPI_REAL8,srcs,Tag,MPI_COMM_WORLD,Mixt%Recv_Req(iConn),ier)
  end do

  ! Set buffer send data
  do iConn = 1, Grid%nConn
    do iCell_List = 1, Grid%ConnSend(iConn)%nCell
      iCell = Grid%ConnSend(iConn)%Cell(iCell_List)
      do iVar = 0, Mixt%nGVar-1
        iLocal = Mixt%nGVar*(iCell_List-1)+iVar
        Mixt%ConnSend(iConn)%limiter(iLocal) = Mixt%limiter(iVar,iCell)
      end do
    end do
  end do

  ! Non-blocking send comm
  do iConn = 1, Grid%nConn
    iDomain = Grid%ConnList(iConn); dest = iDomain-1; Tag = rank+1
    nBuffer = Grid%ConnSend(iConn)%nCell * Mixt%nGVar
    call MPI_Isend(Mixt%ConnSend(iConn)%limiter,nBuffer,MPI_REAL8,dest,Tag,MPI_COMM_WORLD,Mixt%Send_Req(iConn),ier)
  end do

  ! Wait for the set of non-blocking comm
  call MPI_Waitall(Grid%nConn,Mixt%Send_Req,Mixt%Send_Stat,ier)
  call MPI_Waitall(Grid%nConn,Mixt%Recv_Req,Mixt%Recv_Stat,ier)

  ! Set primitive variables from buffer recv
  do iConn = 1, Grid%nConn
    do iCell_List = 1, Grid%ConnRecv(iConn)%nCell
      iCell = Grid%ConnRecv(iConn)%Cell(iCell_List)
      do iVar = 0, Mixt%nGVar-1
        iLocal = Mixt%nGVar*(iCell_List-1)+iVar
        Mixt%limiter(iVar,iCell) = Mixt%ConnRecv(iConn)%limiter(iLocal)
      end do
    end do
  end do

end subroutine

subroutine SetBoundaryCondition(Mixt,Grid,Conf)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf

  integer :: iBC

  ! Check membrane broking before BC
  call CheckMembrane(Mixt,Grid,Conf)

  ! Set primitive in dummy cell
  do iBC = 1, Conf%nBC
    select case(Conf%BCtype(iBC))
    case(1) ! inviscid wall
      call BC_Inviscid_Wall(Mixt,Grid,Conf,iBC)
    case(2) ! viscous wall
      call BC_Viscous_Wall(Mixt,Grid,Conf,iBC)
    case(3) ! pressure inlet
      call BC_Pressure_Inlet(Mixt,Grid,Conf,iBC)
    case(4) ! injection inlet
      call BC_Injection_Inlet(Mixt,Grid,Conf,iBC)
    case(5) ! pressure outlet
      call BC_Pressure_Outlet(Mixt,Grid,Conf,iBC)
    case(6) ! pressure farfield
      call BC_Pressure_Farfield(Mixt,Grid,Conf,iBC)
    case(7) ! Symmetric Wall
      call BC_Inviscid_Wall(Mixt,Grid,Conf,iBC)
    case(8) ! Axisymmetric Wall
      call BC_Inviscid_Wall(Mixt,Grid,Conf,iBC)
    case(9) ! membrane nozzle
      call BC_Membrane_Nozzle(Mixt,Grid,Conf,iBC)
    case(10) ! propellant
      call BC_Propellant(Mixt,Grid,Conf,iBC)
    case(11) ! ignitor
      call BC_Ignitor(Mixt,Grid,Conf,iBC)
    case(12) ! pressure sensor
      call BC_Pressure_Sensor(Mixt,Grid,Conf,iBC)
    case(13) ! nozzle throat
      call BC_Inviscid_Wall(Mixt,Grid,Conf,iBC)
    case(14) ! propellant APN model
      call BC_APN_Model(Mixt,Grid,Conf,iBC)
    end select
  end do

end subroutine

subroutine CheckMembrane(Mixt,Grid,Conf)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf

  integer :: iList, List1, List2
  integer :: iBC, iFace, iCell

  integer :: rank, ier
  real(8) :: p_avg, p_sum, a_sum, Area
  real(8) :: p_sum_total, a_sum_total

  ! if membrane is already broken, return
  if( Mixt%MembraneBroken ) return

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  do iBC = 1, Conf%nBC
    if( Conf%BCType(iBC) == 12 .and. Conf%bv(1,iBC) > 0 ) then
      List1 = Grid%ListIndex(iBC)
      List2 = Grid%ListIndex(iBC+1) - 1
      p_sum = 0d0; a_sum = 0d0; p_avg = 0
      do iList = List1, List2
        iFace = Grid%List(iList)
        iCell = Grid%f2c(1,iFace)
        if( Conf%Axisymmetric ) then
          Area = Grid%fa(iFace) * 2d0 * pi * Grid%fc(2,iFace)
        else ! not axisymmetric case
          Area = Grid%fa(iFace)
        end if
        p_sum = p_sum + Mixt%pv(Grid%nDim+1,iCell) * Area
        a_sum = a_sum + Area
      end do
      call MPI_Allreduce(p_sum,p_sum_total,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call MPI_Allreduce(a_sum,a_sum_total,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      if( a_sum_total > 0 ) p_avg = p_sum_total / a_sum_total
      if( p_avg >= Conf%bv(1,iBC) ) then
        Mixt%MembraneBroken = .true.
        if( rank == 0 ) write(*,*) 'FLUS >> Membrane Nozzle Broken'
      end if
    end if
  end do

  ! set MPI barrier
  call MPI_Barrier(MPI_COMM_WORLD,ier)

end subroutine

subroutine BC_Inviscid_Wall(Mixt,Grid,Conf,iBC)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf
  integer, intent(in) :: iBC

  integer :: iList, List1, List2
  integer :: iBFace, iCell1, iCell2

  real(8) :: RelVelocity

  List1 = Grid%ListIndex(iBC)
  List2 = Grid%ListIndex(iBC+1) - 1

  ! Inviscid Wall Implementation( Moving )
  do iList = List1, List2
    iBFace = Grid%List(iList)
    iCell1 = Grid%f2c(1,iBFace)
    iCell2 = Grid%f2c(2,iBFace)

    ! Mean flow boundary condition
    RelVelocity= sum( Mixt%pv(1:nDim,iCell1) * Grid%fn(:,iBFace) ) - Grid%fv(iBFace)
    Mixt%pv(0,iCell2) = Mixt%pv(0,iCell1)
    Mixt%pv(1:nDim,iCell2) = Mixt%pv(1:nDim,iCell1) - 2d0 * RelVelocity * Grid%fn(:,iBFace)
    Mixt%pv(nDimP1,iCell2) = Mixt%pv(nDimP1,iCell1)

    ! Turbulence boundary condition
    if( Conf%FlowModel == 3 ) then
      Mixt%pv(nDimP2,iCell2) = Mixt%pv(nDimP2,iCell1)
      Mixt%pv(nDimP3,iCell2) = Mixt%pv(nDimP3,iCell1)
    end if

    ! Viscosity boundary condition
    select case(Conf%FlowModel)
    case(2) ! N-S equation
      Mixt%pv(iLamMu,iCell2) = Mixt%pv(iLamMu,iCell1)
      Mixt%pv(iEdyMu,iCell2) = 0d0
    case(3) ! RANS equation
      Mixt%pv(iLamMu,iCell2) = Mixt%pv(iLamMu,iCell1)
      Mixt%pv(iEdyMu,iCell2) = Mixt%pv(iEdyMu,iCell1)
      Mixt%pv(iBlend,iCell2) = Mixt%pv(iBlend,iCell1)
    end select

  end do

end subroutine

subroutine BC_Viscous_Wall(Mixt,Grid,Conf,iBC)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf
  integer, intent(in) :: iBC

  integer :: iList, List1, List2
  integer :: iBFace, iCell1, iCell2

  real(8) :: Temperature, Pressure, Density
  real(8) :: Dist, LamMu, Omega

  List1 = Grid%ListIndex(iBC)
  List2 = Grid%ListIndex(iBC+1) - 1

  ! Viscous Wall Implementation( Not Moving )
  do iList = List1, List2
    iBFace = Grid%List(iList)
    iCell1 = Grid%f2c(1,iBFace)
    iCell2 = Grid%f2c(2,iBFace)

    ! Mean flow boundary condition
    Mixt%pv(0,iCell2) = Mixt%pv(0,iCell1)
    Mixt%pv(1:nDim,iCell2) = -Mixt%pv(1:nDim,iCell1)
    Mixt%pv(nDimP1,iCell2) = Mixt%pv(nDimP1,iCell1)

    ! Turbulence boundary condition
    if( Conf%FlowModel == 3 ) then
      Dist = Grid%WallDist(iCell1)
      Temperature = Mixt%pv(0,iCell1)
      Pressure = Mixt%pv(nDimP1,iCell1)
      Density = Pressure / ( Mixt%GasConstant * Temperature )
      LamMu = Mixt%pv(iLamMu,iCell1)
      Omega = 60d0*LamMu/(Density*0.075d0*Dist**2)
      Mixt%pv(nDimP2,iCell2) = -Mixt%pv(nDimP2,iCell1)
      Mixt%pv(nDimP3,iCell2) = 2d0 * Omega - Mixt%pv(nDimP3,iCell1)
    end if

    ! Viscosity boundary condition
    select case(Conf%FlowModel)
    case(2) ! N-S equation
      Mixt%pv(iLamMu,iCell2) = Mixt%pv(iLamMu,iCell1)
      Mixt%pv(iEdyMu,iCell2) = 0d0
    case(3) ! RANS equation
      Mixt%pv(iLamMu,iCell2) = Mixt%pv(iLamMu,iCell1)
      Mixt%pv(iEdyMu,iCell2) = -Mixt%pv(iEdyMu,iCell1)
      Mixt%pv(iBlend,iCell2) = Mixt%pv(iBlend,iCell1)
    end select

  end do

end subroutine

subroutine BC_Pressure_Inlet(Mixt,Grid,Conf,iBC)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf
  integer, intent(in) :: iBC

  integer :: iList, List1, List2
  integer :: iBFace, iCell1, iCell2

  real(8) :: Total_Pressure, Total_Temperature
  real(8) :: Temperature, Velocity(nDim), Pressure
  real(8) :: SoundSpeed, Mach
  real(8) :: Density, Kine, Omega, LamMu, EddyMu

  List1 = Grid%ListIndex(iBC)
  List2 = Grid%ListIndex(iBC+1) - 1

  ! Set total pressure / temperature
  Total_Pressure = Conf%bv(1,iBC)
  Total_Temperature = Conf%bv(2,iBC)

  if( Conf%bv(3,iBC) >= 1d0 ) then ! Supersonic Inlet Implementation( Not Moving )

    ! Precompute primitive variables
    Mach = Conf%bv(3,iBC)
    Pressure = Total_Pressure / ( 1d0 + 0.5d0 * Mixt%Gamm1 * Mach**2 )**( Mixt%Gamma / Mixt%Gamm1 )
    Temperature = Total_Temperature / ( 1d0 + 0.5d0 * Mixt%Gamm1 * Mach**2 )
    Density = Pressure / ( Mixt%GasConstant * Temperature )
    SoundSpeed = dsqrt( Mixt%Gamma * Mixt%GasConstant * Temperature )
    Velocity(1) = Mach * SoundSpeed * dcos( pi * Conf%bv(4,iBC) / 180d0 )
    Velocity(2) = Mach * SoundSpeed * dsin( pi * Conf%bv(4,iBC) / 180d0 )
    if( nDim == 3 ) Velocity(3) = 0d0

    do iList = List1, List2
      iBFace = Grid%List(iList)
      iCell1 = Grid%f2c(1,iBFace)
      iCell2 = Grid%f2c(2,iBFace)

      ! Mean flow boundary condition
      Mixt%pv(0,iCell2) = 2d0 * Temperature - Mixt%pv(0,iCell1)
      Mixt%pv(1:nDim,iCell2) = 2d0 * Velocity(:) - Mixt%pv(1:nDim,iCell1)
      Mixt%pv(nDimP1,iCell2) = 2d0 * Pressure - Mixt%pv(nDimP1,iCell1)

      ! Turbulence boundary condition
      if( Conf%FlowModel == 3 ) then
        Omega = C1 * dsqrt( sum( Velocity(1:nDim)**2 ) ) / Conf%Length_Ref
        call SetSutherlandLaw(Conf,Temperature,LamMu)
        EddyMu = LamMu * 10d0**(-C2)
        Kine = EddyMu * Omega / Density
        Mixt%pv(nDimP2,iCell2) = 2d0 * Kine - Mixt%pv(nDimP2,iCell1)
        Mixt%pv(nDimP3,iCell2) = 2d0 * Omega - Mixt%pv(nDimP3,iCell1)
      end if

      ! Viscosity boundary condition
      select case(Conf%FlowModel)
      case(2) ! N-S equation
        Temperature = Mixt%pv(0,iCell2)
        call SetSutherlandLaw(Conf,Temperature,LamMu)
        Mixt%pv(iLamMu,iCell2) = LamMu
        Mixt%pv(iEdyMu,iCell2) = 0d0
      case(3) ! RANS equation
        Temperature = Mixt%pv(0,iCell2)
        call SetSutherlandLaw(Conf,Temperature,LamMu)
        Mixt%pv(iLamMu,iCell2) = LamMu
        Mixt%pv(iEdyMu,iCell2) = LamMu * 10d0**(-C2)
        Mixt%pv(iBlend,iCell2) = Mixt%pv(iBlend,iCell1)
      end select

    end do

  else ! Subsonic Inlet Implementaion ( Not Moving )

    do iList = List1, List2
      iBFace = Grid%List(iList)
      iCell1 = Grid%f2c(1,iBFace)
      iCell2 = Grid%f2c(2,iBFace)

      ! Mean flow boundary condition
      Pressure = Mixt%pv(nDimP1,iCell1)
      Mach = dsqrt( 2d0 / Mixt%Gamm1 * ( ( Total_Pressure / Pressure )**( Mixt%Gamm1 / Mixt%Gamma ) - 1d0 ) )
      Temperature = Total_Temperature / ( 1d0 + 0.5d0 * Mixt%Gamm1 * Mach**2 )
      Density = Pressure / ( Mixt%GasConstant * Temperature )
      SoundSpeed = dsqrt( Mixt%Gamma * Mixt%GasConstant * Temperature )
      Velocity(1) = Mach * SoundSpeed * dcos( pi * Conf%bv(4,iBC) / 180d0 )
      Velocity(2) = Mach * SoundSpeed * dsin( pi * Conf%bv(4,iBC) / 180d0 )
      if( nDim == 3 ) Velocity(3) = 0d0
      Mixt%pv(0,iCell2) = 2d0 * Temperature - Mixt%pv(0,iCell1)
      Mixt%pv(1:nDim,iCell2) = 2d0 * Velocity(:) - Mixt%pv(1:nDim,iCell1)
      Mixt%pv(nDimP1,iCell2) = Pressure

      ! Turbulence boundary condition
      if( Conf%FlowModel == 3 ) then
        Omega = C1 * dsqrt( sum( Velocity(1:nDim)**2 ) ) / Conf%Length_Ref
        call SetSutherlandLaw(Conf,Temperature,LamMu)
        EddyMu = LamMu * 10d0**(-C2)
        Kine = EddyMu * Omega / Density
        Mixt%pv(nDimP2,iCell2) = 2d0 * Kine - Mixt%pv(nDimP2,iCell1)
        Mixt%pv(nDimP3,iCell2) = 2d0 * Omega - Mixt%pv(nDimP3,iCell1)
      end if

      ! Viscosity boundary condition
      select case(Conf%FlowModel)
      case(2) ! N-S equation
        Temperature = Mixt%pv(0,iCell2)
        call SetSutherlandLaw(Conf,Temperature,LamMu)
        Mixt%pv(iLamMu,iCell2) = LamMu
        Mixt%pv(iEdyMu,iCell2) = 0d0
      case(3) ! RANS equation
        Temperature = Mixt%pv(0,iCell2)
        call SetSutherlandLaw(Conf,Temperature,LamMu)
        Mixt%pv(iLamMu,iCell2) = LamMu
        Mixt%pv(iEdyMu,iCell2) = LamMu * 10d0**(-C2)
        Mixt%pv(iBlend,iCell2) = Mixt%pv(iBlend,iCell1)
      end select

    end do
  end if

end subroutine

subroutine BC_Injection_Inlet(Mixt,Grid,Conf,iBC)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf
  integer, intent(in) :: iBC

  integer :: iList, List1, List2
  integer :: iBFace, iCell1, iCell2

  real(8) :: Temperature, Velocity(nDim), Pressure, Density
  real(8) :: ProjVelocity, Kine, Omega, LamMu, EddyMu

  List1 = Grid%ListIndex(iBC)
  List2 = Grid%ListIndex(iBC+1) - 1

  ! Injection Inlet Implementation ( Not Moving )
  do iList = List1, List2
    iBFace = Grid%List(iList)
    iCell1 = Grid%f2c(1,iBFace)
    iCell2 = Grid%f2c(2,iBFace)

    ! Mean flow boundary condition
    Temperature = Conf%bv(2,iBC)
    Pressure = Mixt%pv(nDimP1,iCell1)
    Density = Pressure / ( Mixt%GasConstant * Temperature )
    ProjVelocity  = Conf%bv(1,iBC) / Density
    Velocity(:) = - ProjVelocity * Grid%fn(:,iBFace)
    Mixt%pv(0,iCell2) = 2d0 * Temperature - Mixt%pv(0,iCell1)
    Mixt%pv(1:nDim,iCell2) = 2d0 * Velocity(:) - Mixt%pv(1:nDim,iCell1)
    Mixt%pv(nDimP1,iCell2) = Pressure

    ! Turbulence boundary condition
    if( Conf%FlowModel == 3 ) then
      Omega = C1 * dsqrt( sum( Velocity(1:nDim)**2 ) ) / Conf%Length_Ref
      call SetSutherlandLaw(Conf,Temperature,LamMu)
      EddyMu = LamMu * 10d0**(-C2)
      Kine = EddyMu * Omega / Density
      Mixt%pv(nDimP2,iCell2) = 2d0 * Kine - Mixt%pv(nDimP2,iCell1)
      Mixt%pv(nDimP3,iCell2) = 2d0 * Omega - Mixt%pv(nDimP3,iCell1)
    end if

    ! Viscosity boundary condition
    select case(Conf%FlowModel)
    case(2) ! N-S equation
      Temperature = Mixt%pv(0,iCell2)
      call SetSutherlandLaw(Conf,Temperature,LamMu)
      Mixt%pv(iLamMu,iCell2) = LamMu
      Mixt%pv(iEdyMu,iCell2) = 0d0
    case(3) ! RANS equation
      Temperature = Mixt%pv(0,iCell2)
      call SetSutherlandLaw(Conf,Temperature,LamMu)
      Mixt%pv(iLamMu,iCell2) = LamMu
      Mixt%pv(iEdyMu,iCell2) = LamMu * 10d0**(-C2)
      Mixt%pv(iBlend,iCell2) = Mixt%pv(iBlend,iCell1)
    end select

  end do

end subroutine

subroutine BC_Pressure_Outlet(Mixt,Grid,Conf,iBC)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf
  integer, intent(in) :: iBC

  integer :: iList, List1, List2
  integer :: iBFace, iCell1, iCell2

  real(8) :: Mach

  List1 = Grid%ListIndex(iBC)
  List2 = Grid%ListIndex(iBC+1) - 1

  ! Pressure Outlet Implementation ( Not Moving )
  do iList = List1, List2
    iBFace = Grid%List(iList)
    iCell1 = Grid%f2c(1,iBFace)
    iCell2 = Grid%f2c(2,iBFace)

    ! Mean flow boundary condition
    Mach = dsqrt( sum( Mixt%pv(1:nDim,iCell1)**2 ) / ( Mixt%Gamma * Mixt%GasConstant * Mixt%pv(0,iCell1) ) )
    Mixt%pv(0,iCell2) = Mixt%pv(0,iCell1)
    Mixt%pv(1:nDim,iCell2) = Mixt%pv(1:nDim,iCell1)
    if( Mach >= 1d0 ) then ! Supersonic case
      Mixt%pv(nDimP1,iCell2) = Mixt%pv(nDimP1,iCell1)
    else ! Subsonic case
      Mixt%pv(nDimP1,iCell2) = 2d0 * Conf%bv(1,iBC) - Mixt%pv(nDimP1,iCell1)
      if( Mixt%pv(nDimP1,iCell2) <= 0d0 ) Mixt%pv(nDimP1,iCell2) = 1d-15
    end if

    ! Turbulence boundary condition
    if( Conf%FlowModel == 3 ) then
      Mixt%pv(nDimP2,iCell2) = Mixt%pv(nDimP2,iCell1)
      Mixt%pv(nDimP3,iCell2) = Mixt%pv(nDimP3,iCell1)
    end if

    ! Viscosity boundary condition
    select case(Conf%FlowModel)
    case(2) ! N-S equation
      Mixt%pv(iLamMu,iCell2) = Mixt%pv(iLamMu,iCell1)
      Mixt%pv(iEdyMu,iCell2) = 0d0
    case(3) ! RANS equation
      Mixt%pv(iLamMu,iCell2) = Mixt%pv(iLamMu,iCell1)
      Mixt%pv(iEdyMu,iCell2) = Mixt%pv(iEdyMu,iCell1)
      Mixt%pv(iBlend,iCell2) = Mixt%pv(iBlend,iCell1)
    end select

  end do

end subroutine

subroutine BC_Pressure_Farfield(Mixt,Grid,Conf,iBC)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf
  integer, intent(in) :: iBC

  integer :: iList, List1, List2
  integer :: iBFace, iCell1, iCell2

  real(8) :: Temperature, Velocity(nDim), Pressure
  real(8) :: Density, SoundSpeed, ProjVelocity
  real(8) :: RefTemperature, RefVelocity(nDim), RefPressure
  real(8) :: RefDensity, RefSoundSpeed, RefProjVelocity
  real(8) :: bTemperature, bVelocity(nDim), bPressure, bDensity
  real(8) :: Kine, Omega, LamMu, EddyMu

  List1 = Grid%ListIndex(iBC)
  List2 = Grid%ListIndex(iBC+1) - 1

  ! Freestream Variables
  RefPressure = Conf%bv(1,iBC)
  RefTemperature = Conf%bv(2,iBC)
  RefDensity = RefPressure / ( Mixt%GasConstant * RefTemperature )
  RefSoundSpeed  = dsqrt( Mixt%Gamma * Mixt%GasConstant * RefTemperature )
  RefVelocity(1) = Conf%bv(3,iBC) * RefSoundSpeed * dcos( pi * Conf%bv(4,iBC) / 180d0 )
  RefVelocity(2) = Conf%bv(3,iBC) * RefSoundSpeed * dsin( pi * Conf%bv(4,iBC) / 180d0 )
  if( nDim == 3 ) RefVelocity(3) = 0d0

  do iList = List1, List2
    iBFace = Grid%List(iList)
    iCell1 = Grid%f2c(1,iBFace)
    iCell2 = Grid%f2c(2,iBFace)

    ! Inner Domain Variables
    Temperature = Mixt%pv(0,iCell1)
    Velocity(:) = Mixt%pv(1:nDim,iCell1)
    Pressure = Mixt%pv(nDimP1,iCell1)
    Density = Pressure / ( Mixt%GasConstant * Temperature )
    SoundSpeed  = dsqrt( Mixt%Gamma * Mixt%GasConstant * Temperature )

    ! Inner / Freestream Normal Velocity
    ProjVelocity = sum( Velocity(:) * Grid%fn(:,iBFace) )
    RefProjVelocity = sum( RefVelocity(:) * Grid%fn(:,iBFace) )

    ! Pressure Farfield Implementation ( Not Moving )
    if( ProjVelocity >= 0d0 ) then ! Outflow

      ! Mean flow boundary condition
      if( ProjVelocity - SoundSpeed >= 0d0 ) then ! Supersonic
        Mixt%pv(0,iCell2) = Temperature
        Mixt%pv(1:nDim,iCell2) = Velocity(:)
        Mixt%pv(nDimP1,iCell2) = Pressure
      else ! Subsonic
        bPressure = RefPressure
        bDensity = Density + ( bPressure - Pressure ) / RefSoundSpeed**2
        bVelocity(:) = Velocity(:) + ( Pressure - bPressure ) / ( RefDensity * RefSoundSpeed ) * Grid%fn(:,iBFace)
        bTemperature = bPressure / ( Mixt%GasConstant * bDensity )
        Mixt%pv(0,iCell2) = 2d0 * bTemperature - Temperature
        Mixt%pv(1:nDim,iCell2) = 2d0 * bVelocity(:) - Velocity(:)
        Mixt%pv(nDimP1,iCell2) = 2d0 * bPressure - Pressure
      end if

      ! Turbulence boundary condition
      if( Conf%FlowModel == 3 ) then
        Mixt%pv(nDimP2,iCell2) = Mixt%pv(nDimP2,iCell1)
        Mixt%pv(nDimP3,iCell2) = Mixt%pv(nDimP3,iCell1)
      end if

      ! Viscosity boundary condition
      select case(Conf%FlowModel)
      case(2) ! N-S equation
        Temperature = Mixt%pv(0,iCell2)
        call SetSutherlandLaw(Conf,Temperature,LamMu)
        Mixt%pv(iLamMu,iCell2) = LamMu
        Mixt%pv(iEdyMu,iCell2) = 0d0
      case(3) ! RANS equation
        Temperature = Mixt%pv(0,iCell2)
        call SetSutherlandLaw(Conf,Temperature,LamMu)
        Mixt%pv(iLamMu,iCell2) = LamMu
        Mixt%pv(iEdyMu,iCell2) = Mixt%pv(iEdyMu,iCell1)
        Mixt%pv(iBlend,iCell2) = Mixt%pv(iBlend,iCell1)
      end select

    else ! Inflow

      if( ProjVelocity + SoundSpeed <= 0d0 ) then ! Supersonic

        ! Mean flow boundary condition
        Mixt%pv(0,iCell2) = 2d0 * RefTemperature - Temperature
        Mixt%pv(1:nDim,iCell2) = 2d0 * RefVelocity(:) - Velocity(:)
        Mixt%pv(nDimP1,iCell2) = 2d0 * RefPressure - Pressure

        ! Turbulence boundary condition
        if( Conf%FlowModel == 3 ) then
          Omega = C1 * dsqrt( sum( RefVelocity(1:nDim)**2 ) ) / Conf%Length_Ref
          call SetSutherlandLaw(Conf,RefTemperature,LamMu)
          EddyMu = LamMu * 10d0**(-C2)
          Kine = EddyMu * Omega / RefDensity
          Mixt%pv(nDimP2,iCell2) = 2d0 * Kine - Mixt%pv(nDimP2,iCell1)
          Mixt%pv(nDimP3,iCell2) = 2d0 * Omega - Mixt%pv(nDimP3,iCell1)
        end if

      else ! Subsonic

        ! Mean flow boundary condition
        bPressure = 0.5d0 * ( RefPressure + Pressure - RefDensity * RefSoundSpeed * ( RefProjVelocity - ProjVelocity ) )
        bDensity = RefDensity + ( bPressure - RefPressure ) / RefSoundSpeed**2
        bVelocity(:) = RefVelocity(:) - ( RefPressure - bPressure ) / ( RefDensity * RefSoundSpeed ) * Grid%fn(:,iBFace)
        bTemperature = bPressure / ( Mixt%GasConstant * bDensity )
        Mixt%pv(0,iCell2) = 2d0 * bTemperature - Temperature
        Mixt%pv(1:nDim,iCell2) = 2d0 * bVelocity(:) - Velocity(:)
        Mixt%pv(nDimP1,iCell2) = 2d0 * bPressure - Pressure

        ! Turbulence boundary condition
        if( Conf%FlowModel == 3 ) then
          Omega = C1 * dsqrt( sum( bVelocity(1:nDim)**2 ) ) / Conf%Length_Ref
          call SetSutherlandLaw(Conf,bTemperature,LamMu)
          EddyMu = LamMu * 10d0**(-C2)
          Kine = EddyMu * Omega / RefDensity
          Mixt%pv(nDimP2,iCell2) = 2d0 * Kine - Mixt%pv(nDimP2,iCell1)
          Mixt%pv(nDimP3,iCell2) = 2d0 * Omega - Mixt%pv(nDimP3,iCell1)
        end if

      end if

      ! Viscosity boundary condition
      select case(Conf%FlowModel)
      case(2) ! N-S equation
        Temperature = Mixt%pv(0,iCell2)
        call SetSutherlandLaw(Conf,Temperature,LamMu)
        Mixt%pv(iLamMu,iCell2) = LamMu
        Mixt%pv(iEdyMu,iCell2) = 0d0
      case(3) ! RANS equation
        Temperature = Mixt%pv(0,iCell2)
        call SetSutherlandLaw(Conf,Temperature,LamMu)
        Mixt%pv(iLamMu,iCell2) = LamMu
        Mixt%pv(iEdyMu,iCell2) = LamMu * 10d0**(-C2)
        Mixt%pv(iBlend,iCell2) = Mixt%pv(iBlend,iCell1)
      end select

    end if
  end do

end subroutine

subroutine BC_Membrane_Nozzle(Mixt,Grid,Conf,iBC)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf
  integer, intent(in) :: iBC

  select case(Mixt%MembraneBroken)
  case(.false.) ! Inviscid / viscous Wall
    select case(Conf%FlowModel)
    case(1); call BC_Inviscid_Wall(Mixt,Grid,Conf,iBC)
    case(2); call BC_Viscous_Wall(Mixt,Grid,Conf,iBC)
    case(3); call BC_Viscous_Wall(Mixt,Grid,Conf,iBC)
    end select
  case(.true.) ! Pressure Outlet
    call BC_Pressure_Outlet(Mixt,Grid,Conf,iBC)
  end select

end subroutine

subroutine BC_Propellant(Mixt,Grid,Conf,iBC)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf
  integer, intent(in) :: iBC

  integer :: iList, List1, List2
  integer :: iBFace, iCell1, iCell2

  real(8) :: Temperature, Velocity(nDim), Pressure
  real(8) :: Density, RelVelocity, Coeff
  real(8) :: Kine, Omega, LamMu, EddyMu, Dist

  List1 = Grid%ListIndex(iBC)
  List2 = Grid%ListIndex(iBC+1) - 1

  do iList = List1, List2
    iBFace = Grid%List(iList)
    iCell1 = Grid%f2c(1,iBFace)
    iCell2 = Grid%f2c(2,iBFace)

    ! Propellant Boundary Condition Implementation
    select case(Mixt%bFlag(iBFace))
    case(0) ! Wall( Moving )

      select case(Conf%FlowModel)
      case(1) ! Inviscid Wall Implementation( Moving )

        ! Mean flow boundary condition
        RelVelocity= sum( Mixt%pv(1:nDim,iCell1) * Grid%fn(:,iBFace) ) - Grid%fv(iBFace)
        Mixt%pv(0,iCell2) = Mixt%pv(0,iCell1)
        Mixt%pv(1:nDim,iCell2) = Mixt%pv(1:nDim,iCell1) - 2d0 * RelVelocity * Grid%fn(:,iBFace)
        Mixt%pv(nDimP1,iCell2) = Mixt%pv(nDimP1,iCell1)

        ! Turbulence boundary condition
        if( Conf%FlowModel == 3 ) then
          Mixt%pv(nDimP2,iCell2) = Mixt%pv(nDimP2,iCell1)
          Mixt%pv(nDimP3,iCell2) = Mixt%pv(nDimP3,iCell1)
        end if

        ! Viscosity boundary condition
        select case(Conf%FlowModel)
        case(2) ! N-S equation
          Mixt%pv(iLamMu,iCell2) = Mixt%pv(iLamMu,iCell1)
          Mixt%pv(iEdyMu,iCell2) = 0d0
        case(3) ! RANS equation
          Mixt%pv(iLamMu,iCell2) = Mixt%pv(iLamMu,iCell1)
          Mixt%pv(iEdyMu,iCell2) = Mixt%pv(iEdyMu,iCell1)
          Mixt%pv(iBlend,iCell2) = Mixt%pv(iBlend,iCell1)
        end select

      case(2:3) ! Viscous Wall ( Not Moving )

        ! Mean flow boundary condition
        Mixt%pv(0,iCell2) = Mixt%pv(0,iCell1)
        Mixt%pv(1:nDim,iCell2) = -Mixt%pv(1:nDim,iCell1)
        Mixt%pv(nDimP1,iCell2) = Mixt%pv(nDimP1,iCell1)

        ! Turbulence boundary condition
        if( Conf%FlowModel == 3 ) then
          Dist = Grid%WallDist(iCell1)
          Temperature = Mixt%pv(0,iCell1) ! same as iCell2
          Pressure = Mixt%pv(nDimP1,iCell1) ! same as iCell2
          Density = Pressure / ( Mixt%GasConstant * Temperature )
          call SetSutherlandLaw(Conf,Temperature,LamMu)
          Omega = 60d0*LamMu/(Density*0.075d0*Dist**2)
          Mixt%pv(nDimP2,iCell2) = -Mixt%pv(nDimP2,iCell1)
          Mixt%pv(nDimP3,iCell2) = 2d0 * Omega - Mixt%pv(nDimP3,iCell1)
        end if

        ! Viscosity boundary condition
        select case(Conf%FlowModel)
        case(2) ! N-S equation
          Mixt%pv(iLamMu,iCell2) = Mixt%pv(iLamMu,iCell1)
          Mixt%pv(iEdyMu,iCell2) = 0d0
        case(3) ! RANS equation
          Mixt%pv(iLamMu,iCell2) = Mixt%pv(iLamMu,iCell1)
          Mixt%pv(iEdyMu,iCell2) = -Mixt%pv(iEdyMu,iCell1)
          Mixt%pv(iBlend,iCell2) = Mixt%pv(iBlend,iCell1)
        end select

      end select

    case(1,2) ! Injection Inlet ( Moving )

      ! Mean flow boundary condition
      Pressure = Mixt%pv(nDimP1,iCell1)
      Coeff = Pressure / ( Mixt%GasConstant * Mixt%fTemp(iBFace) )
      Density = 0.5d0 * ( Coeff + dsqrt( Coeff**2 + 2d0 * Mixt%Gamm1 * Mixt%mFlux(iBFace)**2 / ( Mixt%Gamma * Mixt%GasConstant * Mixt%fTemp(iBFace) ) ) )
      Temperature = Pressure / ( Mixt%GasConstant * Density )
      RelVelocity = Mixt%mFlux(iBFace) / Density - Grid%fv(iBFace)
      Velocity(:) = - RelVelocity * Grid%fn(:,iBFace)
      Mixt%pv(0,iCell2) = 2d0 * Temperature - Mixt%pv(0,iCell1)
      Mixt%pv(1:nDim,iCell2) = 2d0 * Velocity(:) - Mixt%pv(1:nDim,iCell1)
      Mixt%pv(nDimP1,iCell2) = Pressure

      ! Turbulence boundary condition
      if( Conf%FlowModel == 3 ) then
        Omega = C1 * dsqrt( sum( Velocity(1:nDim)**2 ) ) / Conf%Length_Ref
        call SetSutherlandLaw(Conf,Temperature,LamMu)
        EddyMu = LamMu * 10d0**(-C2)
        Kine = EddyMu * Omega / Density
        Mixt%pv(nDimP2,iCell2) = 2d0 * Kine - Mixt%pv(nDimP2,iCell1)
        Mixt%pv(nDimP3,iCell2) = 2d0 * Omega - Mixt%pv(nDimP3,iCell1)
      end if

      ! Viscosity boundary condition
      select case(Conf%FlowModel)
      case(2) ! N-S equation
        Temperature = Mixt%pv(0,iCell2)
        call SetSutherlandLaw(Conf,Temperature,LamMu)
        Mixt%pv(iLamMu,iCell2) = LamMu
        Mixt%pv(iEdyMu,iCell2) = 0d0
      case(3) ! RANS equation
        Temperature = Mixt%pv(0,iCell2)
        call SetSutherlandLaw(Conf,Temperature,LamMu)
        Mixt%pv(iLamMu,iCell2) = LamMu
        Mixt%pv(iEdyMu,iCell2) = LamMu * 10d0**(-C2)
        Mixt%pv(iBlend,iCell2) = Mixt%pv(iBlend,iCell1)
      end select

    end select
  end do

end subroutine

subroutine BC_Ignitor(Mixt,Grid,Conf,iBC)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf
  integer, intent(in) :: iBC

  integer :: iList, List1, List2
  integer :: iBFace, iCell1, iCell2

  integer :: ier, iIgni
  real(8) :: Total_Temperature
  real(8) :: Temperature, Velocity(nDim), Pressure
  real(8) :: Density, ProjVelocity
  real(8) :: a_sum, a_sum_total
  real(8) :: Pcmax, K1, K2, tb, tm
  real(8) :: Time1, Time2, Coeff, mFlux
  real(8) :: Kine, Omega, LamMu, EddyMu

  real(8), save :: inv_Area = 0d0
  logical, save :: ComputeArea = .true.
  integer, save :: Data_Flag = 1

  List1 = Grid%ListIndex(iBC)
  List2 = Grid%ListIndex(iBC+1) - 1

  ! Compute ignitor area
  if( ComputeArea ) then
    ComputeArea = .false.
    a_sum = 0d0
    do iList = List1, List2
      iBFace = Grid%List(iList)
      if( Conf%Axisymmetric ) then
        a_sum = a_sum + Grid%fa(iBFace) * 2d0 * pi * Grid%fc(2,iBFace)
      else ! not axisymmetric case
        a_sum = a_sum + Grid%fa(iBFace)
      end if
    end do
    call MPI_Allreduce(a_sum,a_sum_total,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
    if( a_sum_total > 0d0 ) inv_Area = 1d0 / a_sum_total
  end if

  ! Compute mass flux of ignitor
  mFlux = 0d0
  if( Conf%indIgnitor ) then ! ignitor model
    Pcmax = Conf%bv(2,iBC)
    K1 = Conf%bv(3,iBC)
    K2 = Conf%bv(4,iBC)
    tb = Conf%bv(5,iBC)
    tm = 1d0 / ( K1 * K2 ) * dlog( K1 * K2 * tb + 1 )
    if( Mixt%Time <= tb ) mFlux = K2 * Pcmax / ( 1d0 - tm / tb ) * ( 1d0 - Mixt%Time / tb )
  else ! ignitor file
    do iIgni = Data_Flag, Conf%nIgni-1
      Time1 = Conf%tIgni(iIgni); Time2 = Conf%tIgni(iIgni+1)
      if( Mixt%Time >= Time1 .and. Mixt%Time <= Time2 ) then
        mFlux = Conf%mIgni(iIgni) + ( Conf%mIgni(iIgni+1) - Conf%mIgni(iIgni) ) / ( Time2 - Time1 ) * ( Mixt%Time - Time1 )
        Data_Flag = iIgni; exit
      end if
    end do
  end if
  mFlux = mFlux * inv_Area

  ! Injection Inlet
  Total_Temperature = Conf%bv(1,iBC)
  do iList = List1, List2
    iBFace = Grid%List(iList)
    iCell1 = Grid%f2c(1,iBFace)
    iCell2 = Grid%f2c(2,iBFace)

    ! Mean flow boundary condition
    Coeff = 1d0 + ( Mixt%Gamma - 1 ) * 0.5d0 * 1d0**2
    Temperature = Total_Temperature / Coeff
    ProjVelocity = dsqrt( Mixt%Gamma * Mixt%GasConstant * Temperature ) * 1d0
    Density = mFlux / ProjVelocity
    Pressure = Density * Mixt%GasConstant * Temperature
    if( Pressure <= Mixt%pv(nDimP1,iCell1) ) then
      Pressure = Mixt%pv(nDimP1,iCell1)
      Coeff = Pressure / ( Mixt%GasConstant * Total_Temperature )
      Density = 0.5d0 * ( Coeff + dsqrt( Coeff**2 + 2d0 * Mixt%Gamm1 * mFlux**2 / ( Mixt%Gamma * Mixt%GasConstant * Total_Temperature ) ) )
      ProjVelocity = mFlux / Density
    end if
    Velocity(:) = - ProjVelocity * Grid%fn(:,iBFace)
    Mixt%pv(0,iCell2) = 2d0 * Temperature - Mixt%pv(0,iCell1)
    Mixt%pv(1:nDim,iCell2) = 2d0 * Velocity(:) - Mixt%pv(1:nDim,iCell1)
    Mixt%pv(nDimP1,iCell2) = 2d0 * Pressure - Mixt%pv(nDimP1,iCell1)

    ! Turbulence boundary condition
    if( Conf%FlowModel == 3 ) then
      Omega = C1 * dsqrt( sum( Velocity(1:nDim)**2 ) ) / Conf%Length_Ref
      call SetSutherlandLaw(Conf,Temperature,LamMu)
      EddyMu = LamMu * 10d0**(-C2)
      Kine = EddyMu * Omega / Density
      Mixt%pv(nDimP2,iCell2) = 2d0 * Kine - Mixt%pv(nDimP2,iCell1)
      Mixt%pv(nDimP3,iCell2) = 2d0 * Omega - Mixt%pv(nDimP3,iCell1)
    end if

    ! Viscosity boundary condition
    select case(Conf%FlowModel)
    case(2) ! N-S equation
      Temperature = Mixt%pv(0,iCell2)
      call SetSutherlandLaw(Conf,Temperature,LamMu)
      Mixt%pv(iLamMu,iCell2) = LamMu
      Mixt%pv(iEdyMu,iCell2) = 0d0
    case(3) ! RANS equation
      Temperature = Mixt%pv(0,iCell2)
      call SetSutherlandLaw(Conf,Temperature,LamMu)
      Mixt%pv(iLamMu,iCell2) = LamMu
      Mixt%pv(iEdyMu,iCell2) = LamMu * 10d0**(-C2)
      Mixt%pv(iBlend,iCell2) = Mixt%pv(iBlend,iCell1)
    end select

  end do

end subroutine

subroutine BC_Pressure_Sensor(Mixt,Grid,Conf,iBC)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf
  integer, intent(in) :: iBC

  ! Wall boundary condition
  select case(Conf%FlowModel)
  case(1); call BC_Inviscid_Wall(Mixt,Grid,Conf,iBC)
  case(2); call BC_Viscous_Wall(Mixt,Grid,Conf,iBC)
  case(3); call BC_Viscous_Wall(Mixt,Grid,Conf,iBC)
  end select

end subroutine

subroutine BC_APN_Model(Mixt,Grid,Conf,iBC)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf
  integer, intent(in) :: iBC

  integer :: iList, List1, List2
  integer :: iBFace, iCell1, iCell2

  real(8) :: Temperature, Velocity(nDim), Pressure
  real(8) :: Density, ProjVelocity
  real(8) :: Coeff, mFlux
  real(8) :: Kine, Omega, LamMu, EddyMu

  List1 = Grid%ListIndex(iBC)
  List2 = Grid%ListIndex(iBC+1) - 1

  ! Injection Inlet Implementation ( Not Moving )
  do iList = List1, List2
    iBFace = Grid%List(iList)
    iCell1 = Grid%f2c(1,iBFace)
    iCell2 = Grid%f2c(2,iBFace)

    ! Mean flow boundary condition
    Pressure = Mixt%pv(nDimP1,iCell1)
    Coeff = Pressure / ( Mixt%GasConstant * Conf%bv(2,iBC) )
    mFlux = Conf%bv(1,iBC) * Conf%bv(3,iBC) * ( Pressure / 6894733.26 )**Conf%bv(4,iBC)
    Density = 0.5d0 * ( Coeff + dsqrt( Coeff**2 + 2d0 * Mixt%Gamm1 * mFlux**2 / ( Mixt%Gamma * Mixt%GasConstant * Conf%bv(2,iBC) ) ) )
    Temperature = Pressure / ( Mixt%GasConstant * Density )
    ProjVelocity = mFlux / Density
    Velocity(:) = - ProjVelocity * Grid%fn(:,iBFace)
    Mixt%pv(0,iCell2) = 2d0 * Temperature - Mixt%pv(0,iCell1)
    Mixt%pv(1:nDim,iCell2) = 2d0 * Velocity(:) - Mixt%pv(1:nDim,iCell1)
    Mixt%pv(nDimP1,iCell2) = Pressure

    ! Turbulence boundary condition
    if( Conf%FlowModel == 3 ) then
      Omega = C1 * dsqrt( sum( Velocity(1:nDim)**2 ) ) / Conf%Length_Ref
      call SetSutherlandLaw(Conf,Temperature,LamMu)
      EddyMu = LamMu * 10d0**(-C2)
      Kine = EddyMu * Omega / Density
      Mixt%pv(nDimP2,iCell2) = 2d0 * Kine - Mixt%pv(nDimP2,iCell1)
      Mixt%pv(nDimP3,iCell2) = 2d0 * Omega - Mixt%pv(nDimP3,iCell1)
    end if

    ! Viscosity boundary condition
    select case(Conf%FlowModel)
    case(2) ! N-S equation
      Temperature = Mixt%pv(0,iCell2)
      call SetSutherlandLaw(Conf,Temperature,LamMu)
      Mixt%pv(iLamMu,iCell2) = LamMu
      Mixt%pv(iEdyMu,iCell2) = 0d0
    case(3) ! RANS equation
      Temperature = Mixt%pv(0,iCell2)
      call SetSutherlandLaw(Conf,Temperature,LamMu)
      Mixt%pv(iLamMu,iCell2) = LamMu
      Mixt%pv(iEdyMu,iCell2) = LamMu * 10d0**(-C2)
      Mixt%pv(iBlend,iCell2) = Mixt%pv(iBlend,iCell1)
    end select

  end do

end subroutine

subroutine SetBoundaryGradient(Mixt,Grid,Conf)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf

  integer :: iBC

  ! Implement boundary gradient / limiter
  do iBC = 1, Conf%nBC
    select case(Conf%BCType(iBC))
    case(1) ! inviscid wall
      call BG_Inviscid_Wall(Mixt,Grid,Conf,iBC)
    case(2) ! viscous wall
      call BG_Viscous_Wall(Mixt,Grid,Conf,iBC)
    case(9) ! membrane nozzle
      call BG_Membrane_Nozzle(Mixt,Grid,Conf,iBC)
    case(10) ! propellant
      call BG_Propellant(Mixt,Grid,Conf,iBC)
    end select
  end do

end subroutine

subroutine BG_Inviscid_Wall(Mixt,Grid,Conf,iBC)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf
  integer, intent(in) :: iBC

  integer :: iList, List1, List2
  integer :: iBFace, iCell1, iCell2, iVar, iDim
  real(8) :: ProjGrad, ProjVelocity, TangVelocity(nDim)
  real(8) :: NormVelGradient(nDim), TangVelGradient(nDim)
  real(8) :: Norm, Tang, Length, Inv_Length, ft(nDim), fb(nDim)
  real(8) :: a(nDim,nDim), Inv_a(nDim,nDim), b(nDim,nDim)
  real(8) :: Inv_Det_a

  List1 = Grid%ListIndex(iBC)
  List2 = Grid%ListIndex(iBC+1) - 1

  do iList = List1, List2
    iBFace = Grid%List(iList)
    iCell1 = Grid%f2c(1,iBFace)
    iCell2 = Grid%f2c(2,iBFace)

    ! Set symmetric gradient ( Temperature, Pressure )
    ProjGrad = sum( Mixt%gpv(:,0,iCell1) * Grid%fn(:,iBFace) )
    Mixt%gpv(:,0,iCell2) = Mixt%gpv(:,0,iCell1) - 2d0 * ProjGrad * Grid%fn(:,iBFace)
    ProjGrad = sum( Mixt%gpv(:,nDimP1,iCell1) * Grid%fn(:,iBFace) )
    Mixt%gpv(:,nDimP1,iCell2) = Mixt%gpv(:,nDimP1,iCell1) - 2d0 * ProjGrad * Grid%fn(:,iBFace)

    ! Compute face tangent vector
    ProjVelocity = sum( Mixt%pv(1:nDim,iCell1) * Grid%fn(:,iBFace) )
    TangVelocity(:) = Mixt%pv(1:nDim,iCell1) - ProjVelocity * Grid%fn(:,iBFace)
    Length = dsqrt( sum( TangVelocity(:)**2 ) )
    if( Length > 1d-12 ) then ! set tangential vector parallel to tagential velocity
      Inv_Length = 1d0 / Length
      ft(:) = TangVelocity(:) * Inv_Length
    else ! set arbitrary tangential vector
      ft(:) = Grid%xyz(:,Grid%f2n(1,iBFace)) - Grid%fc(:,iBFace)
      Inv_Length = 1d0 / dsqrt( sum( ft(:)**2 ) )
      ft(:) = ft(:) * Inv_Length
    end if

    ! Compute physical cell tangent / normal velocity gradient
    do iDim = 1, nDim
      TangVelGradient(iDim) = sum( Mixt%gpv(iDim,1:nDim,iCell1) * ft(:) )
      NormVelGradient(iDim) = sum( Mixt%gpv(iDim,1:nDim,iCell1) * Grid%fn(:,iBFace) )
    end do

    ! Compute ghost cell tangent / normal velocity gradient
    Norm = sum( TangVelGradient(:) * Grid%fn(:,iBFace) )
    Tang = sum( NormVelGradient(:) * ft(:) )
    TangVelGradient(:) = TangVelGradient(:) - 2d0 * Norm * Grid%fn(:,iBFace)
    NormVelGradient(:) = NormVelGradient(:) - 2d0 * Tang * ft(:)

    ! Transform coordinate ( x y z )
    a(:,1) = Grid%fn(:,iBFace)
    a(:,2) = ft(:)
    b(:,1) = NormVelGradient(:)
    b(:,2) = TangVelGradient(:)
    select case(nDim)
    case(2)
      Inv_Det_a = 1d0 / ( a(1,1) * a(2,2) - a(2,1) * a(1,2) )
      Inv_a(1,1) = a(2,2) * Inv_Det_a
      Inv_a(2,1) = -a(2,1) * Inv_Det_a
      Inv_a(1,2) = -a(1,2) * Inv_Det_a
      Inv_a(2,2) = a(1,1) * Inv_Det_a
      do iVar = 1, nDim
        Mixt%gpv(1,iVar,iCell2) = sum( Inv_a(:,iVar) * b(1,:) )
        Mixt%gpv(2,iVar,iCell2) = sum( Inv_a(:,iVar) * b(2,:) )
      end do
    case(3) ! compute binormal vector
      fb(1) = ft(2) * Grid%fn(3,iBFace) - ft(3) * Grid%fn(2,iBFace)
      fb(2) = ft(3) * Grid%fn(1,iBFace) - ft(1) * Grid%fn(3,iBFace)
      fb(3) = ft(1) * Grid%fn(2,iBFace) - ft(2) * Grid%fn(1,iBFace)
      Inv_Length = 1d0 / dsqrt( sum( fb(:)**2 ) )
      a(:,3) = fb(:) * Inv_Length
      b(:,3) = 0d0
      Inv_Det_a = 1d0 / ( a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) - a(2,1)*(a(1,2)*a(3,3) - a(3,2)*a(1,3)) + a(3,1)*(a(1,2)*a(2,3) - a(2,2)*a(1,3)) )
      Inv_a(1,1) = +( a(2,2)*a(3,3) - a(2,3)*a(3,2) ) * Inv_Det_a
      Inv_a(2,1) = -( a(2,1)*a(3,3) - a(2,3)*a(3,1) ) * Inv_Det_a
      Inv_a(3,1) = +( a(2,1)*a(3,2) - a(2,2)*a(3,1) ) * Inv_Det_a
      Inv_a(1,2) = -( a(1,2)*a(3,3) - a(1,3)*a(3,2) ) * Inv_Det_a
      Inv_a(2,2) = +( a(1,1)*a(3,3) - a(1,3)*a(3,1) ) * Inv_Det_a
      Inv_a(3,2) = -( a(1,1)*a(3,2) - a(1,2)*a(3,1) ) * Inv_Det_a
      Inv_a(1,3) = +( a(1,2)*a(2,3) - a(1,3)*a(2,2) ) * Inv_Det_a
      Inv_a(2,3) = -( a(1,1)*a(2,3) - a(1,3)*a(2,1) ) * Inv_Det_a
      Inv_a(3,3) = +( a(1,1)*a(2,2) - a(1,2)*a(2,1) ) * Inv_Det_a
      do iVar = 1, nDim
        Mixt%gpv(1,iVar,iCell2) = sum( Inv_a(:,iVar) * b(1,:) )
        Mixt%gpv(2,iVar,iCell2) = sum( Inv_a(:,iVar) * b(2,:) )
        Mixt%gpv(3,iVar,iCell2) = sum( Inv_a(:,iVar) * b(3,:) )
      end do
    end select

    ! Turbulence boundary gradient
    if( Conf%FlowModel == 3 ) then
      ProjGrad = sum( Mixt%gpv(:,nDimP2,iCell1) * Grid%fn(:,iBFace) )
      Mixt%gpv(:,nDimP2,iCell2) = Mixt%gpv(:,nDimP2,iCell1) - 2d0 * ProjGrad * Grid%fn(:,iBFace)
      ProjGrad = sum( Mixt%gpv(:,nDimP3,iCell1) * Grid%fn(:,iBFace) )
      Mixt%gpv(:,nDimP3,iCell2) = Mixt%gpv(:,nDimP3,iCell1) - 2d0 * ProjGrad * Grid%fn(:,iBFace)
    end if

    ! Copy limiter
    Mixt%limiter(:,iCell2) = Mixt%limiter(:,iCell1)

  end do

end subroutine

subroutine BG_Viscous_Wall(Mixt,Grid,Conf,iBC)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf
  integer, intent(in) :: iBC

  integer :: iList, List1, List2
  integer :: iBFace, iCell1, iCell2, iVar
  real(8) :: ProjGrad

  List1 = Grid%ListIndex(iBC)
  List2 = Grid%ListIndex(iBC+1) - 1

  do iList = List1, List2
    iBFace = Grid%List(iList)
    iCell1 = Grid%f2c(1,iBFace)
    iCell2 = Grid%f2c(2,iBFace)

    ! Set symmetric gradient ( Temperature, Pressure )
    ProjGrad = sum( Mixt%gpv(:,0,iCell1) * Grid%fn(:,iBFace) )
    Mixt%gpv(:,0,iCell2) = Mixt%gpv(:,0,iCell1) - 2d0 * ProjGrad * Grid%fn(:,iBFace)
    ProjGrad = sum( Mixt%gpv(:,nDimP1,iCell1) * Grid%fn(:,iBFace) )
    Mixt%gpv(:,nDimP1,iCell2) = Mixt%gpv(:,nDimP1,iCell1) - 2d0 * ProjGrad * Grid%fn(:,iBFace)
    ! Set reverse tangential gradient ( Velocity )
    do iVar = 1, nDim
      ProjGrad = sum( Mixt%gpv(:,iVar,iCell1) * Grid%fn(:,iBFace) )
      Mixt%gpv(:,iVar,iCell2) = 2d0 * ProjGrad * Grid%fn(:,iBFace) - Mixt%gpv(:,iVar,iCell1)
    end do

    ! Turbulence boundary gradient
    if( Conf%FlowModel == 3 ) then
      ProjGrad = sum( Mixt%gpv(:,nDimP2,iCell1) * Grid%fn(:,iBFace) )
      Mixt%gpv(:,nDimP2,iCell2) = 2d0 * ProjGrad * Grid%fn(:,iBFace) - Mixt%gpv(:,nDimP2,iCell1)
      ProjGrad = sum( Mixt%gpv(:,nDimP3,iCell1) * Grid%fn(:,iBFace) )
      Mixt%gpv(:,nDimP3,iCell2) = Mixt%gpv(:,nDimP3,iCell1) - 2d0 * ProjGrad * Grid%fn(:,iBFace)
    end if

    ! Copy limiter
    Mixt%limiter(:,iCell2) = Mixt%limiter(:,iCell1)

  end do

end subroutine

subroutine BG_Membrane_Nozzle(Mixt,Grid,Conf,iBC)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf
  integer, intent(in) :: iBC

  if( .not. Mixt%MembraneBroken ) then
    select case(Conf%FlowModel)
    case(1); call BG_Inviscid_Wall(Mixt,Grid,Conf,iBC)
    case(2); call BG_Viscous_Wall(Mixt,Grid,Conf,iBC)
    case(3); call BG_Viscous_Wall(Mixt,Grid,Conf,iBC)
    end select
  end if

end subroutine

subroutine BG_Propellant(Mixt,Grid,Conf,iBC)
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Grid), intent(in) :: Grid
  type(t_Conf), intent(in) :: Conf
  integer, intent(in) :: iBC

  integer :: iList, List1, List2
  integer :: iBFace, iCell1, iCell2, iVar, iDim
  real(8) :: ProjGrad, ProjVelocity, TangVelocity(nDim)
  real(8) :: NormVelGradient(nDim), TangVelGradient(nDim)
  real(8) :: Norm, Tang, Length, Inv_Length, ft(nDim), fb(nDim)
  real(8) :: a(nDim,nDim), Inv_a(nDim,nDim), b(nDim,nDim)
  real(8) :: Inv_Det_a

  List1 = Grid%ListIndex(iBC)
  List2 = Grid%ListIndex(iBC+1) - 1

  do iList = List1, List2
    iBFace = Grid%List(iList)
    iCell1 = Grid%f2c(1,iBFace)
    iCell2 = Grid%f2c(2,iBFace)

    ! Propellant Boundary Condition Implementation
    if( Mixt%bFlag(iBFace) == 0 ) then
      select case(Conf%FlowModel)
      case(1) ! inviscid wall

        ! Set symmetric gradient ( Temperature, Pressure )
        ProjGrad = sum( Mixt%gpv(:,0,iCell1) * Grid%fn(:,iBFace) )
        Mixt%gpv(:,0,iCell2) = Mixt%gpv(:,0,iCell1) - 2d0 * ProjGrad * Grid%fn(:,iBFace)
        ProjGrad = sum( Mixt%gpv(:,nDimP1,iCell1) * Grid%fn(:,iBFace) )
        Mixt%gpv(:,nDimP1,iCell2) = Mixt%gpv(:,nDimP1,iCell1) - 2d0 * ProjGrad * Grid%fn(:,iBFace)

        ! Compute face tangent vector
        ProjVelocity = sum( Mixt%pv(1:nDim,iCell1) * Grid%fn(:,iBFace) )
        TangVelocity(:) = Mixt%pv(1:nDim,iCell1) - ProjVelocity * Grid%fn(:,iBFace)
        Length = dsqrt( sum( TangVelocity(:)**2 ) )
        if( Length > 1d-12 ) then
          Inv_Length =  1d0 / Length
          ft(:) = TangVelocity(:) * Inv_Length
        else
          ft(1) = -Grid%fn(2,iBFace)
          ft(2) = Grid%fn(1,iBFace)
          if( nDim == 3 ) ft(3) = 0d0
        end if

        ! Compute physical cell tangent / normal velocity gradient
        do iDim = 1, nDim
          TangVelGradient(iDim) = sum( Mixt%gpv(iDim,1:nDim,iCell1) * ft(:) )
          NormVelGradient(iDim) = sum( Mixt%gpv(iDim,1:nDim,iCell1) * Grid%fn(:,iBFace) )
        end do

        ! Compute ghost cell tangent / normal velocity gradient
        Norm = sum( TangVelGradient(:) * Grid%fn(:,iBFace) )
        Tang = sum( NormVelGradient(:) * ft(:) )
        TangVelGradient(:) = TangVelGradient(:) - 2d0 * Norm * Grid%fn(:,iBFace)
        NormVelGradient(:) = NormVelGradient(:) - 2d0 * Tang * ft(:)

        ! Transform coordinate ( x y z )
        a(:,1) = Grid%fn(:,iBFace)
        a(:,2) = ft(:)
        b(:,1) = NormVelGradient(:)
        b(:,2) = TangVelGradient(:)
        select case(nDim)
        case(2)
          Inv_Det_a = 1d0 / ( a(1,1) * a(2,2) - a(2,1) * a(1,2) )
          Inv_a(1,1) = a(2,2) * Inv_Det_a
          Inv_a(2,1) = -a(2,1) * Inv_Det_a
          Inv_a(1,2) = -a(1,2) * Inv_Det_a
          Inv_a(2,2) = a(1,1) * Inv_Det_a
          do iVar = 1, nDim
            Mixt%gpv(1,iVar,iCell2) = sum( Inv_a(:,iVar) * b(1,:) )
            Mixt%gpv(2,iVar,iCell2) = sum( Inv_a(:,iVar) * b(2,:) )
          end do
        case(3) ! compute binormal vector
          fb(1) = ft(2) * Grid%fn(3,iBFace) - ft(3) * Grid%fn(2,iBFace)
          fb(2) = ft(3) * Grid%fn(1,iBFace) - ft(1) * Grid%fn(3,iBFace)
          fb(3) = ft(1) * Grid%fn(2,iBFace) - ft(2) * Grid%fn(1,iBFace)
          Inv_Length = 1d0 / dsqrt( sum( fb(:)**2 ) )
          a(:,3) = fb(:) * Inv_Length
          b(:,3) = 0d0
          Inv_Det_a = 1d0 / ( a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) - a(2,1)*(a(1,2)*a(3,3) - a(3,2)*a(1,3)) + a(3,1)*(a(2,1)*a(2,3) - a(2,2)*a(1,3)) )
          Inv_a(1,1) = +( a(2,2)*a(3,3) - a(2,3)*a(3,2) ) * Inv_Det_a
          Inv_a(2,1) = -( a(2,1)*a(3,3) - a(2,3)*a(3,1) ) * Inv_Det_a
          Inv_a(3,1) = +( a(2,1)*a(3,2) - a(2,2)*a(3,1) ) * Inv_Det_a
          Inv_a(1,2) = -( a(1,2)*a(3,3) - a(1,3)*a(3,2) ) * Inv_Det_a
          Inv_a(2,2) = +( a(1,1)*a(3,3) - a(1,3)*a(3,1) ) * Inv_Det_a
          Inv_a(3,2) = -( a(1,1)*a(3,2) - a(1,2)*a(3,1) ) * Inv_Det_a
          Inv_a(1,3) = +( a(1,2)*a(2,3) - a(1,3)*a(2,2) ) * Inv_Det_a
          Inv_a(2,3) = -( a(1,1)*a(2,3) - a(1,3)*a(2,1) ) * Inv_Det_a
          Inv_a(3,3) = +( a(1,1)*a(2,2) - a(1,2)*a(2,1) ) * Inv_Det_a
          do iVar = 1, nDim
            Mixt%gpv(1,iVar,iCell2) = sum( Inv_a(:,iVar) * b(1,:) )
            Mixt%gpv(2,iVar,iCell2) = sum( Inv_a(:,iVar) * b(2,:) )
            Mixt%gpv(3,iVar,iCell2) = sum( Inv_a(:,iVar) * b(3,:) )
          end do
        end select

        ! Turbulence boundary gradient
        if( Conf%FlowModel == 3 ) then
          ProjGrad = sum( Mixt%gpv(:,nDimP2,iCell1) * Grid%fn(:,iBFace) )
          Mixt%gpv(:,nDimP2,iCell2) = Mixt%gpv(:,nDimP2,iCell1) - 2d0 * ProjGrad * Grid%fn(:,iBFace)
          ProjGrad = sum( Mixt%gpv(:,nDimP3,iCell1) * Grid%fn(:,iBFace) )
          Mixt%gpv(:,nDimP3,iCell2) = Mixt%gpv(:,nDimP3,iCell1) - 2d0 * ProjGrad * Grid%fn(:,iBFace)
        end if

      case(2:3) ! viscous wall

        ! Set symmetric gradient ( Temperature, Pressure )
        ProjGrad = sum( Mixt%gpv(:,0,iCell1) * Grid%fn(:,iBFace) )
        Mixt%gpv(:,0,iCell2) = Mixt%gpv(:,0,iCell1) - 2d0 * ProjGrad * Grid%fn(:,iBFace)
        ProjGrad = sum( Mixt%gpv(:,nDimP1,iCell1) * Grid%fn(:,iBFace) )
        Mixt%gpv(:,nDimP1,iCell2) = Mixt%gpv(:,nDimP1,iCell1) - 2d0 * ProjGrad * Grid%fn(:,iBFace)
        ! Set reverse tangential gradient ( Velocity )
        do iVar = 1, nDim
          ProjGrad = sum( Mixt%gpv(:,iVar,iCell1) * Grid%fn(:,iBFace) )
          Mixt%gpv(:,iVar,iCell2) = 2d0 * ProjGrad * Grid%fn(:,iBFace) - Mixt%gpv(:,iVar,iCell1)
        end do

        ! Turbulence boundary gradient
        if( Conf%FlowModel == 3 ) then
          ProjGrad = sum( Mixt%gpv(:,nDimP2,iCell1) * Grid%fn(:,iBFace) )
          Mixt%gpv(:,nDimP2,iCell2) = 2d0 * ProjGrad * Grid%fn(:,iBFace) - Mixt%gpv(:,nDimP2,iCell1)
          ProjGrad = sum( Mixt%gpv(:,nDimP3,iCell1) * Grid%fn(:,iBFace) )
          Mixt%gpv(:,nDimP3,iCell2) = Mixt%gpv(:,nDimP3,iCell1) - 2d0 * ProjGrad * Grid%fn(:,iBFace)
        end if

      end select

      ! Copy limiter
      Mixt%limiter(:,iCell2) = Mixt%limiter(:,iCell1)
    end if
  end do

end subroutine

subroutine SumConservative(Mixt,Grid,cvSum)
  type(t_Mixt), intent(in) :: Mixt
  type(t_Grid), intent(in) :: Grid
  real(8), intent(out) :: cvSum(:)

  integer :: ier, iCell
  real(8) :: cvSumLocal(0:Mixt%nCVar-1)

  ! Compute sum of conservative variable
  ! note: conserved variable cv is already multiplied by volume
  cvSumLocal(:) = 0d0
  do iCell = 1, Grid%nCell
    cvSumLocal(:) = cvSumLocal(:) + Mixt%cv(:,iCell)
  end do

  ! MPI all reduce for summation
  call MPI_Allreduce(cvSumLocal,cvSum,Mixt%nCVar,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)

end subroutine

end module
