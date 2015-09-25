module process_module

use mpi
use config_module
use grid_module
use mixture_module
use numerics_module
use time_module
use post_module
use interface_module
implicit none; private
public InitialProcess
public SolverProcess
public RemeshProcess
public PostProcess
public FinishProcess

contains

subroutine InitialProcess(Grid,Mixt,Conf)
  type(t_Grid), intent(out) :: Grid
  type(t_Mixt), intent(out) :: Mixt
  type(t_Conf), intent(out) :: Conf

  integer :: rank, ier
  character(17) :: ConfFileName = './input/fluid.inp'
  character(20) :: BdryFileName = './input/fluid.inp.bc'
  type(t_Grid) :: Global_Grid

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Set configuration info
  call SetConf(Conf,ConfFileName)
  call SetBdryConf(Conf,BdryFileName)

  ! Set global grid from file and set partition info
  if( rank == 0 ) then
    call SetGridFromFile(Global_Grid,Conf)
    call SetPartGrid(Global_Grid)
  end if

  ! Set local grid and delete global grid
  call SetLocalGrid(Grid,Global_Grid)
  call DelGrid(Global_Grid)

  ! Set node / face / cell data structure
  call SetNodeDataStruc(Grid)
  call SetFaceDataStruc(Grid)
  call SetCellDataStruc(Grid)
  call SetListDataStruc(Grid,Conf)

  ! Re-number cell data structure
  call SetGridReNumber(Grid)

  ! Set global to local mapping
  call SetGlobalToLocal(Grid)

  ! Set connected domain send / recv list
  ! and boundary domain recv list
  call SetConnSendRecv(Grid)
  call SetBoundRecv(Grid)

  ! Set face / cell metric
  call SetFaceMetric(Grid)
  call SetCellMetric(Grid)

  ! Set wall distance for kw-SST
  if( Conf%FlowModel == 3 ) then
    call SetWallDistance(Grid,Conf)
  end if

  ! Preprocess for grid velocity
  if( Conf%MoveGrid ) then
    call PreGridVelocity(Grid)
  end if

  ! Set mixture from file or initial condition
  if( Conf%Restart ) then ! set mixture from restart file
    call SetMixtFromFile(Mixt,Grid,Conf)
  else ! set mixture from initial condition
    call SetMixtFromInitial(Mixt,Grid,Conf)
  end if

  ! Set convection flux
  call SetConv(Grid,Conf)

  ! Set post map
  call SetPostMap(Grid)

  ! Post pressure at monitoring boundary
  call Post_Pressure(Grid,Mixt,Conf)

end subroutine

subroutine SolverProcess(Grid,Mixt,Conf,DeltaTime,Bxyz,BCFlag,InVars)
  type(t_Grid), intent(inout) :: Grid
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Conf), intent(inout) :: Conf
  real(8), intent(in) :: DeltaTime
  real(8), intent(in) :: Bxyz(:,:)
  integer, intent(in) :: BCFlag(:)
  real(8), intent(in) :: InVars(:,:)

  integer :: rank, ier

  integer :: Iter
  real(8) :: TargetTime
  real(8) :: Time_Start, Time_End
  real(8) :: L2Norm_Ref, L2Norm_Abs, L2Norm_Rel
  real(8), save :: TimeStepRef
  logical, save :: TimeStepRef_Flag = .false.
  real(8), parameter :: Factor = 0.01d0

  ! Temporary boundary grid data
  integer :: nBNode_temp
  integer :: nBFace_temp
  real(8), allocatable :: Bxyz_temp(:,:)
  integer, allocatable :: Bf2n_temp(:,:)
  integer, allocatable :: BPatch_temp(:)
  integer, allocatable :: BCFlag_temp(:)
  real(8), allocatable :: InVars_temp(:,:)
  real(8), allocatable :: OutVars_temp(:,:)

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Set target time
  TargetTime = Mixt%Time + DeltaTime

  ! Set reference time step
  if( Conf%MoveGrid .and. .not. TimeStepRef_Flag ) then
    call ComputeTimeStep(Mixt,Grid,Conf)
    TimeStepRef = Mixt%TimeStepMin * Factor
    TimeStepRef_Flag = .true.
  end if

  ! Compute grid velocity
  if( Conf%MoveGrid ) then
    call SetGridVelocity(Grid)
  end if

  ! Write brief info
  if( rank == 0 ) then
    if( Conf%Unsteady ) then
      Time_Start = MPI_Wtime()
      write(*,*)
      write(*,*) 'FLUS >> Start Fluid Region Step'
      write(*,*) 'FLUS >> ------------------------------------------------------'
      write(*,*) 'FLUS >>    Iteration         Time                Dt           '
      write(*,*) 'FLUS >> ------------------------------------------------------'
    else
      Time_Start = MPI_Wtime()
      write(*,*)
      write(*,*) 'FLUS >> Start Fluid Region Step'
      write(*,*) 'FLUS >> ------------------------------------------------------'
      write(*,*) 'FLUS >>    Iteration         Error(abs)          Error(rel)   '
      write(*,*) 'FLUS >> ------------------------------------------------------'
    end if
  end if

  ! Main Loop
  if( .not. Conf%Unsteady ) then
    Mixt%Time = TargetTime
    Mixt%L2Norm = 0d0
  end if

  Iter = 0
  do; Iter = Iter + 1

    ! Compute time step
    call ComputeTimeStep(Mixt,Grid,Conf,TargetTime)

    ! Time integration method
    select case(Conf%iTime)
    case(1)
      call GridMovingTest(Grid,Mixt)
    case(2)
      call ExplicitEuler(Grid,Mixt,Conf)
    case(3)
      call TVDRK3(Grid,Mixt,Conf)
    case(4)
      call LUSGS(Grid,Mixt,Conf)
    end select

    ! Compute error
    if( .not. Conf%Unsteady ) then
      call ComputeError(Mixt,Grid)
      if( Iter == 1 ) L2Norm_Ref = Mixt%L2Norm
      L2Norm_Abs = Mixt%L2Norm
      L2Norm_Rel = L2Norm_Abs / L2Norm_Ref
    end if

    ! Print subiteration
    if( rank == 0 .and. mod(Iter,Conf%nPrint) == 0 ) then
      if( Conf%Unsteady .and. .not. Conf%DualTimeStep ) then
        write(*,'(x,a8,i11,ES20.4,ES20.4)') 'FLUS >> ', Iter, Mixt%Time, Mixt%TimeStepMin
      else
        write(*,'(x,a8,i11,ES20.4,ES20.4)') 'FLUS >> ', Iter, log10(L2Norm_Abs), log10(L2Norm_Rel)
      end if
    end if

    if( Conf%Unsteady ) then
      if( dabs( Mixt%Time - TargetTime ) < 1d-12 ) exit
    else
      if( L2Norm_Abs < 1d-12 ) exit
      if( L2Norm_Rel < 10d0**Conf%Error ) exit
    end if

    ! Check remeshing during solver processing
    if( Conf%MoveGrid .and. Mixt%TimeStepMin < TimeStepRef ) then
      call Post_FLUS_Format(Grid,Mixt,Conf,'./output/before.flu')
      ! Create boundary grid data
      call InterfaceSet(Grid,Conf,nBNode_temp,nBFace_temp,Bxyz_temp,Bf2n_temp,BPatch_temp,BCFlag_temp,InVars_temp,OutVars_temp)
      ! Remeshing from boundary grid
      call RemeshProcess(Grid,Mixt,Conf,nBNode_temp,nBFace_temp,Bxyz_temp,Bf2n_temp,BPatch_temp)
      ! Delete boundary grid data
      call InterfaceDel(Bxyz_temp,Bf2n_temp,BPatch_temp,BCFlag_temp,InVars_temp,OutVars_temp)
      ! Set interface data
      call InterfaceIn(Grid,Mixt,TargetTime-Mixt%Time,Bxyz,BCFlag,InVars)
      ! Set Ref Time Step
      call ComputeTimeStep(Mixt,Grid,Conf)
      TimeStepRef = Mixt%TimeStepMin * Factor
      ! Set Grid Velocity
      call SetGridVelocity(Grid)
      call Post_FLUS_Format(Grid,Mixt,Conf,'./output/after.flu')
    end if

  end do

  if( rank == 0 ) then
    Time_End = MPI_Wtime()
    write(*,*) 'FLUS >> ------------------------------------------------------'
    write(*,*) 'FLUS >>    Computational Time = ', Time_End - Time_Start
    write(*,*) 'FLUS >> ------------------------------------------------------'
    write(*,*) 'FLUS >> End Fluid Region Step'
    write(*,*)
  end if

  ! Post pressure at monitor boundary
  call Post_Pressure(Grid,Mixt,Conf)
end subroutine

subroutine RemeshProcess(Grid,Mixt,Conf,nBNode,nBFace,Bxyz,Bf2n,BPatch)
  type(t_Grid), intent(inout) :: Grid
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Conf), intent(inout) :: Conf
  integer, intent(in) :: nBNode
  integer, intent(in) :: nBFace
  real(8), intent(in) :: Bxyz(:,:)
  integer, intent(in) :: Bf2n(:,:)
  integer, intent(in) :: BPatch(:)

  integer :: rank, io, ier
  real(8) :: cvSumOld(0:Mixt%nCVar-1)
  real(8) :: cvSumNew(0:Mixt%nCVar-1)
  real(8) :: cvSumErr(0:Mixt%nCVar-1)
  logical, save :: Flag_New_File = .true.

  type(t_Grid) :: Global_Grid
  type(t_Grid) :: GridOld
  type(t_Mixt) :: MixtOld

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Compute sum of the old conservative variable
  call SumConservative(Mixt,Grid,cvSumOld)

  ! Create backup and delete grid / mixture for remesh
  ! (note: check whether grid/mixt type have pointer variable)
  GridOld = Grid; call DelGrid(Grid)
  MixtOld = Mixt; call DelMixt(Mixt)

  ! Reset flow model (Euler equation)
  ! note: After remeshing, we use Euler equation
  Conf%FlowModel = 1

  ! Set global grid from surface and set partition info
  if( rank == 0 ) then
    call SetGridFromSurface(Global_Grid,Grid%nDim,nBNode,nBFace,Bxyz,Bf2n,BPatch)
    call SetPartGrid(Global_Grid)
  end if

  ! Set loal grid and delete global grid
  call SetLocalGrid(Grid,Global_Grid)
  call DelGrid(Global_Grid)

  ! Set node / face / cell data structure
  call SetNodeDataStruc(Grid)
  call SetFaceDataStruc(Grid)
  call SetCellDataStruc(Grid)
  call SetListDataStruc(Grid,Conf)

  ! Set global to local mapping
  call SetGlobalToLocal(Grid)

  ! Set connected domain send / recv list
  ! and boundary domain recv list
  call SetConnSendRecv(Grid)
  call SetBoundRecv(Grid)

  ! Set face / cell metric
  call SetFaceMetric(Grid)
  call SetCellMetric(Grid)

  ! Preprocess for grid velocity
  if( Conf%MoveGrid ) then
    call PreGridVelocity(Grid)
  end if

  ! Set mixture from backup using data transfer
  call SetMixtFromBackup(Mixt,Grid,Conf,MixtOld,GridOld)

  ! Delete backup grid / mixture
  call DelGrid(GridOld)
  call DelMixt(MixtOld)

  ! Delete old post map and reset post map
  call DelPostMap()
  call SetPostMap(Grid)

  ! Compute sum of the new conservative variable
  call SumConservative(Mixt,Grid,cvSumNew)

  ! Compute error during remesh process
  cvSumErr(:) = ( cvSumNew(:) - cvSumOld(:) ) / ( cvSumOld(:) + 1d-15 )
  if( rank == 0 ) then
    ! Print error during remesh process
    write(*,*) "FLUS >> Conservation error during remesh process"
    write(*,*) "FLUS >> Mass       =", cvSumErr(0)
    write(*,*) "FLUS >> X-Momentum =", cvSumErr(1)
    write(*,*) "FLUS >> Y-Momentum =", cvSumErr(2)
    select case(Grid%nDim)
    case(2)
      wriTE(*,*) "FLUS >> Energy     =", cvSumErr(3)
    case(3)
      write(*,*) "FLUS >> Z-Momentum =", cvSumErr(3)
      write(*,*) "FLUS >> Energy     =", cvSumErr(4)
    end select
    ! Write error log
    if( Flag_New_File ) then
      Flag_New_File = .false.
      open(newunit=io,file='./output/remesh_error.csv',status='REPLACE',action='WRITE')
      select case(Grid%nDim)
      case(2); write(io,*) "Mass, X-Momentum, Y-Momentum, Energy"
      case(3); write(io,*) "Mass, X-Momentum, Y-Momentum, Z-Momentum, Energy"
      end select
    else
      open(newunit=io,file='./output/remesh_error.csv',status='OLD',action='WRITE',access='APPEND')
    end if
    select case(Grid%nDim)
    case(2); write(io,*) cvSumErr(0),',',cvSumErr(1),',',cvSumErr(2),',',cvSumErr(3)
    case(3); write(io,*) cvSumErr(0),',',cvSumErr(1),',',cvSumErr(2),',',cvSumErr(3),',',cvSumErr(4)
    end select
    close(io)
  end if

end subroutine

subroutine PostProcess(Grid,Mixt,Conf,Iter)
  type(t_Grid), intent(in) :: Grid
  type(t_Mixt), intent(in) :: Mixt
  type(t_Conf), intent(in) :: Conf
  integer, intent(in) :: Iter

  integer :: rank, ier
  character(32) :: char_Iter, FileName

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Set post file name
  write(char_Iter,"(i6.6)") Iter
  FileName = "./output/fluid_"//trim(char_Iter)//".flu"

  ! Post FLUS file
  call Post_FLUS_Format(Grid,Mixt,Conf,trim(FileName))
  if( rank == 0 ) write(*,*) 'FLUS >> Fluid region is postprocessed at Time', Mixt%Time

end subroutine

subroutine FinishProcess(Grid,Mixt,Conf)
  type(t_Grid), intent(inout) :: Grid
  type(t_Mixt), intent(inout) :: Mixt
  type(t_Conf), intent(inout) :: Conf

  ! Delete grid variable
  call DelGrid(Grid)

  ! Delete mixture variable
  call DelMixt(Mixt)

  ! Delete configuration
  call DelConf(Conf)

  ! Delete post map
  call DelPostMap()

end subroutine

end module
