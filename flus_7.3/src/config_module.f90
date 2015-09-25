module config_module

use mpi
implicit none; private
public t_Conf
public DelConf
public SetConf
public SetBdryConf

type t_Conf
  ! Grid configuration
  integer :: GridFileFormat ! grid file format
  character(32) :: GridFileName ! grid file name
  ! Simulation configuration
  logical :: Restart ! restart index
  logical :: MoveGrid ! grid moving index
  integer :: nPrint ! simulation status print interval
  ! Solver configuration
  integer :: FlowModel ! flow model index
  integer :: iFlux ! flux scheme index
  integer :: iLimit ! limiter scheme index
  real(8) :: vkt ! limiter coefficient
  ! ODE configuration
  integer :: iTime ! time scheme index
  logical :: Unsteady ! unsteady simulation index
  logical :: DualTimeStep ! dual-time stepping index
  logical :: LocalTimeStep ! local time stepping index
  real(8) :: CFL ! Courant-Friedrichs-Lewy condition
  real(8) :: CFL_DualTimeStep ! CFL for dual-time stepping
  real(8) :: Error ! error for convergence
  ! Boundary condition configuration
  integer :: nBC ! number of boundary condition
  character(32), allocatable :: BCName(:) ! BC name
  integer, allocatable :: BCType(:) ! BC type
  integer, allocatable :: iPatch(:) ! Patch index for BC type
  real(8), allocatable :: bv(:,:) ! values for BC
  logical :: Axisymmetric ! axisymmetric index
  logical :: indIgnitor ! ignitor model index
  integer :: nIgni ! number of data in ignitor file
  real(8), allocatable :: tIgni(:) ! time data in ignitor file
  real(8), allocatable :: mIgni(:) ! mass flux data in ignitor file
  ! Fluid initial condition configuration
  real(8) :: Gamma ! heat capacity ratio
  real(8) :: GasConstant ! specific gas constant
  real(8) :: Prandtl_Lam ! laminar Prandtl number
  real(8) :: Prandtl_Turb ! turbulent Prandtl number
  real(8) :: Viscosity_Ref ! reference viscosity
  real(8) :: Temperature_Ref ! reference temperature
  real(8) :: Temperature_Eff ! effective temperature
  real(8) :: Length_Ref ! reference length
  real(8) :: iniTemp ! initial temperature
  real(8) :: iniPres ! initial pressure
  real(8) :: iniMach ! initial mach number
  real(8) :: iniAngl ! initial angle of attach
end type

contains

subroutine DelConf(Conf)
  type(t_Conf), intent(inout) :: Conf

  if( allocated( Conf%BCName ) ) deallocate( Conf%BCName )
  if( allocated( Conf%BCType ) ) deallocate( Conf%BCType )
  if( allocated( Conf%iPatch ) ) deallocate( Conf%iPatch )
  if( allocated( Conf%bv ) ) deallocate( Conf%bv )
  if( allocated( Conf%tIgni ) ) deallocate( Conf%tIgni )
  if( allocated( Conf%mIgni ) ) deallocate( Conf%mIgni )

end subroutine

subroutine SetConf(Conf,FileName)
  type(t_Conf), intent(out) :: Conf
  character(*), intent(in) :: FileName

  integer :: rank, size, io, error, ier
  character(32) :: Dummy_Char

  ! Get rank of each core / number of cores
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)
  call MPI_Comm_size(MPI_COMM_WORLD,size,ier)

  ! Open configuration file
  open(newunit=io,file=trim(FileName),status='OLD',action='READ',iostat=ier)

  ! Check file existence and print progress
  if( ier > 0 ) then
    if( rank == 0 ) then
      write(*,*) 'FLUS >> Cannot find input file: ', trim(FileName)
      call MPI_Abort(MPI_COMM_WORLD,error,ier)
    end if
  else
    if( rank == 0 ) write(*,*) 'FLUS >> Reading input file: ', trim(FileName)
  end if

  ! Set MPI barrier
  call MPI_Barrier(MPI_COMM_WORLD,ier)

  ! Read configuration file
  read(io,*)
  read(io,*)
  read(io,*) Dummy_Char, Dummy_Char, Conf%GridFileFormat
  read(io,*) Dummy_Char, Dummy_Char, Conf%GridFileName
  read(io,*)
  read(io,*) Dummy_Char, Dummy_Char, Conf%Restart
  read(io,*) Dummy_Char, Dummy_Char, Conf%MoveGrid
  read(io,*) Dummy_Char, Dummy_Char, Conf%nPrint
  read(io,*)
  read(io,*) Dummy_Char, Dummy_Char, Conf%FlowModel
  read(io,*) Dummy_Char, Dummy_Char, Conf%iFlux
  read(io,*) Dummy_Char, Dummy_Char, Conf%iLimit
  read(io,*) Dummy_Char, Dummy_Char, Conf%vkt
  read(io,*)
  read(io,*) Dummy_Char, Dummy_Char, Conf%iTime
  read(io,*) Dummy_Char, Dummy_Char, Conf%Unsteady
  read(io,*) Dummy_Char, Dummy_Char, Conf%DualTimeStep
  read(io,*) Dummy_Char, Dummy_Char, Conf%LocalTimeStep
  read(io,*) Dummy_Char, Dummy_Char, Conf%CFL
  read(io,*) Dummy_Char, Dummy_Char, Conf%CFL_DualTimeStep
  read(io,*) Dummy_Char, Dummy_Char, Conf%Error

  ! Close Conf file
  close(io)

  ! Print configuration info
  if( rank == 0 ) then

    ! Restart information
    if( Conf%Restart ) write(*,*) 'FLUS >> Starting from restart file'

    ! Move grid information
    if( Conf%MoveGrid ) write(*,*) 'FLUS >> Moving grid based on ALE'

    ! Flow model information
    select case(Conf%FlowModel)
    case(1); write(*,*) 'FLUS >> Euler simulation'
    case(2); write(*,*) 'FLUS >> N-S simulation'
    case(3); write(*,*) 'FLUS >> RANS simulation'
    case default
      write(*,*) 'FLUS >> Invalid flow model flag'
      call MPI_Abort(MPI_COMM_WORLD,error,ier)
    end select

    ! Flux scheme information
    select case(Conf%iFlux)
    case(1); write(*,*) 'FLUS >> Flux method: Roe FDS'
    case(2); write(*,*) 'FLUS >> Flux method: RoeM2'
    case(3); write(*,*) 'FLUS >> Flux method: AUSM+'
    case(4); write(*,*) 'FLUS >> Flux method: AUSMPW+'
    case default
      write(*,*) 'FLUS >> Invalid flux method'
      call MPI_Abort(MPI_COMM_WORLD,error,ier)
    end select

    ! Limiter information
    select case(Conf%iLimit)
    case(1); write(*,*) 'FLUS >> Limit method: 1st order'
    case(2); write(*,*) 'FLUS >> Limit method: Barth limiter'
    case(3); write(*,*) 'FLUS >> Limit method: Venkatakrishnan limiter, K=', sngl(Conf%vkt)
    case(4); write(*,*) 'FLUS >> Limit method: MLP-u1'
    case(5); write(*,*) 'FLUS >> Limit method: MLP-u2, K=', sngl(Conf%vkt)
    case default
      write(*,*) 'FLUS >> Invalid limiter method'
      call MPI_Abort(MPI_COMM_WORLD,error,ier)
    end select

    ! Time scheme information
    select case(Conf%iTime)
    case(1); write(*,*) 'FLUS >> Time method: Grid moving test'
    case(2); write(*,*) 'FLUS >> Time method: Euler-Explicit'
    case(3); write(*,*) 'FLUS >> Time method: TVD-RK 3rd'
    case(4); write(*,*) 'FLUS >> Time method: LUSGS'
    case default
      write(*,*) 'FLUS >> Invalid time method'
      call MPI_Abort(MPI_COMM_WORLD,error,ier)
    end select

    ! Time scheme information
    if( Conf%Unsteady ) then
      if( Conf%DualTimeStep ) then
        if( Conf%LocalTimeStep ) then
          write(*,*) 'FLUS >> Unsteady dual time stepping with local time step'
        else
          write(*,*) 'FLUS >> Unsteady dual time stepping with global time step'
        end if
      else
        if( Conf%LocalTimeStep ) then
          write(*,*) 'FLUS >> Invalid option: Unsteady with local time step'
          call MPI_Abort(MPI_COMM_WORLD,error,ier)
        else
          write(*,*) 'FLUS >> Unsteady with global time step'
        end if
      end if
    else
      if( Conf%LocalTimeStep ) then
        write(*,*) 'FLUS >> Steady with local time step'
      else
        write(*,*) 'FLUS >> Steady with global time step'
      end if
    end if

    ! Mesh moving method / data transfer method / post format info
    write(*,*) 'FLUS >> Mesh smooth: IDW smoothing'
    write(*,*) 'FLUS >> Data transfer: Piecewise linear interpolation'
    write(*,*) 'FLUS >> Post format: FLUS format'

    ! Print number of cores
    write(*,*) 'FLUS >> Number of cores: ', size

    write(*,*) 'FLUS >>'

  end if

end subroutine

subroutine SetBdryConf(Conf,FileName)
  type(t_Conf), intent(inout) :: Conf
  character(*), intent(in) :: FileName

  integer :: rank, io, error, ier
  integer :: iBC, BCFlag, BCType, iPatch, iIgni
  real(8) :: Dummy_Real
  character(32) :: String, Dummy_Char
  character(32) :: IgnitorFileName

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Open boundary condition file
  open(newunit=io,file=trim(FileName),status='OLD',action='READ',iostat=ier)

  ! Check file existence and print progress
  if( ier > 0 ) then
    if( rank == 0 ) then
      write(*,*) 'FLUS >> Cannot find boundary input file: ', trim(FileName)
      call MPI_Abort(MPI_COMM_WORLD,error,ier)
    end if
  else
    if( rank == 0 ) write(*,*) 'FLUS >> Reading boundary input file: ', trim(FileName)
  end if

  ! Set MPI barrier
  call MPI_Barrier(MPI_COMM_WORLD,ier)

  ! Count number of boundary condition
  Conf%nBC = 0
  do ! Read one string
    read(io,*,iostat=ier) String
    if( ier < 0 ) exit
    ! Check condition section
    BCFlag = scan(String,'#',back=.false.)
    if( BCFlag /= 0 ) Conf%nBC = Conf%nBC + 1
  end do
  ! Exclude initial condition section
  Conf%nBC = Conf%nBC - 1

  ! Rewind BC file
  rewind(io)

  ! Allocate boundary information
  allocate( Conf%BCName(Conf%nBC) )
  allocate( Conf%BCType(Conf%nBC) )
  allocate( Conf%iPatch(Conf%nBC) )
  allocate( Conf%bv(5,Conf%nBC) )

  iBC = 0
  do ! Read one string
    read(io,*,iostat=ier) String
    if( ier < 0 ) exit
    ! Check condition section
    BCFlag = scan(String,'#',back=.false.)
    if( BCFlag /= 0 ) then
      ! Read patch info and BC type
      read(io,*) Dummy_Char, iPatch
      read(io,*) Dummy_Char, BCType
      if( BCType == 0 ) then ! Initial condition when BCType is equal to 0
        read(io,*) Dummy_Char, Conf%Gamma   ! gamma
        read(io,*) Dummy_Char, Conf%GasConstant ! specific gas constant
        read(io,*) Dummy_Char, Conf%Prandtl_Lam ! laminar Prandtl number
        read(io,*) Dummy_Char, Conf%Prandtl_Turb ! turbulent Prandtl number
        read(io,*) Dummy_Char, Conf%Viscosity_Ref ! reference viscosity
        read(io,*) Dummy_Char, Conf%Temperature_Ref ! reference temperature
        read(io,*) Dummy_Char, Conf%Temperature_Eff ! effective temperature
        read(io,*) Dummy_Char, Conf%Length_Ref ! reference length
        read(io,*) Dummy_Char, Conf%iniPres ! static pressure
        read(io,*) Dummy_Char, Conf%iniTemp ! static temperature
        read(io,*) Dummy_Char, Conf%iniMach ! mach number
        read(io,*) Dummy_Char, Conf%iniAngl ! flow angle
      else ! Boudary condition when BCType is not equal to 0
        iBC = iBC + 1
        Conf%BCName(iBC) = trim(String)
        Conf%BCType(iBC) = BCType
        Conf%iPatch(iBC) = iPatch
        ! bv(:,*) = variables for boundary calculation
        select case(BCType)
        case(1) ! Inviscid wall
        case(2) ! Viscous wall
        case(3) ! Pressure inlet
          read(io,*) Dummy_Char, Conf%bv(1,iBC) ! total pressure
          read(io,*) Dummy_Char, Conf%bv(2,iBC) ! total temperature
          read(io,*) Dummy_Char, Conf%bv(3,iBC) ! mach number
          read(io,*) Dummy_Char, Conf%bv(4,iBC) ! flow angle
        case(4) ! Injection inlet
          read(io,*) Dummy_Char, Conf%bv(1,iBC) ! mass flux
          read(io,*) Dummy_Char, Conf%bv(2,iBC) ! temperature
        case(5) ! Pressure outlet
          read(io,*) Dummy_Char, Conf%bv(1,iBC) ! static pressure
        case(6) ! Pressure farfield
          read(io,*) Dummy_Char, Conf%bv(1,iBC) ! pressure
          read(io,*) Dummy_Char, Conf%bv(2,iBC) ! temperature
          read(io,*) Dummy_Char, Conf%bv(3,iBC) ! mach number
          read(io,*) Dummy_Char, Conf%bv(4,iBC) ! flow angle
        case(7) ! Symmetric
        case(8) ! Axisymmetric
          Conf%Axisymmetric = .true.
        case(9) ! Membrane nozzle
          read(io,*) Dummy_Char, Conf%bv(1,iBC) ! static pressure
        case(10) ! Propellant
        case(11) ! Ignitor model
          read(io,*) Dummy_Char, Conf%bv(1,iBC) ! total temperature
          read(io,*) Dummy_Char, Conf%indIgnitor ! ignitor model flag
          if( Conf%indIgnitor ) then ! ignitor model
            read(io,*) Dummy_Char, Conf%bv(2,iBC) ! Pcmax
            read(io,*) Dummy_Char, Conf%bv(3,iBC) ! K1
            read(io,*) Dummy_Char, Conf%bv(4,iBC) ! K2
            read(io,*) Dummy_Char, Conf%bv(5,iBC) ! tb
          else ! ignitor file
            read(io,*) Dummy_Char, IgnitorFileName ! ignitor file name
          end if
        case(12) ! Pressure sensor
          read(io,*) Dummy_Char, Conf%bv(1,iBC) ! membrane pressure
        case(13) ! Nozzle throat
        case(14) ! APN model
          read(io,*) Dummy_Char, Conf%bv(1,iBC) ! propellant desnity
          read(io,*) Dummy_Char, Conf%bv(2,iBC) ! total temperature
          read(io,*) Dummy_Char, Conf%bv(3,iBC) ! a in APN
          read(io,*) Dummy_Char, Conf%bv(4,iBC) ! n in APN
        end select
      end if
    end if
  end do

  ! Close BC file
  close(io)

  ! Print boundary condition
  if( rank == 0 ) then
    do iBC = 1, Conf%nBC
      select case(Conf%BCType(iBC))
      case(1);  write(*,*) 'FLUS >> ', trim(Conf%BCName(iBC)), ' ( Inviscid Wall )'
      case(2);  write(*,*) 'FLUS >> ', trim(Conf%BCName(iBC)), ' ( Viscous Wall )'
      case(3);  write(*,*) 'FLUS >> ', trim(Conf%BCName(iBC)), ' ( Pressure Inlet )'
      case(4);  write(*,*) 'FLUS >> ', trim(Conf%BCName(iBC)), ' ( Injection Inlet )'
      case(5);  write(*,*) 'FLUS >> ', trim(Conf%BCName(iBC)), ' ( Pressure Outlet )'
      case(6);  write(*,*) 'FLUS >> ', trim(Conf%BCName(iBC)), ' ( Pressure Farfield )'
      case(7);  write(*,*) 'FLUS >> ', trim(Conf%BCName(iBC)), ' ( Symmetric Wall )'
      case(8);  write(*,*) 'FLUS >> ', trim(Conf%BCName(iBC)), ' ( Axisymmetric Wall )'
      case(9);  write(*,*) 'FLUS >> ', trim(Conf%BCName(iBC)), ' ( Membrane Wall )'
      case(10); write(*,*) 'FLUS >> ', trim(Conf%BCName(iBC)), ' ( Propellant Wall )'
      case(11); write(*,*) 'FLUS >> ', trim(Conf%BCName(iBC)), ' ( Ignitor Model )'
      case(12); write(*,*) 'FLUS >> ', trim(Conf%BCName(iBC)), ' ( Pressure Sensor )'
      case(13); write(*,*) 'FLUS >> ', trim(Conf%BCName(iBC)), ' ( Nozzle Throat )'
      case(14); write(*,*) 'FLUS >> ', trim(Conf%BCName(iBC)), ' ( APN Model )'
      end select
      write(*,*) 'FLUS >> Zone ', Conf%iPatch(iBC)
    end do
    write(*,*) 'FLUS >>'
  end if

  ! Read ignitor boundary condition file
  do iBC = 1, Conf%nBC
    if( Conf%BCType(iBC) == 11 .and. .not. Conf%indIgnitor ) then
    
      ! Open ignitor file
      open(newunit=io,file=trim(IgnitorFileName),status='OLD',action='READ',iostat=ier)
  
      ! Check file existence and print progress
      if( ier > 0 ) then
        write(*,*) 'FLUS >> Error : cannot find input file: '//trim(IgnitorFileName)
        call MPI_Abort(MPI_COMM_WORLD,error,ier)
      end if
  
      ! Count Number of Data
      Conf%nIgni = 0
      do while(.true.)
        read(io,*,end=20) Dummy_Real, Dummy_Real
        Conf%nIgni = Conf%nIgni + 1
      end do
20    rewind(io)
  
      ! Read Ignitor Data
      allocate( Conf%tIgni(Conf%nIgni) )
      allocate( Conf%mIgni(Conf%nIgni) )
      do iIgni = 1, Conf%nIgni
        read(io,*) Conf%tIgni(iIgni), Conf%mIgni(iIgni)
      end do

      ! Close ignitor file
      close(io)
    
    end if
  end do

end subroutine

end module
