program flus_exec

  use mpi
  use flus_module
  implicit none
  integer :: rank, io, ier
  integer :: nDim, TestFlag
  real(8) :: Time, TargetTime, DeltaTime
  integer :: Iter, MaxIter, PostIter
  character(99) :: Dummy_Char

  ! fluid surface variables
  integer :: nNode
  integer :: nFace
  real(8), allocatable :: xyz(:,:)
  integer, allocatable :: f2n(:,:)
  integer, allocatable :: Patch(:)
  integer, allocatable :: BCFlag(:)
  real(8), allocatable :: InVars(:,:)
  real(8), allocatable :: OutVars(:,:)
  real(8) :: Ref_Pressure

  ! Initialize MPI process
  call MPI_Init(ier)

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Read test input file
  open(newunit=io,file='./input/test.inp')
  read(io,*)
  read(io,*)
  read(io,*) Dummy_Char, Dummy_Char, nDim
  read(io,*)
  read(io,*) Dummy_Char, Dummy_Char, Time
  read(io,*) Dummy_Char, Dummy_Char, TargetTime
  read(io,*) Dummy_Char, Dummy_Char, DeltaTime
  read(io,*)
  read(io,*) Dummy_Char, Dummy_Char, MaxIter
  read(io,*) Dummy_Char, Dummy_Char, PostIter
  read(io,*)
  read(io,*) Dummy_Char, Dummy_Char, TestFlag
  close(io)

  ! Initialize iteration
  Iter = 0

  ! Initialize fluid region
  call flus(1,Iter,DeltaTime,nNode,nFace,xyz,f2n,Patch,BCFlag,InVars,OutVars,Ref_Pressure)

  ! Post surface grid
  call PostSurfaceGrid('test.plt')

  ! Post initial fluid region
  call flus(4,Iter,DeltaTime,nNode,nFace,xyz,f2n,Patch,BCFlag,InVars,OutVars,Ref_Pressure)

  ! Solve fluid region
  do Iter = 1, MaxIter

    ! Check time step
    if( Time + DeltaTime >= TargetTime ) DeltaTime = TargetTime - Time

    ! Print current system step
    if( rank == 0 ) then
      write(*,*)
      write(*,*) 'TEST >> ======================================================'
      write(*,*) 'TEST >> Start system iteration step'
      write(*,*) 'TEST >> System iteration   =', Iter
      write(*,*) 'TEST >> System time        =', Time
      write(*,*) 'TEST >> System time step   =', DeltaTime
      write(*,*) 'TEST >> ======================================================'
      write(*,*)
    end if

    ! Grid moving test
    select case(TestFlag)
    case(1); call CircleTest
    case(2); call TubeTest
    case(3); call PropTest ! For Time Test Using APN
    end select

    ! Solve fluid region
    call flus(2,Iter,DeltaTime,nNode,nFace,xyz,f2n,Patch,BCFlag,InVars,OutVars,Ref_Pressure)

    ! Remesh fluid region
!    if( Iter == 5 ) then
!      call flus(3,Iter,DeltaTime,nNode,nFace,xyz,f2n,Patch,BCFlag,InVars,OutVars,Ref_Pressure)
!      call PostSurfaceGrid('test2.plt')
!    end if

    ! Update time
    Time = Time + DeltaTime

    ! Check target time
    if( dabs(Time-TargetTime) < 1d-12 ) exit

    ! Post fluid region
    if( mod(Iter,PostIter) == 0 ) then
      call flus(4,Iter,DeltaTime,nNode,nFace,xyz,f2n,Patch,BCFlag,InVars,OutVars,Ref_Pressure)
    end if

  end do

  ! Finalize fluid region
  call flus(4,Iter,DeltaTime,nNode,nFace,xyz,f2n,Patch,BCFlag,InVars,OutVars,Ref_Pressure)
  call flus(5,Iter,DeltaTime,nNode,nFace,xyz,f2n,Patch,BCFlag,InVars,OutVars,Ref_Pressure)

  ! Print finish status
  if( Iter > MaxIter ) Iter = MaxIter
  if( rank == 0 ) then
    write(*,*)
    write(*,*) 'TEST >> ======================================================'
    write(*,*) 'TEST >> Iteration          =', Iter
    write(*,*) 'TEST >> System Time        =', Time
    write(*,*) 'TEST >> ======================================================'
    write(*,*)
  end if

  ! Finalize MPI process
  call MPI_Finalize(ier)

  contains

  subroutine PostSurfaceGrid(FileName)
    character(*), intent(in) :: FileName

    integer :: iDim, iVar, iNode, iFace

    if( rank == 0 ) then
      open(newunit=io,file=trim(FileName))
      write(io,*) 'title="surface"'
      if( nDim == 2 ) then
        write(io,*) 'variables="x","y","var1","var2","var3"'
        write(io,*) 'zone NODES=', nNode, ', ELEMENTS=',nFace
        write(io,*) 'ZONETYPE=FELINESEG, DATAPACKING=BLOCK, VARLOCATION=([3,4,5]=CELLCENTERED)'
      end if
      if( nDim == 3 ) then
        write(io,*) 'variables="x","y","z","var1","var2","var3"'
        write(io,*) 'zone NODES=', nNode, ', ELEMENTS=',nFace
        write(io,*) 'ZONETYPE=FETRIANGLE, DATAPACKING=BLOCK, VARLOCATION=([4,5,6]=CELLCENTERED)'
      end if
      write(io,*) 'SOLUTIONTIME=',Time
      do iDim = 1, nDim
        do iNode = 1, nNode
          write(io,*) xyz(iDim,iNode)
        end do
      end do
      do iVar = 1, 3
        do iFace = 1, nFace
          write(io,*) OutVars(iVar,iFace)
        end do
      end do
      do iFace =1 , nFace
        write(io,*) f2n(:,iface)
      end do
      close(io)
    end if

  end subroutine

  subroutine CircleTest
    integer :: iDim, iNode, iFace
    integer :: ipoint(nNode)

    ipoint(:) = 0
    do iFace = 1, nFace
      if( Patch(iFace) == 13 ) then
        do iDim = 1, nDim
          iNode = f2n(iDim,iFace)
          if( ipoint(iNode) == 0 ) then
            xyz(1,iNode) = xyz(1,iNode) - 1d0 * DeltaTime
            ipoint(iNode) = 1
          end if
        end do
      end if
    end do

  end subroutine

  subroutine TubeTest
    integer :: iLocal, iNode, iFace
    integer :: ipoint(nNode)
    real(8), parameter :: vel0 = 1d0

    ipoint(:) = 0
    do iFace = 1, nFace
      if( Patch(iFace) == 4 ) then
        do iLocal = 1, nDim
          iNode = f2n(iLocal,iFace)
          if( ipoint(iNode) == 0 ) then
            if( xyz(1,iNode) < 0.5d0 ) then
              xyz(1,iNode) = xyz(1,iNode) + DeltaTime*vel0*2d0*xyz(1,iNode)
            else
              xyz(1,iNode) = xyz(1,iNode) + DeltaTime*vel0*2d0*( 1d0 - xyz(1,iNode) )
            end if
            ipoint(iNode) = 1
          end if
        end do
      end if
    end do

  end subroutine

  subroutine PropTest
    integer :: iFace
    real(8), parameter :: rhop = 1800d0
    real(8), parameter :: a_APN = 0.0080325d0
    real(8), parameter :: n_APN = 0.3345d0
    real(8), parameter :: pref = 6894733.26d0
    real(8), parameter :: T_Star = 3533d0

    do iFace = 1, nFace ! APN Model Propellant Properties
      ! Mass Flow Rate / Adiabatic Flame Temperature
      InVars(1,iFace) = rhop * a_APN * ( OutVars(1,iFace) / pref )**n_APN
      InVars(2,iFace) = T_Star
      if( BCFlag(iFace) == 0 ) then
        BCFlag(iFace) = 1
      end if
    end do

  end subroutine

end program
