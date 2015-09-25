module converter_module

use mpi
use ifport
use lib_vtk_io
implicit none; private
include "tecio.f90"
public DelData
public ReadFromFLUS
public PostToPLT
public PostToVTK

! Global variable data
integer :: nDim ! grid dimension
integer :: nNodeCell ! number of nodes per cell
integer :: nNodeFace ! number of nodes per face
integer :: nNode ! number of nodes
integer :: nCell ! number of cells
integer :: nFace ! number of faces
integer :: nVar ! number variables
character(32) :: VarName(20) ! variable name
real(8) :: SolutionTime ! solution time
real(8), allocatable :: xyz(:,:) ! node coordinates
integer, allocatable :: c2n(:,:) ! node connectivity of the cell
integer, allocatable :: f2n(:,:) ! node connectivity of the face
integer, allocatable :: Patch(:) ! patch of the face
real(8), allocatable :: Var(:,:) ! variable of the cell
integer, allocatable :: Part(:) ! partition of the cell

contains

subroutine DelData

  ! deallocate global variables
  if( allocated( xyz ) ) deallocate( xyz )
  if( allocated( c2n ) ) deallocate( c2n )
  if( allocated( f2n ) ) deallocate( f2n )
  if( allocated( Patch ) ) deallocate( Patch )
  if( allocated( Var ) ) deallocate( Var )
  if( allocated( Part ) ) deallocate( Part )

end subroutine

subroutine ReadFromFLUS(InputName)
  character(*), intent(in) :: InputName

  integer :: rank, io, ier

  ! FLUS(*.flu) file format layout
  ! Grid dimension
  ! Number of nodes
  ! Number of cells
  ! Number of faces
  ! Number of variables
  ! Variable name
  ! Time data
  ! Node coordinates
  ! Node connectivity of the cell
  ! Node connectivity of the face
  ! Patch of the face
  ! Variable of the cell
  ! Partition of the cell (temporary)

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Open FLUS file using MPI-IO
  call MPI_File_open(MPI_COMM_WORLD,trim(InputName),MPI_MODE_RDONLY,MPI_INFO_NULL,io,ier)

  ! Read nDim and set element info
  call MPI_File_read_all(io,nDim,1,MPI_INTEGER,MPI_STATUS_IGNORE,ier)
  select case(nDim)
  case(2) ! triangle / line for 2-D
    nNodeCell = 3; nNodeFace = 2
  case(3) ! tetrahedron / triangle for 3-D
    nNodeCell = 4; nNodeFace = 3
  end select

  ! Read number of nodes / cells / faces / variables
  call MPI_File_read_all(io,nNode,1,MPI_INTEGER,MPI_STATUS_IGNORE,ier)
  call MPI_File_read_all(io,nCell,1,MPI_INTEGER,MPI_STATUS_IGNORE,ier)
  call MPI_File_read_all(io,nFace,1,MPI_INTEGER,MPI_STATUS_IGNORE,ier)
  call MPI_File_read_all(io,nVar, 1,MPI_INTEGER,MPI_STATUS_IGNORE,ier)

  ! Read variable name
  call MPI_File_read_all(io,VarName,32*nVar,MPI_CHARACTER,MPI_STATUS_IGNORE,ier)

  ! Read time data
  call MPI_File_read_all(io,SolutionTime,1,MPI_REAL8,MPI_STATUS_IGNORE,ier)

  ! Read node coordinates
  allocate( xyz(nDim,nNode) )
  call MPI_File_read_all(io,xyz,nDim*nNode,MPI_REAL8,MPI_STATUS_IGNORE,ier)

  ! Read node connectivity of the cell
  allocate( c2n(nNodeCell,nCell) )
  call MPI_File_read_all(io,c2n,nNodeCell*nCell,MPI_INTEGER,MPI_STATUS_IGNORE,ier)

  ! Read node connectivity of the face
  allocate( f2n(nNodeFace,nFace) )
  call MPI_File_read_all(io,f2n,nNodeFace*nFace,MPI_INTEGER,MPI_STATUS_IGNORE,ier)

  ! Read patch of the face
  allocate( Patch(nFace) )
  call MPI_File_read_all(io,Patch,nFace,MPI_INTEGER,MPI_STATUS_IGNORE,ier)

  ! Read variable of the cell
  allocate( Var(nCell,nVar) )
  call MPI_File_read_all(io,Var,nCell*nVar,MPI_REAL8,MPI_STATUS_IGNORE,ier)

  ! Read partition of the cell
  allocate( Part(nCell) )
  call MPI_File_read_all(io,Part,nCell,MPI_INTEGER,MPI_STATUS_IGNORE,ier)

end subroutine

subroutine PostToPLT(MidName,OutputName)
  character(*), intent(in) :: MidName
  character(*), intent(out) :: OutputName

  integer :: iDim, iVar, ier
  integer :: FileType, ZoneType
  integer :: Debug, Visdouble, Disdouble
  integer :: VarLocation(nDim+nVar+1)
  character(256) :: VarsName
  character(1) :: NulChar = Char(0)
  integer :: Null(*)
  pointer    (nullptr,Null)
  logical, save :: MakeDir = .true.

  ! Make directory for plt file
  if( MakeDir ) then
    ier = MAKEDIRQQ('plt_format')
    write(*,*) '  Make directory: plt_format'
    MakeDir = .false.
  end if

  ! Set output file name
  OutputName = './plt_format/'//trim(MidName)//'.plt'

  ! Set variables name
  VarsName = 'x y'
  if( nDim == 3 ) VarsName = trim(VarsName)//' z'
  do iVar = 1, nVar
    VarsName = trim(VarsName)//' '//trim(VarName(iVar))
  end do
  VarsName = trim(VarsName)//' Partition'

  ! Set variable location
  VarLocation(:) = 0
  VarLocation(1:nDim) = 1

  ! Set tecio variables
  FileType = 0
  select case(nDim)
  case(2); ZoneType = 2
  case(3); ZoneType = 4
  end select
  Debug = 0
  Visdouble = 0
  Disdouble = 0
  nullptr = 0

  ! Open plt file
  ier = tecini112("Result"//NulChar,trim(VarsName)//NulChar,trim(OutputName)//NulChar,'.',FileType,Debug,Visdouble)

  ! Write zone header
  ier = teczne112('Fluid',2*(nDim-1),nNode,nCell,1,0,0,0,SolutionTime,1,0,1,0,0,0,0,0,Null,VarLocation,Null,0)

  ! Write node coordinates
  do iDim = 1, nDim
    ier = tecdat112(nNode,sngl(xyz(iDim,:)),Disdouble)
  end do

  ! Write variables of the cell
  do iVar = 1, nVar
    ier = tecdat112(nCell,sngl(Var(:,iVar)),Disdouble)
  end do

  ! Write parition of the cell
  ier = tecdat112(nCell,sngl(Part(:)),Disdouble)

  ! Write node connectivity of the cell
  ier = tecnod112(c2n(:,:))

  ! Close plt file
  ier = tecend112()

end subroutine

subroutine PostToVTK(MidName,OutputName)
  character(*), intent(in) :: MidName
  character(*), intent(out) :: OutputName

  integer :: iVar, iCell, io, ier
  integer :: Con, Cell_Type, Flag1, Flag2
  real(8), allocatable :: Zero(:)
  integer, allocatable :: Connect(:)
  integer, allocatable :: Offset(:)
  integer(1), allocatable :: ElemType(:)
  logical, save :: MakeDir = .true.

  ! Make directory and master file for vtk file
  if( MakeDir ) then
    ier = MAKEDIRQQ('vtk_format')
    write(*,*) '  Make directory: vtk_format'
    MakeDir = .false.
    ! Create master pvd file
    open(newunit=io,file='./vtk_format/master.pvd')
    write(io,*) '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
    write(io,*) '<Collection>'
    write(io,*) '<DataSet timestep="', SolutionTime, '" file="', trim(MidName),'.vtu"/>'
    write(io,*) '</Collection>'
    write(io,*) '</VTKFile>'
    close(io)
  else
    open(newunit=io,file='./vtk_format/master.pvd',access='APPEND')
    backspace(io)
    backspace(io)
    write(io,*) '<DataSet timestep="', SolutionTime, '" file="', trim(MidName),'.vtu"/>'
    write(io,*) '</Collection>'
    write(io,*) '</VTKFile>'
    close(io)
  end if

  ! Set output file name
  OutputName = './vtk_format/'//trim(MidName)//'.vtu'

  ! Open vtk xml file
  io = VTK_INI_XML('Binary',trim(OutputName),'UnstructuredGrid')

  ! Write piece info
  select case(nDim)
  case(2)
    allocate( Zero(nNode) ); Zero(:) = 0d0
    io = VTK_GEO_XML(nNode,nCell,xyz(1,:),xyz(2,:),Zero(:))
    deallocate( Zero )
  case(3)
    io = VTK_GEO_XML(nNode,nCell,xyz(1,:),xyz(2,:),xyz(3,:))
  end select

  ! Set element info
  select case(nDim)
  case(2) ! triangle
    Con = 3; Cell_Type = 5
  case(3) ! tetrahedron
    Con = 4; Cell_Type = 10
  end select

  ! Write node connectivity of the cell
  allocate( Connect(Con*nCell) )
  allocate( Offset(nCell) )
  allocate( ElemType(nCell) )
  do iCell = 1, nCell
    Flag1 = Con * (iCell - 1) + 1
    Flag2 = Flag1 + nDim
    Connect(Flag1:Flag2) = c2n(:,iCell) - 1
    Offset(iCell) = Con * iCell
    ElemType(iCell) = Cell_Type
  end do
  io = VTK_CON_XML(nCell, Connect(:), Offset(:), ElemType(:))
  deallocate( Connect )
  deallocate( Offset )
  deallocate( ElemType )

  ! Start data writing
  io = VTK_DAT_XML('CELL', 'OPEN')

  ! Write variable of the cell
  do iVar = 1, nVar
    io = VTK_VAR_XML(nCell,trim(VarName(iVar)),Var(:,iVar))
  end do

  ! Write partition of the cell
  io = VTK_VAR_XML(nCell,'Partition',Part(:))

  ! End data writing
  io = VTK_DAT_XML('Cell', 'CLOSE')

  ! End geometry writing
  io = VTK_GEO_XML()

  ! Close vtk xml file
  io = VTK_END_XML()

end subroutine

end module

program converter

  use mpi
  use converter_module
  implicit none
  integer :: rank, argc, ier
  integer :: Flag1, Flag2
  character(32) :: Option
  character(256) :: InputName
  character(256) :: MidName
  character(256) :: OutputName

  ! Initialize MPI
  call MPI_Init(ier)

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Get option (first argument)
  call GETARG(1,Option)

  ! Check option
  if( rank == 0 ) write(*,*)
  select case(trim(Option))
  case('-plt') ! plt format for tecplot
    if( rank == 0 ) write(*,*) 'File conversion start(plt)'
  case('-vtk') ! vtk format for paraview
    if( rank == 0 ) write(*,*) 'File conversion start(vtk)'
  case default ! invalid option
    if( rank == 0 ) write(*,*) 'Unknown option "'//trim(Option)//'"'
    if( rank == 0 ) write(*,*) 'Accepted values: -plt, -vtk'
    call MPI_Finalize(ier); stop
  end select

  do argc = 2, IARGC()

    ! Get argument
    call GETARG(argc,InputName)

    ! Read FLUS file
    call ReadFromFLUS(InputName)

    ! Set mid name (only file name)
    Flag1 = scan(InputName,'/', BACK=.TRUE.)
    Flag2 = scan(InputName,'.', BACK=.TRUE.)
    MidName = trim(InputName(Flag1+1:Flag2-1))

    ! Convert to other format
    select case(trim(Option))
    case('-plt'); if( rank == 0 ) call PostToPLT(MidName,OutputName)
    case('-vtk'); if( rank == 0 ) call PostToVTK(MidName,OutputName)
    end select

    ! Delete data
    call DelData()

    ! Print progress
    if( rank == 0 ) write(*,*) '  ', trim(InputName), ' -> ', trim(OutputName)

  end do

  ! Print progress
  if( rank == 0 ) then
    write(*,*) 'File conversion complete'
    write(*,*)
  end if

  ! Finalize MPI
  call MPI_Finalize(ier)

end program
