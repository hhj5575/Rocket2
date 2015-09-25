module post_module

use mpi
use config_module
use grid_module
use mixture_module
implicit none; private
public DelPostMap
public SetPostMap
public Post_FLUS_Format
public Post_Pressure

integer, allocatable :: NodeFileMap(:) ! node map in file
integer, allocatable :: CellFileMap(:) ! cell map in file
integer, allocatable :: FaceFileMap(:) ! face map in file
integer, allocatable :: NodeProcMap(:) ! node map in process
integer, allocatable :: CellProcMap(:) ! node map in process
integer, allocatable :: FaceProcMap(:) ! node map in process
integer :: intsize ! MPI_INTEGER size
integer :: real8size ! MPI_REAL8 size
integer :: CoordType ! MPI coordinate type
integer :: NodeCellType ! MPI nodecell type
integer :: NodeFaceType ! MPI nodeface type
integer :: NodeType ! MPI type for writing node coordinates
integer :: CellType ! MPI type for writing node connectivity of the cell
integer :: FaceType ! MPI type for writing node connectivity of the face
integer :: FaceIntType ! MPI type for writing integer data in face
integer :: CellIntType ! MPI type for writing integer data in cell
integer :: CellReal8Type ! MPI type for writing real8 data in cell

contains

subroutine DelPostMap

  integer :: ier

  if( allocated( NodeFileMap ) ) deallocate( NodeFileMap )
  if( allocated( CellFileMap ) ) deallocate( CellFileMap )
  if( allocated( FaceFileMap ) ) deallocate( FaceFileMap )
  if( allocated( NodeProcMap ) ) deallocate( NodeProcMap )
  if( allocated( CellProcMap ) ) deallocate( CellProcMap )
  if( allocated( FaceProcMap ) ) deallocate( FaceProcMap )

  call MPI_Type_free(CoordType,ier)
  call MPI_Type_free(NodeCellType,ier)
  call MPI_Type_free(NodeFaceType,ier)
  call MPI_Type_free(NodeType,ier)
  call MPI_Type_free(CellType,ier)
  call MPI_Type_free(FaceType,ier)
  call MPI_Type_free(FaceIntType,ier)
  call MPI_Type_free(CellIntType,ier)
  call MPI_Type_free(CellReal8Type,ier)

end subroutine

subroutine SetPostMap(Grid)
  type(t_Grid), intent(in) :: Grid

  integer :: ier
  integer :: iNode, iNodeLocal, iNodeGlobal
  integer :: iCell, iCellLocal, iCellGlobal
  integer :: iFace, iFaceLocal, iFaceGlobal

  ! Create node file/proc map
  allocate( NodeFileMap(Grid%nNode) )
  allocate( NodeProcMap(Grid%nNode) )
  iNode = 0
  do iNodeGlobal = 1, size(Grid%Global_to_Local_Node)
    iNodeLocal = Grid%Global_to_Local_Node(iNodeGlobal)
    if( iNodeLocal == 0 .or. iNodeLocal > Grid%nNode ) cycle
    iNode = iNode + 1
    NodeFileMap(iNode) = iNodeGlobal-1
    NodeProcMap(iNode) = iNodeLocal
  end do

  ! Create cell file/proc map
  allocate( CellFileMap(Grid%nCell) )
  allocate( CellProcMap(Grid%nCell) )
  iCell = 0
  do iCellGlobal = 1, size(Grid%Global_to_Local_Cell)
    iCellLocal = Grid%Global_to_Local_Cell(iCellGlobal)
    if( iCellLocal == 0 .or. iCellLocal > Grid%nCell ) cycle
    iCell = iCell + 1
    CellFileMap(iCell) = iCellGlobal-1
    CellProcMap(iCell) = iCellLocal
  end do

  ! Create face file / proc map
  allocate( FaceFileMap(Grid%nFaceBound) )
  allocate( FaceProcMap(Grid%nFaceBound) )
  iFace = 0
  do iFaceGlobal = 1, size(Grid%Global_to_Local_Face)
    iFaceLocal = Grid%Global_to_Local_Face(iFaceGlobal)
    if( iFaceLocal == 0 ) cycle
    iFace = iFace + 1
    FaceFileMap(iFace) = iFaceGlobal-1
    FaceProcMap(iFace) = iFaceLocal
  end do

  ! Set type size
  call MPI_Type_size(MPI_INTEGER,intsize,ier)
  call MPI_Type_size(MPI_REAL8,real8size,ier)

  ! Create coordinate type
  call MPI_Type_contiguous(Grid%nDim,MPI_REAL8,CoordType,ier)
  call MPI_Type_commit(CoordType,ier)

  ! Create nodecell type
  call MPI_Type_contiguous(Grid%nNodeCell,MPI_INTEGER,NodeCellType,ier)
  call MPI_Type_commit(NodeCellType,ier)

  ! Create nodeface type
  call MPI_Type_contiguous(Grid%nNodeFace,MPI_INTEGER,NodeFaceType,ier)
  call MPI_Type_commit(NodeFaceType,ier)

  ! Create node type for writing coordinates
  call MPI_Type_create_indexed_block(Grid%nNode,1,NodeFileMap,CoordType,NodeType,ier)
  call MPI_Type_commit(NodeType,ier)

  ! Create cell type for writing node connectivity of the cell
  call MPI_Type_create_indexed_block(Grid%nCell,1,CellFileMap,NodeCellType,CellType,ier)
  call MPI_Type_commit(CellType,ier)

  ! Create face type for writing node connectivity of the face
  call MPI_Type_create_indexed_block(Grid%nFaceBound,1,FaceFileMap,NodeFaceType,FaceType,ier)
  call MPI_Type_commit(FaceType,ier)

  ! Create face int type for writing integer data in face
  call MPI_Type_create_indexed_block(Grid%nFaceBound,1,FaceFileMap,MPI_INTEGER,FaceIntType,ier)
  call MPI_Type_commit(FaceIntType,ier)

  ! Create cell int type for writing integer data in cell
  call MPI_Type_create_indexed_block(Grid%nCell,1,CellFileMap,MPI_INTEGER,CellIntType,ier)
  call MPI_Type_commit(CellIntType,ier)

  ! Create cell real8 type for writing real8 data in cell
  call MPI_Type_create_indexed_block(Grid%nCell,1,CellFileMap,MPI_REAL8,CellReal8Type,ier)
  call MPI_Type_commit(CellReal8Type,ier)

end subroutine

subroutine Post_FLUS_Format(Grid,Mixt,Conf,FileName)
  type(t_Grid), intent(in) :: Grid
  type(t_Mixt), intent(in) :: Mixt
  type(t_Conf), intent(in) :: Conf
  character(*), intent(in) :: FileName

  integer :: iLocal, iVar, io, ier
  integer :: iNode, iNodeLocal
  integer :: iCell, iCellLocal
  integer :: iFace, iFaceLocal
  integer :: Header(5)
  character(32) :: VarName(0:11)
  real(8), allocatable :: Coord(:,:)
  integer, allocatable :: Nodes(:,:)
  integer, allocatable :: DataInt(:)
  real(8), allocatable :: DataReal8(:)
  integer(kind=MPI_OFFSET_KIND) :: disp

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

  ! Open FLUS file using MPI-IO
  call MPI_File_open(MPI_COMM_WORLD,trim(FileName),MPI_MODE_WRONLY+MPI_MODE_CREATE,MPI_INFO_NULL,io,ier)

  ! Set file header
  Header(1) = Grid%nDim
  Header(2) = Grid%nNodeGlobal
  Header(3) = Grid%nCellGlobal
  Header(4) = Grid%nFaceGlobal
  Header(5) = Mixt%nPVar+2

  ! Write header
  call MPI_File_write_all(io,Header,5,MPI_INTEGER,MPI_STATUS_IGNORE,ier)

  ! Set variable name
  select case(Conf%FlowModel)
  case(1) ! Euler
    select case(Grid%nDim)
    case(2)
      VarName(0) = 'Temperature'; VarName(1) = 'VelocityX'; VarName(2) = 'VelocityY'
      VarName(3) = 'Pressure'   ; VarName(4) = 'Density'  ; VarName(5) = 'Mach'
    case(3)
      VarName(0) = 'Temperature'; VarName(1) = 'VelocityX'; VarName(2) = 'VelocityY'
      VarName(3) = 'VelocityZ'  ; VarName(4) = 'Pressure' ; VarName(5) = 'Density'
      VarName(6) = 'Mach'
    end select
  case(2) ! N-S
    select case(Grid%nDim)
    case(2)
      VarName(0) = 'Temperature'; VarName(1) = 'VelocityX'           ; VarName(2) = 'VelocityY'
      VarName(3) = 'Pressure'   ; VarName(4) = 'ViscosityMolecular'  ; VarName(5) = 'ViscosityEddy'
      VarName(6) = 'Density'    ; VarName(7) = 'Mach'
    case(3)
      VarName(0) = 'Temperature'  ; VarName(1) = 'VelocityX'; VarName(2) = 'VelocityY'
      VarName(3) = 'VelocityZ'    ; VarName(4) = 'Pressure' ; VarName(5) = 'ViscosityMolecular'
      VarName(6) = 'ViscosityEddy'; VarName(7) = 'Density'  ; VarName(8) = 'Mach'
    end select
  case(3) ! RANS
    select case(Grid%nDim)
    case(2)
      VarName(0) = 'Temperature'       ; VarName(1) = 'VelocityX'             ; VarName(2) = 'VelocityY'
      VarName(3) = 'Pressure'          ; VarName(4) = 'TurbulentEnergyKinetic'; VarName(5) = 'TurbulentDissipationRate'
      VarName(6) = 'ViscosityMolecular'; VarName(7) = 'ViscosityEddy'         ; VarName(8) = 'BlendingFunction'
      VarName(9) = 'Density'           ; VarName(10) = 'Mach'
    case(3)
      VarName(0) = 'Temperature'             ; VarName(1) = 'VelocityX'         ; VarName(2) = 'VelocityY'
      VarName(3) = 'VelocityZ'               ; VarName(4) = 'Pressure'          ; VarName(5) = 'TurbulentEnergyKinetic'
      VarName(6) = 'TurbulentDissipationRate'; VarName(7) = 'ViscosityMolecular'; VarName(8) = 'ViscosityEddy'
      VarName(9) = 'BlendingFunction'        ; VarName(10) = 'Density'          ; VarName(11) = 'Mach'
    end select
  end select

  ! Write variable name
  call MPI_File_write_all(io,VarName,32*Header(5),MPI_CHARACTER,MPI_STATUS_IGNORE,ier)

  ! Write time data
  call MPI_File_write_all(io,Mixt%Time,1,MPI_REAL8,MPI_STATUS_IGNORE,ier)

  ! Set displacement & view for writing node coordinates
  disp = 5 * intsize + 32 * Header(5) + real8size
  call MPI_File_set_view(io,disp,CoordType,NodeType,'native',MPI_INFO_NULL,ier)

  ! Write node coordinates
  allocate( Coord(Grid%nDim,Grid%nNode) )
  do iNode = 1, Grid%nNode
    iNodeLocal = NodeProcMap(iNode)
    Coord(:,iNode) = Grid%xyz(:,iNodeLocal)
  end do
  call MPI_File_write_all(io,Coord,Grid%nNode,CoordType,MPI_STATUS_IGNORE,ier)
  deallocate( Coord )

  ! Set displacement & view for writing node connectivity of the cell
  disp = disp + Grid%nDim * Grid%nNodeGlobal * real8size
  call MPI_File_set_view(io,disp,NodeCellType,CellType,'native',MPI_INFO_NULL,ier)

  ! Write node connectivity of the cell
  allocate( Nodes(Grid%nNodeCell,Grid%nCell) )
  do iCell = 1, Grid%nCell
    iCellLocal = CellProcMap(iCell)
    do iLocal = 1, Grid%nNodeCell
      iNodeLocal = Grid%c2n(iLocal,iCellLocal)
      Nodes(iLocal,iCell) = Grid%NodeGID(iNodeLocal)
    end do
  end do
  call MPI_File_write_all(io,Nodes,Grid%nCell,NodeCellType,MPI_STATUS_IGNORE,ier)
  deallocate( Nodes )

  ! Set displacement & view for writing node connectivity of the face
  disp = disp + Grid%nNodeCell * Grid%nCellGlobal * intsize
  call MPI_File_set_view(io,disp,NodeFaceType,FaceType,'native',MPI_INFO_NULL,ier)

  ! Write node connectivity of the face
  allocate( Nodes(Grid%nNodeFace,Grid%nFaceBound) )
  do iFace = 1, Grid%nFaceBound
    iFaceLocal = FaceProcMap(iFace)
    do iLocal = 1, Grid%nNodeFace
      iNodeLocal = Grid%f2n(iLocal,iFaceLocal)
      Nodes(iLocal,iFace) = Grid%NodeGID(iNodeLocal)
    end do
  end do
  call MPI_File_write_all(io,Nodes,Grid%nFaceBound,NodeFaceType,MPI_STATUS_IGNORE,ier)
  deallocate( Nodes )

  ! Set displacement & view for writing patch of the face
  disp = disp + Grid%nNodeFace * Grid%nFaceGlobal * intsize
  call MPI_File_set_view(io,disp,MPI_INTEGER,FaceIntType,'native',MPI_INFO_NULL,ier)

  ! Writing patch of the face
  allocate( DataInt(Grid%nFaceBound) )
  do iFace = 1, Grid%nFaceBound
    iFaceLocal = FaceProcMap(iFace)
    DataInt(iFace) = Grid%Patch(iFaceLocal)
  end do
  call MPI_File_write_all(io,DataInt,Grid%nFaceBound,MPI_INTEGER,MPI_STATUS_IGNORE,ier)
  deallocate( DataInt )

  ! Writing variable of the cell
  allocate( DataReal8(Grid%nCell) )

  do iVar = 0, Mixt%nPVar-1

    ! Set displacement & view for writing primitive variable of the cell
    if( iVar == 0 ) disp = disp + Grid%nFaceGlobal * intsize
    if( iVar /= 0 ) disp = disp + Grid%nCellGlobal * real8size
    call MPI_File_set_view(io,disp,MPI_REAL8,CellReal8Type,'native',MPI_INFO_NULL,ier)

    ! Writing primitive variable of the cell
    do iCell = 1, Grid%nCell
      iCellLocal = CellProcMap(iCell)
      DataReal8(iCell) = Mixt%pv(iVar,iCellLocal)
    end do
    call MPI_File_write_all(io,DataReal8,Grid%nCell,MPI_REAL8,MPI_STATUS_IGNORE,ier)

  end do

  ! Set displacement & view for writing density of the cell
  disp = disp + Grid%nCellGlobal * real8size
  call MPI_File_set_view(io,disp,MPI_REAL8,CellReal8Type,'native',MPI_INFO_NULL,ier)

  ! Writing density of the cell
  do iCell = 1, Grid%nCell
    iCellLocal = CellProcMap(iCell)
    DataReal8(iCell) = Mixt%pv(Grid%nDim+1,iCellLocal) / ( Mixt%GasConstant * Mixt%pv(0,iCellLocal) )
  end do
  call MPI_File_write_all(io,DataReal8,Grid%nCell,MPI_REAL8,MPI_STATUS_IGNORE,ier)

  ! Set displacement & view for writing mach of the cell
  disp = disp + Grid%nCellGlobal * real8size
  call MPI_File_set_view(io,disp,MPI_REAL8,CellReal8Type,'native',MPI_INFO_NULL,ier)

  ! Writing mach of the cell
  do iCell = 1, Grid%nCell
    iCellLocal = CellProcMap(iCell)
    DataReal8(iCell) = dsqrt( sum( Mixt%pv(1:Grid%nDim,iCellLocal)**2 ) / ( Mixt%Gamma * Mixt%GasConstant * Mixt%pv(0,iCellLocal) ) )
  end do
  call MPI_File_write_all(io,DataReal8,Grid%nCell,MPI_REAL8,MPI_STATUS_IGNORE,ier)

  deallocate( DataReal8 )

  ! Set displacement & view for writing partition of the cell
  disp = disp + Grid%nCellGlobal * real8size
  call MPI_File_set_view(io,disp,MPI_INTEGER,CellIntType,'native',MPI_INFO_NULL,ier)

  ! Writing partition of the cell
  allocate( DataInt(Grid%nCell) )
  do iCell = 1, Grid%nCell
    iCellLocal = CellProcMap(iCell)
    DataInt(iCell) = Grid%CellPart(iCellLocal)
  end do
  call MPI_File_write_all(io,DataInt,Grid%nCell,MPI_INTEGER,MPI_STATUS_IGNORE,ier)
  deallocate( DataInt )

  ! Close FLUS file
  call MPI_File_close(io,ier)

end subroutine

subroutine Post_Pressure(Grid,Mixt,Conf)
  TYPE(t_Grid), intent(in) :: Grid
  TYPE(t_Mixt), intent(in) :: Mixt
  TYPE(t_Conf), intent(in) :: Conf

  integer :: rank, io, ier
  integer :: iBC, iFace, iCell
  integer :: iList, List1, List2
  real(8) :: p_avg, Area
  real(8) :: p_sum, p_sum_total
  real(8) :: a_sum, a_sum_total
  real(8), parameter :: pi = 4d0*datan(1d0)
  logical, save :: CreateFile = .true.

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Create pressure post file
  if( rank == 0 .and. CreateFile ) then
    CreateFile = .false.
    open(newunit=io,file='./output/pressure_history.csv',status='REPLACE',action='WRITE')
    write(io,'(A)',advance='NO') 'Time(s)'
    do iBC = 1, Conf%nBC
      select case(Conf%BCType(iBC))
      case(12); write(io,'(A)',advance='NO') ', '//trim(Conf%BCName(iBC)(2:))//'(Pa)'
      end select
    end do
    write(io,'(A)') ' '
    close(io)
  end if

  ! Set MPI barrier
  call MPI_Barrier(MPI_COMM_WORLD,ier)

  ! Write pressure at monitor
  if( rank == 0 ) open(newunit=io,file='./output/pressure_history.csv',status='OLD',action='WRITE',access='APPEND')
  if( rank == 0 ) write(io,'(ES)',advance='NO') Mixt%Time
  do iBC = 1, Conf%nBC
    select case(Conf%BCType(iBC))
    case(12) ! Pressure sensor boundary
      ! Set start / end boundary number
      List1 = Grid%ListIndex(iBC)
      List2 = Grid%ListIndex(iBC+1) - 1
      ! Calculate average pressure
      p_sum = 0d0; a_sum = 0d0; p_avg = 0d0
      do iList = List1, List2
        iFace = Grid%List(iList)
        iCell = Grid%f2c(1,iFace)
        if( Conf%Axisymmetric ) then
          Area = Grid%fa(iFace) * 2d0 * pi * Grid%fc(2,iFace)
          p_sum = p_sum + Mixt%pv(Grid%nDim+1,iCell) * Area
          a_sum = a_sum + Area
        else ! not axisymmetric case
          p_sum = p_sum + Mixt%pv(Grid%nDim+1,iCell) * Grid%fa(iFace)
          a_sum = a_sum + Grid%fa(iFace)
        end if
      end do
      call MPI_Allreduce(p_sum,p_sum_total,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call MPI_Allreduce(a_sum,a_sum_total,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      if( a_sum_total > 0 ) p_avg = p_sum_total / a_sum_total
      if( rank == 0 ) write(io,'(A,ES)',advance='NO') ',', p_avg
    end select
  end do
  if( rank == 0 ) write(io,'(A)') ' '
  if( rank == 0 ) close(io)

end subroutine

end module
