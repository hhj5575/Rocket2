module grid_module

use mpi
use config_module
implicit none; private
public t_Grid
public DelGrid
public SetElementInfo
public AllocNodeData
public AllocFaceData
public AllocCellData
public AllocListData
public SetGridFromFile
public SetGridFromSurface
public SetPartGrid
public SetLocalGrid
public SetNodeDataStruc
public SetFaceDataStruc
public SetCellDataStruc
public SetListDataStruc
public SetGridReNumber
public SetGlobalToLocal
public SetConnSendRecv
public SetBoundRecv
public SetFaceMetric
public SetCellMetric
public SetWallDistance
public PreGridVelocity
public SetGridVelocity
public SetFaceVelocity
public UpdateGrid

! SendRecv data type
type t_SendRecv
  ! SendRecv node / face / cell info
  integer :: nNode ! number of nodes for send / recv
  integer, allocatable :: Node(:) ! node list for send / recv
  integer :: nFace ! number of faces for send / recv
  integer, allocatable :: Face(:) ! face list for send / recv
  integer :: nCell ! number of cells for send / recv
  integer, allocatable :: Cell(:) ! cell list for send / recv
  ! SendRecv data related to grid
  real(8), allocatable :: gv(:) ! grid velocity for send / recv
end type

! Grid data type
type t_Grid
  ! Global data
  integer :: nNodeGlobal ! number of global nodes
  integer :: nFaceGlobal ! number of global faces (boundary)
  integer :: nCellGlobal ! number of global cells
  integer :: nNodeGlobalBound ! number of global boundary nodes
  integer, allocatable :: Global_to_Local_Node(:) ! global to local node index
  integer, allocatable :: Global_to_Local_Face(:) ! global to local face index (boundary)
  integer, allocatable :: Global_to_Local_Cell(:) ! global to local cell index
  ! Dimension / element data
  integer :: nDim ! grid dimension
  integer :: nNodeFace ! number of nodes on face
  integer :: nCellFace ! number of cells on face
  integer :: nNodeCell ! number of nodes on cell
  integer :: nFaceCell ! number of faces on cell
  integer, allocatable :: LocalNode(:,:) ! local node index of the face on cell
  ! Node data (order: boundary ~ interior ~ ghost)
  integer :: nNode      ! number of domain nodes
  integer :: nNodeBound ! number of boundary nodes
  integer :: nNodeGhost ! number of ghost nodes
  integer :: nNodeTotal ! number of total nodes
  integer, allocatable :: NodeGID(:) ! node global ID
  integer, allocatable :: NodePart(:) ! node partition ID
  real(8), allocatable :: xyz(:,:) ! node coordinates
  real(8), allocatable :: gv(:,:) ! node velocity
  real(8), allocatable :: IDW(:,:) ! IDW weighting
  real(8), allocatable :: IDW_SUM(:) ! sum of the IDW weighting
  ! Face data (order: boundary ~ interior ~ ghost)
  integer :: nFace      ! number of domain faces
  integer :: nFaceBound ! number of boundary faces
  integer :: nFaceGhost ! number of ghost faces
  integer :: nFaceTotal ! number of total faces
  integer, allocatable :: FaceGID(:) ! face global ID (boundary)
  integer, allocatable :: Patch(:) ! face patch information (boundary)
  integer, allocatable :: f2n(:,:) ! node connectivity of the face
  integer, allocatable :: f2c(:,:) ! neighboring cell of the face
  real(8), allocatable :: fc(:,:) ! face center coordinate
  real(8), allocatable :: fn(:,:) ! face unit normal vector
  real(8), allocatable :: fa(:) ! face area
  real(8), allocatable :: fv(:) ! face velocity
  ! Cell data (order: domain ~ ghost ~ dummy)
  integer :: nCell      ! number of domain cells
  integer :: nCellGhost ! number of ghost cells
  integer :: nCellDummy ! number of dummy cells
  integer :: nCellTotal ! number of total cells
  integer :: nCellWhole ! number of whole cells
  integer, allocatable :: CellGID(:) ! cell global ID
  integer, allocatable :: CellPart(:) ! cell partition ID
  integer, allocatable :: c2n(:,:) ! node connectivity of the cell
  integer, allocatable :: c2c(:,:) ! neighboring cell of the cell
  integer, allocatable :: c2f(:,:) ! neighboring face of the cell
  real(8), allocatable :: cen(:,:) ! cell center coordinate
  real(8), allocatable :: vol(:) ! cell volume
  real(8), allocatable :: WallDist(:) ! wall distance for turbulence model
  ! List data (face for boundary implementation)
  integer, allocatable :: ListIndex(:) ! list index of the face
  integer, allocatable :: List(:) ! face list of the face
  ! Connected domain send/recv grid data
  integer :: nConn ! number of connected domains
  integer, allocatable :: ConnList(:) ! connected domain list
  type(t_SendRecv), allocatable :: ConnSend(:) ! connected domain send grid data
  type(t_SendRecv), allocatable :: ConnRecv(:) ! connected domain recv grid data
  ! Boundary domain recv data
  integer :: nBound ! number of boundary domains
  integer, allocatable :: BoundList(:) ! boundary domain list
  type(t_SendRecv), allocatable :: BoundRecv(:) ! boundary domain recv data
end type

contains

subroutine DelGrid(Grid)
  type(t_Grid), intent(inout) :: Grid

  integer :: iConn, iBound

  ! Global data
  if( allocated( Grid%Global_to_Local_Node ) ) deallocate( Grid%Global_to_Local_Node )
  if( allocated( Grid%Global_to_Local_Face ) ) deallocate( Grid%Global_to_Local_Face )
  if( allocated( Grid%Global_to_Local_Cell ) ) deallocate( Grid%Global_to_Local_Cell )
  ! Dimension / element data
  if( allocated( Grid%LocalNode ) ) deallocate( Grid%LocalNode )
  ! Node data
  if( allocated( Grid%NodeGID ) ) deallocate( Grid%NodeGID )
  if( allocated( Grid%NodePart ) ) deallocate( Grid%NodePart )
  if( allocated( Grid%xyz ) ) deallocate( Grid%xyz )
  if( allocated( Grid%gv ) ) deallocate( Grid%gv )
  if( allocated( Grid%IDW ) ) deallocate( Grid%IDW )
  if( allocated( Grid%IDW_SUM ) ) deallocate( Grid%IDW_SUM )
  ! Face data
  if( allocated( Grid%FaceGID ) ) deallocate( Grid%FaceGID )
  if( allocated( Grid%Patch ) ) deallocate( Grid%Patch )
  if( allocated( Grid%f2n ) ) deallocate( Grid%f2n )
  if( allocated( Grid%f2c ) ) deallocate( Grid%f2c )
  if( allocated( Grid%fc ) ) deallocate( Grid%fc )
  if( allocated( Grid%fn ) ) deallocate( Grid%fn )
  if( allocated( Grid%fa ) ) deallocate( Grid%fa )
  if( allocated( Grid%fv ) ) deallocate( Grid%fv )
  ! Cell data
  if( allocated( Grid%CellGID ) ) deallocate( Grid%CellGID )
  if( allocated( Grid%CellPart ) ) deallocate( Grid%CellPart )
  if( allocated( Grid%c2n ) ) deallocate( Grid%c2n )
  if( allocated( Grid%c2c ) ) deallocate( Grid%c2c )
  if( allocated( Grid%c2f ) ) deallocate( Grid%c2f )
  if( allocated( Grid%cen ) ) deallocate( Grid%cen )
  if( allocated( Grid%vol ) ) deallocate( Grid%vol )
  if( allocated( Grid%WallDist ) ) deallocate( Grid%WallDist )
  ! List data
  if( allocated( Grid%ListIndex ) ) deallocate( Grid%ListIndex )
  if( allocated( Grid%List ) ) deallocate( Grid%List )
  ! Connected domain send/recv data
  if( allocated( Grid%ConnList ) ) deallocate( Grid%ConnList )
  if( allocated( Grid%ConnSend ) ) then
    do iConn = 1, Grid%nConn
      if( allocated( Grid%ConnSend(iConn)%Node ) ) deallocate( Grid%ConnSend(iConn)%Node )
      if( allocated( Grid%ConnSend(iConn)%Face ) ) deallocate( Grid%ConnSend(iConn)%Face )
      if( allocated( Grid%ConnSend(iConn)%Cell ) ) deallocate( Grid%ConnSend(iConn)%Cell )
      if( allocated( Grid%ConnSend(iConn)%gv   ) ) deallocate( Grid%ConnSend(iConn)%gv   )
    end do
    deallocate( Grid%ConnSend )
  end if
  if( allocated( Grid%ConnRecv ) ) then
    do iConn = 1, Grid%nConn
      if( allocated( Grid%ConnRecv(iConn)%Node ) ) deallocate( Grid%ConnRecv(iConn)%Node )
      if( allocated( Grid%ConnRecv(iConn)%Face ) ) deallocate( Grid%ConnRecv(iConn)%Face )
      if( allocated( Grid%ConnRecv(iConn)%Cell ) ) deallocate( Grid%ConnRecv(iConn)%Cell )
      if( allocated( Grid%ConnRecv(iConn)%gv   ) ) deallocate( Grid%ConnRecv(iConn)%gv   )
    end do
    deallocate( Grid%ConnRecv )
  end if
  ! Boundary domain recv data
  if( allocated( Grid%BoundList ) ) deallocate( Grid%BoundList )
  if( allocated( Grid%BoundRecv ) ) then
    do iBound = 1, Grid%nBound
      if( allocated( Grid%BoundRecv(iBound)%Node ) ) deallocate( Grid%BoundRecv(iBound)%Node )
      if( allocated( Grid%BoundRecv(iBound)%Face ) ) deallocate( Grid%BoundRecv(iBound)%Face )
      if( allocated( Grid%BoundRecv(iBound)%Cell ) ) deallocate( Grid%BoundRecv(iBound)%Cell )
      if( allocated( Grid%BoundRecv(iBound)%gv   ) ) deallocate( Grid%BoundRecv(iBound)%gv   )
    end do
    deallocate( Grid%BoundRecv )
  end if

end subroutine

subroutine SetElementInfo(Grid)
  type(t_Grid), intent(inout) :: Grid

  ! We support triangle / tetrahedron for volume element
  ! And line / triangle for surface element
  ! Based on the CGNS node ordering convention
  select case(Grid%nDim)
  case(2) ! triangle
    Grid%nNodeCell = 3; Grid%nFaceCell = 3
    Grid%nNodeFace = 2; Grid%nCellFace = 2
    allocate( Grid%LocalNode(2,3) )
    Grid%LocalNode(1,1) = 1; Grid%LocalNode(2,1) = 2
    Grid%LocalNode(1,2) = 2; Grid%LocalNode(2,2) = 3
    Grid%LocalNode(1,3) = 3; Grid%LocalNode(2,3) = 1
  case(3) ! tetrahedron
    Grid%nNodeCell = 4; Grid%nFaceCell = 4
    Grid%nNodeFace = 3; Grid%nCellFace = 2
    allocate( Grid%LocalNode(3,4) )
    Grid%LocalNode(1,1) = 1; Grid%LocalNode(2,1) = 3; Grid%LocalNode(3,1) = 2
    Grid%LocalNode(1,2) = 1; Grid%LocalNode(2,2) = 2; Grid%LocalNode(3,2) = 4
    Grid%LocalNode(1,3) = 1; Grid%LocalNode(2,3) = 4; Grid%LocalNode(3,3) = 3
    Grid%LocalNode(1,4) = 2; Grid%LocalNode(2,4) = 3; Grid%LocalNode(3,4) = 4
  end select

end subroutine

subroutine AllocNodeData(Grid)
  type(t_Grid), intent(inout) :: Grid

  ! Set number of total nodes
  Grid%nNodeTotal = Grid%nNode + Grid%nNodeGhost

  ! Allocate / initialize node data
  allocate( Grid%NodeGID(Grid%nNodeTotal) )
  allocate( Grid%NodePart(Grid%nNodeTotal) )
  allocate( Grid%xyz(Grid%nDim,Grid%nNodeTotal) )
  allocate( Grid%gv(Grid%nDim,Grid%nNodeTotal) )
  Grid%NodeGID(:) = 0
  Grid%NodePart(:) = 0
  Grid%xyz(:,:) = 0d0
  Grid%gv(:,:) = 0d0

end subroutine

subroutine AllocFaceData(Grid)
  type(t_Grid), intent(inout) :: Grid

  ! Set number of total faces
  Grid%nFaceTotal = Grid%nFace + Grid%nFaceGhost

  ! Allocate / initialize face data
  allocate( Grid%FaceGID(Grid%nFaceBound) )
  allocate( Grid%Patch(Grid%nFaceBound) )
  allocate( Grid%f2n(Grid%nNodeFace,Grid%nFaceTotal) )
  allocate( Grid%f2c(Grid%nCellFace,Grid%nFaceTotal) )
  allocate( Grid%fc(Grid%nDim,Grid%nFaceTotal) )
  allocate( Grid%fn(Grid%nDim,Grid%nFaceTotal) )
  allocate( Grid%fa(Grid%nFaceTotal) )
  allocate( Grid%fv(Grid%nFaceTotal) )
  Grid%FaceGID(:) = 0
  Grid%Patch(:) = 0
  Grid%f2n(:,:) = 0
  Grid%f2c(:,:) = 0
  Grid%fc(:,:) = 0d0
  Grid%fn(:,:) = 0d0
  Grid%fa(:) = 0d0
  Grid%fv(:) = 0d0

end subroutine

subroutine AllocCellData(Grid)
  type(t_Grid), intent(inout) :: Grid

  ! Set number of total/whole cells
  Grid%nCellTotal = Grid%nCell + Grid%nCellGhost
  Grid%nCellWhole = Grid%nCell + Grid%nCellGhost + Grid%nCellDummy

  ! Allocate / initialize cell data
  allocate( Grid%CellGID(Grid%nCellWhole) )
  allocate( Grid%CellPart(Grid%nCellWhole) )
  allocate( Grid%c2n(Grid%nNodeCell,Grid%nCellWhole) )
  allocate( Grid%c2c(Grid%nFaceCell,Grid%nCellWhole) )
  allocate( Grid%c2f(Grid%nFaceCell,Grid%nCellWhole) )
  allocate( Grid%cen(Grid%nDim,Grid%nCellWhole) )
  allocate( Grid%vol(Grid%nCellWhole) )
  allocate( Grid%WallDist(Grid%nCellWhole) )
  Grid%CellGID(:) = 0
  Grid%CellPart(:) = 0
  Grid%c2n(:,:) = 0
  Grid%c2c(:,:) = 0
  Grid%c2f(:,:) = 0
  Grid%cen(:,:) = 0d0
  Grid%vol(:) = 0d0
  Grid%WallDist(:) = 0d0

end subroutine

subroutine AllocListData(Grid,Conf)
  type(t_Grid), intent(inout) :: Grid
  type(t_Conf), intent(in) :: Conf

  allocate( Grid%ListIndex(Conf%nBC+1) )
  allocate( Grid%List(Grid%nFaceBound) )
  Grid%ListIndex(:) = 0
  Grid%List(:) = 0

end subroutine

subroutine SetGridFromFile(Grid,Conf)
  type(t_Grid), intent(out) :: Grid
  type(t_Conf), intent(in) :: Conf

  integer :: error, ier
  integer :: GridFileFormat
  character(32) :: GridFileName
  integer :: iLocal, iNode, iCell, iFace
  integer :: iNewID, iBoundary, iInterior
  logical :: Sort_Flag
  integer, allocatable :: NewID(:)
  real(8), allocatable :: nxyz(:,:)
  logical, allocatable :: isBoundary(:)

  ! Set grid file format / file name
  GridFileFormat = Conf%GridFileFormat
  GridFileName = trim(Conf%GridFileName)

  ! Read grid file
  select case(GridFileFormat)
  case(0); call Read_FLUS_Format(Grid,trim(GridFileName))
  case(1); call Read_GMSH_Format(Grid,trim(GridFileName))
  case(2); call Read_FIDAP_Format(Grid,trim(GridFileName))
  case default
    write(*,*) 'FLUS >> Invalid grid file format'
    call MPI_Abort(MPI_COMM_WORLD,error,ier)
  end select

  ! Print grid info
  write(*,*) 'FLUS >> Grid  dimension:', Grid%nDim
  write(*,*) 'FLUS >> Number of nodes:', Grid%nNode
  write(*,*) 'FLUS >> Number of faces:', Grid%nFace
  write(*,*) 'FLUS >> Number of cells:', Grid%nCell
  write(*,*) 'FLUS >>'

  ! Sort global boundary nodes for FSI interface
  ! Due do this process, we can construct boundary interface
  ! simply collecting boundary nodes / faces of each partition

  ! Allocate boundary checking flag
  allocate( isBoundary(Grid%nNode) )
  isBoundary(:) = .false.

  ! Count number of boundary nodes
  Grid%nNodeBound = 0
  do iFace = 1, Grid%nFace
    do iLocal = 1, Grid%nNodeFace
      iNode = Grid%f2n(iLocal,iFace)
      if( isBoundary(iNode) ) cycle
      Grid%nNodeBound = Grid%nNodeBound + 1
      isBoundary(iNode) = .true.
    end do
  end do

  ! Check boundary node whether node ID
  ! is greater than number of boundary nodes
  Sort_Flag = .false.
  do iNode = 1, Grid%nNodeBound
    if( .not. isBoundary(iNode) ) then
      Sort_Flag = .true.; exit
    end if
  end do

  ! Sort boundary nodes
  if( Sort_Flag ) then

    ! Print progress
    write(*,*) 'FLUS >> Sorting global boundary nodes for FSI interface'
    write(*,*) 'FLUS >>'

    ! Allocate new node ID
    allocate( NewID(Grid%nNode) )
    NewID(:) = 0

    ! Create array of new node ID
    iBoundary = 0; iInterior = Grid%nNodeBound
    do iNode = 1, Grid%nNode
      if( isBoundary(iNode) ) then ! boundary node
        iBoundary = iBoundary + 1
        NewID(iNode) = iBoundary
      else ! interior node
        iInterior = iInterior + 1
        NewID(iNode) = iInterior
      end if
    end do

    ! Change node coordinates based on new ID
    allocate( nxyz(Grid%nDim,Grid%nNode) )
    do iNode = 1, Grid%nNode
      iNewID = NewID(iNode)
      nxyz(:,iNewID) = Grid%xyz(:,iNode)
    end do
    Grid%xyz(:,:) = nxyz(:,:)
    deallocate( nxyz )

    ! Change node connectivity of the cell based on new ID
    do iCell = 1, Grid%nCell
      do iLocal = 1, Grid%nNodeCell
        iNode = Grid%c2n(iLocal,iCell)
        iNewID = NewID(iNode)
        Grid%c2n(iLocal,iCell) = iNewID
      end do
    end do

    ! Change node connectivity of the face based on new ID (boundary)
    do iFace = 1, Grid%nFace
      do iLocal = 1, Grid%nNodeFace
        iNode = Grid%f2n(iLocal,iFace)
        iNewID = NewID(iNode)
        Grid%f2n(iLocal,iFace) = iNewID
      end do
    end do

    ! Deallocate array of the new node ID
    deallocate( NewID )

  end if

  ! Deallocate boundary checking flag
  deallocate( isBoundary )

end subroutine

subroutine Read_FLUS_Format(Grid,FileName)
  type(t_Grid), intent(out) :: Grid
  character(*), intent(in) :: FileName

  integer :: io, error, ier, counter
  integer :: nDim, nNode, nCell, nFace, nVar
  integer :: intsize, real8size
  integer(kind=MPI_OFFSET_KIND) :: disp

  ! Open FLUS file using MPI-IO
  call MPI_File_open(MPI_COMM_SELF,trim(FileName),MPI_MODE_RDONLY,MPI_INFO_NULL,io,ier)

  ! Check file existence and print progress
  if( ier > 0 ) then
    write(*,*) 'FLUS >> Cannot find FLUS grid file: ', trim(FileName)
    call MPI_Abort(MPI_COMM_WORLD,error,ier)
  else
    write(*,*) 'FLUS >> Reading FLUS grid file: ', trim(FileName)
  end if

  ! Read header (grid dimension & number of nodes / cells / faces / variables)
  call MPI_File_read_all(io,nDim, 1,MPI_INTEGER,MPI_STATUS_IGNORE,ier)
  call MPI_File_read_all(io,nNode,1,MPI_INTEGER,MPI_STATUS_IGNORE,ier)
  call MPI_File_read_all(io,nCell,1,MPI_INTEGER,MPI_STATUS_IGNORE,ier)
  call MPI_File_read_all(io,nFace,1,MPI_INTEGER,MPI_STATUS_IGNORE,ier)
  call MPI_File_read_all(io,nVar, 1,MPI_INTEGER,MPI_STATUS_IGNORE,ier)

  ! Set grid dimension and element info
  Grid%nDim = nDim
  call SetElementInfo(Grid)

  ! Set number of nodes and allocate node data
  Grid%nNode      = nNode
  Grid%nNodeBound = 0
  Grid%nNodeGhost = 0
  call AllocNodeData(Grid)

  ! Set number of faces and allocate face data
  Grid%nFace      = nFace
  Grid%nFaceBound = nFace
  Grid%nFaceGhost = 0
  call AllocFaceData(Grid)

  ! Set number of cells and allocate cell data
  Grid%nCell      = nCell
  Grid%nCellGhost = 0
  Grid%nCellDummy = nFace
  call AllocCellData(Grid)

  ! Set type size
  call MPI_Type_size(MPI_INTEGER,intsize,ier)
  call MPI_Type_size(MPI_REAL8,real8size,ier)

  ! Set displacement & view for reading grid information
  disp = 5 * intsize + 32 * nVar + real8size
  call MPI_File_set_view(io,disp,MPI_REAL8,MPI_REAL8,'native',MPI_INFO_NULL,ier)

  ! Read node coordinates
  counter = Grid%nDim * Grid%nNode
  call MPI_File_read_all(io,Grid%xyz,counter,MPI_REAL8,MPI_STATUS_IGNORE,ier)

  ! Read node connectivity of the cell
  counter = Grid%nNodeCell * Grid%nCell
  call MPI_File_read_all(io,Grid%c2n,counter,MPI_INTEGER,MPI_STATUS_IGNORE,ier)

  ! Read node connectivity of the face
  counter = Grid%nNodeFace * Grid%nFace
  call MPI_File_read_all(io,Grid%f2n,counter,MPI_INTEGER,MPI_STATUS_IGNORE,ier)

  ! Read patch of the face
  counter = Grid%nFace
  call MPI_File_read_all(io,Grid%Patch,counter,MPI_INTEGER,MPI_STATUS_IGNORE,ier)

  ! Close FLUS file
  call MPI_File_close(io,ier)

end subroutine

subroutine Read_GMSH_Format(Grid,FileName)
  type(t_Grid), intent(out) :: Grid
  character(*), intent(in) :: FileName

  integer :: error, io, ier
  integer :: nDim, nNode, nCell, nFace
  integer :: nLine, nTris, nTets, nElem
  integer :: iNode, iCell, iFace, iElem
  integer :: NodeNum, ElemNum, ElemType
  integer :: nTags, Tag(5), Nodes(4)
  character(32) :: String
  logical :: isCell

  integer :: iPhysical, nPhysical
  integer, allocatable :: PhysicalDim(:)
  integer, allocatable :: PhysicalTag(:)
  character(32), allocatable :: PhysicalName(:)

  integer :: iBC, nBC
  integer, allocatable :: BCTag(:)
  character(32), allocatable :: BCName(:)

  ! Open GMSH file
  open(newunit=io,file=trim(FileName),status='OLD',action='READ',iostat=ier)

  ! Check file existence and print progress
  if( ier > 0 ) then
    write(*,*) 'FLUS >> Cannot find GMSH grid file: ', trim(FileName)
    call MPI_Abort(MPI_COMM_WORLD,error,ier)
  else
    write(*,*) 'FLUS >> Reading GMSH grid file: ', trim(FileName)
  end if

  ! Read file (first phase)
  do ! Read one string
    read(io,*,iostat=ier) String
    if( ier < 0 ) exit

    ! Read physics info
    if( trim(String) == '$PhysicalNames' ) then
      read(io,*) nPhysical
      allocate( PhysicalDim(nPhysical) )
      allocate( PhysicalTag(nPhysical) )
      allocate( PhysicalName(nPhysical) )
      do iPhysical = 1, nPhysical
        read(io,*) PhysicalDim(iPhysical), PhysicalTag(iPhysical), PhysicalName(iPhysical)
      end do
    end if

    ! read number of node
    if( trim(String) == '$Nodes' ) then
      read(io,*) nNode
    end if

    ! read number of element
    if( trim(String) == '$Elements' ) then
      nLine = 0; nTris = 0; nTets = 0
      read(io,*) nElem
      do iElem = 1, nElem
        read(io,*) ElemNum, ElemType
        select case(ElemType)
        case(1); nLine = nLine + 1
        case(2); nTris = nTris + 1
        case(4); nTets = nTets + 1
        case default
          write(*,*) 'FLUS >> There are unsuppored elements in GMSH grid file'
          call MPI_Abort(MPI_COMM_WORLD,error,ier)
        end select
      end do
    end if
  end do

  ! Set grid dimension
  nDim = maxval(PhysicalDim(:))

  ! Count number of boundary
  nBC = 0
  do iPhysical = 1, nPhysical
    if( PhysicalDim(iPhysical) == nDim-1 ) nBC = nBC + 1
  end do

  ! Sort boundary section
  allocate( BCTag(nBC) )
  allocate( BCName(nBC) )
  iBC = 0
  do iPhysical = 1, nPhysical
    if( PhysicalDim(iPhysical) == nDim-1 ) then
      iBC = iBC + 1
      BCTag(iBC) = PhysicalTag(iPhysical)
      BCName(iBC) = trim(PhysicalName(iPhysical))
    end if
  end do

  ! deallocate physical info
  deallocate( PhysicalDim )
  deallocate( PhysicalTag )
  deallocate( PhysicalName )

  ! Set number of face / cell
  select case(nDim)
  case(2)
    nFace = nLine
    nCell = nTris
  case(3)
    nFace = nTris
    nCell = nTets
  end select

  ! Rewind file
  rewind(io)

  ! Set grid dimension and element info
  Grid%nDim = nDim
  call SetElementInfo(Grid)

  ! Set number of nodes and allocate node data
  Grid%nNode      = nNode
  Grid%nNodeBound = 0
  Grid%nNodeGhost = 0
  call AllocNodeData(Grid)

  ! Set number of faces and allocate face data
  Grid%nFace      = nFace
  Grid%nFaceBound = nFace
  Grid%nFaceGhost = 0
  call AllocFaceData(Grid)

  ! Set number of cells and allocate cell data
  Grid%nCell      = nCell
  Grid%nCellGhost = 0
  Grid%nCellDummy = nFace
  call AllocCellData(Grid)

  ! Read file (second phase)
  do ! read one string
    read(io,*,iostat=ier) String
    if( ier < 0 ) exit

    ! Set node coordinates
    if( trim(String) == '$Nodes' ) then
      read(io,*)
      do iNode = 1, Grid%nNode
        read(io,*) NodeNum, Grid%xyz(:,iNode)
      end do
    end if

    ! read number of element
    if( trim(String) == '$Elements' ) then
      iCell = 0
      iFace = 0
      read(io,*)
      do iElem = 1, nElem

        ! read one string
        read(io,*) ElemNum, ElemType; backspace(io)
        select case(ElemType)
        case(1); read(io,*) ElemNum, ElemType, nTags, Tag(1:nTags), Nodes(1:2) ! line
        case(2); read(io,*) ElemNum, ElemType, nTags, Tag(1:nTags), Nodes(1:3) ! triangle
        case(4); read(io,*) ElemNum, ElemType, nTags, Tag(1:nTags), Nodes(1:4) ! tetrahedron
        end select

        ! Check element type
        isCell = .true.
        do iBC = 1, nBC
          if( Tag(1) == BCTag(iBC) ) then
            ! Add face element
            isCell = .false.
            iFace = iFace + 1
            ! Set patch info / node connectivity of the face
            Grid%Patch(iFace) = Tag(1)
            do iNode = 1, Grid%nNodeFace
              Grid%f2n(iNode,iFace) = Nodes(iNode)
            end do
            exit
          end if
        end do

        if( isCell ) then
          ! Add cell element
          iCell = iCell + 1
          ! Set node connectivity of the cell
          do iNode = 1, Grid%nNodeCell
            Grid%c2n(iNode,iCell) = Nodes(iNode)
          end do
        end if

      end do
    end if

  end do

  ! Close GMSH file
  close(io)

  ! deallocate bc info
  deallocate( BCTag )
  deallocate( BCName )

end subroutine

subroutine Read_FIDAP_Format(Grid,FileName)
  type(t_Grid), intent(out) :: Grid
  character(*), intent(in) :: FileName

  integer :: error, io, ier
  integer :: nDim, NodeNum, ElemNum, ElemType
  integer :: iNode, iFace, iCell, iElem, iGroup
  integer :: nNode, nFace, nCell, nElem, nGroup
  integer :: Patch, Dummy_Int
  character(32) :: Dummy_Char

  ! Open FIDAP file
  open(newunit=io,file=trim(FileName),status='OLD',action='READ',iostat=ier)

  ! Check file existence and print progress
  if( ier > 0 ) then
    write(*,*) 'FLUS >> Cannot find FIDAP grid file: ', trim(FileName)
    call MPI_Abort(MPI_COMM_WORLD,error,ier)
  else
    write(*,*) 'FLUS >> Reading FIDAP grid file: ', trim(FileName)
  end if

  ! Read number of nodes / faces / cells
  nFace = 0; nCell = 0
  read(io,'(4/)')
  read(io,*) nNode, nElem, nGroup, nDim
  read(io,'(<nNode+9>/)')
  do iGroup = 1, nGroup
    read(io,*) Dummy_Char, Dummy_Int, Dummy_Char, nElem, Dummy_Char, Dummy_Char, Dummy_Char, ElemType
    select case(ElemType)
    case(0); nFace = nFace + nElem ! boundary element
    case(2); nCell = nCell + nElem ! triangle element (2-D)
    case(5); nCell = nCell + nElem ! tetrahedron element (3-D)
    end select
    read(io,'(<nElem>/)')
  end do

  ! Rewind file
  rewind(io)

  ! Set grid dimension and element info
  Grid%nDim = nDim
  call SetElementInfo(Grid)

  ! Set number of nodes and allocate node data
  Grid%nNode      = nNode
  Grid%nNodeBound = 0
  Grid%nNodeGhost = 0
  call AllocNodeData(Grid)

  ! Set number of faces and allocate face data
  Grid%nFace      = nFace
  Grid%nFaceBound = nFace
  Grid%nFaceGhost = 0
  call AllocFaceData(Grid)

  ! Set number of cells and allocate cell data
  Grid%nCell      = nCell
  Grid%nCellGhost = 0
  Grid%nCellDummy = nFace
  call AllocCellData(Grid)

  ! Read node coordinates
  read(io,'(12/)')
  do iNode = 1, Grid%nNode
    read(io,*) NodeNum, Grid%xyz(:,iNode)
  end do

  ! Read node connectivity of the face / cell and patch
  iFace = 0; iCell = 0
  read(io,'(2/)')
  do iGroup = 1, nGroup
    read(io,*) Dummy_Char, Dummy_Int, Dummy_Char, nElem, Dummy_Char, Dummy_Char, Dummy_Char, ElemType, Dummy_Char, Patch
    read(io,*)
    select case(ElemType)
    case(0)
      do iElem = 1, nElem
        iFace = iFace + 1
        Grid%Patch(iFace) = Patch
        read(io,*) ElemNum, Grid%f2n(:,iFace)
      end do
    case(2,5)
      do iElem = 1, nElem
        iCell = iCell + 1
        read(io,*) ElemNum, Grid%c2n(:,iCell)
      end do
    end select
  end do

  ! Close FIDAP file
  close(io)

end subroutine

!subroutine Read_CGNS_Format(Grid,FileName)
!  type(t_Grid), intent(out) :: Grid
!  character(*), intent(in) :: FileName
!
!  integer :: counter, iIndex, iMarker
!  integer :: iElem, iNode, iCell
!  logical :: isBoundary
!
!  ! CGNS variables
!  integer :: fn, B, Z, C, S, ier, file_type
!  integer :: nbases, nzones, ngrids, ncoords, nsections
!  character(32) :: basename, zonename, coordname
!  integer :: cell_dim, phys_dim, vertices, cells, boundverts
!  integer :: nbndry, parent_flag, npe
!  cgsize_t :: cgsize(3), zonetype, datatype
!  cgsize_t :: range_min, range_max, start, end
!  real(8), allocatable :: coord_array(:), gridCoords(:,:)
!  cgsize_t :: type, ElementDataSize
!  integer, allocatable :: elemTypeCGNS(:)
!  integer, allocatable :: elemIndex(:)
!  integer, allocatable :: nElems(:)
!  integer, allocatable :: dataSize(:)
!  character(32), allocatable :: sectionName(:)
!  logical, allocatable :: isInternal(:)
!  integer :: interiorElems, boundaryElems, nMarkers
!  character(32) :: ElementSectionName
!  integer :: indexMax, elemMax
!  cgsize_t, allocatable :: connElems(:,:,:)
!  cgsize_t, allocatable :: Elements(:)
!
!  ! Check CGNS file format
!  call cg_is_cgns_f(trim(FileName),file_type,ier)
!  if( ier /= CG_OK ) call cg_error_exit_f()
!
!  ! Open CGNS file
!  call cg_open_f(trim(FileName),CG_MODE_READ,fn,ier)
!  if( ier /= CG_OK ) call cg_error_exit_f()
!
!  ! Read number of bases
!  call cg_nbases_f(fn,nbases,ier)
!  if( ier /= CG_OK ) call cg_error_exit_f()
!
!  ! Check number of bases and set base index
!  if( nbases > 1 ) then ! error if number of bases is greater than 1
!    write(*,*) 'Multi-base CGNS file not supported'; stop
!  else ! set base index to 1
!    B = 1
!  end if
!
!  ! Read base name and dimension info
!  call cg_base_read_f(fn,B,basename,cell_dim,phys_dim,ier)
!  if( ier /= CG_OK ) call cg_error_exit_f()
!
!  ! Read number of zones
!  call cg_nzones_f(fn,B,nzones,ier)
!  if( ier /= CG_OK ) call cg_error_exit_f()
!
!  ! Check number of zones and set zone index
!  if( nzones > 1 ) then ! error if number of zones is greater than 1
!    write(*,*) 'Multi-zone CGNS file not supported'; stop
!  else ! set zone index to 1
!    Z = 1
!  end if
!
!  ! Read zone name and get isize info
!  call cg_zone_read_f(fn,B,Z,zonename,cgsize,ier)
!  if( ier /= CG_OK ) call cg_error_exit_f()
!
!  ! Rename zone size info
!  vertices   = cgsize(1)
!  cells      = cgsize(2)
!  boundverts = cgsize(3)
!
!  ! Read zone type
!  call cg_zone_type_f(fn,B,Z,zonetype,ier)
!  if( ier /= CG_OK ) call cg_error_exit_f()
!
!  ! Check zone type
!  if( zonetype /= Unstructured ) then
!    write(*,*) 'Structured CGNS file not supported'; stop
!  end if
!
!  ! Read number of grids
!  call cg_ngrids_f(fn,B,Z,ngrids,ier)
!  if( ier /= CG_OK ) call cg_error_exit_f()
!
!  ! Check number of grids
!  if( ngrids > 1) then ! error if number of grids is greater than 1
!    write(*,*) 'Multi-grids per zone CGNS file not supported'; stop
!  end if
!
!  ! Read number of grid coordinates dimension
!  call cg_ncoords_f(fn,B,Z,ncoords,ier)
!  if( ier /= CG_OK ) call cg_error_exit_f()
!
!  ! Set vertices range and allocate coordinates array
!  range_min = 1; range_max = vertices
!  allocate( coord_array(range_min:range_max) )
!  allocate( gridCoords(range_min:range_max,ncoords) )
!
!  ! Loop over number of coordinates dimension
!  do C = 1, ncoords
!
!    ! Read coordinates info
!    call cg_coord_info_f(fn,B,Z,C,datatype,coordname,ier)
!    if( ier /= CG_OK ) call cg_error_exit_f()
!
!    ! Read coordinates ( always in double precision )
!    datatype = RealDouble
!    call cg_coord_read_f(fn,B,Z,coordname,datatype,range_min,range_max,coord_array,ier)
!    if( ier /= CG_OK ) call cg_error_exit_f()
!
!    ! Set coordinates
!    gridCoords(:,C) = coord_array(:)
!
!  end do
!
!  ! Read number of sections
!  call cg_nsections_f(fn,B,Z,nsections,ier)
!  if( ier /= CG_OK ) call cg_error_exit_f()
!
!  ! Allocate variables describing each section
!  allocate( elemTypeCGNS(nsections) )
!  allocate( elemIndex(nsections) )
!  allocate( nElems(nsections) )
!  allocate( dataSize(nsections) )
!  allocate( sectionName(nsections) )
!  allocate( isInternal(nsections) )
!
!  ! Loop over number of sections and get section information
!  interiorElems = 0; boundaryElems = 0; nMarkers = 0
!  indexMax = 0; elemMax = 0
!  do S = 1, nsections
!
!    ! Read section info
!    call cg_section_read_f(fn,B,Z,S,ElementSectionName,type,start,end,nbndry,parent_flag,ier)
!    if( ier /= CG_OK ) call cg_error_exit_f()
!
!    ! Read element data size
!    call cg_ElementDataSize_f(fn,B,Z,S,ElementDataSize,ier)
!    if( ier /= CG_OK ) call cg_error_exit_f()
!
!    ! Get number of nodes per element
!    call cg_npe_f(type,npe,ier)
!    if( ier /= CG_OK ) call cg_error_exit_f()
!
!    ! Set section describing variables
!    elemTypeCGNS(S) = type
!    elemIndex(S) = npe
!    nElems(S) = end - start + 1
!    dataSize(S) = ElementDataSize
!    sectionName(S) = trim(ElementSectionName)
!
!    ! In MIXED type, element index for 8 nodes (HEXA_8), plus 1 element type
!    if( elemTypeCGNS(S) == MIXED ) elemIndex(S) = 9
!
!    ! Check interior / boundary section
!    if( cell_dim == 2 ) then
!      select case(elemTypeCGNS(S))
!      case(BAR_2)
!        isInternal(S) = .false.
!        nMarkers = nMarkers + 1
!        boundaryElems = boundaryElems + nElems(S)
!      case(TRI_3) !case(TRI_3,QUAD_4)
!        isInternal(S) = .true.
!        interiorElems = interiorElems + nElems(S)
!      case(MIXED)
!        ! MIXED element treated later
!      case default
!        write(*,*) 'Not supported element type in CGNS file (2-D)'; stop
!      end select
!    else if( cell_dim == 3 ) then
!      select case(elemTypeCGNS(S))
!      case(TRI_3) !case(TRI_3,QUAD_4)
!        isInternal(S) = .false.
!        nMarkers = nMarkers + 1
!        boundaryElems = boundaryElems + nElems(S)
!      case(TETRA_4) !case(TETRA_4,PYRA_5,PENTA_6,HEXA_8)
!        isInternal(S) = .true.
!        interiorElems = interiorElems + nElems(S)
!      case(MIXED)
!        ! MIXED element treated later
!      case default
!        write(*,*) 'Not supported element type in CGNS file (3-D)'; stop
!      end select
!    end if
!
!    ! Set maximum element index and number of elements
!    if( elemIndex(S) > indexMax .or. S == 1 ) indexMax = elemIndex(S)
!    if( nElems(S) > elemMax .or. S == 1 ) elemMax = nElems(S)
!
!  end do
!
!  ! Allocate variables for connectivity info
!  allocate( connElems(indexMax,elemMax,nsections) )
!
!  ! Loop over number of section and get connectivity info
!  do S = 1, nsections
!
!    ! Allocate elements connectivity array
!    allocate( Elements(dataSize(S)) )
!
!    ! Read elements connectivity
!    call cg_elements_read_f(fn,B,Z,S,Elements,Null,ier)
!    if( ier /= CG_OK ) call cg_error_exit_f()
!
!    if( elemTypeCGNS(S) == MIXED ) then
!
!      counter = 0
!      do iElem = 1, nElems(S)
!        counter = counter + 1
!
!        ! Get element type
!        type = Elements(counter)
!
!        ! Check interior / boundary section
!        if( cell_dim == 2 ) then
!          select case(type)
!          case(BAR_2)
!            isBoundary = .true.
!          case(TRI_3) !case(TRI_3,QUAD_4)
!            isBoundary = .false.
!          case default
!            write(*,*) 'Not supported element type in CGNS file (2-D)'; stop
!          end select
!        else if( cell_dim == 3 ) then
!          select case(type)
!          case(TRI_3) !case(TRI_3,QUAD_4)
!            isBoundary = .true.
!          case(TETRA_4) !case(TETRA_4,PYRA_5,PENTA_6,HEXA_8)
!            isBoundary = .false.
!          case default
!            write(*,*) 'Not supported element type in CGNS file (3-D)'; stop
!          end select
!        end if
!
!        ! Get number of nodes per element
!        call cg_npe_f(type,npe)
!        if( ier /= CG_OK ) call cg_error_exit_f()
!
!        connElems(1,iElem,S) = type
!        do iIndex = 1, npe
!          counter = counter + 1
!          connElems(iIndex+1,iElem,S) = Elements(counter)
!        end do
!
!      end do
!
!      ! Check interior / boundary section
!      if( isBoundary ) then
!        isInternal(S) = .false.
!        nMarkers = nMarkers + 1
!        boundaryElems = boundaryElems + nElems(S)
!      else
!        isInternal(S) = .true.
!        interiorElems = interiorElems + nElems(S)
!      end if
!
!    else ! non-mixed type element
!
!      counter = 0
!      do iElem = 1, nElems(S)
!        do iIndex = 1, elemIndex(S)
!          counter = counter + 1
!          connElems(iIndex,iElem,S) = Elements(counter)
!        end do
!      end do
!
!    end if
!
!    deallocate( Elements )
!
!  end do
!
!  ! Close CGNS file
!  call cg_close_f(fn,ier)
!  if( ier /= CG_OK ) call cg_error_exit_f()
!
!  ! Set grid dimension and volume element info according to dimension
!  Grid%nDim = cell_dim
!  call SetElementInfo(Grid)
!
!  ! Set number of global nodes
!  Grid%Global_nNode = vertices
!
!  ! Set number of nodes and allocate node data
!  Grid%nNode = Grid%Global_nNode
!  Grid%nNodeGhost = 0
!  Grid%nNodeTotal = Grid%nNode + Grid%nNodeGhost
!  call AllocNodeData(Grid)
!
!  ! Set node coordinates and global ID
!  do iNode = 1, Grid%nNode
!    Grid%xyz(:,iNode) = gridCoords(iNode,:)
!    Grid%NodeGID(iNode) = iNode
!  end do
!
!  ! Set number of global cells
!  Grid%Global_nCell = interiorElems
!
!  ! Set number of cells and allocate cell data
!  Grid%nCell = Grid%Global_nCell
!  Grid%nCellGhost = 0
!  Grid%nCellTotal = Grid%nCell + Grid%nCellGhost
!  call AllocCellData(Grid)
!
!  ! Set cell connectivity and global ID
!  iCell = 0
!  do S = 1, nsections
!    if( .not. isInternal(S) ) cycle
!    do iElem = 1, nElems(S)
!      iCell = iCell + 1
!      if( elemTypeCGNS(S) == MIXED ) then
!        type = connElems(1,iElem,S)
!        call cg_npe_f(type,npe,ier)
!        if( ier /= CG_OK ) call cg_error_exit_f()
!        do iIndex = 2, npe+1
!          Grid%c2n(iIndex,iCell) = connElems(iIndex,iElem,S)
!        end do
!      else
!        do iIndex = 1, elemIndex(S)
!          Grid%c2n(iIndex,iCell) = connElems(iIndex,iElem,S)
!        end do
!      end if
!      Grid%CellGID(iCell) = iCell
!    end do
!  end do
!
!  ! Set number of markers / maximum number of boundary face
!  ! Allocate boundary face data
!  Grid%nMarker = nMarkers
!  iMarker = 0
!  do S = 1, nsections
!    if( isInternal(S) ) cycle
!    iMarker = iMarker + 1
!    if( nElems(S) > Grid%Max_nFace_Bound .or. iMarker == 1 ) Grid%Max_nFace_Bound = nElems(S)
!  end do
!  call AllocBFaceData(Grid)
!
!  ! Loop over section and set boundary info
!  iMarker = 0
!  do S = 1, nsections
!    if( isInternal(S) ) cycle
!    iMarker = iMarker + 1
!
!    ! Set number of boundary face
!    Grid%nFace_Bound(iMarker) = nElems(S)
!
!    ! Set boundary tag
!    Grid%Tag_Bound(iMarker) = trim(sectionName(S))
!
!    ! Set boundary face connectivity
!    do iElem = 1, nElems(S)
!      if( elemTypeCGNS(S) == MIXED ) then
!        type = connElems(1,iElem,S)
!        call cg_npe_f(type,npe,ier)
!        if( ier /= CG_OK ) call cg_error_exit_f()
!        do iIndex = 2, npe+1
!          Grid%f2n_Bound(iIndex,iElem,iMarker) = connElems(iIndex,iElem,S)
!        end do
!      else
!        do iIndex = 1, elemIndex(S)
!          Grid%f2n_Bound(iIndex,iElem,iMarker) = connElems(iIndex,iElem,S)
!        end do
!      end if
!    end do
!  end do
!
!  ! Deallocate local variables
!  deallocate( coord_array )
!  deallocate( gridCoords )
!  deallocate( elemTypeCGNS )
!  deallocate( elemIndex )
!  deallocate( nElems )
!  deallocate( dataSize )
!  deallocate( sectionName )
!  deallocate( isInternal )
!  deallocate( connElems )
!
!end subroutine

subroutine SetGridFromSurface(Grid,nDim,nBNode,nBFace,Bxyz,Bf2n,BPatch)
  type(t_Grid), intent(out) :: Grid
  integer, intent(in) :: nDim
  integer, intent(in) :: nBNode
  integer, intent(in) :: nBFace
  real(8), intent(in) :: Bxyz(:,:)
  integer, intent(in) :: Bf2n(:,:)
  integer, intent(in) :: BPatch(:)

  integer :: error, ier
  integer :: iLocal, iNode, iFace
  logical, allocatable :: isBoundary(:)

  ! Create grid using external program
  select case(nDim)
  case(2); call MakeGrid_Trigen(Grid,nBNode,nBFace,Bxyz,Bf2n,BPatch)
  case(3); call MakeGrid_Tetgen(Grid,nBNode,nBFace,Bxyz,Bf2n,BPatch)
  end select

  ! Print grid info
  write(*,*) 'FLUS >> Grid  dimension:', Grid%nDim
  write(*,*) 'FLUS >> Number of nodes:', Grid%nNode
  write(*,*) 'FLUS >> Number of faces:', Grid%nFace
  write(*,*) 'FLUS >> Number of cells:', Grid%nCell
  write(*,*) 'FLUS >>'

  ! Just check global boundary nodes
  ! note: if error occur, modify this routine

  ! Allocate boundary checking flag
  allocate( isBoundary(Grid%nNode) )
  isBoundary(:) = .false.

  ! Count number of boundary nodes
  Grid%nNodeBound = 0
  do iFace = 1, Grid%nFace
    do iLocal = 1, Grid%nNodeFace
      iNode = Grid%f2n(iLocal,iFace)
      if( isBoundary(iNode) ) cycle
      Grid%nNodeBound = Grid%nNodeBound + 1
      isBoundary(iNode) = .true.
    end do
  end do

  ! Check boundary node whether node ID
  ! is greater than number of boundary nodes
  do iNode = 1, Grid%nNodeBound
    if( .not. isBoundary(iNode) ) then
      write(*,*) 'FLUS >> Sorting flag is true after remesh, check!!!'
      call MPI_Abort(MPI_COMM_WORLD,error,ier)
    end if
  end do

  ! Deallocate boundary checking flag
  deallocate( isBoundary )

end subroutine

subroutine MakeGrid_Trigen(Grid,nBNode,nBFace,Bxyz,Bf2n,BPatch)
  type(t_Grid), intent(out) :: Grid
  integer, intent(in) :: nBNode
  integer, intent(in) :: nBFace
  real(8), intent(in) :: Bxyz(:,:)
  integer, intent(in) :: Bf2n(:,:)
  integer, intent(in) :: BPatch(:)

  integer :: iNode, nNode, nNodeEst
  integer :: iFace, nFace, nFaceEst
  integer :: iCell, nCell, nCellEst
  real(8), allocatable :: xyz(:,:)
  integer, allocatable :: Patch(:)
  integer, allocatable :: f2n(:,:)
  integer, allocatable :: c2n(:,:)

  ! Print progress
  write(*,*) 'FLUS >> Creating grid from surface using TRIGEN'

  ! Estimate new grid node / cell
  nNodeEst = nBNode * 50
  nFaceEst = nBFace * 50
  nCellEst = nBNode * 50
  allocate( xyz(2,nNodeEst) )
  allocate( Patch(nFaceEst) )
  allocate( f2n(2,nFaceEst) )
  allocate( c2n(3,nCellEst) )

  ! Call trigen library
  call Trigen(nBNode,nBFace,Bxyz,Bf2n,BPatch,nNode,nCell,nFace,xyz,c2n,f2n,Patch)

  ! Set grid dimension and element info
  Grid%nDim = 2
  call SetElementInfo(Grid)

  ! Set number of nodes and allocate node data
  Grid%nNode      = nNode
  Grid%nNodeBound = 0
  Grid%nNodeGhost = 0
  call AllocNodeData(Grid)

  ! Set number of faces and allocate face data
  Grid%nFace      = nFace
  Grid%nFaceBound = nFace
  Grid%nFaceGhost = 0
  call AllocFaceData(Grid)

  ! Set number of cells and allocate cell data
  Grid%nCell      = nCell
  Grid%nCellGhost = 0
  Grid%nCellDummy = nFace
  call AllocCellData(Grid)

  ! Set node info
  do iNode = 1, Grid%nNode
    Grid%xyz(:,iNode) = xyz(:,iNode)
  end do

  ! Set face info
  do iFace = 1, Grid%nFace
    Grid%Patch(iFace) = Patch(iFace)
    Grid%f2n(:,iFace) = f2n(:,iFace)
  end do

  ! Set cell info
  do iCell = 1, Grid%nCell
    Grid%c2n(:,iCell) = c2n(:,iCell)
  end do

  ! deallocate local variable
  deallocate( xyz )
  deallocate( Patch )
  deallocate( f2n )
  deallocate( c2n )

end subroutine

subroutine MakeGrid_Tetgen(Grid,nBNode,nBFace,Bxyz,Bf2n,BPatch)
  use, intrinsic :: iso_c_binding
  type(t_Grid), intent(out) :: Grid
  integer(kind=c_int), intent(in) :: nBNode
  integer(kind=c_int), intent(in) :: nBFace
  real(kind=c_double), intent(in) :: Bxyz(:,:)
  integer(kind=c_int), intent(in) :: Bf2n(:,:)
  integer(kind=c_int), intent(in) :: BPatch(:)

  integer :: iNode, nNode, nNodeEst
  integer :: iFace, nFace, nFaceEst
  integer :: iCell, nCell, nCellEst
  real(8), allocatable :: xyz(:,:)
  integer, allocatable :: Patch(:)
  integer, allocatable :: f2n(:,:)
  integer, allocatable :: c2n(:,:)

  ! Tetgen option
  character(kind=c_char, len=32) :: opt='pYYq1.01O9/7'

  ! Interface for Fortran / C++ communication
  interface
    subroutine tetgen(opt,nBNode,nBFace,Bxyz,Bf2n,BPatch,nNode,nCell,nFace,xyz,c2n,f2n,Patch) bind(c)
      import :: c_int, c_double, c_char
      character(kind=c_char) :: opt(*)
      integer(kind=c_int) :: nBNode
      integer(kind=c_int) :: nBFace
      real(kind=c_double) :: Bxyz(3,nBNode)
      integer(kind=c_int) :: Bf2n(3,nBFace)
      integer(kind=c_int) :: BPatch(nBFace)
      integer(kind=c_int) :: nNode
      integer(kind=c_int) :: nCell
      integer(kind=c_int) :: nFace
      real(kind=c_double) :: xyz(3,*)
      integer(kind=c_int) :: c2n(4,*)
      integer(kind=c_int) :: f2n(3,*)
      integer(kind=c_int) :: Patch(*)
    end subroutine
  end interface

  ! Print progress
  write(*,*) 'FLUS >> Creating grid from surface using TETGEN'

  ! Estimate new grid node / cell
  nNodeEst = nBNode * 50
  nFaceEst = nBFace * 50
  nCellEst = nBNode * 50
  allocate( xyz(3,nNodeEst) )
  allocate( Patch(nFaceEst) )
  allocate( f2n(3,nFaceEst) )
  allocate( c2n(4,nCellEst) )

  ! Call tetgen library
  call Tetgen(trim(opt)//c_null_char,nBNode,nBFace,Bxyz,Bf2n,BPatch,nNode,nCell,nFace,xyz,c2n,f2n,Patch)

  ! Set grid dimension and element info
  Grid%nDim = 3
  call SetElementInfo(Grid)

  ! Set number of nodes and allocate node data
  Grid%nNode      = nNode
  Grid%nNodeBound = 0
  Grid%nNodeGhost = 0
  call AllocNodeData(Grid)

  ! Set number of faces and allocate face data
  Grid%nFace      = nFace
  Grid%nFaceBound = nFace
  Grid%nFaceGhost = 0
  call AllocFaceData(Grid)

  ! Set number of cells and allocate cell data
  Grid%nCell      = nCell
  Grid%nCellGhost = 0
  Grid%nCellDummy = nFace
  call AllocCellData(Grid)

  ! Set node info
  do iNode = 1, Grid%nNode
    Grid%xyz(:,iNode) = xyz(:,iNode)
  end do

  ! Set face info
  do iFace = 1, Grid%nFace
    Grid%Patch(iFace) = Patch(iFace)
    Grid%f2n(:,iFace) = f2n(:,iFace)
  end do

  ! Set cell info
  do iCell = 1, Grid%nCell
    Grid%c2n(:,iCell) = c2n(:,iCell)
  end do

  ! deallocate local variable
  deallocate( xyz )
  deallocate( Patch )
  deallocate( f2n )
  deallocate( c2n )

end subroutine

subroutine SetPartGrid(Grid)
  type(t_Grid), intent(inout) :: Grid

  integer :: size, ier
  integer :: iLocal, iNode, iCell

  ! Metis variables
  integer :: ne, nn, etype, numflag, nparts, edgecut
  integer, allocatable :: elmnts(:), epart(:), npart(:)

  ! Get number of cores
  call MPI_Comm_size(MPI_COMM_WORLD,size,ier)

  ! Set METIS variables
  nn = Grid%nNode
  ne = Grid%nCell
  allocate( elmnts(Grid%nNodeCell*ne) )
  if( Grid%nDim == 2 ) etype = 1 ! triangle
  if( Grid%nDim == 3 ) etype = 2 ! tetrahedron
  numflag = 1
  nparts = size
  allocate( epart(ne) )
  allocate( npart(nn) )

  ! Set partition info if number of parts is greater than 1
  if( nparts > 1 ) then

    ! Print progress
    write(*,*) 'FLUS >> Paritioning global grid using METIS 4.0.3'

    ! Set node connectivity of cell
    do iCell = 1, Grid%nCell
      do iLocal = 1, Grid%nNodeCell
        iNode = Grid%nNodeCell*(iCell-1)+iLocal
        elmnts(iNode) = Grid%c2n(iLocal,iCell)
      end do
    end do

    ! Grid partitioning using METIS 4.0.3
    call METIS_PartMeshDual(ne,nn,elmnts,etype,numflag,nparts,edgecut,epart,npart)

    ! Set node / cell partition info
    do iNode = 1, Grid%nNode
      Grid%NodePart(iNode) = npart(iNode)
    end do
    do iCell = 1, Grid%nCell
      Grid%CellPart(iCell) = epart(iCell)
    end do

  else

    ! Set node / cell partition to 1
    do iNode = 1, Grid%nNode
      Grid%NodePart(iNode) = 1
    end do
    do iCell = 1, Grid%nCell
      Grid%CellPart(iCell) = 1
    end do

  end if

  deallocate( elmnts )
  deallocate( epart )
  deallocate( npart )

end subroutine

subroutine SetLocalGrid(Grid,Global)
  type(t_Grid), intent(out) :: Grid
  type(t_Grid), intent(inout) :: Global

  integer :: rank, size, srcs, dest, ier, counter
  integer :: iLocal, iLocal_Part, iBuffer, iDim
  integer :: iNode, iNodeTotal, iNodeDomain, iNodeGhost
  integer :: iCell, iCellTotal, iCellDomain, iCellGhost
  integer :: iFace, iFaceTotal
  integer :: iadj, xadj_max
  integer, allocatable :: xadj(:), adjncy(:,:), adjncy_temp(:,:)
  integer :: iDomain, jDomain, nDomain, DomainListSize
  integer, allocatable :: DomainList(:)
  logical :: CheckDomain
  integer :: iNode_Part, Max_nNode_Part
  integer :: iFace_Part, Max_nFace_Part
  integer, allocatable :: nNode_Part(:), Node_Part(:,:)
  integer, allocatable :: nFace_Part(:), Face_Part(:,:)
  integer, allocatable :: Node_Flag(:)

  ! MPI comm const / array
  integer :: Buffer_Send_nDim
  integer :: Buffer_Recv_nDim
  integer :: Buffer_Send_nNodeTotal, Buffer_Send_nNodeDomain, Buffer_Send_nNodeGhost
  integer :: Buffer_Recv_nNodeTotal, Buffer_Recv_nNodeDomain, Buffer_Recv_nNodeGhost
  integer :: Buffer_Send_nCellTotal, Buffer_Send_nCellDomain, Buffer_Send_nCellGhost
  integer :: Buffer_Recv_nCellTotal, Buffer_Recv_nCellDomain, Buffer_Recv_nCellGhost
  integer :: Buffer_Send_nFaceTotal
  integer :: Buffer_Recv_nFaceTotal
  integer, allocatable :: Buffer_Send_NodeGID(:), Buffer_Send_NodePart(:)
  integer, allocatable :: Buffer_Recv_NodeGID(:), Buffer_Recv_NodePart(:)
  real(8), allocatable :: Buffer_Send_xyz(:)
  real(8), allocatable :: Buffer_Recv_xyz(:)
  integer, allocatable :: Buffer_Send_CellGID(:), Buffer_Send_CellPart(:), Buffer_Send_c2n(:)
  integer, allocatable :: Buffer_Recv_CellGID(:), Buffer_Recv_CellPart(:), Buffer_Recv_c2n(:)
  integer, allocatable :: Buffer_Send_FaceGID(:), Buffer_Send_Patch(:), Buffer_Send_f2n(:)
  integer, allocatable :: Buffer_Recv_FaceGID(:), Buffer_Recv_Patch(:), Buffer_Recv_f2n(:)

  ! MPI comm request / status
  integer :: Send_Req(9), Send_Stat(MPI_STATUS_SIZE,9)
  integer :: Recv_Req(9), Recv_Stat(MPI_STATUS_SIZE,9)

  ! Get rank of each core / number of cores
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)
  call MPI_Comm_size(MPI_COMM_WORLD,size,ier)

  ! Set number of domains
  nDomain = size

  ! Create data structure needed for partitioning grid
  if( rank == 0 ) then

    ! 1) Set adjacent cells of the node

    ! Allocate / initialize xadj / adjncy
    ! Approximate initial maximum xadj = 100
    xadj_max = 100
    allocate( xadj(Global%nNode) )
    allocate( adjncy(xadj_max,Global%nNode) )
    xadj(:) = 0; adjncy(:,:) = 0

    ! Set xadj / adjncy data using domain cell
    do iCell = 1, Global%nCell
      do iLocal = 1, Global%nNodeCell
        iNode = Global%c2n(iLocal,iCell)
        xadj(iNode) = xadj(iNode) + 1
        if( xadj(iNode) > xadj_max ) then
          write(*,*) 'FLUS >> Number of adjacent cells is greater than upper bound'
          write(*,*) 'FLUS >> Increase upper bound of the number of adjacent cells (+100)'
          ! Create temporary adjncy
          allocate( adjncy_temp(xadj_max,Global%nNode) )
          adjncy_temp(:,:) = adjncy(:,:)
          ! Reallocate / initialize adjncy
          xadj_max = xadj_max + 100 ! Increase maximum xadj
          deallocate( adjncy )
          allocate( adjncy(xadj_max,Global%nNode) )
          adjncy(:,:) = 0
          ! Copy data and deallocate temporary adjncy
          adjncy(1:xadj_max-100,:) = adjncy_temp(:,:)
          deallocate( adjncy_temp )
        end if
        adjncy(xadj(iNode),iNode) = iCell
      end do
    end do

    ! 2) Set node data in partition list for speed-up
    ! Partition include node if one of the adjacent cell is included in partition

    ! Allocate domain list
    allocate( DomainList(nDomain) )

    ! Divide the nodes in partition list
    allocate( nNode_Part(nDomain) )
    nNode_Part(:) = 0

    ! Count node list based on adjacent cell parition ID
    do iNode = 1, Global%nNode
      DomainListSize = 0
      DomainList(:) = 0
      do iadj = 1, xadj(iNode)
        iCell = adjncy(iadj,iNode)
        iDomain = Global%CellPart(iCell)

        CheckDomain = .true.
        do jDomain = 1, DomainListSize
          if( DomainList(jDomain) == iDomain ) then
            CheckDomain = .false.; exit
          end if
        end do

        if( CheckDomain ) then
          DomainListSize = DomainListSize + 1
          DomainList(DomainListSize) = iDomain
          nNode_Part(iDomain) = nNode_Part(iDomain) + 1
        end if

      end do
    end do

    ! Find maximum number of nodes per partition
    Max_nNode_Part = 0
    do iDomain = 1, nDomain
      if( nNode_Part(iDomain) > Max_nNode_Part ) Max_nNode_Part = nNode_Part(iDomain)
    end do

    ! Allocate the node partition array
    allocate( Node_Part(Max_nNode_Part,nDomain) )
    nNode_Part(:) = 0; Node_Part(:,:) = 0

    ! Create node list based on partition ID
    do iNode = 1, Global%nNode
      DomainListSize = 0
      DomainList(:) = 0
      do iadj = 1, xadj(iNode)
        iCell = adjncy(iadj,iNode)
        iDomain = Global%CellPart(iCell)

        CheckDomain = .true.
        do jDomain = 1, DomainListSize
          if( DomainList(jDomain) == iDomain ) then
            CheckDomain = .false.; exit
          end if
        end do

        if( CheckDomain ) then
          DomainListSize = DomainListSize + 1
          DomainList(DomainListSize) = iDomain
          nNode_Part(iDomain) = nNode_Part(iDomain) + 1
          Node_Part(nNode_Part(iDomain),iDomain) = iNode
        end if

      end do
    end do

    ! Deallocate domain list
    deallocate( DomainList )

    ! 3) Set face data in partition list for speed-up

    ! Create flag array for checking node
    allocate( Node_Flag(Global%nNode) )
    Node_Flag(:) = 0

    ! Divide the faces in partition list
    allocate( nFace_Part(nDomain) )
    nFace_Part(:) = 0

    ! Count face list based on parent cell partition ID
    do iFace = 1, Global%nFace
      ! Select one arbitrary node and find parent cell
      iNode = Global%f2n(1,iFace)
      do iadj = 1, xadj(iNode)
        iCell = adjncy(iadj,iNode)
        Node_Flag(Global%c2n(:,iCell)) = 1
        counter = sum(Node_Flag(Global%f2n(:,iFace)))
        Node_Flag(Global%c2n(:,iCell)) = 0
        if( counter == Global%nNodeFace ) then
          iDomain = Global%CellPart(iCell)
          nFace_Part(iDomain) = nFace_Part(iDomain) + 1
          exit
        end if
      end do
    end do

    ! Find maximum number of faces per partition
    Max_nFace_Part = 0
    do iDomain = 1, nDomain
      if( nFace_Part(iDomain) > Max_nFace_Part ) Max_nFace_Part = nFace_Part(iDomain)
    end do

    ! Allocate the face partition array
    allocate( Face_Part(Max_nFace_Part,nDomain) )
    nFace_Part(:) = 0; Face_Part(:,:) = 0

    ! Set face list based on parent cell partition ID
    do iFace = 1, Global%nFace
      ! Select one arbitrary node and find parent cell
      iNode = Global%f2n(1,iFace)
      do iadj = 1, xadj(iNode)
        iCell = adjncy(iadj,iNode)
        Node_Flag(Global%c2n(:,iCell)) = 1
        counter = sum(Node_Flag(Global%f2n(:,iFace)))
        Node_Flag(Global%c2n(:,iCell)) = 0
        if( counter == Global%nNodeFace ) then
          iDomain = Global%CellPart(iCell)
          nFace_Part(iDomain) = nFace_Part(iDomain) + 1
          Face_Part(nFace_Part(iDomain),iDomain) = iFace
          exit
        end if
      end do
    end do

    ! Deallocate node flag
    deallocate( Node_Flag )

    ! 4) Allocate global to local node / cell / face mapping array

    ! Allocate global to local node / cell / face mapping
    allocate( Global%Global_to_Local_Node(Global%nNode) )
    allocate( Global%Global_to_Local_Cell(Global%nCell) )
    allocate( Global%Global_to_Local_Face(Global%nFace) )

  end if

  ! Divide global grid and distribute local grid to each core
  do iDomain = 1, nDomain

    if( rank == 0 ) then ! master core

      ! Set buffer send dimension
      Buffer_Send_nDim = Global%nDim

      ! Initialize global to local mapping
      Global%Global_to_Local_Cell(:) = 0
      Global%Global_to_Local_Node(:) = 0

      ! Count number of buffer send nodes / cells
      Buffer_Send_nCellTotal = 0; Buffer_Send_nCellDomain = 0; Buffer_Send_nCellGhost = 0
      Buffer_Send_nNodeTotal = 0; Buffer_Send_nNodeDomain = 0; Buffer_Send_nNodeGhost = 0
      do iLocal_Part = 1, nNode_Part(iDomain)
        iNode_Part = Node_Part(iLocal_Part,iDomain)

        do iadj = 1, xadj(iNode_Part)
          iCell = adjncy(iadj,iNode_Part)

          if( Global%Global_to_Local_Cell(iCell) == 0 ) then
            Global%Global_to_Local_Cell(iCell) = 1
            Buffer_Send_nCellTotal = Buffer_Send_nCellTotal + 1
            if( Global%CellPart(iCell) /= iDomain ) then
              Buffer_Send_nCellGhost = Buffer_Send_nCellGhost + 1
            else
              Buffer_Send_nCellDomain = Buffer_Send_nCellDomain + 1
            end if
          end if

          do iLocal = 1, Global%nNodeCell
            iNode = Global%c2n(iLocal,iCell)

            if( Global%Global_to_Local_Node(iNode) == 0 ) then
              Global%Global_to_Local_Node(iNode) = 1
              Buffer_Send_nNodeTotal = Buffer_Send_nNodeTotal + 1
              if( Global%NodePart(iNode) /= iDomain ) then
                Buffer_Send_nNodeGhost = Buffer_Send_nNodeGhost + 1
              else
                Buffer_Send_nNodeDomain = Buffer_Send_nNodeDomain + 1
              end if
            end if

          end do

        end do

      end do

      ! Count number of buffer send faces
      ! (Ghost faces does not exist, so partition faces are same as total faces)
      Buffer_Send_nFaceTotal = nFace_Part(iDomain)

      ! Allocate buffer send node info
      allocate( Buffer_Send_NodeGID(Buffer_Send_nNodeTotal) )
      allocate( Buffer_Send_NodePart(Buffer_Send_nNodeTotal) )
      allocate( Buffer_Send_xyz(Buffer_Send_nNodeTotal*Buffer_Send_nDim) )

      ! Allocate buffer send cell info
      allocate( Buffer_Send_CellGID(Buffer_Send_nCellTotal) )
      allocate( Buffer_Send_CellPart(Buffer_Send_nCellTotal) )
      allocate( Buffer_Send_c2n(Buffer_Send_nCellTotal*Global%nNodeCell) )

      ! Allocate buffer send face info
      allocate( Buffer_Send_FaceGID(Buffer_Send_nFaceTotal) )
      allocate( Buffer_Send_Patch(Buffer_Send_nFaceTotal) )
      allocate( Buffer_Send_f2n(Buffer_Send_nFaceTotal*Global%nNodeFace) )

      if( iDomain-1 /= 0 ) then ! iDomain is different from the master rank 0
        dest = iDomain-1
        ! Non-blocking send comm
        call MPI_Isend(Buffer_Send_nDim,       1,MPI_INTEGER,dest,1,MPI_COMM_WORLD,Send_Req(1),ier)
        call MPI_Isend(Buffer_Send_nNodeTotal, 1,MPI_INTEGER,dest,2,MPI_COMM_WORLD,Send_Req(2),ier)
        call MPI_Isend(Buffer_Send_nNodeDomain,1,MPI_INTEGER,dest,3,MPI_COMM_WORLD,Send_Req(3),ier)
        call MPI_Isend(Buffer_Send_nNodeGhost, 1,MPI_INTEGER,dest,4,MPI_COMM_WORLD,Send_Req(4),ier)
        call MPI_Isend(Buffer_Send_nCellTotal, 1,MPI_INTEGER,dest,5,MPI_COMM_WORLD,Send_Req(5),ier)
        call MPI_Isend(Buffer_Send_nCellDomain,1,MPI_INTEGER,dest,6,MPI_COMM_WORLD,Send_Req(6),ier)
        call MPI_Isend(Buffer_Send_nCellGhost, 1,MPI_INTEGER,dest,7,MPI_COMM_WORLD,Send_Req(7),ier)
        call MPI_Isend(Buffer_Send_nFaceTotal, 1,MPI_INTEGER,dest,8,MPI_COMM_WORLD,Send_Req(8),ier)
        ! Wait for the set of non-blocking comm
        call MPI_Waitall(8,Send_Req,Send_Stat,ier)
      else ! iDomain is master rank. so, simply copy buffer send
        Buffer_Recv_nDim        = Buffer_Send_nDim
        Buffer_Recv_nNodeTotal  = Buffer_Send_nNodeTotal
        Buffer_Recv_nNodeDomain = Buffer_Send_nNodeDomain
        Buffer_Recv_nNodeGhost  = Buffer_Send_nNodeGhost
        Buffer_Recv_nCellTotal  = Buffer_Send_nCellTotal
        Buffer_Recv_nCellDomain = Buffer_Send_nCellDomain
        Buffer_Recv_nCellGhost  = Buffer_Send_nCellGhost
        Buffer_Recv_nFaceTotal  = Buffer_Send_nFaceTotal
      end if

    end if

    if( rank == iDomain-1 ) then

      if( rank /= 0 ) then ! receive the size of buffers
        srcs = 0 ! master core
        ! Non-blocking receive comm
        call MPI_Irecv(Buffer_Recv_nDim,       1,MPI_INTEGER,srcs,1,MPI_COMM_WORLD,Recv_Req(1),ier)
        call MPI_Irecv(Buffer_Recv_nNodeTotal, 1,MPI_INTEGER,srcs,2,MPI_COMM_WORLD,Recv_Req(2),ier)
        call MPI_Irecv(Buffer_Recv_nNodeDomain,1,MPI_INTEGER,srcs,3,MPI_COMM_WORLD,Recv_Req(3),ier)
        call MPI_Irecv(Buffer_Recv_nNodeGhost, 1,MPI_INTEGER,srcs,4,MPI_COMM_WORLD,Recv_Req(4),ier)
        call MPI_Irecv(Buffer_Recv_nCellTotal, 1,MPI_INTEGER,srcs,5,MPI_COMM_WORLD,Recv_Req(5),ier)
        call MPI_Irecv(Buffer_Recv_nCellDomain,1,MPI_INTEGER,srcs,6,MPI_COMM_WORLD,Recv_Req(6),ier)
        call MPI_Irecv(Buffer_Recv_nCellGhost, 1,MPI_INTEGER,srcs,7,MPI_COMM_WORLD,Recv_Req(7),ier)
        call MPI_Irecv(Buffer_Recv_nFaceTotal, 1,MPI_INTEGER,srcs,8,MPI_COMM_WORLD,Recv_Req(8),ier)
        ! Wait for thie set of non-blocking comm
        call MPI_Waitall(8,Recv_Req,Recv_Stat,ier)
      end if

      ! Set grid dimension and element info
      ! This routine must be called here for allocate buffer recv array
      Grid%nDim = Buffer_Recv_nDim
      call SetElementInfo(Grid)

      ! Allocate buffer recv node info
      allocate( Buffer_Recv_NodeGID(Buffer_Recv_nNodeTotal) )
      allocate( Buffer_Recv_NodePart(Buffer_Recv_nNodeTotal) )
      allocate( Buffer_Recv_xyz(Buffer_Recv_nNodeTotal*Buffer_Recv_nDim) )

      ! Allocate buffer recv cell info
      allocate( Buffer_Recv_CellGID(Buffer_Recv_nCellTotal) )
      allocate( Buffer_Recv_CellPart(Buffer_Recv_nCellTotal) )
      allocate( Buffer_Recv_c2n(Buffer_Recv_nCellTotal*Grid%nNodeCell) )

      ! Allocate buffer recv boundary face info
      allocate( Buffer_Recv_FaceGID(Buffer_Recv_nFaceTotal) )
      allocate( Buffer_Recv_Patch(Buffer_Recv_nFaceTotal) )
      allocate( Buffer_Recv_f2n(Buffer_Recv_nFaceTotal*Grid%nNodeFace) )

    end if

    if( rank == 0 ) then ! master core

      ! Initialize global to local mapping
      Global%Global_to_Local_Node(:) = 0
      Global%Global_to_Local_Cell(:) = 0
      Global%Global_to_Local_Face(:) = 0

      ! Set buffer send node / cell info
      iCellDomain = 0; iCellGhost = Buffer_Send_nCellDomain
      iNodeDomain = 0; iNodeGhost = Buffer_Send_nNodeDomain
      do iLocal_Part = 1, nNode_Part(iDomain)
        iNode_Part = Node_Part(iLocal_Part,iDomain)

        do iadj = 1, xadj(iNode_Part)
          iCell = adjncy(iadj,iNode_Part)

          ! Check and set buffer send cell info
          if( Global%Global_to_Local_Cell(iCell) == 0 ) then
            if( Global%CellPart(iCell) == iDomain ) then
              iCellDomain = iCellDomain + 1
              iCellTotal = iCellDomain
            else
              iCellGhost = iCellGhost + 1
              iCellTotal = iCellGhost
            end if
            Global%Global_to_Local_Cell(iCell) = iCellTotal
            Buffer_Send_CellGID(iCellTotal) = iCell
            Buffer_Send_CellPart(iCellTotal) = Global%CellPart(iCell)

            ! Check and set buffer send node info
            do iLocal = 1, Global%nNodeCell
              iNode = Global%c2n(iLocal,iCell)
              if( Global%Global_to_Local_Node(iNode) == 0 ) then
                if( Global%NodePart(iNode) == iDomain ) then
                  iNodeDomain = iNodeDomain + 1
                  iNodeTotal = iNodeDomain
                else
                  iNodeGhost = iNodeGhost + 1
                  iNodeTotal = iNodeGhost
                end if
                Global%Global_to_Local_Node(iNode) = iNodeTotal
                Buffer_Send_NodeGID(iNodeTotal) = iNode
                Buffer_Send_NodePart(iNodeTotal) = Global%NodePart(iNode)
                do iDim = 1, Buffer_Send_nDim
                  iBuffer = Buffer_Send_nDim*(iNodeTotal-1)+iDim
                  Buffer_Send_xyz(iBuffer) = Global%xyz(iDim,iNode)
                end do
              end if
            end do

            ! Set buffer send node connectivity of the cell
            do iLocal = 1, Global%nNodeCell
              iNode = Global%c2n(iLocal,iCell)
              iBuffer = Global%nNodeCell*(iCellTotal-1)+iLocal
              Buffer_Send_c2n(iBuffer) = Global%Global_to_Local_Node(iNode)
            end do
          end if

        end do

      end do

      ! Set buffer send face info
      iFaceTotal = 0
      do iFace_Part = 1, nFace_Part(iDomain)
        iFace = Face_Part(iFace_Part,iDomain)
        iFaceTotal = iFaceTotal + 1
        Global%Global_to_Local_Face(iFace) = iFaceTotal
        Buffer_Send_FaceGID(iFaceTotal) = iFace
        Buffer_Send_Patch(iFaceTotal) = Global%Patch(iFace)
        do iLocal = 1, Global%nNodeFace
          iNode = Global%f2n(iLocal,iFace)
          iBuffer = Global%nNodeFace*(iFaceTotal-1)+iLocal
          Buffer_Send_f2n(iBuffer) = Global%Global_to_Local_Node(iNode)
        end do
      end do

      if( iDomain-1 /= 0 ) then ! iDomain is different from the master rank 0
        dest = iDomain-1
        ! Non-blocking send comm
        call MPI_Isend(Buffer_Send_NodeGID, Buffer_Send_nNodeTotal,                 MPI_INTEGER,dest,1,MPI_COMM_WORLD,Send_Req(1),ier)
        call MPI_Isend(Buffer_Send_NodePart,Buffer_Send_nNodeTotal,                 MPI_INTEGER,dest,2,MPI_COMM_WORLD,Send_Req(2),ier)
        call MPI_Isend(Buffer_Send_xyz,     Buffer_Send_nNodeTotal*Buffer_Send_nDim,MPI_REAL8,  dest,3,MPI_COMM_WORLD,Send_Req(3),ier)
        call MPI_Isend(Buffer_Send_CellGID, Buffer_Send_nCellTotal,                 MPI_INTEGER,dest,4,MPI_COMM_WORLD,Send_Req(4),ier)
        call MPI_Isend(Buffer_Send_CellPart,Buffer_Send_nCellTotal,                 MPI_INTEGER,dest,5,MPI_COMM_WORLD,Send_Req(5),ier)
        call MPI_Isend(Buffer_Send_c2n,     Buffer_Send_nCellTotal*Global%nNodeCell,MPI_INTEGER,dest,6,MPI_COMM_WORLD,Send_Req(6),ier)
        call MPI_Isend(Buffer_Send_FaceGID, Buffer_Send_nFaceTotal,                 MPI_INTEGER,dest,7,MPI_COMM_WORLD,Send_Req(7),ier)
        call MPI_Isend(Buffer_Send_Patch,   Buffer_Send_nFaceTotal,                 MPI_INTEGER,dest,8,MPI_COMM_WORLD,Send_Req(8),ier)
        call MPI_Isend(Buffer_Send_f2n,     Buffer_Send_nFaceTotal*Global%nNodeFace,MPI_INTEGER,dest,9,MPI_COMM_WORLD,Send_Req(9),ier)
        ! Wait for thie set of non-blocking comm
        call MPI_Waitall(9,Send_Req,Send_Stat,ier)
      else ! iDomain is master rank. so, simply copy buffer send
        Buffer_Recv_NodeGID(:)  = Buffer_Send_NodeGID(:)
        Buffer_Recv_NodePart(:) = Buffer_Send_NodePart(:)
        Buffer_Recv_xyz(:)      = Buffer_Send_xyz(:)
        Buffer_Recv_CellGID(:)  = Buffer_Send_CellGID(:)
        Buffer_Recv_CellPart(:) = Buffer_Send_CellPart(:)
        Buffer_Recv_c2n(:)      = Buffer_Send_c2n(:)
        Buffer_Recv_FaceGID(:)  = Buffer_Send_FaceGID(:)
        Buffer_Recv_Patch(:)    = Buffer_Send_Patch(:)
        Buffer_Recv_f2n(:)      = Buffer_Send_f2n(:)
      end if

      ! Deallocate buffer send node / cell info
      deallocate( Buffer_Send_NodeGID )
      deallocate( Buffer_Send_NodePart )
      deallocate( Buffer_Send_xyz )
      deallocate( Buffer_Send_CellGID )
      deallocate( Buffer_Send_CellPart )
      deallocate( Buffer_Send_c2n )
      deallocate( Buffer_Send_FaceGID )
      deallocate( Buffer_Send_Patch )
      deallocate( Buffer_Send_f2n )

    end if

    if( rank == iDomain-1 ) then

      if( rank /= 0 ) then ! receive buffer array
        srcs = 0 ! master core
        ! Non-blocking recv comm
        call MPI_Irecv(Buffer_Recv_NodeGID, Buffer_Recv_nNodeTotal,                 MPI_INTEGER,srcs,1,MPI_COMM_WORLD,Recv_Req(1),ier)
        call MPI_Irecv(Buffer_Recv_NodePart,Buffer_Recv_nNodeTotal,                 MPI_INTEGER,srcs,2,MPI_COMM_WORLD,Recv_Req(2),ier)
        call MPI_Irecv(Buffer_Recv_xyz,     Buffer_Recv_nNodeTotal*Buffer_Recv_nDim,MPI_REAL8,  srcs,3,MPI_COMM_WORLD,Recv_Req(3),ier)
        call MPI_Irecv(Buffer_Recv_CellGID, Buffer_Recv_nCellTotal,                 MPI_INTEGER,srcs,4,MPI_COMM_WORLD,Recv_Req(4),ier)
        call MPI_Irecv(Buffer_Recv_CellPart,Buffer_Recv_nCellTotal,                 MPI_INTEGER,srcs,5,MPI_COMM_WORLD,Recv_Req(5),ier)
        call MPI_Irecv(Buffer_Recv_c2n,     Buffer_Recv_nCellTotal*Grid%nNodeCell,  MPI_INTEGER,srcs,6,MPI_COMM_WORLD,Recv_Req(6),ier)
        call MPI_Irecv(Buffer_Recv_FaceGID, Buffer_Recv_nFaceTotal,                 MPI_INTEGER,srcs,7,MPI_COMM_WORLD,Recv_Req(7),ier)
        call MPI_Irecv(Buffer_Recv_Patch,   Buffer_Recv_nFaceTotal,                 MPI_INTEGER,srcs,8,MPI_COMM_WORLD,Recv_Req(8),ier)
        call MPI_Irecv(Buffer_Recv_f2n,     Buffer_Recv_nFaceTotal*Grid%nNodeFace,  MPI_INTEGER,srcs,9,MPI_COMM_WORLD,Recv_Req(9),ier)
        ! Wait for thie set of non-blocking comm
        call MPI_Waitall(9,Recv_Req,Recv_Stat,ier)
      end if

      ! Set number of nodes and allocate node data
      Grid%nNode      = Buffer_Recv_nNodeDomain
      Grid%nNodeBound = 0
      Grid%nNodeGhost = Buffer_Recv_nNodeGhost
      call AllocNodeData(Grid)

      ! Set number of cells and allocate cell data
      Grid%nCell      = Buffer_Recv_nCellDomain
      Grid%nCellGhost = Buffer_Recv_nCellGhost
      Grid%nCellDummy = Buffer_Recv_nFaceTotal
      call AllocCellData(Grid)

      ! Set number of faces and allocate face data
      Grid%nFace      = Buffer_Recv_nFaceTotal
      Grid%nFaceBound = Buffer_Recv_nFaceTotal
      Grid%nFaceGhost = 0
      call AllocFaceData(Grid)

      ! Set node data
      do iNode = 1, Grid%nNodeTotal
        Grid%NodeGID(iNode)  = Buffer_Recv_NodeGID(iNode)
        Grid%NodePart(iNode) = Buffer_Recv_NodePart(iNode)
        do iDim = 1, Grid%nDim
          iBuffer = Grid%nDim*(iNode-1)+iDim
          Grid%xyz(iDim,iNode) = Buffer_Recv_xyz(iBuffer)
        end do
      end do

      ! Set cell data
      do iCell = 1, Grid%nCellTotal
        Grid%CellGID(iCell)  = Buffer_Recv_CellGID(iCell)
        Grid%CellPart(iCell) = Buffer_Recv_CellPart(iCell)
        do iLocal = 1, Grid%nNodeCell
          iBuffer = Grid%nNodeCell*(iCell-1)+iLocal
          Grid%c2n(iLocal,iCell) = Buffer_Recv_c2n(iBuffer)
        end do
      end do

      ! Set face data
      do iFace = 1, Grid%nFaceTotal
        Grid%FaceGID(iFace) = Buffer_Recv_FaceGID(iFace)
        Grid%Patch(iFace)   = Buffer_Recv_Patch(iFace)
        do iLocal = 1, Grid%nNodeFace
          iBuffer = Grid%nNodeFace*(iFace-1)+iLocal
          Grid%f2n(iLocal,iFace) = Buffer_Recv_f2n(iBuffer)
        end do
      end do

      ! Deallocate receive buffers
      deallocate( Buffer_Recv_NodeGID )
      deallocate( Buffer_Recv_NodePart )
      deallocate( Buffer_Recv_xyz )
      deallocate( Buffer_Recv_CellGID )
      deallocate( Buffer_Recv_CellPart )
      deallocate( Buffer_Recv_c2n )
      deallocate( Buffer_Recv_FaceGID )
      deallocate( Buffer_Recv_Patch )
      deallocate( Buffer_Recv_f2n )

    end if

    ! Print local grid info
    if( rank == 0 ) then
      write(*,*) 'FLUS >> Setting local grid at rank', iDomain-1
      write(*,*) 'FLUS >> Number of nodes : ', Buffer_Send_nNodeDomain, Buffer_Send_nNodeGhost
      write(*,*) 'FLUS >> Number of cells : ', Buffer_Send_nCellDomain, Buffer_Send_nCellGhost
      write(*,*) 'FLUS >> Number of faces : ', Buffer_Send_nFaceTotal
    end if

  end do

  ! Deallocate data structure
  if( rank == 0 ) then
    deallocate( xadj, adjncy )
    deallocate( nNode_Part, Node_Part )
    deallocate( nFace_Part, Face_Part )
  end if

  ! Set number of global nodes / cells / faces
  call MPI_Allreduce(Grid%nNode,Grid%nNodeGlobal,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ier)
  call MPI_Allreduce(Grid%nFace,Grid%nFaceGlobal,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ier)
  call MPI_Allreduce(Grid%nCell,Grid%nCellGlobal,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ier)

  ! Set number of global boundary nodes
  if( rank == 0 ) Grid%nNodeGlobalBound = Global%nNodeBound
  call MPI_Bcast(Grid%nNodeGlobalBound,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  ! Set MPI barrier
  call MPI_Barrier(MPI_COMM_WORLD,ier)

  ! Print local grid info
  if( rank == 0 ) then
    write(*,*) 'FLUS >> Setting local grid complete'
    write(*,*) 'FLUS >>'
  end if

end subroutine

subroutine SetNodeDataStruc(Grid)
  type(t_Grid), intent(inout) :: Grid

  integer :: rank, ier
  integer :: iLocal, iNode, iFace, iCell
  integer :: iNewID, iBoundary, iInterior
  logical :: Sort_Flag
  integer, allocatable :: NewID(:)
  integer, allocatable :: NewNodeGID(:)
  integer, allocatable :: NewNodePart(:)
  real(8), allocatable :: Newxyz(:,:)
  logical, allocatable :: isBoundary(:)

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Print progress
  if( rank == 0 ) write(*,*) 'FLUS >> Creating node data structure'

  ! Allocate boundary checking flag
  allocate( isBoundary(Grid%nNode) )
  isBoundary(:) = .false.

  ! Count number of boundary nodes
  Grid%nNodeBound = 0
  do iNode = 1, Grid%nNode
    if( Grid%NodeGID(iNode) > Grid%nNodeGlobalBound ) cycle
    Grid%nNodeBound = Grid%nNodeBound + 1
    isBoundary(iNode) = .true.
  end do

  ! Check boundary node whether node ID is greater than number of boundary nodes
  Sort_Flag = .false.
  do iNode = 1, Grid%nNodeBound
    if( .not. isBoundary(iNode) ) then
      Sort_Flag = .true.; exit
    end if
  end do

  ! Sort boundary nodes
  if( Sort_Flag ) then

    ! Allocate new node ID
    allocate( NewID(Grid%nNode) )
    NewID(:) = 0

    ! Create array of new node ID
    iBoundary = 0; iInterior = Grid%nNodeBound
    do iNode = 1, Grid%nNode
      if( isBoundary(iNode) ) then ! boundary node
        iBoundary = iBoundary + 1
        NewID(iNode) = iBoundary
      else ! interior node
        iInterior = iInterior + 1
        NewID(iNode) = iInterior
      end if
    end do

    ! Change node global ID / partition based on new ID
    allocate( NewNodeGID(Grid%nNode) )
    allocate( NewNodePart(Grid%nNode) )
    do iNode = 1, Grid%nNode
      iNewID = NewID(iNode)
      NewNodeGID(iNewID) = Grid%NodeGID(iNode)
      NewNodePart(iNewID) = Grid%NodePart(iNode)
    end do
    Grid%NodeGID(1:Grid%nNode) = NewNodeGID(:)
    Grid%NodePart(1:Grid%nNode) = NewNodePart(:)
    deallocate( NewNodeGID )
    deallocate( NewNodePart )

    ! Change node coordinates based on new ID
    allocate( Newxyz(Grid%nDim,Grid%nNode) )
    do iNode = 1, Grid%nNode
      iNewID = NewID(iNode)
      Newxyz(:,iNewID) = Grid%xyz(:,iNode)
    end do
    Grid%xyz(:,1:Grid%nNode) = Newxyz(:,:)
    deallocate( Newxyz )

    ! Change node connectivity of the cell based on new ID
    do iCell = 1, Grid%nCellTotal
      do iLocal = 1, Grid%nNodeCell
        iNode = Grid%c2n(iLocal,iCell)
        if( iNode > Grid%nNode ) cycle
        iNewID = NewID(iNode)
        Grid%c2n(iLocal,iCell) = iNewID
      end do
    end do

    ! Change node connectivity of the face based on new ID
    do iFace = 1, Grid%nFaceTotal
      do iLocal = 1, Grid%nNodeFace
        iNode = Grid%f2n(iLocal,iFace)
        if( iNode > Grid%nNode ) cycle
        iNewID = NewID(iNode)
        Grid%f2n(iLocal,iFace) = iNewID
      end do
    end do

    ! Deallocate array of the new node ID
    deallocate( NewID )

  end if

  ! Deallocate boundary checking flag
  deallocate( isBoundary )

end subroutine

subroutine SetFaceDataStruc(Grid)
  type(t_Grid), intent(inout) :: Grid

  integer :: rank, error, ier
  integer :: iLocal, iCell, iFace, iFaceBound, iCount
  integer :: nFaceTotal, nFaceDomain, nFaceBound, nFaceEst
  integer :: Min_Node, Current_Face, Previous_Face, counter
  logical :: Face_Checking
  integer, allocatable :: Node_Flag(:)
  integer, allocatable :: f2n_Temp(:,:), Nodes(:)
  integer, allocatable :: f2c_Temp(:,:), Cells(:)
  integer, allocatable :: First_Face(:), Next_Face(:)
  integer, allocatable :: FaceGID_backup(:), Patch_backup(:)

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Print progress
  if( rank == 0 ) write(*,*) 'FLUS >> Creating face data structure'

  ! Estimate number of face
  nFaceEst = 6 * Grid%nCellWhole

  ! Allocate temporary face data structure
  allocate(f2n_Temp(Grid%nNodeFace,nFaceEst))
  allocate(f2c_Temp(Grid%nCellFace,nFaceEst))
  f2n_Temp(:,:) = 0d0
  f2c_Temp(:,:) = 0d0

  ! Allocate temporary data
  allocate( Nodes(Grid%nNodeFace) )
  allocate( Cells(Grid%nCellFace) )

  ! Node flag for checking existing face
  allocate( Node_Flag(Grid%nNodeTotal) )
  Node_Flag(:) = 0

  ! 1) First face connected to node(small number)
  ! 2) Next face connected to previous face
  allocate( First_Face(Grid%nNodeTotal) )
  allocate( Next_Face(nFaceEst) )
  First_Face(:) = 0
  Next_Face(:) = 0

  ! Construct face structure
  nFaceTotal = 0
  do iCell = 1, Grid%nCellTotal
    do iFace = 1, Grid%nFaceCell
      Nodes(:) = Grid%c2n(Grid%LocalNode(:,iFace),iCell)
      Min_Node = minval(Nodes(:))
      Current_Face = First_Face(Min_Node)
      Previous_Face = 0
      Face_Checking = .false.
      do while( .not. Face_Checking )
        if( Current_Face == 0 ) then
          ! Define new face
          nFaceTotal = nFaceTotal + 1
          if( nFaceTotal > nFaceEst ) then
            write(*,*) 'FLUS >> Number of esimated faces too small'
            call MPI_Abort(MPI_COMM_WORLD,error,ier)
          end if
          f2n_Temp(:,nFaceTotal) = Nodes(:)
          f2c_Temp(1,nFaceTotal) = iCell
          if( Previous_Face == 0 ) First_Face(Min_Node) = nFaceTotal
          if( Previous_Face /= 0 ) Next_Face(Previous_Face) = nFaceTotal
          Face_Checking = .true.
        else
          ! Checking existing face
          Node_Flag(Nodes(:)) = 1
          iCount = sum( Node_Flag(f2n_Temp(:,Current_Face)) )
          Node_Flag(Nodes(:)) = 0
          if( iCount == Grid%nNodeFace ) then
            f2c_Temp(2,Current_Face) = iCell
            Face_Checking = .true.
          end if
        end if
        ! Check next face
        if( Face_Checking == .false. ) then
          Previous_Face = Current_Face
          Current_Face = Next_Face(Current_Face)
        end if
      end do
    end do
  end do

  ! Deallocate first / next face
  deallocate( First_Face )
  deallocate( Next_Face )

  ! Count number of domain face
  do iFace = 1, nFaceTotal
    ! Face data structure created in domain ~ ghost cell order
    ! So, just check whether first adjacent cell is greater than domain cell
    if( f2c_Temp(1,iFace) > Grid%nCell ) exit
  end do
  nFaceDomain = iFace-1

  ! Reorder face data (boundary ~ interior)
  nFaceBound = 0
  do iFace = 1, nFaceDomain
    if( f2c_Temp(2,iFace) == 0 ) then
      nFaceBound = nFaceBound + 1
      Nodes(:) = f2n_Temp(:,iFace)
      f2n_Temp(:,iFace) = f2n_Temp(:,nFaceBound)
      f2n_Temp(:,nFaceBound) = Nodes(:)
      Cells(:) = f2c_Temp(:,iFace)
      f2c_Temp(:,iFace) = f2c_Temp(:,nFaceBound)
      f2c_Temp(:,nFaceBound) = Cells(:)
    end if
  end do

  ! Check number of boundary faces
  if( nFaceBound /= Grid%nFace ) then
    write(*,*) 'Given grid is something wrong. need to check boundary face'
    write(*,*) 'From grid file (or surface grid), number of boundary face =', Grid%nFace
    write(*,*) 'From face data structure, number of boundary face =', nFaceBound
    call MPI_Abort(MPI_COMM_WORLD,error,ier)
  end if

  ! Reorder boundary face data (based on existing face)
  do iFace = 1, Grid%nFace
    Node_Flag(Grid%f2n(:,iFace)) = 1
    do iFaceBound = iFace, nFaceBound
      counter = sum( Node_Flag(f2n_Temp(:,iFaceBound)) )
      if( counter == Grid%nNodeFace ) then
        Node_Flag(Grid%f2n(:,iFace)) = 0
        Nodes(:) = f2n_Temp(:,iFace)
        f2n_Temp(:,iFace) = f2n_Temp(:,iFaceBound)
        f2n_Temp(:,iFaceBound) = Nodes(:)
        Cells(:) = f2c_Temp(:,iFace)
        f2c_Temp(:,iFace) = f2c_Temp(:,iFaceBound)
        f2c_Temp(:,iFaceBound) = Cells(:)
        exit
      end if
    end do
  end do

  ! Deallocate temporary data
  deallocate( Nodes )
  deallocate( Cells )

  ! Deallocate node checking flag
  deallocate( Node_Flag )

  ! Backup existing face data (global ID / patch info)
  allocate( FaceGID_backup(Grid%nFace) )
  allocate( Patch_backup(Grid%nFace) )
  FaceGID_backup(:) = Grid%FaceGID(:)
  Patch_backup(:) = Grid%Patch(:)

  ! Deallocate existing face data
  if( allocated( Grid%FaceGID ) ) deallocate( Grid%FaceGID )
  if( allocated( Grid%Patch ) ) deallocate( Grid%Patch )
  if( allocated( Grid%f2n ) ) deallocate( Grid%f2n )
  if( allocated( Grid%f2c ) ) deallocate( Grid%f2c )
  if( allocated( Grid%fc ) ) deallocate( Grid%fc )
  if( allocated( Grid%fn ) ) deallocate( Grid%fn )
  if( allocated( Grid%fa ) ) deallocate( Grid%fa )
  if( allocated( Grid%fv ) ) deallocate( Grid%fv )

  ! Reset number of faces and allocate face data
  ! And now, Grid%nFace is different from Grid%nFaceBound
  ! So, do not confuse number of domain / boundary faces
  Grid%nFace      = nFaceDomain
  Grid%nFaceBound = nFaceBound
  Grid%nFaceGhost = nFaceTotal-nFaceDomain
  Grid%nFaceTotal = nFaceTotal
  call AllocFaceData(Grid)

  ! Set boundary face global ID / patch info using backup data
  do iFaceBound = 1, Grid%nFaceBound
    Grid%FaceGID(iFaceBound) = FaceGID_backup(iFaceBound)
    Grid%Patch(iFaceBound) = Patch_backup(iFaceBound)
  end do

  ! Deallocate backup data
  deallocate( FaceGID_backup )
  deallocate( Patch_backup )

  ! Set node connectivity / adjacent cell data
  ! Due to this process, boundary face orient fixed in outward
  do iFace = 1, Grid%nFaceTotal
    Grid%f2n(:,iFace) = f2n_Temp(:,iFace)
    Grid%f2c(:,iFace) = f2c_Temp(:,iFace)
  end do

  ! Deallocate temporary face data structure
  deallocate( f2n_Temp )
  deallocate( f2c_Temp )

  ! Set adjacent dummy cell data structure
  iCell = Grid%nCellTotal
  do iFaceBound = 1, Grid%nFaceBound
    iCell = iCell + 1
    Grid%f2c(2,iFaceBound) = iCell
    do iLocal = 1, Grid%nNodeFace
      Grid%c2n(iLocal,iCell) = Grid%f2n(iLocal,iFaceBound)
    end do
  end do

end subroutine

subroutine SetCellDataStruc(Grid)
  type(t_Grid), intent(inout) :: Grid

  integer :: rank ,ier
  integer :: iFace, iCell1, iCell2
  integer, allocatable :: nConn(:)

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Print progress
  if( rank == 0 ) write(*,*) 'FLUS >> Creating cell data structure'

  allocate( nConn(Grid%nCellWhole) )
  nConn(:) = 0
  do iFace = 1, Grid%nFaceTotal
    iCell1 = Grid%f2c(1,iFace)
    iCell2 = Grid%f2c(2,iFace)
    nConn(iCell1) = nConn(iCell1) + 1
    Grid%c2c(nConn(iCell1),iCell1) = iCell2
    Grid%c2f(nConn(iCell1),iCell1) = iFace
    if( iCell2 < 1 ) cycle
    nConn(iCell2) = nConn(iCell2) + 1
    Grid%c2c(nConn(iCell2),iCell2) = iCell1
    Grid%c2f(nConn(iCell2),iCell2) = iFace
  end do
  deallocate(nConn)

end subroutine

subroutine SetListDataStruc(Grid,Conf)
  type(t_Grid), intent(inout) :: Grid
  type(t_Conf), intent(in) :: Conf

  integer :: rank, error, ier
  integer :: iLocal, iBC, iFaceBound
  integer :: iPatch, MaxPatch, counter
  integer, allocatable :: Patch2BC(:)

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Print progress
  if( rank == 0 ) write(*,*) 'FLUS >> Creating list data structure'

  ! Allocate list data
  call AllocListData(Grid,Conf)

  ! Set maximum patch number
  MaxPatch = Conf%iPatch(1)
  do iBC = 2, Conf%nBC
    if( Conf%iPatch(iBC) > MaxPatch ) MaxPatch = Conf%iPatch(iBC)
  end do

  ! Allocate patch to BC
  allocate( Patch2BC(0:MaxPatch) )
  Patch2BC(:) = 0

  ! Set patch to BC data
  do iBC = 1, Conf%nBC
    Patch2BC(Conf%iPatch(iBC)) = iBC
  end do

  ! Count number of faces per boundary condition
  do iFaceBound = 1, Grid%nFaceBound
    iPatch = Grid%Patch(iFaceBound)
    iBC = Patch2BC(iPatch)
    Grid%ListIndex(iBC) = Grid%ListIndex(iBC) + 1
  end do

  ! Check all boundary face and boundary condition
  counter = sum( Grid%ListIndex(1:Conf%nBC) )
  if( counter /= Grid%nFaceBound ) then
    write(*,*) 'FLUS >> Grid file patchs are not fully described in boundary input file'
    call MPI_Abort(MPI_COMM_WORLD,error,ier)
  end if

  ! Set list array
  do iBC = Conf%nBC, 2, -1
    Grid%ListIndex(iBC) = sum( Grid%ListIndex(1:iBC-1) )
  end do
  Grid%ListIndex(1) = 0
  do iFaceBound = 1, Grid%nFaceBound
    iPatch = Grid%Patch(iFaceBound)
    iBC = Patch2BC(iPatch)
    Grid%ListIndex(iBC) = Grid%ListIndex(iBC) + 1
    iLocal = Grid%ListIndex(iBC)
    Grid%List(iLocal) = iFaceBound
  end do

  ! Set list index (starting face)
  do iBC = Conf%nBC+1, 2, -1
    Grid%ListIndex(iBC) = Grid%ListIndex(iBC-1)+1
  end do
  Grid%ListIndex(1) = 1

  ! Deallocate patch to BC
  deallocate( Patch2BC )

end subroutine

subroutine SetGridReNumber(Grid)
  type(t_Grid), intent(inout) :: Grid

  integer :: rank, ier, counter
  integer :: iFace, iCell, iLocal, iNeiCell, iNewID
  integer, allocatable :: NewID(:), OldID(:)
  integer, allocatable :: c2n_Temp(:,:)
  integer, allocatable :: c2f_Temp(:,:)
  integer, allocatable :: c2c_Temp(:,:)
  integer, allocatable :: f2c_Temp(:,:)
  integer, allocatable :: NewCellGID(:)
  integer, allocatable :: NewCellPart(:)
  logical, allocatable :: Cell_Flag(:)

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Print progress
  if( rank == 0 ) write(*,*) 'FLUS >> Re-numbering cell data structure'

  ! Allocate new id for cell re-numbering (Id after Re-numbering)
  allocate(NewID(Grid%nCell))
  NewID(:) = 0

  ! Allocate old id for cell re-numbering (ID before Re-numbering)
  allocate(OldID(Grid%nCell))
  OldID(:) = 0

  ! Allocate flag array for checking cell
  allocate( Cell_Flag(Grid%nCell) )
  Cell_Flag(:) = .false.

  ! Allocate temporary cell data structure
  allocate(c2n_Temp(Grid%nNodeCell,Grid%nCell))
  allocate(c2f_Temp(Grid%nFaceCell,Grid%nCell))
  allocate(f2c_Temp(Grid%nCellFace,Grid%nFace))
  allocate(c2c_Temp(Grid%nFaceCell,Grid%nCell))
  c2n_Temp(:,:) = 0d0
  c2f_Temp(:,:) = 0d0
  f2c_Temp(:,:) = 0d0
  c2c_Temp(:,:) = 0d0

  ! Set starting cell with domain cell of the iFace = 1
  iNewID = 1
  iCell = Grid%f2c(1,1)
  OldID(iNewID) = iCell
  Cell_Flag(iCell) = .true.

  ! Set re-numbered cell information on iNewID
  do counter = 1, Grid%nCell
    do iLocal = 1, Grid%nFaceCell
      iCell = OldID(counter)
      iNeiCell = Grid%c2c(iLocal,iCell)
      if( iNeiCell < 1 .or. iNeiCell > Grid%nCell ) cycle
      if( Cell_Flag(iNeiCell) ) cycle
      iNewID = iNewID + 1
      OldID(iNewID) = iNeiCell
      Cell_Flag(iNeiCell) = .true.
    end do
    ! if partition separated
    if( iNewID == counter ) then
      do iCell = 1, Grid%nCell
        if( Cell_Flag(iCell) ) cycle
        iNewID = iNewID + 1
        OldID(iNewID) = iCell
        Cell_Flag(iCell) = .true.
        exit
      end do
    end if
  end do

  ! Set interactive relation between NewID & OldID
  do counter = 1, Grid%nCell
    NewID( OldID(counter) ) = counter
  end do

  ! Sort connecting information (c2n, c2f, c2c)
  do iNewID = 1, Grid%nCell
    c2n_Temp(:,iNewID) = Grid%c2n(:,OldID(iNewID))
    c2f_Temp(:,iNewID) = Grid%c2f(:,OldID(iNewID))
    c2c_Temp(:,iNewID) = Grid%c2c(:,OldID(iNewID))
  end do
  do iNewID = 1, Grid%nCell
    Grid%c2n(:,iNewID) = c2n_Temp(:,iNewID)
    Grid%c2f(:,iNewID) = c2f_Temp(:,iNewID)
    Grid%c2c(:,iNewID) = c2c_Temp(:,iNewID)
  end do

  ! Change connecting information(c2c)
  ! Check whole cell because ghost/dummy c2c include domain cell info
  do iCell = 1, Grid%nCellWhole
    do iLocal = 1, Grid%nFaceCell
      iNeiCell = Grid%c2c(iLocal,iCell)
      if( iNeiCell < 1 .or. iNeiCell > Grid%nCell ) cycle
      Grid%c2c(iLocal,iCell) = NewID(iNeiCell)
    end do
  end do

  ! Change connecting information(f2c)
  do iFace = 1, Grid%nFace
    do iLocal = 1, Grid%nCellFace
      iCell = Grid%f2c(iLocal,iFace)
      if( iCell > Grid%nCell ) cycle
      Grid%f2c(iLocal,iFace) = NewID(iCell)
    end do
  end do

  ! Change cell global ID / Partition based on re-numbered ID
  allocate( NewCellGID(Grid%nCell) )
  allocate( NewCellPart(Grid%nCell) )
  do iCell = 1, Grid%nCell
    iNewID = NewID(iCell)
    NewCellGID(iNewID) = Grid%CellGID(iCell)
    NewCellPart(iNewID) = Grid%CellPart(iCell)
  end do
  Grid%CellGID(1:Grid%nCell) = NewCellGID(:)
  Grid%CellPart(1:Grid%nCell) = NewCellPart(:)
  deallocate( NewCellGID )
  deallocate( NewCellPart )

end subroutine

subroutine SetGlobalToLocal(Grid)
  type(t_Grid), intent(inout) :: Grid

  integer :: rank ,ier
  integer :: iNode, iCell, iFace
  integer :: GlobalID, MaxGlobalID

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Print progress
  if( rank == 0 ) write(*,*) 'FLUS >> Creating global to local mapping data'

  ! Find maximum global node ID
  MaxGlobalID = 0
  do iNode = 1, Grid%nNodeTotal
    if( Grid%NodeGID(iNode) > MaxGlobalID ) then
      MaxGlobalID = Grid%NodeGID(iNode)
    end if
  end do

  ! Allocate / initialize global to local node mapping array
  allocate( Grid%Global_to_Local_Node(MaxGlobalID) )
  Grid%Global_to_Local_Node(:) = 0

  ! Set global to local node mapping
  do iNode = 1, Grid%nNodeTotal
    GlobalID = Grid%NodeGID(iNode)
    Grid%Global_to_Local_Node(GlobalID) = iNode
  end do

  ! Find maximum global face ID (boundary)
  MaxGlobalID = 0
  do iFace = 1, Grid%nFaceBound
    if( Grid%FaceGID(iFace) > MaxGlobalID ) then
      MaxGlobalID = Grid%FaceGID(iFace)
    end if
  end do

  ! Allocate / initialize global to local face mapping array (boundary)
  allocate( Grid%Global_to_Local_Face(MaxGlobalID) )
  Grid%Global_to_Local_Face(:) = 0

  ! Set global to local face mapping (boundary)
  do iFace = 1, Grid%nFaceBound
    GlobalID = Grid%FaceGID(iFace)
    Grid%Global_to_Local_Face(GlobalID) = iFace
  end do

  ! Find maximum global cell ID
  MaxGlobalID = 0
  do iCell = 1, Grid%nCellTotal
    if( Grid%CellGID(iCell) > MaxGlobalID ) then
      MaxGlobalID = Grid%CellGID(iCell)
    end if
  end do

  ! Allocate / initialize global to local cell mapping array
  allocate( Grid%Global_to_Local_Cell(MaxGlobalID) )
  Grid%Global_to_Local_Cell(:) = 0

  ! Set global to local cell mapping
  do iCell = 1, Grid%nCellTotal
    GlobalID = Grid%CellGID(iCell)
    Grid%Global_to_Local_Cell(GlobalID) = iCell
  end do

end subroutine

subroutine SetConnSendRecv(Grid)
  type(t_Grid), intent(inout) :: Grid

  integer :: rank, size, ier
  integer :: srcs, dest, Tag1, Tag2, Req1, Req2
  integer :: iNodeGhost, iCellGhost, GlobalID, LocalID
  integer :: iSendNode, iSendCell, iRecvNode, iRecvCell
  integer :: nSendNode, nSendCell, nRecvNode, nRecvCell
  integer :: iConn, Connected, iDomain, nDomain
  integer, allocatable :: DomainList(:)

  ! MPI comm variables array
  integer, allocatable :: Send_Req(:), Send_Stat(:,:)
  integer, allocatable :: Recv_Req(:), Recv_Stat(:,:)

  ! Get rank of each core / number of cores
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)
  call MPI_Comm_size(MPI_COMM_WORLD,size,ier)

  ! Print progress
  if( rank == 0 ) write(*,*) 'FLUS >> Creating connected domain send/recv data'

  ! Set number of domains
  nDomain = size

  ! Allocate domain flag
  allocate( DomainList(nDomain) )
  DomainList(:) = 0

  ! Count number of connection and set domain list
  Grid%nConn = 0
  do iNodeGhost = Grid%nNode+1, Grid%nNodeTotal
    iDomain = Grid%NodePart(iNodeGhost)
    if( DomainList(iDomain) /= 0 ) cycle
    Grid%nConn = Grid%nConn + 1
    DomainList(iDomain) = Grid%nConn
  end do
  do iCellGhost = Grid%nCell+1, Grid%nCellTotal
    iDomain = Grid%CellPart(iCellGhost)
    if( DomainList(iDomain) /= 0 ) cycle
    Grid%nConn = Grid%nConn + 1
    DomainList(iDomain) = Grid%nConn
  end do

  ! Update connection list
  do iDomain = 1, nDomain
    srcs = iDomain-1
    call MPI_Scatter(DomainList,1,MPI_INTEGER,Connected,1,MPI_INTEGER,srcs,MPI_COMM_WORLD,ier)
    if( Connected /= 0 .and. DomainList(iDomain) == 0 ) then
      Grid%nConn = Grid%nConn + 1
      DomainList(iDomain) = Grid%nConn
    end if
  end do

  ! Allocate connection list
  allocate( Grid%ConnList(Grid%nConn) )
  Grid%ConnList(:) = 0

  ! Set connection list
  do iDomain = 1, nDomain
    if( DomainList(iDomain) /= 0 ) then
      Grid%ConnList(DomainList(iDomain)) = iDomain
    end if
  end do

  ! Allocate / initialize receive list
  allocate( Grid%ConnRecv(Grid%nConn) )
  Grid%ConnRecv(:)%nNode = 0
  Grid%ConnRecv(:)%nCell = 0

  ! Count number of receive data
  do iNodeGhost = Grid%nNode+1, Grid%nNodeTotal
    iDomain = Grid%NodePart(iNodeGhost)
    iConn = DomainList(iDomain)
    Grid%ConnRecv(iConn)%nNode = Grid%ConnRecv(iConn)%nNode + 1
  end do
  do iCellGhost = Grid%nCell+1, Grid%nCellTotal
    iDomain = Grid%CellPart(iCellGhost)
    iConn = DomainList(iDomain)
    Grid%ConnRecv(iConn)%nCell = Grid%ConnRecv(iConn)%nCell + 1
  end do

  ! Allocate receive list data array / initialize recv list data
  do iConn = 1, Grid%nConn
    nRecvNode = Grid%ConnRecv(iConn)%nNode
    nRecvCell = Grid%ConnRecv(iConn)%nCell
    allocate( Grid%ConnRecv(iConn)%Node(nRecvNode) )
    allocate( Grid%ConnRecv(iConn)%Cell(nRecvCell) )
    Grid%ConnRecv(iConn)%nNode = 0; Grid%ConnRecv(iConn)%Node(:) = 0
    Grid%ConnRecv(iConn)%nCell = 0; Grid%ConnRecv(iConn)%Cell(:) = 0
  end do

  ! Set receive data array with global ID
  do iNodeGhost = Grid%nNode+1, Grid%nNodeTotal
    iDomain = Grid%NodePart(iNodeGhost)
    iConn = DomainList(iDomain)
    Grid%ConnRecv(iConn)%nNode = Grid%ConnRecv(iConn)%nNode + 1
    iRecvNode = Grid%ConnRecv(iConn)%nNode
    Grid%ConnRecv(iConn)%Node(iRecvNode) = Grid%NodeGID(iNodeGhost)
  end do
  do iCellGhost = Grid%nCell+1, Grid%nCellTotal
    iDomain = Grid%CellPart(iCellGhost)
    iConn = DomainList(iDomain)
    Grid%ConnRecv(iConn)%nCell = Grid%ConnRecv(iConn)%nCell + 1
    iRecvCell = Grid%ConnRecv(iConn)%nCell
    Grid%ConnRecv(iConn)%Cell(iRecvCell) = Grid%CellGID(iCellGhost)
  end do

  ! Deallocate domain flag
  deallocate( DomainList )

  ! Allocate / initialize send list
  allocate( Grid%ConnSend(Grid%nConn) )
  Grid%ConnSend(:)%nNode = 0
  Grid%ConnSend(:)%nCell = 0

  ! Allocate send/recv request and status
  allocate( Send_Req(2*Grid%nConn), Send_Stat(MPI_STATUS_SIZE,2*Grid%nConn) )
  allocate( Recv_Req(2*Grid%nConn), Recv_Stat(MPI_STATUS_SIZE,2*Grid%nConn) )

  ! Non-blocking recv comm
  do iConn = 1, Grid%nConn
    iDomain = Grid%ConnList(iConn); srcs = iDomain-1
    Tag1 = iDomain; Tag2 = iDomain+nDomain
    Req1 = iConn;   Req2 = iConn+Grid%nConn
    call MPI_Irecv(Grid%ConnSend(iConn)%nNode,1,MPI_INTEGER,srcs,Tag1,MPI_COMM_WORLD,Recv_Req(Req1),ier)
    call MPI_Irecv(Grid%ConnSend(iConn)%nCell,1,MPI_INTEGER,srcs,Tag2,MPI_COMM_WORLD,Recv_Req(Req2),ier)
  end do

  ! Non-blocking send comm
  do iConn = 1, Grid%nConn
    iDomain = Grid%ConnList(iConn); dest = iDomain-1
    Tag1 = rank+1; Tag2 = rank+1+nDomain
    Req1 = iConn;  Req2 = iConn+Grid%nConn
    call MPI_Isend(Grid%ConnRecv(iConn)%nNode,1,MPI_INTEGER,dest,Tag1,MPI_COMM_WORLD,Send_Req(Req1),ier)
    call MPI_Isend(Grid%ConnRecv(iConn)%nCell,1,MPI_INTEGER,dest,Tag2,MPI_COMM_WORLD,Send_Req(Req2),ier)
  end do

  ! Wait until all MPI process finish
  call MPI_Waitall(2*Grid%nConn,Send_Req,Send_Stat,ier)
  call MPI_Waitall(2*Grid%nConn,Recv_Req,Recv_Stat,ier)

  ! Allocate / initialize send list data array
  do iConn = 1, Grid%nConn
    nSendNode = Grid%ConnSend(iConn)%nNode
    nSendCell = Grid%ConnSend(iConn)%nCell
    allocate( Grid%ConnSend(iConn)%Node(nSendNode) )
    allocate( Grid%ConnSend(iConn)%Cell(nSendCell) )
    Grid%ConnSend(iConn)%Node(:) = 0
    Grid%ConnSend(iConn)%Cell(:) = 0
  end do

  ! Non-blocking recv comm
  do iConn = 1, Grid%nConn
    iDomain = Grid%ConnList(iConn); srcs = iDomain-1
    nSendNode = Grid%ConnSend(iConn)%nNode
    nSendCell = Grid%ConnSend(iConn)%nCell
    Tag1 = iDomain; Tag2 = iDomain+nDomain
    Req1 = iConn;   Req2 = iConn+Grid%nConn
    call MPI_Irecv(Grid%ConnSend(iConn)%Node,nSendNode,MPI_INTEGER,srcs,Tag1,MPI_COMM_WORLD,Recv_Req(Req1),ier)
    call MPI_Irecv(Grid%ConnSend(iConn)%Cell,nSendCell,MPI_INTEGER,srcs,Tag2,MPI_COMM_WORLD,Recv_Req(Req2),ier)
  end do

  ! Non-blocking send comm
  do iConn = 1, Grid%nConn
    iDomain = Grid%ConnList(iConn); dest = iDomain-1
    nRecvNode = Grid%ConnRecv(iConn)%nNode
    nRecvCell = Grid%ConnRecv(iConn)%nCell
    Tag1 = rank+1; Tag2 = rank+1+nDomain
    Req1 = iConn;  Req2 = iConn+Grid%nConn
    call MPI_Isend(Grid%ConnRecv(iConn)%Node,nRecvNode,MPI_INTEGER,dest,Tag1,MPI_COMM_WORLD,Send_Req(Req1),ier)
    call MPI_Isend(Grid%ConnRecv(iConn)%Cell,nRecvCell,MPI_INTEGER,dest,Tag2,MPI_COMM_WORLD,Send_Req(Req2),ier)
  end do

  ! Wait until all MPI process finish
  call MPI_Waitall(2*Grid%nConn,Send_Req,Send_Stat,ier)
  call MPI_Waitall(2*Grid%nConn,Recv_Req,Recv_Stat,ier)

  ! Deallocate send/recv request and status
  deallocate( Send_Req, Send_Stat )
  deallocate( Recv_Req, Recv_Stat )

  ! Transfer global ID to local ID in send/recv list
  do iConn = 1, Grid%nConn
    ! Transfer global ID to local ID in send list
    do iSendNode = 1, Grid%ConnSend(iConn)%nNode
      GlobalID = Grid%ConnSend(iConn)%Node(iSendNode)
      LocalID = Grid%Global_to_Local_Node(GlobalID)
      Grid%ConnSend(iConn)%Node(iSendNode) = LocalID
    end do
    do iSendCell = 1, Grid%ConnSend(iConn)%nCell
      GlobalID = Grid%ConnSend(iConn)%Cell(iSendCell)
      LocalID = Grid%Global_to_Local_Cell(GlobalID)
      Grid%ConnSend(iConn)%Cell(iSendCell) = LocalID
    end do
    ! Transfer global ID to local ID in recv list
    do iRecvNode = 1, Grid%ConnRecv(iConn)%nNode
      GlobalID = Grid%ConnRecv(iConn)%Node(iRecvNode)
      LocalID = Grid%Global_to_Local_Node(GlobalID)
      Grid%ConnRecv(iConn)%Node(iRecvNode) = LocalID
    end do
    do iRecvCell = 1, Grid%ConnRecv(iConn)%nCell
      GlobalID = Grid%ConnRecv(iConn)%Cell(iRecvCell)
      LocalID = Grid%Global_to_Local_Cell(GlobalID)
      Grid%ConnRecv(iConn)%Cell(iRecvCell) = LocalID
    end do
  end do

  ! Allocate send/recv variable data
  do iConn = 1, Grid%nConn
    nSendNode = Grid%ConnSend(iConn)%nNode
    allocate( Grid%ConnSend(iConn)%gv(nSendNode*Grid%nDim) )
    Grid%ConnSend(iConn)%gv(:) = 0d0
    nRecvNode = Grid%ConnRecv(iConn)%nNode
    allocate( Grid%ConnRecv(iConn)%gv(nRecvNode*Grid%nDim) )
    Grid%ConnRecv(iConn)%gv(:) = 0d0
  end do

  ! MPI barrier
  call MPI_Barrier(MPI_COMM_WORLD,ier)

end subroutine

subroutine SetBoundRecv(Grid)
  type(t_Grid), intent(inout) :: Grid

  integer :: rank, size, srcs, ier
  integer :: iBound, iDomain, nDomain
  integer, allocatable :: nNodeBound_Part(:)
  integer, allocatable :: nFaceBound_Part(:)

  ! Get rank of each core / number of cores
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)
  call MPI_Comm_size(MPI_COMM_WORLD,size,ier)

  ! Print progress
  if( rank == 0 ) write(*,*) 'FLUS >> Creating boundary domain recv data'

  ! Set number of domains
  nDomain = size

  ! Allocate / gather number of boundary nodes / faces per partition
  allocate( nNodeBound_Part(nDomain) ); nNodeBound_Part(:) = 0
  allocate( nFaceBound_Part(nDomain) ); nFaceBound_Part(:) = 0
  call MPI_Allgather(Grid%nNodeBound,1,MPI_INTEGER,nNodeBound_Part,1,MPI_INTEGER,MPI_COMM_WORLD,ier)
  call MPI_Allgather(Grid%nFaceBound,1,MPI_INTEGER,nFaceBound_Part,1,MPI_INTEGER,MPI_COMM_WORLD,ier)

  ! Count number of boundary domain
  Grid%nBound = 0
  do iDomain = 1, nDomain
    if( nNodeBound_Part(iDomain) + nFaceBound_Part(iDomain) > 0 ) then
      Grid%nBound = Grid%nBound + 1
    end if
  end do

  ! Allocate / initialize boundary domain list
  allocate( Grid%BoundList(Grid%nBound) )
  Grid%BoundList(:) = 0

  ! Set boundary domain list
  Grid%nBound = 0
  do iDomain = 1, nDomain
    if( nNodeBound_Part(iDomain) + nFaceBound_Part(iDomain) > 0 ) then
      Grid%nBound = Grid%nBound + 1
      Grid%BoundList(Grid%nBound) = iDomain
    end if
  end do

  ! Set master recv list
  allocate( Grid%BoundRecv(Grid%nBound) )
  do iBound = 1, Grid%nBound
    iDomain = Grid%BoundList(iBound)
    Grid%BoundRecv(iBound)%nNode = nNodeBound_Part(iDomain)
    Grid%BoundRecv(iBound)%nFace = nFaceBound_Part(iDomain)
    allocate( Grid%BoundRecv(iBound)%Node(Grid%BoundRecv(iBound)%nNode) )
    allocate( Grid%BoundRecv(iBound)%Face(Grid%BoundRecv(iBound)%nFace) )
  end do

  ! Deallocate number of boundary nodes / faces per partition
  deallocate( nNodeBound_Part )
  deallocate( nFaceBound_Part )

  ! Set boundary domain node / face list
  do iBound = 1, Grid%nBound
    iDomain = Grid%BoundList(iBound); srcs = iDomain-1
    if( rank == iDomain-1 ) then
      Grid%BoundRecv(iBound)%Node(:) = Grid%NodeGID(1:Grid%nNodeBound)
      Grid%BoundRecv(iBound)%Face(:) = Grid%FaceGID(1:Grid%nFaceBound)
    end if
    call MPI_Bcast(Grid%BoundRecv(iBound)%Node(:),Grid%BoundRecv(iBound)%nNode,MPI_INTEGER,srcs,MPI_COMM_WORLD,ier)
    call MPI_Bcast(Grid%BoundRecv(iBound)%Face(:),Grid%BoundRecv(iBound)%nFace,MPI_INTEGER,srcs,MPI_COMM_WORLD,ier)
  end do

  ! MPI barrier
  call MPI_Barrier(MPI_COMM_WORLD,ier)

end subroutine

subroutine SetFaceMetric(Grid)
  type(t_Grid), intent(inout) :: Grid

  integer :: iFace, iLocal, iNode
  integer :: iNode1, iNode2, iNode3
  real(8) :: a(3), b(3), Coeff

  ! Compute face normal
  select case(Grid%nDim)
  case(2)
    do iFace = 1,Grid%nFaceTotal
      iNode1 = Grid%f2n(1,iFace)
      iNode2 = Grid%f2n(2,iFace)
      Grid%fn(1,iFace) =   ( Grid%xyz(2,iNode2) - Grid%xyz(2,iNode1) )
      Grid%fn(2,iFace) = - ( Grid%xyz(1,iNode2) - Grid%xyz(1,iNode1) )
    end do
  case(3)
    do iFace = 1,Grid%nFaceTotal
      iNode1 = Grid%f2n(1,iFace)
      iNode2 = Grid%f2n(2,iFace)
      iNode3 = Grid%f2n(3,iFace)
      a(:) = Grid%xyz(:,iNode2) - Grid%xyz(:,iNode1)
      b(:) = Grid%xyz(:,iNode3) - Grid%xyz(:,iNode1)
      Grid%fn(1,iFace) = 0.5d0 * ( a(2)*b(3) - a(3)*b(2) )
      Grid%fn(2,iFace) = 0.5d0 * ( a(3)*b(1) - a(1)*b(3) )
      Grid%fn(3,iFace) = 0.5d0 * ( a(1)*b(2) - a(2)*b(1) )
    end do
  end select

  ! Compute face unit normal / face area
  do iFace = 1, Grid%nFaceTotal
    Grid%fa(iFace) = dsqrt( sum( Grid%fn(:,iFace)**2 ) )
    Grid%fn(:,iFace) = Grid%fn(:,iFace) / Grid%fa(iFace)
  end do

  ! Compute face center
  Coeff = 1d0 / dble(Grid%nNodeFace)
  do iFace = 1,Grid%nFaceTotal
    Grid%fc(:,iFace) = 0d0
    do iLocal = 1, Grid%nNodeFace
      iNode = Grid%f2n(iLocal,iFace)
      Grid%fc(:,iFace) = Grid%fc(:,iFace) + Grid%xyz(:,iNode)
    end do
    Grid%fc(:,iFace) = Grid%fc(:,iFace) * Coeff
  end do

end subroutine

subroutine SetCellMetric(Grid)
  type(t_Grid), intent(inout) :: Grid

  integer :: iLocal
  integer :: iCell, iNode
  integer :: iFace, iFaceBound
  integer :: iCell1, iCell2
  real(8) :: Coeff, Distance

  ! Compute cell center
  Coeff = 1d0 / dble(Grid%nNodeCell)
  do iCell = 1, Grid%nCellTotal
    Grid%cen(:,iCell) = 0d0
    do iLocal = 1, Grid%nNodeCell
      iNode = Grid%c2n(iLocal,iCell)
      Grid%cen(:,iCell) = Grid%cen(:,iCell) + Grid%xyz(:,iNode)
    end do
    Grid%cen(:,iCell) = Grid%cen(:,iCell) * Coeff
  enddo

  ! Compute dummy cell center
  do iFaceBound = 1, Grid%nFaceBound
    iCell1 = Grid%f2c(1,iFaceBound)
    iCell2 = Grid%f2c(2,iFaceBound)
    ! normal distance icell2 center to ibface
    Distance = sum( (Grid%fc(:,iFaceBound) - Grid%cen(:,iCell1)) * Grid%fn(:,iFaceBound) )
    ! Set dummy cell center
    Grid%cen(:,iCell2) = 2d0 * Distance * Grid%fn(:,iFaceBound) + Grid%cen(:,iCell1)
  end do

  ! Compute cell volume
  Coeff = 1d0 / Grid%nDim
  do iCell = 1,Grid%nCellTotal
    Grid%vol(iCell) = 0d0
    do iLocal = 1, Grid%nFaceCell
      iFace = Grid%c2f(iLocal,iCell)
      if( Grid%f2c(1,iFace) == iCell ) then
        Grid%vol(iCell) = Grid%vol(iCell) + sum( Grid%fc(:,iFace) * Grid%fn(:,iFace) ) * Grid%fa(iFace)
      else
        Grid%vol(iCell) = Grid%vol(iCell) - sum( Grid%fc(:,iFace) * Grid%fn(:,iFace) ) * Grid%fa(iFace)
      end if
    end do
    Grid%vol(iCell) = Grid%vol(iCell) * Coeff
  end do

end subroutine

subroutine SetWallDistance(Grid,Conf)
  type(t_Grid), intent(inout) :: Grid
  type(t_Conf), intent(in) :: Conf

  integer :: rank, srcs, ier
  integer :: iBound, iDomain, iBuffer, iLocal
  integer :: iDim, iCell, iBFace, iBC
  integer :: iNodeBound, iFaceBound
  integer :: iFacePartBound, nFacePartBound
  integer :: iList, List1, List2
  real(8) :: Distance

  ! MPI comm array
  real(8), allocatable :: Buffer_fc(:)
  logical, allocatable :: Buffer_WallFlag(:)

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Initialize wall distance
  Grid%WallDist(:) = 2d0**1023

  do iBound = 1, Grid%nBound
    iDomain = Grid%BoundList(iBound); srcs = iDomain-1

    ! Set number of partition boundary faces
    nFacePartBound = Grid%BoundRecv(iBound)%nFace

    ! Allocate boundary face center coordinates
    allocate( Buffer_fc(nFacePartBound*Grid%nDim) )
    allocate( Buffer_WallFlag(nFacePartBound) )

    ! Set buffer data
    if( rank == iDomain-1 ) then
      do iFaceBound = 1, Grid%nFaceBound
        do iDim = 1, Grid%nDim
          iBuffer = Grid%nDim*(iFaceBound-1)+iDim
          Buffer_fc(iBuffer) = 0d0
          do iLocal = 1, Grid%nNodeFace
            iNodeBound = Grid%f2n(iLocal,iFaceBound)
            Buffer_fc(iBuffer) = Buffer_fc(iBuffer) + Grid%xyz(iDim,iNodeBound)
          end do
          Buffer_fc(iBuffer) = Buffer_fc(iBuffer) / Grid%nNodeFace
        end do
      end do
      do iBC = 1, Conf%nBC
        List1 = Grid%ListIndex(iBC)
        List2 = Grid%ListIndex(iBC+1) - 1
        select case(Conf%BCType(iBC))
        case(1,2,12,13)
          do iList = List1, List2
            iBFace = Grid%List(iList)
            Buffer_WallFlag(iBFace) = .true.
          end do
        case default
          do iList = List1, List2
            iBFace = Grid%List(iList)
            Buffer_WallFlag(iBFace) = .false.
          end do
        end select
      end do
    end if

    ! Broadcast buffer data
    call MPI_Bcast(Buffer_fc,nFacePartBound*Grid%nDim,MPI_REAL8,srcs,MPI_COMM_WORLD,ier)
    call MPI_Bcast(Buffer_WallFlag,nFacePartBound,MPI_LOGICAL,srcs,MPI_COMM_WORLD,ier)

    ! Compute wall distance
    do iCell = 1, Grid%nCell
      do iFacePartBound = 1, nFacePartBound
        if( Buffer_WallFlag(iFacePartBound) ) then
          iBuffer = Grid%nDim*(iFacePartBound-1)
          Distance = dsqrt( sum( ( Grid%cen(:,iCell) - Buffer_fc(iBuffer+1:iBuffer+Grid%nDim) )**2 ) )
          if( Distance < Grid%WallDist(iCell) ) Grid%WallDist(iCell) = Distance
        end if
      end do
    end do

    ! Deallocate buffer data
    deallocate( Buffer_fc )
    deallocate( Buffer_WallFlag )

  end do

  ! MPI barrier
  call MPI_Barrier(MPI_COMM_WORLD,ier)

end subroutine

subroutine PreGridVelocity(Grid)
  type(t_Grid), intent(inout) :: Grid

  integer :: rank, srcs, ier
  integer :: iBound, iBuffer, iLocal1, iLocal2
  integer :: iDim, iNode, iNodeBound, iDomain
  integer :: iNodePartBound, iNodeGlobalBound
  integer :: nNodePartBound, nNodeGlobalBound

  ! MPI comm array
  real(8), allocatable :: Buffer_xyz(:)

  ! Get rank of each core / number of cores
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Print progress
  if( rank == 0 ) write(*,*) 'FLUS >> Preprocess for grid velocity computation'

  ! Get number of global boundary nodes
  nNodeGlobalBound = sum( Grid%BoundRecv(:)%nNode )

  ! Allocate / initialize IDW weighting / summation of that
  if( allocated( Grid%IDW ) ) deallocate( Grid%IDW )
  if( allocated( Grid%IDW_SUM ) ) deallocate( Grid%IDW_SUM )
  allocate( Grid%IDW(nNodeGlobalBound,Grid%nNode) )
  allocate( Grid%IDW_SUM(Grid%nNode) )
  Grid%IDW(:,:) = 0d0
  Grid%IDW_SUM(:) = 0d0

  do iBound = 1, Grid%nBound
    iDomain = Grid%BoundList(iBound); srcs = iDomain-1

    ! Set number of partition boundary nodes
    nNodePartBound = Grid%BoundRecv(iBound)%nNode

    ! Allocate boundary coordinates
    allocate( Buffer_xyz(nNodePartBound*Grid%nDim) )

    ! Set buffer data
    if( rank == iDomain-1 ) then
      do iNodeBound = 1, Grid%nNodeBound
        do iDim = 1, Grid%nDim
          iBuffer = Grid%nDim*(iNodeBound-1)+iDim
          Buffer_xyz(iBuffer) = Grid%xyz(iDim,iNodeBound)
        end do
      end do
    end if

    ! Broadcast buffer data
    call MPI_Bcast(Buffer_xyz,nNodePartBound*Grid%nDim,MPI_REAL8,srcs,MPI_COMM_WORLD,ier)

    ! Compute IDW weighting / summation of the weighting in interior nodes
    do iNode = Grid%nNodeBound+1, Grid%nNode
      iNodeGlobalBound = sum( Grid%BoundRecv(1:iBound-1)%nNode )
      do iNodePartBound = 1, nNodePartBound
        iLocal1 = Grid%nDim*(iNodePartBound-1)+1
        iLocal2 = Grid%nDim*(iNodePartBound-1)+Grid%nDim
        iNodeGlobalBound = iNodeGlobalBound + 1
        Grid%IDW(iNodeGlobalBound,iNode) = 1d0 / sum( ( Grid%xyz(:,iNode) - Buffer_xyz(iLocal1:iLocal2) )**2 )
        Grid%IDW_SUM(iNode) = Grid%IDW_SUM(iNode) + Grid%IDW(iNodeGlobalBound,iNode)
      end do
    end do

    ! Deallocate buffer data
    deallocate( Buffer_xyz )

  end do

  ! MPI barrier
  call MPI_Barrier(MPI_COMM_WORLD,ier)

end subroutine

subroutine SetGridVelocity(Grid)
  type(t_Grid), intent(inout) :: Grid

  integer :: rank, srcs, ier
  integer :: iBound, iBuffer, iLocal1, iLocal2
  integer :: iDim, iNode, iNodeBound, iDomain
  integer :: iNodePartBound, iNodeGlobalBound
  integer :: nNodePartBound

  ! MPI comm array
  real(8), allocatable :: Buffer_gv(:)

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Print progress
  if( rank == 0 ) write(*,*) 'FLUS >> Computing grid velocity using IDW smoothing'

  ! Initialize grid velocity of the interior node
  Grid%gv(:,Grid%nNodeBound+1:Grid%nNode) = 0d0

  do iBound = 1, Grid%nBound
    iDomain = Grid%BoundList(iBound); srcs = iDomain-1

    ! Set number of partition boundary nodes
    nNodePartBound = Grid%BoundRecv(iBound)%nNode

    ! Allocate / initialize buffer send / recv data
    allocate( Buffer_gv(nNodePartBound*Grid%nDim) )

    ! Set buffer data
    if( rank == iDomain-1 ) then
      do iNodeBound = 1, Grid%nNodeBound
        do iDim = 1, Grid%nDim
          iBuffer = Grid%nDim*(iNodeBound-1)+iDim
          Buffer_gv(iBuffer) = Grid%gv(iDim,iNodeBound)
        end do
      end do
    end if

    ! Broadcast buffer data
    call MPI_Bcast(Buffer_gv,nNodePartBound*Grid%nDim,MPI_REAL8,srcs,MPI_COMM_WORLD,ier)

    ! Compute grid velocity (IDW weighting)
    do iNode = Grid%nNodeBound+1, Grid%nNode
      iNodeGlobalBound = sum( Grid%BoundRecv(1:iBound-1)%nNode )
      do iNodePartBound = 1, nNodePartBound
        iLocal1 = Grid%nDim*(iNodePartBound-1)+1
        iLocal2 = Grid%nDim*(iNodePartBound-1)+Grid%nDim
        iNodeGlobalBound = iNodeGlobalBound + 1
        Grid%gv(:,iNode) = Grid%gv(:,iNode) + Grid%IDW(iNodeGlobalBound,iNode) * Buffer_gv(iLocal1:iLocal2)
      end do
    end do

    ! Deallocate buffer data
    deallocate( Buffer_gv )

  end do

  ! Compute grid velocity (normalize)
  do iNode = Grid%nNodeBound+1, Grid%nNode
    Grid%gv(:,iNode) = Grid%gv(:,iNode) / Grid%IDW_SUM(iNode)
  end do

  ! Set grid velocity in ghost node using MPI
  call SetMPIGridVelocity(Grid)

end subroutine

subroutine SetMPIGridVelocity(Grid)
  type(t_Grid), intent(inout) :: Grid

  integer :: rank, dest, srcs, Tag, ier
  integer :: iDim, iNode, iNode_List
  integer :: iLocal, iConn, iDomain, nBuffer

  ! MPI comm request / status
  integer, allocatable :: Send_Req(:), Send_Stat(:,:)
  integer, allocatable :: Recv_Req(:), Recv_Stat(:,:)

  ! Get rank of each core
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)

  ! Allocate MPI comm request / status
  allocate( Send_Req(Grid%nConn), Send_Stat(MPI_STATUS_SIZE,Grid%nConn) )
  allocate( Recv_Req(Grid%nConn), Recv_Stat(MPI_STATUS_SIZE,Grid%nConn) )

  ! Non-blocking recv comm
  do iConn = 1, Grid%nConn
    iDomain = Grid%ConnList(iConn); srcs = iDomain-1; Tag = iDomain
    nBuffer = Grid%ConnRecv(iConn)%nNode * Grid%nDim
    call MPI_Irecv(Grid%ConnRecv(iConn)%gv,nBuffer,MPI_REAL8,srcs,Tag,MPI_COMM_WORLD,Recv_Req(iConn),ier)
  end do

  ! Set buffer send data
  do iConn = 1, Grid%nConn
    do iNode_List = 1, Grid%ConnSend(iConn)%nNode
      iNode = Grid%ConnSend(iConn)%Node(iNode_List)
      do iDim = 1, Grid%nDim
        iLocal = Grid%nDim*(iNode_List-1)+iDim
        Grid%ConnSend(iConn)%gv(iLocal) = Grid%gv(iDim,iNode)
      end do
    end do
  end do

  ! Non-blocking send comm
  do iConn = 1, Grid%nConn
    iDomain = Grid%ConnList(iConn); dest = iDomain-1; Tag = rank+1
    nBuffer = Grid%ConnSend(iConn)%nNode * Grid%nDim
    call MPI_Isend(Grid%ConnSend(iConn)%gv,nBuffer,MPI_REAL8,dest,Tag,MPI_COMM_WORLD,Send_Req(iConn),ier)
  end do

  ! Wait for the set of non-blocking comm
  call MPI_Waitall(Grid%nConn,Send_Req,Send_Stat,ier)
  call MPI_Waitall(Grid%nConn,Recv_Req,Recv_Stat,ier)

  ! Deallocate MPI comm request / status
  deallocate( Send_Req, Send_Stat )
  deallocate( Recv_Req, Recv_Stat )

  ! Set grid velocity from buffer recv
  do iConn = 1, Grid%nConn
    do iNode_List = 1, Grid%ConnRecv(iConn)%nNode
      iNode = Grid%ConnRecv(iConn)%Node(iNode_List)
      do iDim = 1, Grid%nDim
        iLocal = Grid%nDim*(iNode_List-1)+iDim
        Grid%gv(iDim,iNode) = Grid%ConnRecv(iConn)%gv(iLocal)
      end do
    end do
  end do

end subroutine

subroutine SetFaceVelocity(Grid,DeltaTime)
  type(t_Grid), intent(inout) :: Grid
  real(8), intent(in) :: DeltaTime

  integer :: iFace, iLocal
  integer :: iNode1, iNode2, iNode3
  real(8) :: dx, dy, xyz(3,3), Sfn(3), gv_avg(3), Sfv
  real(8) :: Coeff(2), vec1(3), vec2(3), Cross(3)
  real(8), parameter :: inv3 = 1d0 / 3d0

  select case(Grid%nDim)
  case(2)
    do iFace = 1, Grid%nFaceTotal
      iNode1 = Grid%f2n(1,iFace)
      iNode2 = Grid%f2n(2,iFace)
      xyz(1:2,1) = Grid%xyz(:,iNode1) + 0.5d0 * DeltaTime * Grid%gv(:,iNode1)
      xyz(1:2,2) = Grid%xyz(:,iNode2) + 0.5d0 * DeltaTime * Grid%gv(:,iNode2)
      dx = xyz(1,2) - xyz(1,1)
      dy = xyz(2,2) - xyz(2,1)
      Sfn(1) =  dy
      Sfn(2) = -dx
      gv_avg(1:2) = 0.5d0 * ( Grid%gv(:,iNode1) + Grid%gv(:,iNode2) )
      Sfv = gv_avg(1) * Sfn(1) + gv_avg(2) * Sfn(2)
      Grid%fv(iFace) = Sfv / Grid%fa(iFace)
    end do

  case(3)
    Coeff(1) = 0.5d0 - 0.5d0 * dsqrt(inv3)
    Coeff(2) = 0.5d0 + 0.5d0 * dsqrt(inv3)
    do iFace = 1, Grid%nFaceTotal
      iNode1 = Grid%f2n(1,iFace)
      iNode2 = Grid%f2n(2,iFace)
      iNode3 = Grid%f2n(3,iFace)
      Sfn(:) = 0d0
      do iLocal = 1, 2
        xyz(:,1) = Grid%xyz(:,iNode1) + DeltaTime * Grid%gv(:,iNode1) * Coeff(iLocal)
        xyz(:,2) = Grid%xyz(:,iNode2) + DeltaTime * Grid%gv(:,iNode2) * Coeff(iLocal)
        xyz(:,3) = Grid%xyz(:,iNode3) + DeltaTime * Grid%gv(:,iNode3) * Coeff(iLocal)
        vec1(:) = xyz(:,1) - xyz(:,3)
        vec2(:) = xyz(:,2) - xyz(:,3)
        Cross(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
        Cross(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
        Cross(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)
        Sfn(:) = Sfn(:) + Cross(:)
      end do
      Sfn(:) = Sfn(:) * 0.25d0
      gv_avg(:) = ( Grid%gv(:,iNode1) + Grid%gv(:,iNode2) + Grid%gv(:,iNode3) ) * inv3
      Sfv = gv_avg(1) * Sfn(1) + gv_avg(2) * Sfn(2) + gv_avg(3) * Sfn(3)
      Grid%fv(iFace) = Sfv / Grid%fa(iFace)
    end do

  end select

end subroutine

subroutine UpdateGrid(Grid,TimeStep)
  type(t_Grid), intent(inout) :: Grid
  real(8) :: TimeStep

  integer :: iNode

  ! Update node coordinates
  do iNode = 1, Grid%nNodeTotal
    Grid%xyz(:,iNode) = Grid%xyz(:,iNode) + TimeStep * Grid%gv(:,iNode)
  end do

end subroutine

subroutine PostTempPLT(Grid)
  type(t_Grid), intent(in) :: Grid

  integer :: iDim, iNode, iCell, io

  open(newunit=io,file='test.plt')
  write(io,*) 'title="surface"'
  write(io,*) 'variables="x","y","nodepart","cellpart"'
  write(io,*) 'zone zonetype=FETRIANGLE, n=', Grid%nNodeTotal, ', E=',Grid%nCellTotal, ',DATAPACKING=BLOCK, varlocation=([4]=CELLCENTERED)'
  do iDim = 1, 2
    do iNode = 1, Grid%nNodeTotal
      write(io,*) Grid%xyz(iDim,iNode)
    end do
  end do
  do iNode = 1, Grid%nNodeTotal
    write(io,*) Grid%NodePart(iNode)
  end do
  do iCell = 1, Grid%nCellTotal
    write(io,*) Grid%CellPart(iCell)
  end do
  do iCell =1 , Grid%nCellTotal
    write(io,*) Grid%c2n(:,iCell)
  end do
  close(io)

end subroutine

end module
