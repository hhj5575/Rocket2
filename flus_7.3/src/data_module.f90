module data_module

use mpi
use grid_module
use quadtree_module
use octree_module
implicit none; private
public DataTransfer

contains

subroutine DataTransfer(sGrid,tGrid,nData,sData,tData)
  type(t_Grid), intent(in) :: sGrid
  type(t_Grid), intent(in) :: tGrid
  integer, intent(in) :: nData
  real(8), intent(in) :: sData(:,:)
  real(8), intent(out) :: tData(:,:)

  integer :: rank, size, srcs, error, ier
  integer :: iData, iBuffer, iDomain, nDomain
  integer :: iDim, iLocal, iNode, jNode, iCell
  integer :: n1, n2, n3, n4, iadj, nApprox
  real(8) :: Vol1, Vol2, Vol3, Vol4, dVol
  real(8), allocatable :: dVolMin(:)
  integer, allocatable :: xadj(:), adjncy(:,:)
  real(8), allocatable :: bNodeData(:,:)
  real(8), allocatable :: sNodeData(:,:)
  real(8), allocatable :: tNodeData(:,:)
  logical, allocatable :: Check_Flag(:)

  type(t_Grid) :: bGrid
  real(8), allocatable :: Buffer_xyz(:)
  integer, allocatable :: Buffer_c2n(:)
  real(8), allocatable :: Buffer_vol(:)
  real(8), allocatable :: Buffer_Dat(:)

  ! Get rank of each core / number of cores
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ier)
  call MPI_Comm_size(MPI_COMM_WORLD,size,ier)

  ! Set number of domains
  nDomain = size

  ! Allocate / initialize node data
  allocate( sNodeData(nData,sGrid%nNodeTotal) )
  allocate( tNodeData(nData,tGrid%nNodeTotal) )
  sNodeData(:,:) = 0d0
  tNodeData(:,:) = 0d0

  ! Set source node data
  call SetNodeData(sGrid,sData,sNodeData)

  ! Allocate / initialize transfer checking flag
  allocate( Check_Flag(tGrid%nNodeTotal) )
  Check_Flag(:) = .false.

  ! Allocate / initialize minimum delta volume
  allocate( dVolMin(tGrid%nNodeTotal) )
  dVolMin(:) = 2d0**1023

  do iDomain = 1, nDomain
    srcs = iDomain-1

    ! Set dimension / number of nodes / cells of the buffer grid
    if( rank == iDomain-1 ) then
      bGrid%nDim       = sGrid%nDim
      bGrid%nNode      = sGrid%nNode
      bGrid%nNodeGhost = sGrid%nNodeGhost
      bGrid%nNodeTotal = sGrid%nNodeTotal
      bGrid%nCell      = sGrid%nCell
      bGrid%nCellGhost = sGrid%nCellGhost
      bGrid%nCellDummy = sGrid%nCellDummy
      bGrid%nCellTotal = sGrid%nCellTotal
      bGrid%nCellWhole = sGrid%nCellWhole
    end if

    ! Broadcast dimension / number of nodes / cells of the buffer grid
    call MPI_Bcast(bGrid%nDim,      1,MPI_INTEGER,srcs,MPI_COMM_WORLD,ier)
    call MPI_Bcast(bGrid%nNode,     1,MPI_INTEGER,srcs,MPI_COMM_WORLD,ier)
    call MPI_Bcast(bGrid%nNodeGhost,1,MPI_INTEGER,srcs,MPI_COMM_WORLD,ier)
    call MPI_Bcast(bGrid%nNodeTotal,1,MPI_INTEGER,srcs,MPI_COMM_WORLD,ier)
    call MPI_Bcast(bGrid%nCell,     1,MPI_INTEGER,srcs,MPI_COMM_WORLD,ier)
    call MPI_Bcast(bGrid%nCellGhost,1,MPI_INTEGER,srcs,MPI_COMM_WORLD,ier)
    call MPI_Bcast(bGrid%nCellDummy,1,MPI_INTEGER,srcs,MPI_COMM_WORLD,ier)
    call MPI_Bcast(bGrid%nCellTotal,1,MPI_INTEGER,srcs,MPI_COMM_WORLD,ier)
    call MPI_Bcast(bGrid%nCellWhole,1,MPI_INTEGER,srcs,MPI_COMM_WORLD,ier)

    ! Set element info and allocate node / cell data
    call SetElementInfo(bGrid)
    call AllocNodeData(bGrid)
    call AllocCellData(bGrid)

    ! Allocate / initialize data
    allocate( bNodeData(nData,bGrid%nNodeTotal) )
    bNodeData(:,:) = 0d0

    ! Allocate buffer grid / data
    allocate( Buffer_xyz(bGrid%nNodeTotal*bGrid%nDim) )
    allocate( Buffer_c2n(bGrid%nCellTotal*bGrid%nNodeCell) )
    allocate( Buffer_vol(bGrid%nCellTotal) )
    allocate( Buffer_Dat(bGrid%nNodeTotal*nData) )

    if( rank == iDomain-1 ) then

      ! Set buffer coordinates
      do iNode = 1, sGrid%nNodeTotal
        do iDim = 1, sGrid%nDim
          iBuffer = sGrid%nDim*(iNode-1)+iDim
          Buffer_xyz(iBuffer) = sGrid%xyz(iDim,iNode)
        end do
      end do

      ! Set buffer node connectivity
      do iCell = 1, sGrid%nCellTotal
        do iLocal = 1, sGrid%nNodeCell
          iBuffer = sGrid%nNodeCell*(iCell-1)+iLocal
          Buffer_c2n(iBuffer) = sGrid%c2n(iLocal,iCell)
        end do
      end do

      ! Set buffer cell volume
      do iCell = 1, sGrid%nCellTotal
        Buffer_vol(iCell) = sGrid%vol(iCell)
      end do

      ! Set buffer data
      do iNode = 1, sGrid%nNodeTotal
        do iData = 1, nData
          iBuffer = nData*(iNode-1)+iData
          Buffer_Dat(iBuffer) = sNodeData(iData,iNode)
        end do
      end do

    end if

    ! Broadcast buffer grid / data
    call MPI_Bcast(Buffer_xyz,bGrid%nNodeTotal*bGrid%nDim,     MPI_REAL8,  srcs,MPI_COMM_WORLD,ier)
    call MPI_Bcast(Buffer_c2n,bGrid%nCellTotal*bGrid%nNodeCell,MPI_INTEGER,srcs,MPI_COMM_WORLD,ier)
    call MPI_Bcast(Buffer_vol,bGrid%nCellTotal,                MPI_REAL8,  srcs,MPI_COMM_WORLD,ier)
    call MPI_Bcast(Buffer_Dat,bGrid%nNodeTotal*nData,          MPI_REAL8,  srcs,MPI_COMM_WORLD,ier)

    ! Set buffer coordinates
    do iNode = 1, bGrid%nNodeTotal
      do iDim = 1, bGrid%nDim
        iBuffer = bGrid%nDim*(iNode-1)+iDim
        bGrid%xyz(iDim,iNode) = Buffer_xyz(iBuffer)
      end do
    end do

    ! Set buffer node connectivity
    do iCell = 1, bGrid%nCellTotal
      do iLocal = 1, bGrid%nNodeCell
        iBuffer = bGrid%nNodeCell*(iCell-1)+iLocal
        bGrid%c2n(iLocal,iCell) = Buffer_c2n(iBuffer)
      end do
    end do

    ! Set buffer cell volume
    do iCell = 1, bGrid%nCellTotal
      bGrid%vol(iCell) = Buffer_vol(iCell)
    end do

    ! Set buffer data
    do iNode = 1, bGrid%nNodeTotal
      do iData = 1, nData
        iBuffer = nData*(iNode-1)+iData
        bNodeData(iData,iNode) = Buffer_Dat(iBuffer)
      end do
    end do

    ! Deallocate buffer array
    deallocate( Buffer_xyz )
    deallocate( Buffer_c2n )
    deallocate( Buffer_vol )
    deallocate( Buffer_Dat )

    ! Allocate / initialize xadj / adjncy
    ! Approximate maximum xadj = 60
    allocate( xadj(bGrid%nNodeTotal) )
    allocate( adjncy(60,bGrid%nNodeTotal) )
    xadj(:) = 0; adjncy(:,:) = 0

    ! Set xadj / adjncy data using domain cell
    do iCell = 1, bGrid%nCell
      do iLocal = 1, bGrid%nNodeCell
        iNode = bGrid%c2n(iLocal,iCell)
        xadj(iNode) = xadj(iNode) + 1
        if( xadj(iNode) > 60 ) then
          write(*,*) 'FLUS >> Number of adjacent cell is greater than upper bound'
          call MPI_Abort(MPI_COMM_WORLD,error,ier)
        end if
        adjncy(xadj(iNode),iNode) = iCell
      end do
    end do

    select case(bGrid%nDim)
    case(2)
      ! Create quadtree with buffer grid
      call MakeQuadTree(bGrid%nNodeTotal,bGrid%xyz(:,1:bGrid%nNodeTotal))
      ! Tranfer data using consistent interpolation
      do iNode = 1, tGrid%nNodeTotal
        if( Check_Flag(iNode) ) cycle
        ! Find closest cell using quadTree
        call FindQuadTree(tGrid%xyz(:,iNode),jNode)
        if( jNode == -1 ) cycle
        ! Piecewise linear interpolation
        do iadj = 1, xadj(jNode)
          iCell = adjncy(iadj,jNode)
          n1 = bGrid%c2n(1,iCell); n2 = bGrid%c2n(2,iCell); n3 = bGrid%c2n(3,iCell)
          call ComputeVolume_Triangle(bGrid%xyz(:,n2),bGrid%xyz(:,n3),tGrid%xyz(:,iNode),Vol1)
          call ComputeVolume_Triangle(bGrid%xyz(:,n1),bGrid%xyz(:,n3),tGrid%xyz(:,iNode),Vol2)
          call ComputeVolume_Triangle(bGrid%xyz(:,n1),bGrid%xyz(:,n2),tGrid%xyz(:,iNode),Vol3)
          dVol = dabs(Vol1+Vol2+Vol3-bGrid%vol(iCell))
          if( dVol < dVolMin(iNode) ) then
            tNodeData(:,iNode) = ( Vol1 * bNodeData(:,n1) + Vol2 * bNodeData(:,n2) + Vol3 * bNodeData(:,n3) ) / bGrid%vol(iCell)
            dVolMin(iNode) = dVol
            if( dVolMin(iNode) < 1d-12 ) then
              Check_Flag(iNode) = .true.;  exit
            end if
          end if
        end do
      end do
      ! Delete quadtree
      call DeleQuadTree
    case(3)
      ! Create octree with buffer grid
      call MakeOcTree(bGrid%nNodeTotal,bGrid%xyz(:,1:bGrid%nNodeTotal))
      ! Tranfer data using consistent interpolation
      do iNode = 1, tGrid%nNodeTotal
        if( Check_Flag(iNode) ) cycle
        ! Find closest cell using octree
        call FindOcTree(tGrid%xyz(:,iNode),jNode)
        if( jNode == -1 ) cycle
        ! Piecewise linear interpolation
        do iadj = 1, xadj(jNode)
          iCell = adjncy(iadj,jNode)
          n1 = bGrid%c2n(1,iCell); n2 = bGrid%c2n(2,iCell); n3 = bGrid%c2n(3,iCell); n4 = bGrid%c2n(4,iCell)
          call ComputeVolume_Tetrahedron(bGrid%xyz(:,n2),bGrid%xyz(:,n3),bGrid%xyz(:,n4),tGrid%xyz(:,iNode),Vol1)
          call ComputeVolume_Tetrahedron(bGrid%xyz(:,n1),bGrid%xyz(:,n3),bGrid%xyz(:,n4),tGrid%xyz(:,iNode),Vol2)
          call ComputeVolume_Tetrahedron(bGrid%xyz(:,n1),bGrid%xyz(:,n2),bGrid%xyz(:,n4),tGrid%xyz(:,iNode),Vol3)
          call ComputeVolume_Tetrahedron(bGrid%xyz(:,n1),bGrid%xyz(:,n2),bGrid%xyz(:,n3),tGrid%xyz(:,iNode),Vol4)
          dVol = dabs(Vol1+Vol2+Vol3+Vol4-bGrid%vol(iCell))
          if( dVol < dVolMin(iNode) ) then
            tNodeData(:,iNode) = ( Vol1 * bNodeData(:,n1) + Vol2 * bNodeData(:,n2) + Vol3 * bNodeData(:,n3) +  Vol4 * bNodeData(:,n4) ) / bGrid%vol(iCell)
            dVolMin(iNode) = dVol
            if( dVolMin(iNode) < 1d-12 ) then
              Check_Flag(iNode) = .true.;  exit
            end if
          end if
        end do
      end do
      ! Delete octree
      call DeleOcTree
    end select

    ! Delete buffer grid
    call DelGrid(bGrid)

    ! Deallocate adjacent cell data
    deallocate( xadj, adjncy )

    ! Deallocate buffer data
    deallocate( bNodeData )

  end do

  ! Check number of approximated nodes
  nApprox = 0
  do iNode = 1, tGrid%nNodeTotal
    if( .not. Check_Flag(iNode) ) nApprox = nApprox + 1
  end do
  if( rank == 0 ) then
    write(*,*) 'FLUS >> Check approximated node during data transfer'
    write(*,*) 'FLUS >> Rank / nApproximation / MaxVolumeDiff'
  end if
  do iDomain = 1, nDomain
    if( rank == iDomain-1 .and. nApprox > 0 ) write(*,*) rank, nApprox, maxval(dVolMin(:))
    call MPI_Barrier(MPI_COMM_WORLD,ier)
  end do

  ! Set target cell data
  call SetCellData(tGrid,tNodeData,tData)

  ! Deallocate source / target node data
  deallocate( sNodeData )
  deallocate( tNodeData )

  ! Deallocate checking flag /  minimum delta volume
  deallocate( Check_Flag )
  deallocate( dVolMin )

  ! MPI barrier
  call MPI_Barrier(MPI_COMM_WORLD,ier)

end subroutine

subroutine SetNodeData(Grid,CellData,NodeData)
  type(t_Grid), intent(in) :: Grid
  real(8), intent(in) :: CellData(:,:)
  real(8), intent(out) :: NodeData(:,:)

  integer :: iLocal, iCell, iNode
  real(8) :: Coeff, PartVol
  real(8), allocatable :: NodeVol(:)

  ! Allocate / initialize node volume
  allocate( NodeVol(Grid%nNodeTotal) )
  NodeVol(:)  = 0d0

  ! Initialize node data
  NodeData(:,:) = 0d0

  ! Sum data with volume weighting
  Coeff = 1d0 / dble(Grid%nNodeCell)
  do iCell = 1, Grid%nCellTotal
    do iLocal = 1, Grid%nNodeCell
      iNode = Grid%c2n(iLocal,iCell)
      PartVol = Grid%vol(iCell) * Coeff
      NodeVol(iNode) = NodeVol(iNode) + PartVol
      NodeData(:,iNode) = NodeData(:,iNode) + CellData(:,iCell) * PartVol
    end do
  end do

  ! Normalize node data
  do iNode = 1, Grid%nNodeTotal
    NodeData(:,iNode) = NodeData(:,iNode) / NodeVol(iNode)
  end do

  ! Deallocate node volume
  deallocate( NodeVol )

end subroutine

subroutine SetCellData(Grid,NodeData,CellData)
  type(t_Grid), intent(in) :: Grid
  real(8), intent(in) :: NodeData(:,:)
  real(8), intent(out) :: CellData(:,:)

  integer :: iCell, iNode, iLocal
  real(8) :: Coeff

  ! Initialize cell data
  CellData(:,:) = 0d0

  ! Sum data with volume weighting and normalize cell data
  Coeff = 1d0 / dble(Grid%nNodeCell)
  do iCell = 1, Grid%nCellTotal
    do iLocal = 1, Grid%nNodeCell
      iNode = Grid%c2n(iLocal,iCell)
      CellData(:,iCell) = CellData(:,iCell) + NodeData(:,iNode)
    end do
    CellData(:,iCell) = CellData(:,iCell) * Coeff
  end do

end subroutine

subroutine ComputeVolume_Triangle(xyz1,xyz2,xyz3,vol)
  real(8), intent(in) :: xyz1(2), xyz2(2), xyz3(2)
  real(8), intent(out) :: vol

  vol = 0.5d0 * dabs( ( xyz3(1) - xyz1(1) ) * ( xyz2(2) - xyz1(2) ) - ( xyz3(2) - xyz1(2) ) * ( xyz2(1)-xyz1(1) ) )

end subroutine

subroutine ComputeVolume_Tetrahedron(xyz1,xyz2,xyz3,xyz4,vol)
  real(8), intent(in) :: xyz1(3), xyz2(3), xyz3(3), xyz4(3)
  real(8), intent(out) :: vol

  real(8) :: vec1(3), vec2(3)
  real(8) :: vec3(3), vec4(3)

  ! Set vec1 / vec2 / vec3
  vec1(:) = xyz1(:) - xyz2(:)
  vec2(:) = xyz1(:) - xyz3(:)
  vec3(:) = xyz1(:) - xyz4(:)

  ! vec4 = vec1 X vec2
  vec4(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
  vec4(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
  vec4(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)

  ! vol = ABS( vec3 * vec4 ) / 6
  vol = dabs( vec3(1) * vec4(1) + vec3(2) * vec4(2) + vec3(3) * vec4(3) ) / 6d0

end subroutine

end module
