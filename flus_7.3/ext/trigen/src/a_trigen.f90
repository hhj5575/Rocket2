SUBROUTINE TRIGEN(nBNode, nBFace, Bxyz, Bf2n, BPatch, nNode1, nCell1, nFace1, xyz1, c2n1, f2n1, Patch1)

IMPLICIT NONE
INTEGER, INTENT(in   ) :: nBNode, nBFace
REAL(8), INTENT(in   ) :: Bxyz(2,nBNode)
INTEGER, INTENT(in   ) :: Bf2n(2,nBFace)
INTEGER, INTENT(in   ) :: BPatch(nBFace)
INTEGER, INTENT(  out) :: nNode1, nCell1, nFace1
REAL(8), INTENT(  out) :: xyz1(2,*)
INTEGER, INTENT(  out) :: c2n1(3,*)
INTEGER, INTENT(  out) :: f2n1(2,*)
INTEGER, INTENT(  out) :: Patch1(*)

INTEGER :: Node, Cell, Local
INTEGER, ALLOCATABLE :: NodeIndex(:)

ALLOCATE( NodeIndex(nBNode) )
  
WRITE(*,*) 'Allocating trigen internal variables'
CALL allocvars

WRITE(*,*) 'Reading the list of nodes defining inner and outer boundaries'
CALL input(nBNode,nBFace,Bxyz,Bf2n,NodeIndex)

WRITE(*,*) 'Setup point list and edge data structure for boundary nodes'
CALL initbd

WRITE(*,*) 'Connecting points to generate the Delaunay triangulation'
CALL connct

WRITE(*,*) 'Initializing quadtree data structure'
CALL initqd

WRITE(*,*) 'Creating points inside the domain using Rebay algorithm'
CALL rebay

WRITE(*,*) 'Wrting grid (nodes/cells/faces)'
CALL grid(nNode1,nCell1,xyz1,c2n1)

! Reset node coordinates
xyz1(:,1:nBNode) = Bxyz(:,:)

! Reset node connectivity of the cell
DO Cell = 1, nCell1
  DO Local = 1, 3
    Node = c2n1(Local,Cell)
    IF( Node <= nBNode ) c2n1(Local,Cell) = NodeIndex(Node)
  END DO
END DO

! Just copy face information
nFace1 = nBFace
f2n1(:,1:nFace1) = Bf2n(:,1:nBFace)
Patch1(1:nFace1) = BPatch(1:nBFace)

WRITE(*,*) 'Deallcating trigen internal variables'
CALL dealloc

DEALLOCATE( NodeIndex )

WRITE(*,*)
WRITE(*,*) 'Statistics:'
WRITE(*,*)
WRITE(*,*) '  Input points:', nBNode
WRITE(*,*) '  Input faces :', nBFace
WRITE(*,*) '  Input holes :', 0
WRITE(*,*)
WRITE(*,*) '  Grid points :', nNode1
WRITE(*,*) '  Grid cells  :', nCell1
WRITE(*,*) '  Grid faces  :', nFace1
WRITE(*,*)

END SUBROUTINE
