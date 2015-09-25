module quadtree_module

! Quadtree data structure
!y
!-------------!-------------!
!             !             !
!             !             !
!      3      !      4      !
!             !             !
!             !             !
!-------------!-------------!
!             !             !
!             !             !
!      1      !      2      !
!             !             !
!             !             !
!-------------!-------------!x

implicit none; private
public MakeQuadTree
public FindQuadTree
public DeleQuadTree

! Quadtree data type
type t_Tree
  integer :: Level
  integer :: TreeNum
  integer :: indTree
  integer, allocatable :: ElemInfo(:)
  type(t_Tree), pointer :: Root => Null()
  type(t_Tree), pointer :: SubTree(:) => Null()
  real(8) :: xBound(2), yBound(2)
end type t_Tree

! Quadtree data
type(t_Tree), target  :: Tree

! Pointer to source info
integer, pointer :: nElem
integer, pointer :: nNode
real(8), pointer :: xyz(:,:)
integer, pointer :: nCell
integer, pointer :: c2n(:,:)
logical :: Elemtype

contains

subroutine MakeQuadTree(snNode,sxyz,snCell,sc2n)
  integer, intent(in), target :: snNode
  real(8), intent(in), target :: sxyz(2,*)
  integer, intent(in), target :: snCell
  integer, intent(in), target :: sc2n(3,*)
  optional :: snCell, sc2n

  type(t_Tree), pointer :: pTree
  type(t_Tree), pointer :: pSubTree
  integer :: Node, Elem, iElem
  integer :: Sub, tSub
  integer :: nxsgn, nysgn
  real(8) :: xshift, yshift
  real(8) :: xsize, ysize, lsize
  real(8) :: Cordi(2), iCordi(2)
  logical :: Flag
  real(8), parameter :: tol   = 1d-12
  real(8), parameter :: tolp1 = 1d0 + tol

  ! Set pointer ( source node info )
  nNode => snNode
  xyz   => sxyz(:,1:snNode)
  nElem => nNode
  Elemtype = .false.
  if( present(snCell) ) then
    nCell => snCell
    c2n   => sc2n(:,1:snCell)
    nElem => nCell
    Elemtype = .true.
  end if

  ! Set initial tree & bound ( with node data )
  Tree%Level   = 1
  Tree%TreeNum = 0
  Tree%indTree = 0
  allocate( Tree%ElemInfo(4) )
  Tree%ElemInfo(:) = 0
  Tree%xBound(:) = xyz(1,1)
  Tree%yBound(:) = xyz(2,1)
  do Node = 2, nNode
    Tree%xBound(1) = dmin1( Tree%xBound(1), xyz(1,Node) )
    Tree%xBound(2) = dmax1( Tree%xBound(2), xyz(1,Node) )
    Tree%yBound(1) = dmin1( Tree%yBound(1), xyz(2,Node) )
    Tree%yBound(2) = dmax1( Tree%yBound(2), xyz(2,Node) )
  end do

  ! Set initial bound 1% larger
  lsize = dsqrt( ( Tree%xBound(2) - Tree%xBound(1) )**2 &
               + ( Tree%yBound(2) - Tree%yBound(1) )**2 ) * 0.01d0
  Tree%xBound(1) = Tree%xBound(1) - lsize
  Tree%xBound(2) = Tree%xBound(2) + lsize
  Tree%yBound(1) = Tree%yBound(1) - lsize
  Tree%yBound(2) = Tree%yBound(2) + lsize

  ! Set tree
  do Elem = 1, nElem

    ! Set element coordinate
    call Coordinate(Elem,Cordi(:))

    ! Set Root Tree
    pTree => Tree

    Flag = .false.
    do while( Flag == .false. )

      ! Find subtree
      do while( pTree%indTree < 0 )
        xshift = Cordi(1) - 0.5d0 * ( pTree%xBound(1) + pTree%xBound(2) )
        xsize  = dmax1(tol, dabs(xshift))
        nxsgn  = ( idint( dble( tolp1 * xshift / xsize ) ) + 1 ) / 2 * 2
        yshift = Cordi(2) - 0.5d0 * ( pTree%yBound(1) + pTree%yBound(2) )
        ysize  = dmax1(tol, dabs(yshift))
        nysgn  = ( idint( dble( tolp1 * yshift / ysize ) ) + 1 ) / 2 * 2
        tSub   = 1 + nxsgn / 2 + nysgn
        pTree  => pTree%SubTree(tSub)
      end do

      ! Set new subtree data
      if( pTree%indTree == 4 ) then
        pTree%indTree = -1
        allocate( pTree%SubTree(4) )
        do Sub = 1, 4
          pSubTree => pTree%SubTree(Sub)
          pSubTree%Level   = pTree%Level + 1
          pSubTree%TreeNum = Sub
          pSubTree%indTree = 0
          allocate( pSubTree%ElemInfo(4) )
          pSubTree%ElemInfo(:) = 0
          pSubTree%Root => pTree
          nxsgn = 2 - mod(Sub,2)
          pSubTree%xBound(  nxsgn) = pTree%xBound(nxsgn)
          pSubTree%xBound(3-nxsgn) = 0.5d0 * ( pTree%xBound(1) + pTree%xBound(2) )
          nysgn = 1 + ( isign(1,2*Sub-5) + 1 ) / 2
          pSubTree%yBound(  nysgn) = pTree%yBound(nysgn)
          pSubTree%yBound(3-nysgn) = 0.5d0 * ( pTree%yBound(1) + pTree%yBound(2) )
        end do
        ! Assign iNode to new subtree
        do Sub = 1, 4
          iElem  = pTree%ElemInfo(Sub)
          call Coordinate(iElem,iCordi(:))
          xshift = iCordi(1) - 0.5d0 * ( pTree%xBound(1) + pTree%xBound(2) )
          xsize  = dmax1(tol, dabs(xshift))
          nxsgn  = ( idint( dble( tolp1 * xshift / xsize ) ) + 1 ) / 2 * 2
          yshift = iCordi(2) - 0.5d0 * ( pTree%yBound(1) + pTree%yBound(2) )
          ysize  = dmax1(tol, dabs(yshift))
          nysgn  = ( idint( dble( tolp1 * yshift / ysize ) ) + 1 ) / 2 * 2
          tSub   = 1 + nxsgn / 2 + nysgn
          pSubTree => pTree%SubTree(tSub)
          pSubTree%indTree = pSubTree%indTree + 1
          pSubTree%ElemInfo(pSubTree%indTree) = iElem
        end do
        deallocate( pTree%ElemInfo )
      else
        Flag = .true.
      end if

    end do

    ! Assign node to tree
    pTree%indTree = pTree%indTree + 1
    pTree%ElemInfo(pTree%indTree) = Elem

  end do

end subroutine

subroutine Coordinate(Elem,Cordi)
  integer, intent(in) :: Elem
  real(8), intent(out) :: Cordi(2)

  integer :: n1, n2, n3
  real(8), parameter :: inv3  = 1d0 / 3d0

  ! Set Elem Coordinate
  select case(Elemtype)
  case(.false.)
    Cordi(:) = xyz(:,Elem)
  case(.true.)
    n1 = c2n(1,Elem)
    n2 = c2n(2,Elem)
    n3 = c2n(3,Elem)
    Cordi(:) = inv3 * ( xyz(:,n1) + xyz(:,n2) + xyz(:,n3) )
  end select

end subroutine

subroutine FindQuadTree(ixyz,iElem)
  real(8), intent(in) :: ixyz(2)
  integer, intent(out) :: iElem

  type(t_Tree), pointer :: pTree
  integer :: SubNum
  integer :: nxsgn, nysgn
  real(8) :: xshift, yshift
  real(8) :: xsize, ysize
  real(8), parameter :: tol   = 1d-12
  real(8), parameter :: tolp1 = 1d0 + tol

  ! Set root tree
  pTree => Tree

  ! Check whether ixyz is located outside the highest tree
  if( ixyz(1) < pTree%xBound(1) .or. ixyz(2) < pTree%yBound(1) .or. &
      ixyz(1) > pTree%xBound(2) .or. ixyz(2) > pTree%yBound(2) ) then
     iElem = -1
     return
  end if

  ! Find lowest level subtree
  do while( pTree%indTree < 0 )
    xshift = ixyz(1) - 0.5d0 * ( pTree%xBound(1) + pTree%xBound(2) )
    xsize  = dmax1(tol, dabs(xshift))
    nxsgn  = ( idint( dble( tolp1 * xshift / xsize ) ) + 1 ) / 2 * 2
    yshift = ixyz(2) - 0.5d0 * ( pTree%yBound(1) + pTree%yBound(2) )
    ysize  = dmax1(tol, dabs(yshift))
    nysgn  = ( idint( dble( tolp1 * yshift / ysize ) ) + 1 ) / 2 * 2
    SubNum = 1 + nxsgn / 2 + nysgn
    pTree  => pTree%SubTree(SubNum)
  end do

  ! Find closest node or over cell
  select case(Elemtype)
  case(.false.)
    call FindNode(pTree,ixyz,iElem)
  case(.true.)
    call FindCell(pTree,ixyz,iElem)
  end select

end subroutine

subroutine FindNode(pTree,ixyz,iNode)
  type(t_Tree), intent(inout), pointer :: pTree
  real(8), intent(in) :: ixyz(2)
  integer, intent(out) :: iNode

  integer :: Node
  integer :: Sub, iSub
  real(8) :: Dist, DistMin
  real(8) :: xyzl(2), xyzu(2)

  ! Find closest node at the lowest subtree
  DistMin = 2d0**1023
  do Sub = 1, pTree%indTree
    Node = pTree%ElemInfo(Sub)
    Dist = dsqrt( sum( ( ixyz(:) - xyz(:,Node) )**2 ) )
    if( Dist < DistMin ) then
      iNode = Node
      DistMin = Dist
    end if
  end do

  ! Check root tree & find closest node
  xyzl(:) = ixyz(:) - DistMin
  xyzu(:) = ixyz(:) + DistMin
  do while( xyzl(1) < pTree%xBound(1) .or. xyzl(2) < pTree%yBound(1) .or. &
            xyzu(1) > pTree%xBound(2) .or. xyzu(2) > pTree%yBound(2) )
    if( pTree%Level == 1 ) exit
    iSub = pTree%TreeNum
    pTree => pTree%Root
    call FindNode_SubTree(pTree,iSub,ixyz(:),DistMin,iNode)
    xyzl(:) = ixyz(:) - DistMin
    xyzu(:) = ixyz(:) + DistMin
  end do

end subroutine

recursive subroutine FindNode_SubTree(pRoot,iSub,ixyz,DistMin,iNode)
  type(t_Tree), intent(in) :: pRoot
  integer, intent(in) :: iSub
  real(8), intent(in) :: ixyz(2)
  real(8), intent(inout) :: DistMin
  integer, intent(out) :: iNode

  type(t_Tree), pointer :: pTree
  integer :: Node
  integer :: Sub1, Sub2
  real(8) :: Dist

  ! Find subtree recursively
  do Sub1 = 1, 4
    if( Sub1 == iSub ) cycle
    pTree => pRoot%SubTree(Sub1)
    select case(pTree%indTree)
    case(-1)
      call FindNode_SubTree(pTree,0,ixyz(:),DistMin,iNode)
    case default
      ! Find closest node at the lowest subtree
      do Sub2 = 1, pTree%indTree
        Node = pTree%ElemInfo(Sub2)
        Dist = dsqrt( sum( ( ixyz(:) - xyz(:,Node) )**2 ) )
        if( Dist < DistMin ) then
          iNode = Node
          DistMin = Dist
        end if
      end do
    end select
  end do

end subroutine

subroutine FindCell(pTree,ixyz,iCell)
  type(t_Tree), intent(inout), pointer :: pTree
  real(8), intent(in) :: ixyz(2)
  integer, intent(out) :: iCell

  integer :: Cell
  integer :: n1, n2, n3
  integer :: Sub, iSub
  real(8) :: Vol, dVol, dVolMin
  real(8) :: Vol1, Vol2, Vol3
  logical :: Flag
  real(8), parameter :: tol = 1d-12

  ! Find over cell at the lowest subtree
  Flag = .false.
  dVolMin = 2d0**1023
  do Sub = 1, pTree%indTree
    Cell = pTree%ElemInfo(Sub)
    n1 = c2n(1,Cell)
    n2 = c2n(2,Cell)
    n3 = c2n(3,Cell)
    call CalVol(xyz(:,n1), xyz(:,n2), xyz(:,n3), Vol )
    call CalVol(xyz(:,n2), xyz(:,n3), ixyz(:)  , Vol1)
    call CalVol(xyz(:,n1), xyz(:,n3), ixyz(:)  , Vol2)
    call CalVol(xyz(:,n1), xyz(:,n2), ixyz(:)  , Vol3)
    dVol = dabs( Vol1 + Vol2 + Vol3 - Vol )
    if( dVol < tol ) then
      iCell = Cell
      Flag = .true.
      exit
    else if( dVol < dVolMin ) then
      iCell = Cell
      dVolMin = dVol
    end if
  end do

  ! Check root tree & find over cell
  do while( Flag == .false. )
    if( pTree%Level == 1 ) exit
    iSub = pTree%TreeNum
    pTree => pTree%Root
    call FindCell_SubTree(pTree,iSub,ixyz(:),dVolMin,iCell,Flag)
  end do

end subroutine

recursive subroutine FindCell_SubTree(pRoot,iSub,ixyz,dVolMin,iCell,Flag)
  type(t_Tree), intent(in) :: pRoot
  integer, intent(in) :: iSub
  real(8), intent(in) :: ixyz(2)
  real(8), intent(inout) :: dVolMin
  integer, intent(out) :: iCell
  logical, intent(inout) :: Flag

  type(t_Tree), pointer :: pTree
  integer :: Cell
  integer :: n1, n2, n3
  integer :: Sub1, Sub2
  real(8) :: Vol, dVol
  real(8) :: Vol1, Vol2, Vol3
  real(8), parameter :: tol = 1d-12

  ! Find subtree recursively
  do Sub1 = 1, 4
    if( Sub1 == iSub ) cycle
    pTree => pRoot%SubTree(Sub1)
    select case(pTree%indTree)
    case(-1)
      call FindCell_SubTree(pTree,0,ixyz(:),dVolMin,iCell,Flag)
    case default
      ! Find closest node at the lowest subtree
      do Sub2 = 1, pTree%indTree
        Cell = pTree%ElemInfo(Sub2)
        n1 = c2n(1,Cell)
        n2 = c2n(2,Cell)
        n3 = c2n(3,Cell)
        call CalVol(xyz(:,n1), xyz(:,n2), xyz(:,n3), Vol )
        call CalVol(xyz(:,n2), xyz(:,n3), ixyz(:)  , Vol1)
        call CalVol(xyz(:,n1), xyz(:,n3), ixyz(:)  , Vol2)
        call CalVol(xyz(:,n1), xyz(:,n2), ixyz(:)  , Vol3)
        dVol = dabs( Vol1 + Vol2 + Vol3 - Vol )
        if( dVol < tol ) then
          iCell = Cell
          Flag = .true.
          exit
        else if( dVol < dVolMin ) then
          iCell = Cell
          dVolMin = dVol
        end if
      end do
    end select
    if( Flag == .true. ) exit
  end do

end subroutine

subroutine CalVol(xyz1,xyz2,xyz3,vol)
  real(8), intent(in) :: xyz1(2)
  real(8), intent(in) :: xyz2(2)
  real(8), intent(in) :: xyz3(2)
  real(8), intent(out) :: vol

  vol = 0.5d0 * dabs( ( xyz3(1) - xyz1(1) ) * ( xyz2(2) - xyz1(2) ) - ( xyz3(2) - xyz1(2) ) * ( xyz2(1)-xyz1(1) ) )

end subroutine

subroutine DeleQuadTree

  call DeleTree_SubTree(Tree)
  nullify( nNode )
  nullify( xyz )
  if( Elemtype == .true. ) then
    nullify( nCell )
    nullify( c2n )
  end if

end subroutine

recursive subroutine DeleTree_SubTree(pRoot)
  type(t_Tree), intent(inout) :: pRoot

  type(t_Tree), pointer :: pTree
  integer :: Sub

  do Sub = 1, 4
    pTree => pRoot%SubTree(Sub)
    select case(pTree%indTree)
    case(-1)
      call DeleTree_SubTree(pTree)
    case default
      ! Deallocate ElemInfo
      deallocate( pTree%ElemInfo )
    end select
  end do
  deallocate( pRoot%SubTree )

end subroutine

end module
