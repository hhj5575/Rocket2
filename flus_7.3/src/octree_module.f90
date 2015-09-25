module octree_module

!Octree data structure
!y
!-------------!-------------!
!             !             !
!             !             !
!      7      !      8      !
!             !             !
!             !             !
!-------------!-------------!
!             !             !
!             !             !
!      5      !      6      !
!             !             !
!             !             !
!z=1----------!-------------!
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
!z=0----------!-------------!x

implicit none; private
public MakeOcTree
public FindOcTree
public DeleOcTree

! Octree data type
type t_Tree
  integer :: Level
  integer :: TreeNum
  integer :: indTree
  integer, allocatable :: ElemInfo(:)
  type(t_Tree), pointer :: Root => Null()
  type(t_Tree), pointer :: SubTree(:) => Null()
  real(8) :: xBound(2), yBound(2), zBound(2)
end type t_Tree

! Octree data
type(t_Tree), target  :: Tree

! Pointer to source info
integer, pointer :: nElem
integer, pointer :: nNode
real(8), pointer :: xyz(:,:)
integer, pointer :: nCell
integer, pointer :: c2n(:,:)
logical :: Elemtype

contains

subroutine MakeOcTree(snNode,sxyz,snCell,sc2n)
  integer, intent(in), target :: snNode
  real(8), intent(in), target :: sxyz(3,*)
  integer, intent(in), target :: snCell
  integer, intent(in), target :: sc2n(4,*)
  optional :: snCell, sc2n
  
  type(t_Tree), pointer :: pTree
  type(t_Tree), pointer :: pSubTree
  integer :: Node, Elem, iElem
  integer :: Sub, tSub
  integer :: nxsgn, nysgn, nzsgn
  real(8) :: xshift, yshift, zshift
  real(8) :: xsize, ysize, zsize, lsize
  real(8) :: Cordi(3), iCordi(3)
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
  allocate( Tree%ElemInfo(8) )
  Tree%ElemInfo(:) = 0
  Tree%xBound(:) = xyz(1,1)
  Tree%yBound(:) = xyz(2,1)
  Tree%zBound(:) = xyz(3,1)
  do Node = 2, nNode
    Tree%xBound(1) = dmin1( Tree%xBound(1), xyz(1,Node) )
    Tree%xBound(2) = dmax1( Tree%xBound(2), xyz(1,Node) )
    Tree%yBound(1) = dmin1( Tree%yBound(1), xyz(2,Node) )
    Tree%yBound(2) = dmax1( Tree%yBound(2), xyz(2,Node) )
    Tree%zBound(1) = dmin1( Tree%zBound(1), xyz(3,Node) )
    Tree%zBound(2) = dmax1( Tree%zBound(2), xyz(3,Node) )
  end do
  
  ! Set initial bound 1% larger
  lsize = dsqrt( ( Tree%xBound(2) - Tree%xBound(1) )**2 &
               + ( Tree%yBound(2) - Tree%yBound(1) )**2 &
               + ( Tree%zBound(2) - Tree%zBound(1) )**2 ) * 0.01d0
  Tree%xBound(1) = Tree%xBound(1) - lsize
  Tree%xBound(2) = Tree%xBound(2) + lsize
  Tree%yBound(1) = Tree%yBound(1) - lsize
  Tree%yBound(2) = Tree%yBound(2) + lsize
  Tree%zBound(1) = Tree%zBound(1) - lsize
  Tree%zBound(2) = Tree%zBound(2) + lsize
  
  ! Set tree
  do Elem = 1, nElem
  
    ! Set element coordinate
    call Coordinate(Elem,Cordi(:))
  
    ! Set root tree
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
        zshift = Cordi(3) - 0.5d0 * ( pTree%zBound(1) + pTree%zBound(2) )
        zsize  = dmax1(tol, dabs(zshift))
        nzsgn  = ( idint( dble( tolp1 * zshift / zsize ) ) + 1 ) / 2 * 2
        tSub   = 1 + nxsgn / 2 + nysgn + nzsgn * 2
        pTree  => pTree%SubTree(tSub)
      end do
  
      ! Set new subtree data
      if( pTree%indTree == 8 ) then
        pTree%indTree = -1
        allocate( pTree%SubTree(8) )
        do Sub = 1, 8
          pSubTree => pTree%SubTree(Sub)
          pSubTree%Level   = pTree%Level + 1
          pSubTree%TreeNum = Sub
          pSubTree%indTree = 0
          allocate( pSubTree%ElemInfo(8) )
          pSubTree%ElemInfo(:) = 0
          pSubTree%Root => pTree
          nxsgn = 2 - mod(Sub,2)
          pSubTree%xBound(  nxsgn) = pTree%xBound(nxsgn)
          pSubTree%xBound(3-nxsgn) = 0.5d0 * ( pTree%xBound(1) + pTree%xBound(2) )
          nysgn = 2 - mod((Sub+1)/2,2)
          pSubTree%yBound(  nysgn) = pTree%yBound(nysgn)
          pSubTree%yBound(3-nysgn) = 0.5d0 * ( pTree%yBound(1) + pTree%yBound(2) )
          nzsgn = 1 + ( isign(1,2*Sub-9) + 1 ) / 2
          pSubTree%zBound(  nzsgn) = pTree%zBound(nzsgn)
          pSubTree%zBound(3-nzsgn) = 0.5d0 * ( pTree%zBound(1) + pTree%zBound(2) )
        end do
        ! Assign inode to new subtree
        do Sub = 1, 8
          iElem = pTree%ElemInfo(Sub)
          call COORDINATE(iElem,iCordi(:))
          xshift = iCordi(1) - 0.5d0 * ( pTree%xBound(1) + pTree%xBound(2) )
          xsize  = dmax1(tol, dabs(xshift))
          nxsgn  = ( idint( dble( tolp1 * xshift / xsize ) ) + 1 ) / 2 * 2
          yshift = iCordi(2) - 0.5d0 * ( pTree%yBound(1) + pTree%yBound(2) )
          ysize  = dmax1(tol, dabs(yshift))
          nysgn  = ( idint( dble( tolp1 * yshift / ysize ) ) + 1 ) / 2 * 2
          zshift = iCordi(3) - 0.5d0 * ( pTree%zBound(1) + pTree%zBound(2) )
          zsize  = dmax1(tol, dabs(zshift))
          nzsgn  = ( idint( dble( tolp1 * zshift / zsize ) ) + 1 ) / 2 * 2
          tSub   = 1 + nxsgn / 2 + nysgn + nzsgn * 2
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
  real(8), intent(out) :: Cordi(3)
  
  integer :: n1, n2, n3, n4
  real(8), parameter :: inv4  = 0.25d0
  
  ! Set elem coordinate
  select case(Elemtype)
  case(.false.)
    Cordi(:) = xyz(:,Elem)
  case(.true.)
    n1 = c2n(1,Elem)
    n2 = c2n(2,Elem)
    n3 = c2n(3,Elem)
    n4 = c2n(4,Elem)
    Cordi(:) = inv4 * ( xyz(:,n1) + xyz(:,n2) + xyz(:,n3) + xyz(:,n4) )
  end select
  
end subroutine
  
subroutine FindOcTree(ixyz,iElem)
  real(8), intent(in) :: ixyz(3)
  integer, intent(out) :: iElem
  
  type(t_Tree), pointer :: pTree
  integer :: SubNum
  integer :: nxsgn, nysgn, nzsgn
  real(8) :: xshift, yshift, zshift
  real(8) :: xsize, ysize, zsize
  real(8), parameter :: tol   = 1d-12
  real(8), parameter :: tolp1 = 1d0 + tol
  
  ! Set root tree
  pTree => Tree
  
  ! Check whether ixyz is located outside the highest tree
  if( ixyz(1) < pTree%xBound(1) .or. ixyz(2) < pTree%yBound(1) .or. ixyz(3) < pTree%zBound(1) .or. &
      ixyz(1) > pTree%xBound(2) .or. ixyz(2) > pTree%yBound(2) .or. ixyz(3) > pTree%zBound(2) ) then
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
    zshift = ixyz(3) - 0.5d0 * ( pTree%zBound(1) + pTree%zBound(2) )
    zsize  = dmax1(tol, dabs(zshift))
    nzsgn  = ( idint( dble( tolp1 * zshift / zsize ) ) + 1 ) / 2 * 2
    SubNum = 1 + nxsgn / 2 + nysgn + nzsgn * 2
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
  real(8), intent(in) :: ixyz(3)
  integer, intent(out) :: iNode
  
  integer :: Node
  integer :: Sub, iSub
  real(8) :: Dist, DistMin
  real(8) :: xyzl(3), xyzu(3)
  
  ! Find closest node at the lowest subtree
  DistMin = 2d0**1023
  do Sub = 1, pTree%indTree
    Node = pTree%ElemInfo(Sub)
    Dist = DSQRT( SUM( ( ixyz(:) - xyz(:,Node) )**2 ) )
    if( Dist < DistMin ) then
      iNode = Node
      DistMin = Dist
    end if
  end do
  
  ! Check root tree & find closest node
  xyzl(:) = ixyz(:) - DistMin
  xyzu(:) = ixyz(:) + DistMin
  do while( xyzl(1) < pTree%xBound(1) .or. xyzl(2) < pTree%yBound(1) .or. xyzl(3) < pTree%zBound(1) .or. &
            xyzu(1) > pTree%xBound(2) .or. xyzu(2) > pTree%yBound(2) .or. xyzu(3) > pTree%zBound(2) )
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
  real(8), intent(in) :: ixyz(3)
  real(8), intent(inout) :: DistMin
  integer, intent(out) :: iNode
  
  type(t_Tree), pointer :: pTree
  integer :: Node
  integer :: Sub1, Sub2
  real(8) :: Dist
  
  ! Find subtree recursively
  do Sub1 = 1, 8
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
  real(8), intent(in) :: ixyz(3)
  integer, intent(out) :: iCell
  
  integer :: Cell
  integer :: n1, n2, n3, n4
  integer :: Sub, iSub
  real(8) :: Vol, dVol, dVolMin
  real(8) :: Vol1, Vol2, Vol3, Vol4
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
    n4 = c2n(4,Cell)
    call CalVol(xyz(:,n1), xyz(:,n2), xyz(:,n3), xyz(:,n4), Vol )
    call CalVol(xyz(:,n2), xyz(:,n3), xyz(:,n4), ixyz(:)  , Vol1)
    call CalVol(xyz(:,n1), xyz(:,n3), xyz(:,n4), ixyz(:)  , Vol2)
    call CalVol(xyz(:,n1), xyz(:,n2), xyz(:,n4), ixyz(:)  , Vol3)
    call CalVol(xyz(:,n1), xyz(:,n2), xyz(:,n3), ixyz(:)  , Vol4)
    dVol = dabs( Vol1 + Vol2 + Vol3 + Vol4 - Vol )
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
  real(8), intent(in) :: ixyz(3)
  real(8), intent(inout) :: dVolMin
  integer, intent(out) :: iCell
  logical, intent(inout) :: Flag
  
  type(t_Tree), pointer :: pTree
  integer :: Cell
  integer :: n1, n2, n3, n4
  integer :: Sub1, Sub2
  real(8) :: Vol, dVol
  real(8) :: Vol1, Vol2, Vol3, Vol4
  real(8), parameter :: tol = 1d-12
  
  ! Find subtree recursively
  do Sub1 = 1, 8
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
        n4 = c2n(4,Cell)
        call CalVol(xyz(:,n1), xyz(:,n2), xyz(:,n3), xyz(:,n4), Vol )
        call CalVol(xyz(:,n2), xyz(:,n3), xyz(:,n4), ixyz(:)  , Vol1)
        call CalVol(xyz(:,n1), xyz(:,n3), xyz(:,n4), ixyz(:)  , Vol2)
        call CalVol(xyz(:,n1), xyz(:,n2), xyz(:,n4), ixyz(:)  , Vol3)
        call CalVol(xyz(:,n1), xyz(:,n2), xyz(:,n3), ixyz(:)  , Vol4)
        dVol = dabs( Vol1 + Vol2 + Vol3 + Vol4 - Vol )
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
  
subroutine CalVol(xyz1,xyz2,xyz3,xyz4,vol)
  real(8), intent(in) :: xyz1(3)
  real(8), intent(in) :: xyz2(3)
  real(8), intent(in) :: xyz3(3)
  real(8), intent(in) :: xyz4(3)
  real(8), intent(out) :: vol
  
  real(8) :: vec1(3), vec2(3)
  real(8) :: vec3(3), vec4(3)
  
  ! Vector
  vec1(:) = xyz1(:) - xyz2(:)
  vec2(:) = xyz1(:) - xyz3(:)
  vec3(:) = xyz1(:) - xyz4(:)
  
  ! vec4 = vec1 X vec2
  vec4(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
  vec4(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
  vec4(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)
  
  ! vol = abs( vec3 * vec4 ) / 6
  vol = dabs( vec3(1) * vec4(1) + vec3(2) * vec4(2) + vec3(3) * vec4(3) ) / 6d0
  
end subroutine
  
subroutine DeleOcTree
  
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
  
  do Sub = 1, 8
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
