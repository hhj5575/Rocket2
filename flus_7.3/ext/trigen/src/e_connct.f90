!---------------------------------------------------------------------!
!                                                                     !
!    CONNECTS POINTS TO GENERATE THE DELAUNAY TRIANGULATION OF THE    !
!    DOMAIN BETWEEN THE INNER AND OUTER BOUNDARIES.                   !
!    THE ALGORITHM PROCEEDS BY CONSTRUCTING THE TRIANGLES ON A        !
!    FRONT DEFINED BY A LIST OF EDGES, STARTING FROM THE BOUNDARY     !
!    EDGES AND MOVING INTO THE DOMAIN UNTIL ALL TRIANGLES ARE IN      !
!    PLACE.                                                           !
!                                                                     !
!---------------------------------------------------------------------!

subroutine connct

    use limitvars
    use pntvars
    use edgvars
    use trivars
    use cndvars
    use lnkvars
    use frtvars

    implicit none
    integer :: idegn(mtest)
    real(8) :: htest(mtest)
    integer :: l, i, k, n, n1, n2, nold
    integer :: ndegn, ncel1, ncel2, l1, l2
    integer :: imin, iedg
    real(8) :: xmid, ymid, dx, dy, ds
    real(8) :: tol, hmin

    tol         = 1d-15
    ncell       = 0
    nrept       = 0

10  newfr       = 0
    nrept       = nrept  +1

    do 100 l=1,nfront

    n = nfrt(l)
    if (ndg(4,n).ne.0) go to 100
    n1          = ndg(1,n)
    n2          = ndg(2,n)
    xmid        = 0.5d0*(x(1,n1) + x(1,n2))
    ymid        = 0.5d0*(x(2,n1) + x(2,n2))

    call findpt (n1,n2,n)
    if (itest.eq.0) go to 190

!   FIND TRIANGLE FOR WHICH THE HEIGHT OF THE CIRCUMCENTER ABOVE THE
!   EDGE MID-POINT IS MINIMIZED

    dx          = x(1,n2)  -xmid
    dy          = x(2,n2)  -ymid
    ds          = 1d0/dsqrt(dx*dx  +dy*dy)
    htest(1)    = (dx*(ycen(1)  -ymid) - dy*(xcen(1)  -xmid))*ds
    hmin        = htest(1)
    imin        = 1
    if (itest.eq.1) go to 22

    do 20 i=2,itest
    htest(i)    = (dx*(ycen(i)  -ymid)  -dy*(xcen(i)  -xmid))*ds
    if (htest(i).gt.hmin) go to 20
    hmin        = htest(i)
    imin        = i
20  continue

!   CHECK FOR DEGENERACY

    ndegn       = 0
    do 21 i=1,itest
    if (htest(i).gt.hmin+tol) go to 21
    ndegn       = ndegn + 1
    idegn(ndegn) = i
21  continue
    if (ndegn.eq.1) go to 22
    if (itest.eq.0) go to 190
    !call degen (n,ndegn,idegn,imin)
    call degen2(n,ndegn,idegn,imin)
    go to 100

!   ADD NEW TRIANGLE TO LIST AND UPDATE DATA STRUCTURE

22  ncell = ncell  +1
    ndc(1,ncell) = n1
    ndc(2,ncell) = n2
    ndc(3,ncell) = ntest(imin)
    xc(ncell)   = xcen(imin)
    yc(ncell)   = ycen(imin)
    rc(ncell)   = dsqrt((x(1,n1)  -xcen(imin))**2 + (x(2,n1)  -ycen(imin))**2)
    ndg(4,n)    = ncell
    ncelpt(n1)  = ncell
    ncelpt(n2)  = ncell
    ncelpt(ntest(imin)) = ncell

    l1 = min0(n1,ntest(imin))
    l2 = max0(n1,ntest(imin))
    call compar (l1,l2,nold,iedg)
    if (nold.eq.0) go to 25
    ndg(4,nold) = ncell
    go to 30
25  nedge = nedge +1
    ndg(1,nedge) = n1
    ndg(2,nedge) = ntest(imin)
    ndg(3,nedge) = ncell
    ndg(4,nedge) = 0
    newfr = newfr  +1
    nfill(newfr) = nedge
    npoint(nedge) = 0
    if (ipoint(l1).ne.0) npoint(iedg) = nedge
    if (ipoint(l1).eq.0) ipoint(l1) = nedge

30  l1          = min0(n2,ntest(imin))
    l2          = max0(n2,ntest(imin))
    call compar (l1,l2,nold,iedg)
    if (nold.eq.0) go to 35
    ndg(4,nold) = ncell
    go to 100
35  nedge       = nedge +1
    ndg(1,nedge) = ntest(imin)
    ndg(2,nedge) = n2
    ndg(3,nedge) = ncell
    ndg(4,nedge) = 0
    newfr        = newfr  +1
    nfill(newfr) = nedge
    npoint(nedge) = 0
    if (ipoint(l1).eq.0) ipoint(l1) = nedge
    if (ipoint(l1).ne.0) npoint(iedg) = nedge

100 continue

    nfront = 0
    if (newfr.eq.0) go to 130
    do 120 i=1,newfr
    if (ndg(4,nfill(i)).ne.0) go to 120
    nfront = nfront + 1
    nfrt(nfront) = nfill(i)
120 continue
    if (nfront.gt.0) go to 10

!   INITIALIZE NBH ARRAY

130 do 140 k=1,ncell
    npoint(k)  = 0
140 continue

    do 150 l=1,nedge
    ncel1 = ndg(3,l)
    ncel2 = ndg(4,l)
    npoint(ncel2) = npoint(ncel2)  +1
    nbh(npoint(ncel2),ncel2) = ncel1
    if (ncel1.eq.-1) go to 145
    npoint(ncel1) = npoint(ncel1)  +1
    nbh(npoint(ncel1),ncel1) = ncel2
    go to 150
145 ncell      = ncell +1
    ndc(1,ncell) = -1
    ndc(2,ncell) = ndg(1,l)
    ndc(3,ncell) = ndg(2,l)
    nbh(1,ncell) = ncel2
    nbh(2,ncell) = -1
    nbh(3,ncell) = -1
    nbh(npoint(ncel2),ncel2) = ncell
    ndg(4,l)     = ncell ! for 160, ndg(4,i)가 다시 edge정보로는 안쓰임
    nfill(ncell) = 1 ! for 165
    rc(ncell) = 0d0
150 continue

    do 155 i=1,nnode
    ipoint(i) = 0
    npoint(i) = 0
155 continue

    do 160 l=1,nedge
    if (ndg(3,l).ge.0) go to 160
    l1 = ndg(1,l)
    l2 = ndg(2,l)
    if (ipoint(l1).eq.0) ipoint(l1) = ndg(4,l)
    if (ipoint(l2).eq.0) ipoint(l2) = ndg(4,l)
    if (ipoint(l1).ne.0) npoint(l1) = ndg(4,l)
    if (ipoint(l2).ne.0) npoint(l2) = ndg(4,l)
160 continue

    do 165 i=1,nnode
    nfill(ipoint(i)) = nfill(ipoint(i))  +1
    nfill(npoint(i)) = nfill(npoint(i))  +1
    nbh(nfill(ipoint(i)),ipoint(i)) = npoint(i)
    nbh(nfill(npoint(i)),npoint(i)) = ipoint(i)
165 continue

!   ARRANGE IPOINT AND NPOINT FOR LATER ( INSERT )

    do 170 k=1,ncell
    npoint(k)  = 0
170 continue
    do 180 i=1,nnode
    ipoint(i) = 0
180 continue

    return

190 write(*,200)
    stop
200 format('itest is zero in connect routine')

end subroutine
