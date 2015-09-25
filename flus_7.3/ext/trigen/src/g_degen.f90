!---------------------------------------------------------------------!
!                                                                     !
!    CONSTRUCTS A PAIR OF TRIANGLES IN THE DEGENERATE CASE OF FOUR    !
!    CO-CIRCULAR POINTS.                                              !
!                                                                     !
!---------------------------------------------------------------------!

subroutine degen( n, ndegn, idegn, imin )

    use limitvars
    use pntvars
    use edgvars
    use trivars
    use cndvars
    use lnkvars
    use frtvars
    
    implicit none
    integer, intent(in) :: n, imin, idegn(mtest)
    integer, intent(inout) :: ndegn
    real(8) :: ang(2)
    real(8) :: dx1, dy1, dx2, dy2, ds1, ds2
    integer :: n1, n2, m1, m2, k1, k2
    integer :: l1, l2, nold, iedg
    
    if (ndegn.gt.2) go to 100
    n1          = ndg(1,n)
    n2          = ndg(2,n)
    m1          = ntest(idegn(1))
    m2          = ntest(idegn(2))
    dx1         = x(1,n2)  -x(1,n1)
    dy1         = x(2,n2)  -x(2,n1)
    ds1         = dsqrt(dx1*dx1  +dy1*dy1)
    dx2         = x(1,m1)  -x(1,n1)
    dy2         = x(2,m1)  -x(2,n1)
    ds2         = dsqrt(dx2*dx2  +dy2*dy2)
    ang(1)      = (dx1*dx2  +dy2*dy1)/ds1/ds2
    dx2         = x(1,m2)  -x(1,n1)
    dy2         = x(2,m2)  -x(2,n1)
    ds2         = dsqrt(dx2*dx2  +dy2*dy2)
    ang(2)      = (dx1*dx2  +dy2*dy1)/ds1/ds2
    k1          = m1
    if (ang(2).lt.ang(1)) k1 = m2
    k2          = m1  +m2  -k1
    
    ncell       = ncell  +1
    ndc(1,ncell) = n1
    ndc(2,ncell) = n2
    ndc(3,ncell) = k2
    xc(ncell)   = xcen(imin)
    yc(ncell)   = ycen(imin)
    rc(ncell)   = dsqrt((x(1,n1)  -xcen(imin))**2 + (x(2,n1)  -ycen(imin))**2)
    ndg(4,n)    = ncell
    ncelpt(n1)  = ncell
    ncelpt(n2)  = ncell
    ncelpt(k2)  = ncell
    
    l1          = min0(n2,k2)
    l2          = max0(n2,k2)
    call compar (l1,l2,nold,iedg)
    if (nold.eq.0) go to 25
    ndg(4,nold) = ncell
    go to 30
25  nedge       = nedge  +1
    ndg(1,nedge) = k2
    ndg(2,nedge) = n2
    ndg(3,nedge) = ncell
    ndg(4,nedge) = 0
    newfr        = newfr  +1
    nfill(newfr) = nedge
    npoint(nedge) = 0
    if (ipoint(l1).eq.0) ipoint(l1) = nedge
    if (ipoint(l1).ne.0) npoint(iedg) = nedge
    
30  ncell       = ncell  +1
    ndc(1,ncell) = k2
    ndc(2,ncell) = k1
    ndc(3,ncell) = n1
    xc(ncell)   = xc(ncell-1)
    yc(ncell)   = yc(ncell-1)
    rc(ncell)   = rc(ncell-1)
    
    nedge       = nedge  +1
    ndg(1,nedge) = k2
    ndg(2,nedge) = n1
    ndg(3,nedge) = ncell
    ndg(4,nedge) = ncell  -1
    ncelpt(k1)  = ncell
    
    l1          = min0(n1,k1)
    l2          = max0(n1,k1)
    call compar (l1,l2,nold,iedg)
    if (nold.eq.0) go to 35
    ndg(4,nold) = ncell
    go to 40
35  nedge       = nedge  +1
    ndg(1,nedge) = n1
    ndg(2,nedge) = k1
    ndg(3,nedge) = ncell
    ndg(4,nedge) = 0
    newfr        = newfr  +1
    nfill(newfr) = nedge
    npoint(nedge) = 0
    if (ipoint(l1).eq.0) ipoint(l1) = nedge
    if (ipoint(l1).ne.0) npoint(iedg) = nedge
    
40  l1          = min0(k1,k2)
    l2          = max0(k1,k2)
    call compar (l1,l2,nold,iedg)
    if (nold.eq.0) go to 45
    ndg(4,nold) = ncell
    go to 50
45  nedge       = nedge  +1
    ndg(1,nedge) = k1
    ndg(2,nedge) = k2
    ndg(3,nedge) = ncell
    ndg(4,nedge) = 0
    newfr        = newfr  +1
    nfill(newfr) = nedge
    npoint(nedge) = 0
    if (ipoint(l1).eq.0) ipoint(l1) = nedge
    if (ipoint(l1).ne.0) npoint(iedg) = nedge
    
50  return
    
100 ndegn        = ndegn  +2
    write(6,*) n
    write(6,*) idegn(1:ndegn-2)
    write (6,600) ndegn
600 format (/5x,'degeneracy with ',i3,' co-circular points')
    stop
    
end subroutine