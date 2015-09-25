!---------------------------------------------------------------------!
!                                                                     !
!          CONSTRUCTS TRIANGLES IN THE DEGENERATE CASE OF ANY         !
!          CO-CIRCULAR POINTS.                                        !
!                                                                     !
!---------------------------------------------------------------------!

subroutine degen2( n, ndegn, idegn, imin )
    
    use limitvars
    use pntvars
    use edgvars
    use trivars
    use cndvars
    use lnkvars
    use frtvars
    
    implicit none
    integer, intent(inout) :: n, imin, idegn(mtest)
    integer, intent(inout) :: ndegn
    real(8) :: ang(ndegn)
    real(8) :: dx1, dy1, dx2, dy2, ds1, ds2
    integer :: n1, n2, m(ndegn), k(ndegn)
    integer :: l1, l2, nold, iedg, i, ii, ia, ib, imid
    
    n1          = ndg(1,n)
    n2          = ndg(2,n)
    dx1         = x(1,n2)  -x(1,n1)
    dy1         = x(2,n2)  -x(2,n1)
    ds1         = dsqrt(dx1*dx1  +dy1*dy1)
    do i =1,ndegn
    m(i)        = ntest(idegn(i))
    dx2         = x(1,m(i))  -x(1,n1)
    dy2         = x(2,m(i))  -x(2,n1)
    ds2         = dsqrt(dx2*dx2  +dy2*dy2)
    ang(i)      = (dx1*dx2  +dy2*dy1)/ds1/ds2
    end do
    
    k(1) = 1
    do 5 i = 2,ndegn  
    ii = 1
    if( ang(i) .gt. ang(k(1)) ) goto 10
    if( ang(i) .lt. ang(k(i-1)) ) goto 15
    ia = 1
    ib = i-1
30  imid = (ia+ib)/2
    if( ang(i) .lt. ang(k(imid)) ) goto 20
    ib = imid
    goto 25
20  ia = imid
25  if( ib .gt. ia+1 ) goto 30
    ii = ia
    if(ia.ne.ib) ii = ib
10  k(ii+1:i) = k(ii:i-1)
    k(ii) = i
    goto 5
15  k(i) = i
5   end do
    
    do i = 1,ndegn
    k(i) = m(k(i))
    end do
    
    ncell       = ncell  +1
    ndc(1,ncell) = n1
    ndc(2,ncell) = n2
    ndc(3,ncell) = k(1)
    xc(ncell)   = xcen(imin)
    yc(ncell)   = ycen(imin)
    rc(ncell)   = dsqrt((x(1,n1)  -xcen(imin))**2 + (x(2,n1)  -ycen(imin))**2)
    ndg(4,n)    = ncell
    ncelpt(n1)  = ncell
    ncelpt(n2)  = ncell
    ncelpt(k(1))  = ncell
    
    l1          = min0(n2,k(1))
    l2          = max0(n2,k(1))
    call compar (l1,l2,nold,iedg)
    if (nold.eq.0) go to 35
    ndg(4,nold) = ncell
    go to 40
35  nedge       = nedge  +1
    ndg(1,nedge) = k(1)
    ndg(2,nedge) = n2
    ndg(3,nedge) = ncell
    ndg(4,nedge) = 0
    newfr        = newfr  +1
    nfill(newfr) = nedge
    npoint(nedge) = 0
    if (ipoint(l1).eq.0) ipoint(l1) = nedge
    if (ipoint(l1).ne.0) npoint(iedg) = nedge
    
40  do 45 i = 1,ndegn-1
    
    ncell       = ncell  +1
    ndc(1,ncell) = k(i)
    ndc(2,ncell) = k(i+1)
    ndc(3,ncell) = n1
    xc(ncell)   = xc(ncell-1)
    yc(ncell)   = yc(ncell-1)
    rc(ncell)   = rc(ncell-1)
    
    nedge       = nedge  +1
    ndg(1,nedge) = k(i)
    ndg(2,nedge) = n1
    ndg(3,nedge) = ncell
    ndg(4,nedge) = ncell  -1
    ncelpt(k(i+1))  = ncell
    
    l1          = min0(k(i),k(i+1))
    l2          = max0(k(i),k(i+1))
    call compar (l1,l2,nold,iedg)
    if (nold.eq.0) go to 50
    ndg(4,nold) = ncell
    go to 45
50  nedge       = nedge  +1
    ndg(1,nedge) = k(i+1)
    ndg(2,nedge) = k(i)
    ndg(3,nedge) = ncell
    ndg(4,nedge) = 0
    newfr        = newfr  +1
    nfill(newfr) = nedge
    npoint(nedge) = 0
    if (ipoint(l1).eq.0) ipoint(l1) = nedge
    if (ipoint(l1).ne.0) npoint(iedg) = nedge
    
45  end do
    
    l1          = min0(n1,k(ndegn))
    l2          = max0(n1,k(ndegn))
    call compar (l1,l2,nold,iedg)
    if (nold.eq.0) go to 55
    ndg(4,nold) = ncell
    go to 60
55  nedge       = nedge  +1
    ndg(1,nedge) = n1
    ndg(2,nedge) = k(ndegn)
    ndg(3,nedge) = ncell
    ndg(4,nedge) = 0
    newfr        = newfr  +1
    nfill(newfr) = nedge
    npoint(nedge) = 0
    if (ipoint(l1).eq.0) ipoint(l1) = nedge
    if (ipoint(l1).ne.0) npoint(iedg) = nedge  
    
60  return

end subroutine
