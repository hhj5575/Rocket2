!---------------------------------------------------------------------!
!                                                                     !
!    FINDS QUADRANT WHICH CONTAINS THE NEW POINT AND THEN             !
!    FINDS CLOSEST POINT IN THAT QUADRANT. IF QUADRANT IS NOW         !
!    FULL, IT IS DIVIDED INTO FOUR NEW QUADRANTS AND THE              !
!    QUADTREE DATA STRUCTURE IS AMENDED.                              !
!                                                                     !
!---------------------------------------------------------------------!

subroutine quad( n, jnpt )

    use limitvars
    use pntvars
    use edgvars
    use trivars
    use lnkvars
    use qudvars
    use ctlvars
    
    implicit none
    integer, intent(out) :: jnpt
    integer :: jj(4), n, i, j, l, k, nxsgn, nysgn
    integer :: iroot, isqu, nsignx, nsigny
    real(8) :: tolp, tol, dmin
    real(8) :: xshift, yshift, xsize, ysize
    
    tolp       = 1d-15
    if (x(1,n).lt.xquad(1,1)-tolp.or.x(1,n).gt.xquad(2,1)+tolp) go to 200
    if (x(2,n).lt.yquad(1,1)-tolp.or.x(2,n).gt.yquad(2,1)+tolp) go to 200
    tol        = 1.000001d0
    
5   i          = 1
    l=nquad(7,i)
    if (l.ge.0) go to 15
10  xshift     = x(1,n)  -0.5d0*(xquad(1,i)  +xquad(2,i))
    xsize      = dmax1(1d-9,abs(xshift))
    nxsgn      = (ifix(real(tol*xshift/xsize))  +1)/2*2
    yshift     = x(2,n)  -0.5d0*(yquad(1,i)  +yquad(2,i))
    ysize      = dmax1(1d-9,abs(yshift))
    nysgn      = (ifix(real(tol*yshift/ysize))  +1)/2*2
    l          = 1  +nxsgn/2  +nysgn
    i          = nquad(l,i)
    l          = nquad(7,i)
    if (l.lt.0) go to 10
15  iroot      = i
    if (ncode.eq.0) go to 20
    isqu       = i
    call nearpt (isqu,n,jnpt,dmin)
    if (dmin.lt.tolp) go to 210
    return
    
20  if (l.eq.4) go to 30
    
25  l          = l  +1
    nquad(l,i) = n
    nquad(7,i) = l
    return

!   QUADRANT NOW CONTAINS FIVE POINTS. FORM FOUR NEW QUADRANTS
!   AND AMEND TREE STRUCTURE.

30  do 35 k=1,4
    iquad      = iquad  +1
    if (iquad.gt.mquad) go to 220
    nquad(7,iquad) = 0
    nquad(6,iquad) = i
    nquad(5,iquad) = k
    nsignx     = 2  -mod(k,2)
    xquad(nsignx,iquad) = xquad(nsignx,i)
    xquad(3-nsignx,iquad) = 0.5d0*(xquad(1,i)  +xquad(2,i))
    nsigny     = 1  +(isign(1,2*k-5)  +1)/2
    yquad(nsigny,iquad) = yquad(nsigny,i)
    yquad(3-nsigny,iquad) = 0.5d0*(yquad(1,i)  +yquad(2,i))
    do 35 l=1,4
35  nquad(l,iquad) = 0

!   ASSIGN POINTS TO NEW QUADRANTS

    do 40 l=1,4
40  jj(l)      = 0
    do 45 k=1,4
    j          = nquad(k,i)
    xshift     = x(1,j)  -0.5d0*(xquad(1,i)  +xquad(2,i))
    xsize      = dmax1(1d-9,abs(xshift))
    nxsgn      = (ifix(real(tol*xshift/xsize))  +1)/2*2
    yshift     = x(2,j)  -0.5d0*(yquad(1,i)  +yquad(2,i))
    ysize      = dmax1(1d-9,abs(yshift))
    nysgn      = (ifix(real(tol*yshift/ysize))  +1)/2*2
    l          = 1  +nxsgn/2  +nysgn
    jj(l)      = nquad(7,iquad-4+l)  +1
    nquad(7,iquad-4+l) = jj(l)
    nquad(jj(l),iquad-4+l) = j
45  continue

!   CHECK WHETHER ALL FOUR POINTS LIE IN ONE QUADRANT

    nquad(7,i) = -1
    do 50 l=1,4
50  nquad(l,i) = iquad  -4  +l
    do 60 l=1,4
    i          = iquad  -4  +l
    if (jj(l).eq.4) go to 30
60  continue
    i          = iroot

!   QUADRANT I IS ALREADY FULL. LOCATE SUB-QUADRANT IN
!   WHICH POINT LIES

80  xshift     = x(1,n)  -0.5d0*(xquad(1,i)  +xquad(2,i))
    xsize      = dmax1(1d-9,abs(xshift))
    nxsgn      = (ifix(real(tol*xshift/xsize))  +1)/2*2
    yshift     = x(2,n)  -0.5d0*(yquad(1,i)  +yquad(2,i))
    ysize      = dmax1(1d-9,abs(yshift))
    nysgn      = (ifix(real(tol*yshift/ysize))  +1)/2*2
    l          = 1  +nxsgn/2  +nysgn
    i          = nquad(l,i)
    l          = nquad(7,i)
    if (l.lt.0) go to 80
    go to 25
    
200 jnpt       = -1
    return
210 jnpt       = -1
    return
220 write (6,620)
    stop
600 format(/5x,'new point outside domain,  x = ',f10.4,' y = ',f10.4/5x,'point number ',i5)
610 format(/5x,'new point, x = ',f10.4,' y = ',f10.4/5x,'is too close to an existing point'/)
620 format(//5x,'parameter mquad exceeded in routine quad')
    
end subroutine