!---------------------------------------------------------------------!
!                                                                     !
!    DETERMINE TRIANGLES WHOSE CIRCUMCIRCLES CONTAIN THE POINT        !
!    (XPT,YPT) AND CONNECT POINT TO BOUNDARY EDGES OF CAVITY TO       !
!    GENERATE THE DELAUNAY TRIANGULATION CONTAINING THAT POINT.       !
!                                                                     !
!---------------------------------------------------------------------!

subroutine insert( xpt, ypt )

    use limitvars
    use pntvars
    use edgvars
    use trivars
    use lnkvars
    use frtvars
    use cndvars
    use typvars
    use ordvars
    use rebvars
    use ctlvars
    use nrmvars
    
    implicit none
    integer :: ldel(mtest), ntri(mtest), newtri(mtest)
    integer :: newc(mtest), ntriag(2,mtest)
    integer :: jnpt, jstrt, jpre, j, jnew, jcnt, ncnt
    integer :: k1, k2, k3, k, lnbh, n1, n2, n3, i, ndel, ichk, nnew
    integer :: l, l1, l2, l3, lend, m1, m2, m3, mend, lmin, lmax, mmin, mmax, m
    real(8) :: tolsiz, tol, xpt, ypt, area, a1, a2, a3, sum
    real(8) :: x1, x2, y1, y2, a, altde, f, alpha, beta
    
    tolsiz = 0.8d0
    tol    = 1d-15
    nnode  = nnode + 1
    if (nnode.gt.mnode) go to 280
    x(1,nnode) = xpt
    x(2,nnode) = ypt
    ipoint(nnode) = 0

!   FIND TRIANGLE WHICH CONTAINS THE NEW POINT

    ncode      = 1
    call quad (nnode,jnpt)
    if (jnpt.eq.-1) go to 270
    jstrt      = ncelpt(jnpt)
    jpre       = jstrt
    j          = jstrt
    jnew       = j
    jcnt       = 1
    ntri(1)    = jnew
    nfill(jnew) = 1 ! 안비우고 그냥 씀( 이상한거 발견되면 고치기 비우고 안비우고 시간 차이 너무 남 )
    k1 = ndc(1,j)
    k2 = ndc(2,j)
    k3 = ndc(3,j)

    if (dsqrt((xpt-xc(j))**2+(ypt-yc(j))**2).lt.rc(j)) goto 25
10  jnew       = 0
    do 15 k=1,3
    lnbh       = nbh(k,j)
    if (lnbh.eq.jpre) go to 15
    if (ndc(1,lnbh).ne.jnpt.and.ndc(2,lnbh).ne.jnpt.and.ndc(3,lnbh).ne.jnpt) go to 15
    jnew       = lnbh
15  continue
    if (jnew.eq.jstrt) go to 35
    if (ndc(1,jnew).eq.-1) go to 20
    jcnt       = jcnt  +1
    if (jcnt.gt.mtest) go to 255
    ntri(jcnt) = jnew
    nfill(jnew) = 1
    if (dsqrt((xpt-xc(jnew))**2 +(ypt-yc(jnew))**2).lt.rc(jnew)) go to 25
20  jpre       = j
    j          = jnew
    go to 10
25  n1         = ndc(1,jnew)
    n2         = ndc(2,jnew)
    n3         = ndc(3,jnew)
    area       = dabs((x(1,n2)  -x(1,n1))*(x(2,n3)  -x(2,n1))-(x(2,n2)  -x(2,n1))*(x(1,n3)  -x(1,n1)))
    a1         = dabs((x(1,n3)  -x(1,n2))*(ypt  -x(2,n2))-(x(2,n3)  -x(2,n2))*(xpt  -x(1,n2)))
    a2         = dabs((x(1,n1)  -x(1,n3))*(ypt  -x(2,n3))-(x(2,n1)  -x(2,n3))*(xpt  -x(1,n3)))
    a3         = dabs((x(1,n2)  -x(1,n1))*(ypt  -x(2,n1))-(x(2,n2)  -x(2,n1))*(xpt  -x(1,n1)))
    sum=a1+a2+a3-tol
    if (area.gt.a1+a2+a3-tol) go to 50
    jpre       = j
    j          = jnew
    go to 10

!   POINT NOT CONTAINED IN FIRST CELL GROUP. EXAMINE ALL CELLS
!   WHICH HAVE AN EDGE IN COMMON WITH THE FIRST GROUP.

35  ncnt      = 0
    do 45 j=1,jcnt
    do 40 k=1,3
    jnew       = nbh(k,ntri(j))
    if (nfill(jnew).eq.1) go to 40
    if (ndc(1,jnew).eq.-1) go to 40
    ncnt       = ncnt  +1
    if (ncnt.gt.mtest) go to 255
    newtri(ncnt) = jnew
    nfill(jnew) = 1
    if (dsqrt((xpt -xc(jnew))**2 +(ypt -yc(jnew))**2).gt.rc(jnew)) go to 40
    n1         = ndc(1,jnew)
    n2         = ndc(2,jnew)
    n3         = ndc(3,jnew)
    area       = dabs((x(1,n2)  -x(1,n1))*(x(2,n3)  -x(2,n1))-(x(2,n2)  -x(2,n1))*(x(1,n3)  -x(1,n1)))
    a1         = dabs((x(1,n3)  -x(1,n2))*(ypt  -x(2,n2))-(x(2,n3)  -x(2,n2))*(xpt  -x(1,n2)))
    a2         = dabs((x(1,n1)  -x(1,n3))*(ypt  -x(2,n3))-(x(2,n1)  -x(2,n3))*(xpt  -x(1,n3)))
    a3         = dabs((x(1,n2)  -x(1,n1))*(ypt  -x(2,n1))-(x(2,n2)  -x(2,n1))*(xpt  -x(1,n1)))
    sum = a1+a2+a3-tol
    if (area.gt.a1+a2+a3-tol) go to 50
40  continue
45  continue
    
    if (ncnt.eq.0) go to 260
    do 47 i=1,ncnt
    ntri(i)   = newtri(i)
47  continue
    jcnt      = ncnt
    go to 35

!   DETERMINE LIST OF CELLS IN CAVITY

50  dens(nnode)=(dens(n1)*a1+dens(n2)*a2+dens(n3)*a3)/area
    
    do 52 i=1,ncell
    nfill(i)  = 0
52  continue
    
    ndel       = 1
    ldel(ndel) = jnew
    npoint(jnew) = 1
    ichk       = 1
    ntri(1)    = jnew
55  ncnt       = 0
    do 65 i=1,ichk
    l          = ntri(i)
    do 60 k=1,3
    lnbh       = nbh(k,l)
    if (npoint(lnbh).eq.1) go to 60
    if (ndc(1,lnbh).eq.-1) go to 60
    if (dsqrt((xpt-xc(lnbh))**2 +(ypt-yc(lnbh))**2).gt.rc(lnbh)) go to 60
    n1         = ndc(1,lnbh)
    n2         = ndc(2,lnbh)
    n3         = ndc(3,lnbh)
    ncnt       = ncnt  +1
    if (ncnt.gt.mtest) go to 255
    newtri(ncnt) = lnbh
    ndel       = ndel  +1
    ldel(ndel) = lnbh
    npoint(lnbh) = 1
60  continue
65  continue
    if (ncnt.eq.0) go to 75
    do 70 i=1,ncnt
    ntri(i)    = newtri(i)
70  continue
    ichk       = ncnt
    go to 55

!   DETERMINE EDGES ON CAVITY BOUNDARY AND GENERATE LIST OF NEW CELLS

75  ncnt       = 0
    do 115 i=1,ndel
    do 110 k=1,3
    lnbh       = nbh(k,ldel(i))
    if (npoint(lnbh).eq.1) go to 110
    l1         = ndc(1,ldel(i))
    l2         = ndc(2,ldel(i))
    l3         = ndc(3,ldel(i))
    lend       = l3
    m1         = ndc(1,lnbh)
    m2         = ndc(2,lnbh)
    m3         = ndc(3,lnbh)
80  mend       = m3
    lmin       = min0(l1,l2)
    lmax       = max0(l1,l2)
85  mmin       = min0(m1,m2)
    mmax       = max0(m1,m2)
    if (mmin.eq.lmin.and.mmax.eq.lmax) go to 95
    if (m1.eq.mend) go to 90
    m          = m1
    m1         = m2
    m2         = m3
    m3         = m
    go to 85
90  if(l1.eq.lend) goto 200
    l          = l1
    l1         = l2
    l2         = l3
    l3         = l
    go to 80
95  ncnt       = ncnt  +1
    if (ncnt.gt.mtest) go to 210
    ntriag(1,ncnt) = l1
    ntriag(2,ncnt) = l2
    ipoint(l1)  = ipoint(l1)  +1
    ipoint(l2)  = ipoint(l2)  +1
    ndg(ipoint(l1),l1) = ncnt
    ndg(ipoint(l2),l2) = ncnt
    nfill(ncnt) = lnbh
    do 100 j=1,3
    if (nbh(j,lnbh).eq.ldel(i)) go to 105
100 continue
    go to 220
105 ntri(ncnt) = j
110 continue
115 continue

!   CHECK ON CONSISTENCY OF NEW DATA STRUCTURE

    if (ncnt.ne.ndel+2) go to 230
    do 120 i=1,ncnt
    l1         = ntriag(1,i)
    l2         = ntriag(2,i)
    if(ipoint(l1).ne.2) go to 240
    if(ipoint(l2).ne.2) go to 240
120 continue

!   DEFINE NEW CELLS

    do 125 i=1,ndel
    newc(i)    = ldel(i)
    npoint(newc(i)) = 0
125 continue
    nnew       = ndel
    do 130 i=1,2
    nnew       = nnew  +1
    ncell      = ncell +1
    newc(nnew) = ncell
    npoint(ncell) = 0
130 continue

!   RESET IPOINT ARRAY     reset ipoint array

    do 132 i=1,ncnt
    l1         = ntriag(1,i)
    l2         = ntriag(2,i)
    ipoint(l1) = 0
    ipoint(l2) = 0
132 continue

!   DETERMINE CIRCUMCENTER AND CIRCUM-RADIUS FOR NEW TRIANGLES


    do 135 i=1,ncnt
    l1         = ntriag(1,i)
    l2         = ntriag(2,i)
    x1         = x(1,l2)  -x(1,l1)
    y1         = x(2,l2)  -x(2,l1)
    x2         = xpt  -x(1,l1)
    y2         = ypt  -x(2,l1)
    a          = x1*y2  -x2*y1
    if (abs(a).lt.1d-10) go to 250
    altde      = dabs(a)/dsqrt(x1*x1  +y1*y1)
    if (altde.lt.0.25d0*(dens(l1)+dens(l2))) goto 250 ! skew
133 f          = 1d0/a
    alpha      = 0.5d0*(x1*(x(1,l2)  +x(1,l1))  +y1*(x(2,l2)  +x(2,l1)))
    beta       = 0.5d0*(x2*(xpt  +x(1,l1))  +y2*(ypt  +x(2,l1)))
    xcen(i)    = f*(alpha*y2  -beta*y1)
    ycen(i)    = f*(beta*x1  -alpha*x2)
    goto 135
135 continue

!   GENERATE NBH ARRAY FOR NEW TRIANGLES

    do 140 i=1,ncnt
    nbh(ntri(i),nfill(i)) = newc(i)
    nbh(1,newc(i)) = nfill(i)
    l1         = ntriag(1,i)
    l2         = ntriag(2,i)
    k          = ndg(1,l1)  +ndg(2,l1)  -i
    nbh(2,newc(i)) = newc(k)
    k          = ndg(1,l2)  +ndg(2,l2)  -i
    nbh(3,newc(i)) = newc(k)
140 continue

!   UPDATE DATA STRUCTURE

    do 145 i=1,ndel
    if(nacpt(ldel(i)).eq.0) goto 145
    call pulout(ldel(i))
145 continue

    do 150 i=1,ncnt
    nfill(i)   = 0
    l1         = ntriag(1,i)
    l2         = ntriag(2,i)
    ndc(1,newc(i)) = nnode
    ndc(2,newc(i)) = l1
    ndc(3,newc(i)) = l2
      
    rc(newc(i)) = 0.
    nstop(newc(i)) = 0
    xc(newc(i)) = xcen(i)
    yc(newc(i)) = ycen(i)
    rc(newc(i)) = dsqrt((xpt  -xc(newc(i)))**2+(ypt  -yc(newc(i)))**2)
    ncelpt(l1)  = newc(i)
    ncelpt(l2)  = newc(i)
150 continue
    ncelpt(nnode) = newc(1)
    ityp(nnode) = 0
 
    do 152 i=1,ncnt
    l = newc(i)
    nacpt(l) = 0
    n1 = ndc(1,l)
    n2 = ndc(2,l)
    n3 = ndc(3,l)
    if(dens(n1)+dens(n2)+dens(n3).lt.3d0*rc(l)*tolsiz) nacpt(l)=1
152 continue
      
    do 153 i=1,ncnt
    l = newc(i)
    if(nacpt(l).eq.1) call putin(l)
153 continue
    do 154 i=1,nwait
    l=listc(i)
    nacpt(l)=1
154 continue
    do 156 i=1,nwait
    l =listc(i)
    l1=nbh(1,l) 
    l2=nbh(2,l) 
    l3=nbh(3,l) 
    if(nacpt(l1).gt.0.and.nacpt(l2).gt.0.and.nacpt(l3).gt.0) nacpt(l)=2
156 continue
    
    ncode      = 0
    call quad(nnode,jnpt)
    return
    
200 write (6,600) nnode,ndel
    stop
210 write (6,610)
    stop
220 write (6,620)
    stop
230 write (6,630)
    stop
240 write (6,640)
    stop
250 continue
    nstop(lcel) = 1
    ncell     = ncell  -2
    nnode     = nnode  -1
    return
255 write (6,655) 
    stop
260 continue
    nstop (lcel) =1
    nnode    = nnode  -1
    return
265 continue
    nstop (lcel) =1
    nnode    = nnode  -1
    return
270 continue
    nstop (lcel) =1
    nnode    = nnode  -1
    return
280 write (6,680)
    stop
600 format(//5x,'unable to find a cavity edge'/5x,'nnode = ',i5,'  ndel = ',i3)
610 format(//5x,'parameter mtest exceeded in routine insert')
620 format(//5x,'unable to find nbh value to give cavity edge'/5x,'contiguity')
630 format(//5x,'new cavity triangles not equal to ndel+2')
640 format(//5x,'incorrect assembly of contiguity information'/5x,'for new triangles')
655 format(//5x,'number of deleted triangles in cavity exceeds'/5x,'parameter limit mtest')
660 format(//5x,'I: unable to find triangle which contains new point'/5x,'x = ',f10.4,'y = ',f10.4,'nnode = ',i5)
680 format(//5x,'nnode exceeds limit mnode in routine insert')
      
end subroutine