!---------------------------------------------------------------------!
!                                                                     !
!    INTRODUCE POINTS INSIDE THE DOMAIN USING REBAY'S ALGORITHM       !
!                                                                     !
!---------------------------------------------------------------------!

subroutine rebay

    use pntvars
    use edgvars
    use trivars
    use lnkvars
    use qudvars
    use typvars
    use bptvars
    use nrmvars
    use ordvars
    use rebvars

    implicit none
    integer :: l, n, n1, n2, n3, l1, l2, l3, i, npre, k
    integer :: knhbr, lnhb, nend, m, m1, m2, m3, mend, nintr, lnbh
    real(8) :: tolsiz, p, s1sq, s2sq, xpt, ypt, xmid, ymid
    real(8) :: q, densty, rad1, rad2, cirrad, arg, dist
    real(8) :: dx, dy, rs, rn1, rn2, xdif, ydif, rsign

    !call midtest

!   INITIAL CLASSIFICATION OF TRIANGLES AS ACCEPTED OR NON-ACCEPTED

    tolsiz = 0.8d0

    do 10 l=1,ncell
    nacpt(l) = 0
    n1 = ndc(1,l)
    n2 = ndc(2,l)
    n3 = ndc(3,l)
    if(n1.eq.-1) goto 10
    if(dens(n1)+dens(n2)+dens(n3).lt.3d0*rc(l)*tolsiz) nacpt(l)=1
10  continue

    nwait = 0

    do 15 l=1,ncell
    if(nacpt(l).eq.0) goto 15
    l1 = nbh(1,l)
    l2 = nbh(2,l)
    l3 = nbh(3,l)
    call putin(l)
    if(nacpt(l1).gt.0.and.nacpt(l2).gt.0.and.nacpt(l3).gt.0) nacpt(l) = 2
15  continue

    if(nwait.eq.0) return

    npre = nnode
    it = 0

!   FIND ACTIVE TRIANGLE WITH LARGEST CIRCUM-RADIUS TO DENSITY RATIO
!   AND LOCATE EDGE THAT IT SHARES WITH AN ACCEPTED TRIANGLE

18  i=0
20  i=i+1
    l=listc(i)
    if(i.eq.nwait+1) goto 200
    if(nacpt(l).ne.1) goto 20
    it = it + 1

    knhbr = 0
    do 25 k=1,3
    lnhb = nbh(k,l)
    if(nacpt(lnhb).gt.0) goto 25
    knhbr = k
25  continue
    if(knhbr.eq.0) goto 210

    lnhb = nbh(knhbr,l)
    n1 = ndc(1,l)
    n2 = ndc(2,l)
    n3 = ndc(3,l)
    nend = n3
30  m1 = ndc(1,lnhb)
    m2 = ndc(2,lnhb)
    m3 = ndc(3,lnhb)
    mend = m3
35  if(m1+m2.eq.n1+n2.and.iabs(m1-m2).eq.iabs(n1-n2)) goto 50
    if(m1.eq.mend) goto 40
    m = m1
    m1 = m2
    m2 = m3
    m3 = m
    goto 35
40  if (n1.eq.nend) goto 220
    n = n1
    n1 = n2
    n2 = n3
    n3 = n
    goto 30

!   DETERMINE LOCATION ON VORONOI SEGMENT WHERE NEW POINT WILL BE INSERTED

50  p = 0.5d0*dsqrt((x(1,n2)-x(1,n1))**2+(x(2,n2)-x(2,n1))**2)
    s1sq = 0.25d0*((x(1,n3)-x(1,n2))**2+(x(2,n3)-x(2,n2))**2)
    s2sq = 0.25d0*((x(1,n1)-x(1,n3))**2+(x(2,n1)-x(2,n3))**2)
    xpt=xc(l)
    ypt=yc(l)

    if(p*p.gt.s1sq+s2sq) goto 70
    xmid = 0.5d0*(x(1,n1) + x(1,n2))
    ymid = 0.5d0*(x(2,n1) + x(2,n2))
    q=dsqrt((xc(l)-xmid)**2+(yc(l)-ymid)**2)

    if(q.lt.p) goto 70
    densty = 0.5d0*(dens(n1)+dens(n2))
    rad1=dmax1(densty,p)
    rad2=0.5d0*(p*p+q*q)/q
    cirrad=dmin1(rad1,rad2)
    arg=dabs(cirrad*cirrad-p*p)
    dist=cirrad+dsqrt(arg)
    dx=x(1,n2)-x(1,n1)
    dy=x(2,n2)-x(2,n1)
    rs=1d0/dsqrt(dx*dx+dy*dy)
    rn1=dy*rs
    rn2=-dx*rs
    xdif=x(1,n3)-xmid
    ydif=x(2,n3)-ymid
    rsign=1d0
    if(xdif*rn1+ydif*rn2.lt.0d0) rsign=-1d0
    xpt = xmid+rsign*dist*rn1
    ypt = ymid+rsign*dist*rn2

70  nstop(l) = 0
    lcel = l

    call insert(xpt,ypt)

    if(nstop(l).eq.1) goto 80

    nintr=nnode-npre

!    if(mod(it,500).eq.0) then
!        write(6,700) it,nintr,nwait
!        call midtest
!    end if

    if(nwait.eq.0.or.it.eq.200000) goto 100
    goto 18

80  nacpt(l) = 0
    do 85 k=1,3
    lnbh = nbh(k,l)
    if(nacpt(lnbh).eq.2) nacpt(lnbh) =1
85  continue

    call pulout(l)

    if(nwait.eq.0.or.it.eq.2000000) goto 100
    goto 18

!100 call qualty
!    return
100 return

200 write(6,600) nwait
    return
210 write(6,610)
    stop
220 write(6,620)
    stop
600   format(/5x,'no active triangles left.'/ 5x,'non-accepted cells remaining',i5)
610   format(/5x,'an active cell without any accepted neighboring'/ 5x,'triangles has been found')
620   format(/5x,'unable to find a common edge between an active and'/ 5x,'an accepted triangle')
700   format(/5x,'iteration ',i5,' accepted points ',i5/ 5x,'cells still too large ',i5)

end subroutine
