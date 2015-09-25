!---------------------------------------------------------------------!
!                                                                     !
!    USES QUADTREE STRUCTURE TO FIND NEAREST POINT TO X(N,-)          !
!                                                                     !
!---------------------------------------------------------------------!

subroutine nearpt( i, n, jnpt, dmin )

    use limitvars
    use pntvars
    use edgvars
    use trivars
    use qudvars
    
    implicit none
    integer, intent(in) :: n
    integer :: isrch(msrch), ksrch(msrch)
    integer :: l, i, j, k, jnpt, iflag, kc, ic, ii, itry, kk, jj
    real(8) :: dmin, dismin, dist
    real(8) :: xl, xu, yl, yu
    
    l           = nquad(7,i)
    jnpt        = 1  
    dismin      = (xquad(2,1)  -xquad(1,1))**2+(yquad(2,1)  -yquad(1,1))**2
    if (l.eq.0) go to 15

    do 10 k=1,l
    j           = nquad(k,i)
    dist        = (x(1,n)  -x(1,j))**2  +(x(2,n)  -x(2,j))**2
    if (dist.ge.dismin) go to 10
    jnpt        = j
    dismin      = dist
10  continue
    
15  dmin        = sqrt(dismin)
    xl          = x(1,n)  -dmin
    xu          = x(1,n)  +dmin
    yl          = x(2,n)  -dmin
    yu          = x(2,n)  +dmin
20  if (xl.gt.xquad(1,i).and.xu.lt.xquad(2,i).and. yl.gt.yquad(1,i).and.yu.lt.yquad(2,i)) return
    if (i.eq.1) return
    iflag       = i
    i           = nquad(6,i)

!   EXAMINE PARENT OF QUADRANT AND ALL ITS OFFSPRING

    kc          = 1
    ksrch(1)    = i
25  ic          = 0

    do 40 j=1,kc
    ii          = ksrch(j)
    do 40 k=1,4
    itry        = nquad(k,ii)
    if (itry.eq.iflag) go to 40
    if (xl.gt.xquad(2,itry).or.xu.lt.xquad(1,itry).or. yl.gt.yquad(2,itry).or.yu.lt.yquad(1,itry)) go to 40
    l           = nquad(7,itry)
    if (l.ge.0) go to 30
    ic          = ic  +1
    isrch(ic)   = itry
    go to 40
30  if (l.eq.0) go to 40
    do 35 kk=1,l
    jj          = nquad(kk,itry)
    dist        = (x(1,n)  -x(1,jj))**2  +(x(2,n)  -x(2,jj))**2
    if (dist.ge.dismin) go to 35
    jnpt        = jj
    dismin      = dist
35  continue
    dmin        = sqrt(dismin)
    xl          = x(1,n)  -dmin
    xu          = x(1,n)  +dmin
    yl          = x(2,n)  -dmin
    yu          = x(2,n)  +dmin
40  continue
    
    if (ic.gt.msrch) go to 200
    if (ic.eq.0) go to 20
    
    kc          = ic
    do 45 j=1,kc
45  ksrch(j)    = isrch(j)
    go to 25
    
200 write (6,600)
    stop
600 format(5x,'dimension of isrch and ksrch arrays exceeded'/5x,'in routine nearpt. increase size of msrch')

end subroutine
