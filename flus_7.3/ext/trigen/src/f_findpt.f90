!---------------------------------------------------------------------!
!                                                                     !
!    FIND LIST OF CANDIDATE POINTS THAT FORM VALID TRIANGLES WITH     !
!    EDGE (N1,N2)                                                     !
!                                                                     !
!---------------------------------------------------------------------!

subroutine findpt( n1, n2, iedg )

    use limitvars
    use pntvars
    use edgvars
    use trivars
    use cndvars
    
    implicit none
    integer, intent(in) :: n1, n2, iedg
    integer :: i, l
    integer :: nex(2), nchk, ntry, npiv
    real(8) :: x1, y1, x2, y2, x3, y3, a, atry
    real(8) :: xpiv, ypiv, dir, achk
    real(8) :: f, alpha, beta
    
    itest       = 0
    x1          = x(1,n2)  -x(1,n1)
    y1          = x(2,n2)  -x(2,n1)
    
    do 100 i=1,nnode
        
    if (i.eq.n1.or.i.eq.n2) go to 100
    x2          = x(1,i)   -x(1,n1)
    y2          = x(2,i)   -x(2,n1)
    a           = x1*y2  -x2*y1
    if ( a.lt.1d-10) go to 100
    if (nrept.gt.1) go to 50
    nchk        = 0
    ntry        = nbh(1,iedg)
    x3          = x(1,ntry)  -x(1,n1)
    y3          = x(2,ntry)  -x(2,n1)
    atry        = x1*y3  -x3*y1
    if (atry.lt.1d-10) go to 10
    nchk        = nchk  +1
    nex(nchk)   = 1
10  ntry        = nbh(2,iedg)
    x3          = x(1,ntry)  -x(1,n2)
    y3          = x(2,ntry)  -x(2,n2)
    atry        = x1*y3  -x3*y1
    if (atry.lt.1d-10) go to 20
    nchk        = nchk  +1
    nex(nchk)   = 2

!   CHECK FOR INTERSECTION WITH ADJACENT EDGES

20  if (nchk.eq.0) go to 50
    do 30 l=1,nchk
    ntry        = nbh(nex(l),iedg)
    npiv        = ndg(nex(l),iedg)
    xpiv        = x(1,i)  -x(1,npiv)
    ypiv        = x(2,i)  -x(2,npiv)
    x3          = x(1,ntry)  -x(1,npiv)
    y3          = x(2,ntry)  -x(2,npiv)
    dir         = 3  -2*nex(l)
    achk        = dir*(xpiv*y3  -x3*ypiv)
    if (achk.lt.-1d-10) go to 100
30  continue
    
50  itest       = itest  +1
    if (itest.gt.mtest) go to 100
    ntest(itest) = i
    f           = 1d0/a
    alpha       = 0.5d0*(x1*(x(1,n2)  +x(1,n1))  +y1*(x(2,n2)  +x(2,n1)))
    beta        = 0.5d0*(x2*(x(1,i)   +x(1,n1))  +y2*(x(2,i)   +x(2,n1)))
    xcen(itest) = f*(alpha*y2  -beta*y1)
    ycen(itest) = f*(beta*x1  -alpha*x2)
    
100 continue
    
    if (itest.le.mtest) return
    
    write (6,600)
    stop
600 format(/5x,'itest exceeds limit mtest in routine findpt')
    
end subroutine