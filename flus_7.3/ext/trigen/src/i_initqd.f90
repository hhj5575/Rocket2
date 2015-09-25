!---------------------------------------------------------------------!
!                                                                     !
!    INITIALIZES QUADTREE DATA STRUCTURE AFTER BOUNDARY POINTS        !
!    HAVE BEEN TRIANGULATED.                                          !
!                                                                     !
!---------------------------------------------------------------------!

subroutine initqd

    use pntvars
    use qudvars
    use ctlvars
    
    implicit none
    integer :: l, n, jnpt
    real(8) :: xmin, ymin, xmax, ymax
    
    do 10 l=2,6
10  nquad(l,1)  = 0
    nquad(1,1)  = 1
    nquad(7,1)  = 1
    xmin        = x(1,1)
    xmax        = x(1,1)
    ymin        = x(2,1)
    ymax        = x(2,1)

    do 20 n=2,nnode
    xmin        = dmin1(xmin,x(1,n))
    xmax        = dmax1(xmax,x(1,n))
    ymin        = dmin1(ymin,x(2,n))
    ymax        = dmax1(ymax,x(2,n))
20  continue

    xquad(1,1)  = xmin 
    xquad(2,1)  = xmax
    yquad(1,1)  = ymin
    yquad(2,1)  = ymax
    iquad       = 1 
    
    ncode       = 0
    
    do 30 n=2,nnode
    call quad (n,jnpt)
30  continue
    
    return
    
end subroutine