!-----------------------------------------------------------------!
!                                                                 !
!    READ IN LIST OF NODES DEFINING INNER AND OUTER BOUNDARIES    !
!                                                                 !
!-----------------------------------------------------------------!

SUBROUTINE input(nBNode, nBFace, Bxyz, Bf2n, NodeIndex)

USE limitvars
USE bdyvars
IMPLICIT NONE
INTEGER, INTENT(in   ) :: nBNode, nBFace
REAL(8), INTENT(in   ) :: Bxyz(2,nBNode)
INTEGER, INTENT(in   ) :: Bf2n(2,nBFace)
INTEGER, INTENT(  out) :: NodeIndex(nBNode)
integer :: l, i, j
integer :: node
integer :: nintot

! READ IN DEFINITION OF INNER BOUNDARY
ncmp = 0
if( ncmp .gt. mcmp ) goto 200
nintot = 0
if( ncmp .eq. 0 ) goto 15

do 10 l=1, ncmp
  nint(l) = 0
  if( nint(l) .gt. mbdy ) goto 210
  nintot = nintot + nint(l)
  do 10 i=1,nint(l)
    read(5,*) xint(l,i), yint(l,i)
10 continue

! READ IN DEFINITION OF OUTER BOUNDARY
15 nout = nBNode
if( nout .gt. mbdy ) goto 220
node = Bf2n(1,1)
xout(1) = Bxyz(1,node)
yout(1) = Bxyz(2,node)
NodeIndex(1) = node
node = Bf2n(2,1)
do 20 i=2, nout
  DO j = 1, nBFace
    IF( Bf2n(1,j) == node ) THEN
      xout(i) = Bxyz(1,node)
      yout(i) = Bxyz(2,node)
      NodeIndex(i) = node
      node = Bf2n(2,j)
      EXIT
    END IF
  END DO
20 continue

! OPEN(10,file='./temp/bound.plt')
! DO i = 1, nout
!   WRITE(10,*) xout(i), yout(i)
! END DO
! CLOSE(10)

if( nintot + nout .gt. mnode ) goto 230

return

200 write(6,600)
stop
210 write(6,610)
stop
220 write(6,620)
stop
230 write(6,630)
stop
500 format(1x)
600 format(/5x,'number of interior components exceeds limit mcmp'/5x,'program stopped in routine input')
610 format(/5x,'number of points on interior components ',i3,/5x,'exceeds limit mbdy. program stopped in routine input')
620 format(/5x,'number of points on exterior boundary exceeds'/5x,'limit mbdy. program stopped in routine input')
630 format(/5x,'total number of boundary points exceeds linit mnode.'/5x,'program stopped in routine input')

end subroutine
