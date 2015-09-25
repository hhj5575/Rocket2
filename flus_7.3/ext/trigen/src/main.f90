!-----------------------------------------------------------------!
!                                                                 !
!    DETERMINES DELAUNAY TRIANGLATION OF A PLANAR SET OF POINTS   !
!    BY A MOVING FRONT METHOD.                                    !
!    THE DATA INPUT CONSISTS OF THE INNER AND OUTER BOUNDARIES.   !
!    EACH DEFINED BY AN ORDERED SEQUENCE OF POINT.                !
!                                                                 !
!-----------------------------------------------------------------!

PROGRAM main

USE TRIGEN_MOD
IMPLICIT NONE
TYPE(TRIGENIO) :: in, out
INTEGER :: i
REAL(8) :: dang, ang
REAL(8), PARAMETER :: pi = 4d0*datan(1d0)

! SAMPLE
in%numberofpoints = 4
ALLOCATE( in%pointlist(2,in%numberofpoints) )
!dang = pi/DBLE(in%numberofpoints)
!DO i = 1, in%numberofpoints
!  ang = dang * DBLE(i-1)
!  in%pointlist(1,i) = DCOS(ang)
!  in%pointlist(2,i) = DSIN(ang)
!END DO

in%pointlist(1,1) = 0d0
in%pointlist(2,1) = 0d0
in%pointlist(1,2) = 1.2d0
in%pointlist(2,2) = 0d0
in%pointlist(1,3) = 1d0
in%pointlist(2,3) = 1d0
in%pointlist(1,4) = 0d0
in%pointlist(2,4) = 1d0

in%numberoffaces = 4
ALLOCATE( in%facelist(2,in%numberoffaces) )
DO i = 1, in%numberoffaces
  in%facelist(1,i) = i
  in%facelist(2,i) = MOD(i,in%numberoffaces) + 1
END DO

in%numberofholes = 0

ALLOCATE( out%pointlist(2,500000) )
ALLOCATE( out%facelist(2,in%numberoffaces) )
ALLOCATE( out%facemarkerlist(in%numberoffaces) )
ALLOCATE( out%celllist(3,1000000) )

CALL TRIANGULATE(in, out)

open(10,file='./mesh.plt')
write(10,*) 'TITLE = "Mesh"'
write(10,*) 'VARIABLES= "X", "Y"'
write(10,*) 'ZONE F=FEPOINT, N=',out%numberofpoints,'E=',out%numberofcells,'ET=TRIANGLE'
do i=1,out%numberofpoints
  write(10,*) out%pointlist(:,i)
END DO
do i=1,out%numberofcells
  write(10,*) out%celllist(:,i)
END DO
close(10)


END PROGRAM
