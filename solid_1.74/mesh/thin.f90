!******************************************************************
!
!  PROGRAM : make the mesh for a thin object.
!
!  PURPOSE : 얇은 판에서의 요소망 생성
!
!******************************************************************
subroutine mesh_thin(anes, nn, node, thick)

implicit none

integer, intent(in) :: anes, thick
integer, intent(inout) :: nn
real, intent(inout) :: node(anes,3)

integer :: i, j, hp, count
real :: seg 
real :: temp_node(nn*0.5,2), vec(nn*0.5,2)

hp = nn * 0.5
count = 0
do i = hp+1, nn
	count = count + 1
	temp_node(count,1:2) = node(i,2:3)
	vec(count,1:2) = temp_node(count,1:2) - node(count,1:2)
enddo

seg = 1/thick
write (*,*) seg
stop
nn = hp
do i = 1, thick
	do j = 1, hp
		nn = nn + 1
		node(nn,1) = nn
		!node(nn,2:3) = vec(1:2)*(	
	enddo
enddo



end subroutine

