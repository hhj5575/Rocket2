!****************************************************************************
!
!  PROGRAM : Calc d
!
!  PURPOSE : 분할 된 영역의 각 절점간의 평균 거리 계산	
!
!****************************************************************************
subroutine calc_d(anes, nn, node, dn, dbs, db, dd)

implicit none

integer, intent(in) :: anes, nn, dn, dbs(1:dn)
integer, intent(in) :: db(1:dn,1:nn)
real, intent(in) :: node(2,anes)
real, intent(inout) :: dd(dn)

integer :: i, j, j_r
real :: len, total_len

do i = 1, dn
	total_len = 0
	do j = 1, dbs(i)
		if ( j == dbs(i) ) then
			j_r = 1
		else
			j_r = j + 1
		endif
		
		call calc_len( node(:,db(i,j)), node(:,db(i,j_r)), len )
		total_len = total_len + len

	enddo
	dd(i) = total_len / dbs(i)
enddo

end subroutine
