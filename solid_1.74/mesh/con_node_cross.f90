!****************************************************************************
!
!  PROGRAM : Con_node_cross
!
!  PURPOSE : 각 절점에서 연결된 점과 경계선들간의 교차 검사
!			
!****************************************************************************
Subroutine con_node_cross(anes, nn, ne, n_b, node, ele, pbn, me)

implicit none

integer, intent(in) :: anes, pbn, nn, ne, n_b, ele(5,anes)
real, intent(inout) :: me(2,pbn), node(2,anes)

integer :: i, j, k, o, count, k_r, k_l, temp_j, num
integer :: con_node(8,pbn), con_count(pbn), fix_num(3,pbn)
real :: fix_node(2,pbn), fix(2), vector(2)
logical :: repetition, cross
!real ::

con_node = 0
num = 0
do i = 1, pbn
	do j = n_b + 1, nn
		if ( (me(1,i) == node(1,j)) .AND. (me(2,i) == node(2,j)) ) then
			temp_j = j 
			exit
		endif
	enddo
	count = 0
	do j = 1, ne
		do k = 2, 5
			if ( temp_j == ele(k,j) ) then
				k_r = k + 1;  k_l = k - 1
				if ( k == 2 ) k_l = 5
				if ( k == 5 ) k_r = 2 
				repetition = .TRUE.
				do o = 1, count
					if ( con_node(o,i) == ele(k_l,j) ) then 
						repetition = .FALSE.
						exit
					endif
				enddo
				if ( repetition ) then
					do o = 1, nn - n_b
						if ( n_b + o == ele(k_l,j) ) then 
							repetition = .FALSE.
							exit
						endif
					enddo
				endif
				if ( repetition ) then
					count = count + 1
					con_node(count,i) = ele(k_l,j)
				endif
				
				repetition = .TRUE.
				do o = 1, count
					if ( con_node(o,i) == ele(k_r,j) ) then 
						repetition = .FALSE.
						exit
					endif
				enddo
				if ( repetition ) then
					do o = 1, nn - n_b
						if ( n_b + o == ele(k_r,j) ) then 
							repetition = .FALSE.
							exit
						endif
					enddo
				endif
				if ( repetition ) then
					count = count + 1
					con_node(count,i) = ele(k_r,j)
				endif
			endif
		enddo
	enddo
	con_count(i) = count
	!write (*,*) nint(node(temp_j,1)), con_node(i,1:con_count(i))

	do j = 1, con_count(i)
		do k = 1, pbn
			k_r = k + 1
			if ( k == pbn ) k_r = 1
			call calc_cross(node(:,temp_j), node(:,con_node(j,i)), me(:,k), me(:,k_r), cross)

			if ( cross .EQV. .FALSE. ) then
				num = num + 1
				call calc_cross_p(node(:,temp_j), node(:,con_node(j,i)), me(:,k), me(:,k_r), fix)
				fix_num(:,num) = (/ i, temp_j, con_node(j,i) /)
				fix_node(:,num) = fix(:)
				!write (*,*) node(temp_j,1), node(con_node(i,j),1)
				exit
			endif
		enddo
	enddo
enddo

do i = 1, num
	vector(:) = ( node(:,fix_num(3,i)) - fix_node(:,i) ) / 4
	node(:,fix_num(2,i)) = fix_node(:,i) + vector(:)
	me(:,fix_num(1,i)) = node(:,fix_num(2,i))
enddo

end subroutine