!****************************************************************************
!
!  PROGRAM : Check element angle
!
!  PURPOSE : Element내의 각도가 180가 되면 수정
!
!****************************************************************************
subroutine check_ele_angle(anes, ne, node, ele, first_nn, pre_ne)

implicit none

integer, intent(in) :: anes, ne, first_nn, pre_ne, ele(5,anes)
real, intent(inout) :: node(2,anes)

integer :: i, j, j_l, j_r, j_c
real :: alpha, vector(2) 

do i = pre_ne, ne
	do j = 2, 5
		j_l = j - 1;  j_r = j + 1
		if ( j == 2 ) j_l = 5
		if ( j == 5 ) j_r = 2
		
		if ( ele(j,i) > first_nn ) then
			
			call calc_angle( node(:,ele(j_l,i)), node(:,ele(j,i)), node(:,ele(j_r,i)), alpha )
			if ( ( alpha > 175. ) .AND. ( alpha < 185. ) ) then
				j_c = j + 2
				if ( j_c > 5 ) j_c = j_c - 4
				
				vector(:) = node(:,ele(j,i)) - node(:,ele(j_c,i))
				
				node(:,ele(j,i)) = node(:,ele(j,i)) + vector(:) * 0.15
			endif

		endif
	enddo
enddo

end subroutine

		

		

