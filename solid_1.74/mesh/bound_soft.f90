!****************************************************************************
!
!  PROGRAM : Angle check
!
!  PURPOSE : 각도 검사 후 봉합 
!		     - 각도가 0 ~ 60, 290 ~ 360 도 일때 봉합
!			 - 봉합 할 수 있는 점이 없을때까지 반복
!
!****************************************************************************
Subroutine bound_soft(anes, bs, b_num, b, d, ne, nn, ele, node)

implicit none

integer, intent(in) :: anes, bs, ne, nn, ele(5,anes), b_num(bs)
real, intent(in) :: d
real, intent(inout) :: b(2,bs), node(2,anes)

integer :: i, j, k, o, i_l, i_r, ele_l, ele_c, ele_r, ele_sum_max, ele_sum
real :: alpha, vector(2), cross_p(2)

real, parameter :: pi = 3.141592653589793238462643383279502884197169399375

do i = 1, bs
	
	if (i == 1) then
		i_l = bs
		i_r = i + 1
	else if (i == bs) then
		i_l = i - 1
		i_r = 1
	else
		i_l = i - 1
		i_r = i + 1
	endif

	call calc_angle(b(:,i_l), b(:,i), b(:,i_r), alpha)
	if (alpha .GT. 180) then

		ele_l = b_num(i_l)
		ele_c = b_num(i)
		ele_r = b_num(i_r)
		
		ele_sum_max = 0
		do j = 1, ne
			ele_sum = 0
			do k = 2, 5
				if (ele(k,j) == ele_l) then
					ele_sum = ele_sum + 1
				elseif (ele(k,j) == ele_c) then
					ele_sum = ele_sum + 1
				elseif (ele(k,j) == ele_r) then
					ele_sum = ele_sum + 1
				endif
			enddo
			if (ele_sum .GT. ele_sum_max) then
				ele_sum_max = ele_sum
			endif
		enddo

		if (ele_sum_max .NE. 3) then

			alpha = 360. - alpha
			vector(1) = b(1,i) + 10 * d * (cos(alpha*pi/360) * (b(1,i_l) - b(1,i)) - sin(alpha*pi/360) * (b(2,i_l) - b(2,i)))
			vector(2) = b(2,i) + 10 * d * (sin(alpha*pi/360) * (b(1,i_l) - b(1,i)) + cos(alpha*pi/360) * (b(2,i_l) - b(2,i)))

			call calc_cross_p(b(:,i), vector, b(:,i_l), b(:,i_r), cross_p)
			
			b(:,i) = (b(:,i) + cross_p(:)) / 2.

			do o = 1, nn 
				if (o == b_num(i)) then
					node(:,o) = b(:,i)
				endif
			enddo
				
		endif

	endif

enddo

end subroutine bound_soft